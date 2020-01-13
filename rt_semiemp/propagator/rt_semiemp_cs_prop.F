C> @file rt_tddft_cs_prop.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C>  Closed-shell propagation driver.
C>
C>  Propagates closed shell density matrix using von Neumann
C>  dynamics.
C--------------------------------------------------------------------
      subroutine rt_semiemp_cs_prop(rtdb,g_densre,SUMELEL,QQ,LABI,LABJ
     1,LABK,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,n,RAD,KAN,NA,CO,
     2 NU,NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgs,g_c,
     3 g_temp)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
C#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"
#include "inp.fh"
      
C     == Inputs ==
c     type(rt_params_t), intent(inout) :: params  !< struct containing parameters
      type (rt_field_t)  this_field

      type(rt_field_t) field
C     == Parameters ==
      character(*), parameter     :: pname       = "rt_semiemp_cs_prop:"

      double precision, external :: rt_semiemp_gaussian_window
      
C     == Variables ==
      integer rtdb
      integer g_densre, g_s, g_zfp, g_c, g_cgs
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_zdens_mo, g_zfock_mo !in MO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt  !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init, g_ztemp3,g_ztemp1, g_ztemp2
      integer g_zcheb, g_scr, g_zscr, g_temp, g_zc
      integer g_dipx, g_dipy, g_dipz
      integer g_zscr_ao, g_zscr_mo, g_zscr2_mo, g_zs
      integer g_ft, g_ft_inv, g_zft, g_zw      
      double precision charge, charge_elec, charge_nuc,fval
      double precision fgauss

      double precision GPARAM
      type(rt_energies_t) energies
      integer nt, nrestarts
      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      integer it, rest_it, ig, iorb, test
      integer lmoocc, imoocc

      character*60 outstr
      logical lrestart
      logical got_orbspace
      double precision emin, t, dt, tmax
      integer i, j
      integer LABI(*), LABJ(*), LABK(*), LABL(*), LL, NIN(*), NG(*)
      double precision QQ(*), SUMELEL(*), Q(*), GAMMA(*)
      double precision RAD(*),Z12(*),D12(*),H(*),BB(*),E(3)
      double precision ALPHA(*),CO(*)
      integer NU(*),NW(*),NP(*),KAN(*)      
      integer NEL,NA,ISW1,nrr
C     (time-dependent system properties)
      double precision ga_trace_diag, pstrace
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc
      integer orbspace(2)
      integer idipx, ldipx, idipy, ldipy, idipz, ldipz
      integer ihomo,ilumo
      character*20 tag
      character*18 polarization, field_type

      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props
      integer g_dipx_mo, g_dipy_mo, g_dipz_mo

C     (estimating time remaining)
      double precision est_cputime_start, max
      double precision frequency, width, center
      integer est_it_start

C     == Variables ==
      character*32 basename
      character*32 entry_name

C     Properties
      logical print_dipole
      logical print_quad
      logical print_field
      logical print_moocc
      logical print_dipcontribs
      logical print_energy
      logical print_cputime
      logical print_charge
      logical print_convergence
      logical print_s2


      
      me = ga_nodeid()


C     == Initialization ==
      print_dipole = .false.
      print_quad = .false.
      print_field = .false.
      print_moocc = .false.
      print_dipcontribs = .false.
      print_energy = .false.
      print_cputime = .false.
      print_charge = .false.
      print_convergence = .false.
      print_s2 = .false.
      frequency = 0.0d0
      center = 0.0d0
      width = 0.0d0
      fval = 0.0d0 

      if (me.eq.0) call util_print_centered (luout,
     $     "Closed shell propagation", 40, .true.)
      

      
       nrr = n*(n+1)/2
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
c     n = params%ns_ao            !local alias for clarity
c     m = params%ns_mo            !  ""       ""
      call ga_zero(g_temp)
      ok = .true.
      got_orbspace = .false.
C     (Real and imaginary parts of dens mat, used for analysis)
CXXX  [KAL]: remove or unify these???
      ok=ok.and.ga_create(mt_dbl,n,n,"Re[P]" , 0, 0, g_densre_ao)
      ok=ok.and.ga_create(mt_dbl,n,n,"Im[P]" , 0, 0, g_densim_ao)
C     (Fock and dens mats in AO and MO basis)
      ok=ok.and.ga_create(mt_dcpl,n,n,"P" , 0,0, g_zdens_ao)
      ok=ok.and.ga_create(mt_dcpl,n,n,"F" , 0,0, g_zfock_ao)
      ok=ok.and.ga_create(mt_dcpl,n,n,"Cheby propagator", 0,0, g_zcheb)
      ok=ok.and.ga_create(mt_dbl,n,n,"D.E", 0,0, g_scr)
      ok=ok.and.ga_create(mt_dcpl,n,n,"CMPLX[D.E]", 0,0, g_zscr)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex s", 0,0, g_zs)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex FT", 0,0, g_zft)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex W start", 0,0, g_zw)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp1)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp2)
c     ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0,g_ztemp3)
      ok=ok.and.ga_create(mt_dcpl,n,n,"[F0,P0]", 0,0, g_zfp)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex inv c", 0,0, g_zc)
      ok=ok.and.ga_create(mt_dcpl,n,n,"projeted MO", 0,0, g_zdens_mo)
      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)

       call ga_zero(g_densre_ao)
       call ga_zero(g_densim_ao)
       call ga_zero(g_zdens_ao)
       call ga_zero(g_zfock_ao)
       call ga_zero(g_zcheb)
       call ga_zero(g_scr)
       call ga_zero(g_zscr)
       call ga_zero(g_zs)
       call ga_zero(g_zft)
       call ga_zero(g_zw)
       call ga_zero(g_ztemp1)
       call ga_zero(g_ztemp2)
       call ga_zero(g_zfp)
       call ga_zero(g_zc)
       call ga_zero(g_zdens_mo) 

            if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPX',ldipx,idipx))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)

             if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPY',ldipy,idipy))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)

             if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPZ',ldipz,idipz))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)
      



c     write(6,*)"test",test
C     Check that canonical orthogonalization is set up correctly.
C
c     call canorg_check (params, params%g_s)

      if (.not.rtdb_get(rtdb, "rt_semiemp:print_dipole",
     $     mt_log, 1, print_dipole)) call errquit (pname//
     $     "failed to read print_dipole from rtdb", 0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:field_max',mt_dbl,1,fval))
     $       call errquit(pname//'failed to read field_max rtdb',0,
     $ RTDB_ERR)
      
      if (.not. rtdb_cget (rtdb, "rt_semiemp:polarization",
     $     1, polarization)) call errquit (pname//
     $     "failed to get 'polarization' target into rtdb",0,0)

      if (.not. rtdb_cget (rtdb, "rt_semiemp:field_type",
     $     1, field_type)) call errquit (pname//
     $     "failed to get 'field_type' target into rtdb",0,0)
      
            if ( (field_type.ne."delta").and.
     $         (field_type.ne."gaussian") )
     $      call errquit (pname//"invalid field type: "//field_type,0,0)


      if (inp_compare(.false.,field_type,'delta').or.
     $ inp_compare(.false.,field_type,'gaussian')) then

      if (inp_compare(.false.,polarization,'x'))  then

       call ga_add(fval,g_dipx,z0,g_scr,g_scr)

       elseif (inp_compare(.false.,polarization,'y'))  then

       call ga_add(fval,g_dipy,z0,g_scr,g_scr)

       elseif (inp_compare(.false.,polarization,'z'))  then
       call ga_add(fval,g_dipz,z0,g_scr,g_scr)
       
       end if
      end if

      if (inp_compare(.false.,field_type,'gaussian'))
     $ call errquit(pname//' gaussian pulse is not appropriate for this
     $propagator',0,
     $ RTDB_ERR)
      
      if (inp_compare(.false.,field_type,'gaussian')) then
       if (.not.rtdb_get(rtdb,'rt_semiemp:frequency',mt_dbl,1,
     $   frequency))
     $       call errquit(pname//'failed to read frequency rtdb',0,
     $ RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:center',mt_dbl,1,
     $   center))
     $       call errquit(pname//'failed to read center rtdb',0,
     $ RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:width',mt_dbl,1,
     $   width))
     $       call errquit(pname//'failed to read width rtdb',0,
     $ RTDB_ERR)
       
c      write(6,*)"frequency, center, width =",frequency, center,width

      end if 

c     call ga_scale(g_scr,-1.0d0)
      if (.not. rtdb_get (rtdb, "rt_semiemp:got_orbspace",
     $    mt_log, 1, got_orbspace)) call errquit (pname//
     $     "failed to get 'got_orbspace' target into rtdb",0,0)

      if (got_orbspace) call rt_activespace(rtdb,n)
      call convert_d2z(z1,g_scr,z0,g_scr,g_zscr)

c     call ga_print(g_dipx)
c     call ga_print(g_dipy)
c     call ga_print(g_dipz)
c     call ga_print(g_zscr)
C     Get starting ("ground state") movecs and density matrix in AO
C     basis (note 1 for closed shell). Then convert to dens mat to MO
C     basis.
C
C      call rt_tddft_movecs_zdens (params, 1, g_movecs_ao_gs, g_zdens_ao)

      if (.not. rtdb_cget (rtdb, "rt_semiemp:tag",
     $     1, tag)) call errquit (pname//
     $     "failed to get 'tag' target into rtdb",0,0)


       if (.not.rtdb_get(rtdb,'rt_semiemp:tmax',mt_dbl,1,tmax))
     $       call errquit(pname//'failed to read tmax rtdb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:dt',mt_dbl,1,dt))
     $       call errquit(pname//'failed to read dt rtdb',0,RTDB_ERR)

      t = 0.0d0
c     dt = 0.1
c     tmax = 0.1  

      call convert_d2z (1d0, g_densre, 0d0, g_densre, g_zdens_ao)

      if (inp_compare(.false.,field_type,'delta')) then
      call ga_zgemm('N','N',n,n,n,z1,g_zscr,g_zdens_ao,z0,
     1 g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_zdens_ao,g_zscr,z0,
     1 g_ztemp2)

      call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zw)

      call ga_scale(g_zw, -zi)
      end if 

      call semiemp_cs_fock_build(g_zdens_ao,SUMELEL,QQ,LABI,LABJ
     1 ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfock_ao,1,GPARAM,emin,g_c)

       
C------------------------------------------------------------------
c     AO TO MO transformation
c     call ga_dgemm('N','N',n,n,n,1d0,g_ft,g_core,0d0,g_inv_c)
c     call ga_print(g_ft)
c     call ga_print(g_core)
c     call ga_print(g_zdens_mo)

c     test = ga_spd_invert(g_core)
c     write(6,*)"test",test
      call convert_d2z(1d0,g_cgs,0d0,g_cgs,g_zc)
c     call convert_d2z(1d0,g_dipz,0d0,g_dipz,g_ztemp1)
c     call ga_print(g_zc)
c     call mo_zproj(g_zdens_ao,g_zc,n,g_zdens_mo) 
      if (got_orbspace) then

      call mo_proj(g_dipx,g_cgs,n,g_scr)      
      call ga_trimat(g_scr,Dbl_MB(idipx),n)

      
      call mo_proj(g_dipy,g_cgs,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(idipy),n)

      call mo_proj(g_dipz,g_cgs,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(idipz),n)
      else
 
      call ga_trimat(g_dipx,Dbl_MB(idipx),n)
      call ga_trimat(g_dipy,Dbl_MB(idipy),n)
      call ga_trimat(g_dipz,Dbl_MB(idipz),n)
      end if
c     call ga_print(g_zdens_mo)
c     call convert_d2z(1d0,g_inv_c,0d0,g_inv_c,g_zinv_c)
c     call convert_d2z(1d0,g_temp,0d0,g_temp,g_ztemp1)
c     call ga_add(z1,g_zfock_ao,z1,g_zscr,g_zfock_ao)
c     call ga_zgemm('T','N',n,n,n,z1,g_zc,g_ztemp1,z0,
c    1 g_ztemp2)

c     call ga_zgemm('N','N',n,n,n,z1,g_ztemp2,g_zc,z0,
c    1 g_ztemp1)

c     call ga_print(g_ztemp1)
C--------------------------------------------------------------------------
      call ga_zgemm('N','N',n,n,n,z1,g_zfock_ao,g_zdens_ao,z0,
     1 g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_zdens_ao,g_zfock_ao,z0,
     1 g_ztemp2)

      call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zfp)

c     call ga_print(g_zfp)
c     call ga_print(g_zw)
c     call ga_print(g_zfock_ao)
      call ga_zero(g_ztemp1)
      call convert_d2z (1d0, g_s, 0d0, g_s, g_zs)

            if (.not.rtdb_get(rtdb,'rt_semiemp:nrestarts',
     $           mt_int,1,nrestarts))
     $           call errquit(pname//'Read failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_get (rtdb, "rt_semiemp:restart",
     $     mt_log, 1, lrestart))
     $     lrestart = .false.

      nt = tmax/dt
      it = 0
      if (lrestart) call  rt_semiemp_restart_load(it,n,t,g_zw)
   10 continue
      call ga_zero(g_ztemp1)
      call ga_zero(g_ztemp2)
c     call ga_zero(g_ztemp3)
c     call ga_print(g_zw)
c     call ga_print(g_zdens_ao)
      call ga_add(z1,g_zdens_ao,z1,g_zw,g_ztemp1)

c     if (inp_compare(.false.,field_type,'gaussian')) then
c      fgauss = 0.0d0       
c      fgauss = rt_semiemp_gaussian_window(frequency,
c    $ center,width,t)
c     write(6,*)"fgauss",fgauss
c     call ga_add(z1*fgauss,g_zscr,z0,g_zscr,g_zscr)

c     call ga_zgemm('N','N',n,n,n,z1,g_zscr,g_ztemp1,z0,
c    1 g_ztemp2)

c     call ga_zgemm('N','N',n,n,n,z1,g_ztemp1,g_zscr,z0,
c    1 g_ztemp3)

c     call ga_add(z1, g_ztemp2, -z1, g_ztemp3, g_zw)

c     call ga_scale(g_zw, -zi)
        
c      end if


      if (print_dipole) then
       if (got_orbspace) then
        call ga_zero(g_zdens_mo)
       call mo_zproj(g_ztemp1,g_zc,n,g_zdens_mo)

      call rt_semiemp_dipole(rtdb,g_zdens_mo,RAD,KAN,NA,CO,NU,NW,NP,NIN,
     1 ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,Dbl_MB(idipx),Dbl_MB(idipy),
     1 Dbl_MB(idipz),tag,t,g_s,g_ft,it,nrestarts,got_orbspace)
       else
c     call ga_print(g_ztemp1)
      call rt_semiemp_dipole(rtdb,g_ztemp1,RAD,KAN,NA,CO,NU,NW,NP,NIN,
     1 ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,Dbl_MB(idipx),Dbl_MB(idipy),
     1 Dbl_MB(idipz),tag,t,g_s,g_ft,it,nrestarts,got_orbspace)
       end if
      end if

      if (mod(it,nrestarts).eq.0) 
     $ call rt_semiemp_restart_save(it,n,t,g_zw)
      
c     call ga_print(g_zw)
c     call ga_print(g_zdens_ao)
      call prop_superoperator(g_zdens_ao,g_zfock_ao,GPARAM,n,emin,g_zw,
     1 dt,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,NIN,Q,NG,GAMMA,g_zfp,
     2 g_zcheb,g_c)
c     call prop_chebyshev(g_zdens_ao,g_zfock_ao,GPARAM,n,emin,g_zcheb,
c    1 dt)
c     call ga_zero (g_ztemp)
c     call ga_print(g_zcheb)
      call ga_copy(g_zcheb,g_zw)

c     call semiemp_cs_fock_build(g_zdens_ao,SUMELEL,QQ,LABI,LABJ
c    1 ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfock_ao,2,GPARAM,emin,g_zfp)
      t = t + dt
      it = it + 1      
      if (t.lt.tmax) GO TO 10



C     == Cleanup GAs for time-dependent quantities and params ==
      if (.not.ma_pop_stack(ldipz))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      if (.not.ma_pop_stack(ldipy))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      if (.not.ma_pop_stack(ldipx))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      ok = .true.
      ok = ok .and. ga_destroy(g_zdens_mo)
      ok = ok .and. ga_destroy(g_zc)
      ok = ok .and. ga_destroy(g_zfp)
c     ok = ok .and. ga_destroy(g_ztemp3)
      ok = ok .and. ga_destroy(g_ztemp2)
      ok = ok .and. ga_destroy(g_ztemp1)
      ok = ok .and. ga_destroy(g_zw)
      ok = ok .and. ga_destroy(g_zft)
      ok = ok .and. ga_destroy(g_zs)
      ok = ok .and. ga_destroy(g_zscr)
      ok = ok .and. ga_destroy(g_scr)
      ok = ok .and. ga_destroy(g_zcheb)
      ok = ok .and. ga_destroy(g_zfock_ao)
      ok = ok .and. ga_destroy(g_zdens_ao)
      ok = ok .and. ga_destroy(g_densim_ao)
      ok = ok .and. ga_destroy(g_densre_ao)



      end subroutine
c $Id: rt_tddft_cs_prop.F 28240 2016-05-22 00:18:37Z klopata $
