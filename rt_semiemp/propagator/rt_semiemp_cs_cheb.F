C> @file rt_tddft_cs_prop.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C>  Open-shell propagation driver.
C>
C>  Propagates open shell density matrix using von Neumann
C>  dynamics.
C--------------------------------------------------------------------
      subroutine rt_semiemp_cs_cheb(rtdb,g_densre,SUMELEL,LL,NIN,Q,
     1  NG,GAMMA,g_dipx,g_dipy,g_dipz,n,RAD,KAN,NA,CO,
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
      character(*), parameter     :: pname       = "rt_tddft_cs_cheb: "
      
      double precision, external :: rt_semiemp_gaussian_window


C     == External ==
      logical, external  :: rt_tddft_atpoint
      logical, external  :: rt_tddft_at_snap      
      logical, external  :: rt_tddft_cs_tdfock
      logical, external  :: rt_tddft_checkpoint


      
C     == Variables ==
      integer rtdb
      integer ihomo, ilumo
      integer g_densre, g_s, g_zfpa, g_zfpb, g_c
      integer g_densa, g_densb
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt  !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init, g_ztemp,g_ztemp1, g_ztemp2
      integer g_zfock
      integer g_zcheb, g_scr, g_zscr, g_temp, g_zc, g_zcb, g_temp1
      integer g_dipx, g_dipy, g_dipz
      integer g_zscr_ao, g_zscr_mo, g_zscr2_mo, g_zs
      integer g_ft, g_cgs, g_zft, g_zwa, g_zwb      
      integer g_zfat, g_zfbt, g_zfat_dt, g_zfbt_dt
      integer g_zpda, g_zpdb
      double precision charge, charge_elec, charge_nuc,fval
      double precision fgauss, frequency, width, center
      integer g_dipx_mo, g_dipy_mo, g_dipz_mo
      double precision GPARAM
      type(rt_energies_t) energies
      integer nt, nrestarts
      integer n, m              !aliases for ns_ao, ns_mo
      integer me
      logical ok
      logical got_orbspace
      integer it, rest_it, ig, iorb, test
      integer lmoocc, imoocc
      integer miter
      character*60 outstr
      logical lrestart
      integer idipx, ldipx, idipy, ldipy, idipz, ldipz

      double precision emin, t, dt, tmax
      double precision GPARAMA, GPARAMB
      integer i, j
      integer  LL, NIN(*), NG(*)
      double precision  SUMELEL(*), Q(*), GAMMA(*)
      double precision RAD(*),Z12(*),D12(*),H(*),BB(*),E(3)
      double precision ALPHA(*),CO(*)
      integer NU(*),NW(*),NP(*),KAN(*)      
      integer NEL,NA,ISW1,nrr
      integer g_inv_c 
C     (time-dependent system properties)
      double precision diff1, diff2, norm1, norm2, norm3, norm4
      double precision ga_trace_diag, pstrace
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc
      double precision ealpha, ebeta, mult
      double precision emina, eminb
      double precision energy
      character*20 tag, spin
      character*18 polarization, field_type

      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start, max, pptrace, avgs2
      integer est_it_start
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

C     == Variables ==
      character*32 basename
      character*32 entry_name


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
      got_orbspace = .false.
      
      me = ga_nodeid()


C     == Initialization ==
      
c     write(6,*)"i am in rt_semiemp_cs_cheb"
      if (me.eq.0) call util_print_centered (luout,
     $     "Closed shell propagation", 40, .true.)
      

      
       nrr = n*(n+1)/2
      if (.not. rtdb_get (rtdb, "rt_semiemp:got_orbspace",
     $    mt_log, 1, got_orbspace)) call errquit (pname//
     $     "failed to get 'got_orbspace' target into rtdb",0,0)

C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
c     n = params%ns_ao            !local alias for clarity
c     m = params%ns_mo            !  ""       ""
      call ga_zero(g_temp)
      ok = .true.

c     ok=ok.and.ga_create(mt_dbl,n,n,"c for alpha" , 0, 0, g_ca)
c     ok=ok.and.ga_create(mt_dbl,n,n,"c for beta" , 0, 0, g_cb)

C     (Fock and dens mats in AO and MO basis)
      ok=ok.and.ga_create(mt_dcpl,n,n,"P" , 0,0, g_zdens_ao)
      ok=ok.and.ga_create(mt_dcpl,n,n,"F" , 0,0, g_zfock_ao)
      ok=ok.and.ga_create(mt_dcpl,n,n,"Cheby propagator", 0,0, g_zcheb)
      ok=ok.and.ga_create(mt_dbl,n,n,"D.E", 0,0, g_scr)
      ok=ok.and.ga_create(mt_dcpl,n,n,"CMPLX[D.E]", 0,0, g_zscr)
      ok=ok.and.ga_create(mt_dbl,n,n,"complex temp", 0,0, g_temp1)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp2)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex c ", 0,0, g_zc)
      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)

      call ga_zero(g_zdens_ao) 
      call ga_zero (g_zfock_ao)
      call ga_zero(g_zcheb)
      call ga_zero(g_scr)
      call ga_zero(g_zscr)
      call ga_zero(g_temp1)
      call ga_zero(g_ztemp2)
      call ga_zero(g_zc)


       if (.not.rtdb_get(rtdb,'rt_semiemp:field_max',mt_dbl,1,fval))
     $       call errquit(pname//'failed to read field_max rtdb',0,
     $ RTDB_ERR)
      
      if (.not. rtdb_cget (rtdb, "rt_semiemp:polarization",
     $     1, polarization)) call errquit (pname//
     $     "failed to get 'polarization' target into rtdb",0,0)

      if (.not. rtdb_cget (rtdb, "rt_semiemp:field_type",
     $     1, field_type)) call errquit (pname//
     $     "failed to get 'field_type' target into rtdb",0,0)
      
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

      call ga_scale(g_scr,-1.0d0)

      call convert_d2z(z1,g_scr,z0,g_scr,g_zscr)

       frequency = 0.0d0
       width = 0.0d0
       center = 0.0d0
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

      t = 0.0d0



      call convert_d2z (1d0, g_densre, 0d0, g_densre, g_zdens_ao)


c     call ga_zgemm('N','N',n,n,n,z1,g_zscr,g_zdensa,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_zdensa,g_zscr,z0,
c    1 g_ztemp2)

c     call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zwa)

c     call ga_zgemm('N','N',n,n,n,z1,g_zscr,g_zdensb,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_zdensb,g_zscr,z0,
c    1 g_ztemp2)

c     call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zwb)

c     write(6,*)"ealpha,ebeta",ealpha,ebeta
c     call ga_scale(g_zwa, -zi*ealpha)

c     call ga_scale(g_zwb, -zi*ebeta)

c     call ga_print(g_ca)
c     call ga_print(g_zwb)
      call semiemp_cs_fock_build(g_zdens_ao,SUMELEL,LL,NIN,n,Q,NG,GAMMA,
     1 g_zfock_ao,1,GPARAM,emin,g_c)

       if (inp_compare(.false.,field_type,'delta')) 
     1  call ga_add(z1, g_zfock_ao, z1, g_zscr, g_zfock_ao)


c     call convert_d2z (1d0, g_ca, 0d0, g_ca, g_zca)
c     call convert_d2z (1d0, g_cb, 0d0, g_cb, g_zcb)
C------------------------------------------------------------------
c     AO TO MO transformation
c     call ga_dgemm('N','N',n,n,n,1.0,g_ft,g_core,0.0d0,g_inv_c)
c     call ga_print(g_ft)
c     call ga_print(g_core)
c     call ga_print(g_inv_c)

c     test = ga_spd_invert(g_inv_c)

c     call ga_print(g_inv_c)
      call  convert_d2z(1d0,g_cgs,0d0,g_cgs,g_zc)
c     call ga_add(z1,g_zfock_ao,z1,g_zscr,g_zfock_ao)
c     call ga_zgemm('N','N',n,n,n,z1,g_zinv_c,g_zdens_ao,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_ztemp1,g_zinv_c,z0,
c    1 g_ztemp2)

c     call ga_print(g_ztemp2)
C--------------------------------------------------------------------------
            if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPX',ldipx,idipx))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)

             if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPY',ldipy,idipy))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)

             if (.not.MA_Push_Get(MT_Dbl,n*n,'DIPZ',ldipz,idipz))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)



      if (got_orbspace) then

      call mo_proj(g_dipx,g_cgs,n,g_temp1)
      call ga_trimat(g_temp1,Dbl_MB(idipx),n)
     

      call mo_proj(g_dipy,g_cgs,n,g_temp1)
      call ga_trimat(g_temp1,Dbl_MB(idipy),n)


      call mo_proj(g_dipz,g_cgs,n,g_temp1)
      call ga_trimat(g_temp1,Dbl_MB(idipz),n)

      else

      call ga_trimat(g_dipx,Dbl_MB(idipx),n)
      call ga_trimat(g_dipy,Dbl_MB(idipy),n)
      call ga_trimat(g_dipz,Dbl_MB(idipz),n)
      end if

c     call ga_zgemm('N','N',n,n,n,z1,g_zfocka,g_zdensa,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_zdensa,g_zfocka,z0,
c    1 g_ztemp2)

c     call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zfpa)

c     call ga_zgemm('N','N',n,n,n,z1,g_zfockb,g_zdensb,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_zdensb,g_zfockb,z0,
c    1 g_ztemp2)

c     call ga_add(z1, g_ztemp1, -z1, g_ztemp2, g_zfpb)

c     call ga_print(g_zfpa)
c     call ga_print(g_zfpb) 
C-----------------------------------------------------------------------------------
c     call ga_zero(g_ztemp1)
c     call convert_d2z (1d0, g_s, 0d0, g_s, g_zs)
C     Get starting ("ground state") movecs and density matrix in AO
C     basis (note 1 for closed shell). Then convert to dens mat to MO
C     basis.
C
C      call rt_tddft_movecs_zdens (params, 1, g_movecs_ao_gs,
c     write(6,*)"got_orbspace",got_orbspace
      if (got_orbspace) call rt_activespace(rtdb,n)


      if (.not.rtdb_get(rtdb, "rt_semiemp:print_moocc",
     $     mt_log, 1, print_moocc))
     $ continue

      if (.not.rtdb_get(rtdb, "rt_semiemp:print_dipole",
     $     mt_log, 1, print_dipole)) call errquit (pname//
     $     "failed to read print_dipole from rtdb", 0, RTDB_ERR)

      if (.not.rtdb_get(rtdb, "rt_semiemp:print_energy",
     $     mt_log, 1, print_energy)) 
     $     continue

      if (.not.rtdb_get(rtdb, "rt_semiemp:print_charge",
     $     mt_log, 1, print_charge))
     $     continue


      if (.not. rtdb_cget (rtdb, "rt_semiemp:tag",
     $     1, tag)) call errquit (pname//
     $     "failed to get 'tag' target into rtdb",0,0)


       if (.not.rtdb_get(rtdb,'rt_semiemp:tmax',mt_dbl,1,tmax))
     $       call errquit(pname//'failed to read tmax rtdb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:dt',mt_dbl,1,dt))
     $       call errquit(pname//'failed to read dt rtdb',0,RTDB_ERR)

            if (.not.rtdb_get(rtdb,'rt_semiemp:nrestarts',
     $           mt_int,1,nrestarts))
     $           call errquit(pname//'Read failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_get (rtdb, "rt_semiemp:restart",
     $     mt_log, 1, lrestart))
     $     lrestart = .false.

      nt = tmax/dt
      it = 0
      if (lrestart) call  rt_semiemp_restart_load(it,n,t,g_zdens_ao)
      ihomo = 0
      ilumo = 0

      if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $       call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)


      energy = 0d0
   50 continue

       if (inp_compare(.false.,field_type,'gaussian')) then
       fgauss = 0.0d0
       fgauss = rt_semiemp_gaussian_window(frequency,
     $ center,width,t)
c      write(6,*)"fgauss=",fgauss
      call ga_add(z1, g_zfock_ao, z1*fgauss, g_zscr, g_zfock_ao)
      end if 
      
c     call ga_add(z1,g_zdensa,z1,g_zwa,g_ztemp1)
c     call ga_add(z1,g_zdensb,z1,g_zwb,g_ztemp2)

c     call rt_semiemp_avgs2(rtdb,g_zdensa,g_zdensb,g_densa,g_densb,
c    1 g_scr,g_zscr,n,nel,t)



      if (print_moocc) then
      call convert_z2d(g_zdens_ao,g_densre,g_scr)
      call rt_semiemp_moocc(rtdb,g_densre,n,g_cgs,t,1)
      end if

      if (print_dipole) then
       if (got_orbspace) then
        call mo_zproj(g_zdens_ao,g_zc,n,g_ztemp2)
 
      call rt_semiemp_dipole(rtdb,g_ztemp2,RAD,KAN,NA,CO,NU,NW,NP,NIN,
     1 ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,Dbl_MB(idipx),Dbl_MB(idipy),
     1 Dbl_MB(idipz),tag,t,g_s,g_ft,it,nrestarts,got_orbspace)
       else
      call rt_semiemp_dipole(rtdb,g_zdens_ao,RAD,KAN,NA,CO,NU,NW,NP,NIN,
     1 ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,Dbl_MB(idipx),Dbl_MB(idipy),
     1 Dbl_MB(idipz),tag,t,g_s,g_ft,it,nrestarts,got_orbspace)
       end if
      end if

      if (print_charge) then
      call rt_semiemp_charge(rtdb,g_zdens_ao,g_s,g_ft,n,t)
      end if
 
      if (mod(it,nrestarts).eq.0) 
     $ call rt_semiemp_restart_save(it,n,t,g_zdens_ao)

      if (print_energy) then
      call convert_z2d(g_zfock_ao,g_temp,g_scr)

      call convert_z2d(g_zdens_ao,g_temp1,g_scr)

      call ga_zero(g_scr)

      call ga_sqmat(g_scr,Q,n,1)

      call ga_add(1.0d0,g_scr,1.0d0,g_temp,g_scr)

      energy = 0.5*ga_ddot(g_temp1,g_scr)

      if (me.eq.0) write(6,*)"t,energy",t,energy   
      end if  
      miter = 0
    

      call prop_chebyshev(g_zdens_ao,g_zfock_ao,GPARAM,n,emin, 
     1 g_zcheb,dt)

      call semiemp_cs_fock_build(g_zdens_ao,SUMELEL,LL,NIN,n,Q,NG,GAMMA,
     1  g_zfock_ao,2,GPARAM,emin,g_c)

      t = t + dt
      it = it + 1      
      if (t.le.tmax) GO TO 50



C     == Cleanup GAs for time-dependent quantities and params ==
      if (.not.ma_pop_stack(ldipz))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      if (.not.ma_pop_stack(ldipy))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      if (.not.ma_pop_stack(ldipx))
     & call errquit('rt_semiemp_dipole:cannot pop stack dens',3, MA_ERR)

      ok = .true.
      ok = ok .and. ga_destroy(g_zc)
      ok = ok .and. ga_destroy(g_ztemp2)
      ok = ok .and. ga_destroy(g_temp1)
      ok = ok .and. ga_destroy(g_zscr)
      ok = ok .and. ga_destroy(g_scr)
      ok = ok .and. ga_destroy(g_zcheb)
      ok = ok .and. ga_destroy(g_zfock_ao)
      ok = ok .and. ga_destroy(g_zdens_ao)



      end subroutine
c $Id: rt_tddft_cs_prop.F 28240 2016-05-22 00:18:37Z klopata $
