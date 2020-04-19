C> @file rt_tddft_cs_prop.F
C> @ingroup rt_tddft

C--------------------------------------------------------------------
C>  Open-shell propagation driver.
C>
C>  Propagates open shell density matrix using von Neumann
C>  dynamics.
C--------------------------------------------------------------------
      subroutine rt_semiemp_os_cheb(rtdb,g_densre,SUMELEL,QQ,LABI,LABJ
     1,LABK,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,n,RAD,KAN,NA,CO,
     2 NU,NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgsa,g_cgsb,
     3 g_core,g_temp,g_densa,g_densb)
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
      character(*), parameter     :: pname       = "rt_tddft_os_prop: "
      


C     == External ==

      double precision, external :: rt_semiemp_gaussian_window
      
C     == Variables ==
      integer rtdb
      integer ihomo, ilumo
      integer g_densre, g_s, g_zfpa, g_zfpb, g_core
      integer g_densa, g_densb
      integer g_cgsa, g_cgsb
      integer g_zdensa, g_zdensb,  g_zfocka, g_zfockb !in AO basis
      integer g_zdens_ao, g_zfock_ao !in AO basis
      integer g_densre_ao, g_densim_ao
      integer g_densre_mo, g_densim_mo
      integer g_zfock_mo_n12dt  !stored previous F'(t-0.5*dt) values
      integer g_zfock_ao_init, g_ztemp,g_ztemp1, g_ztemp2
      integer g_zcheba, g_zchebb, g_scr, g_zscr, g_temp
      integer g_dipx, g_dipy, g_dipz
      integer g_zscr_ao, g_zscr_mo, g_zscr2_mo, g_zs
      integer g_ft, g_ft_inv, g_zft, g_zwa, g_zwb      
      integer g_zfat, g_zfbt, g_zfat_dt, g_zfbt_dt
      integer g_zpda, g_zpdb
      double precision charge, charge_elec, charge_nuc,fval
      double precision fgauss, frequency, width, center

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
      integer iadipx, ladipx, iadipy, ladipy, iadipz, ladipz
      integer ibdipx, lbdipx, ibdipy, lbdipy, ibdipz, lbdipz
      double precision emin, t, dt, tmax
      double precision GPARAMA, GPARAMB
      integer i, j
      integer LABI(*), LABJ(*), LABK(*), LABL(*), LL, NIN(*), NG(*)
      double precision QQ(*), SUMELEL(*), Q(*), GAMMA(*)
      double precision RAD(*),Z12(*),D12(*),H(*),BB(*),E(3)
      double precision ALPHA(*),CO(*)
      integer NU(*),NW(*),NP(*),KAN(*)      
      integer NEL,NA,ISW1,nrr
      integer g_inv_c
      integer g_zdensa_mo, g_zdensb_mo, g_zca, g_zcb 
C     (time-dependent system properties)
      double precision diff1, diff2, norm1, norm2, norm3, norm4
      double precision ga_trace_diag, pstrace
      double precision tt       ! current time
      double precision cpusec
      double precision Exc(2)   !at most we have two Vxc
      double precision ealpha, ebeta, mult
      double precision emina, eminb
      character*20 tag, spin
      character*18 polarization, field_type

      type(rt_vector_t) dipmom 
      type(rt_quad_t) quadmom

      
C     (profiling)
      double precision elapsed, elapsed_checks, elapsed_props


C     (estimating time remaining)
      double precision est_cputime_start, max, pptrace, avgs2
      integer est_it_start
      logical print_dipole, print_moocc
C     == Variables ==
      character*32 basename
      character*32 entry_name



      
      me = ga_nodeid()

      print_dipole=.false.
      print_moocc=.false.
C     == Initialization ==
      
      write(6,*)"i am in rt_semiemp_os_cheb"
      if (me.eq.0) call util_print_centered (luout,
     $     "Open shell propagation", 40, .true.)
      
       got_orbspace = .false.
       if (.not. rtdb_get (rtdb, "rt_semiemp:got_orbspace",
     $  mt_log, 1, got_orbspace)) call errquit (pname//
     $     "failed to get 'got_orbspace' target into rtdb",0,0)

      
       nrr = n*(n+1)/2
C
C     Memory allocation of time-varying properties.  Recall that all
C     quantites in AO basis are ns_ao x ns_ao, while those in MO basis are
C     ns_mo x ns_mo.
C
c     n = params%ns_ao            !local alias for clarity
c     m = params%ns_mo            !  ""       ""
      ok = .true.


C     (Fock and dens mats in AO and MO basis)
      ok=ok.and.ga_create(mt_dcpl,n,n,"P" , 0,0, g_zdens_ao)
      ok=ok.and.ga_create(mt_dcpl,n,n,"alpha P" , 0,0, g_zdensa)
      ok=ok.and.ga_create(mt_dcpl,n,n,"beta P" , 0,0, g_zdensb)
      ok=ok.and.ga_create(mt_dcpl,n,n,"alpha F" , 0,0, g_zfocka)
      ok=ok.and.ga_create(mt_dcpl,n,n,"beta F" , 0,0, g_zfockb)
      ok=ok.and.ga_create(mt_dcpl,n,n,"Cheby propagator", 0,0, g_zcheba)
      ok=ok.and.ga_create(mt_dbl,n,n,"D.E", 0,0, g_scr)
      ok=ok.and.ga_create(mt_dcpl,n,n,"CMPLX[D.E]", 0,0, g_zscr)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp1)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp2)
      if (got_orbspace) then
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex alpha Coeff", 0,0, g_zca)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex beta Coeff", 0,0, g_zcb)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex mo RA", 0,0, g_zdensa_mo)
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex mo RB", 0,0, g_zdensb_mo)
      end if
      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)

      call ga_zero(g_zdens_ao) 
      call ga_zero (g_zdensa)
      call ga_zero (g_zdensb)
      call ga_zero (g_zfocka)
      call ga_zero (g_zfockb)
      call ga_zero (g_zcheba)
      call ga_zero(g_scr)
      call ga_zero(g_zscr)
      call ga_zero (g_temp)
      call ga_zero (g_ztemp1)
      call ga_zero (g_ztemp2)


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
c     call ga_print(g_zscr)
c      write(6,*)"t,dt",tmax,dt
      t = 0.0d0
      frequency = 0.0
      width = 0.0
      center = 0.0
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



      call convert_d2z (1d0, g_densre, 0d0, g_densre, g_zdens_ao)
      call convert_d2z (1d0, g_densa, 0d0, g_densa, g_zdensa)
      call convert_d2z (1d0, g_densb, 0d0, g_densb, g_zdensb)



      if (.not. rtdb_cget (rtdb, "rt_semiemp:spin",
     $     1, spin)) call errquit (pname//
     $     "failed to get 'spin' target from rtdb",0,0)

      if (inp_compare(.false.,spin,'total'))  then
         ealpha = 1.0
         ebeta  = 1.0
       elseif (inp_compare(.false.,spin,'alpha'))  then
         ealpha = 1.0
         ebeta = 0.0001
       elseif (inp_compare(.false.,spin,'beta'))  then
         ealpha = 0.0001
         ebeta = 1.0
       end if


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

      call semiemp_os_fock_build(g_zdens_ao,g_zdensa,SUMELEL,QQ,LABI,
     1 LABJ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfocka,g_zfockb,1,GPARAMA,
     2 GPARAMB)

      if (inp_compare(.false.,field_type,'delta')) then
      
       call ga_add(z1, g_zfocka, z1, g_zscr, g_zfocka)

       call ga_add(z1, g_zfockb, z1, g_zscr, g_zfockb)
      end if

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
c     call  convert_d2z(1d0,g_inv_c,0d0,g_inv_c,g_zinv_c)
c     call ga_add(z1,g_zfock_ao,z1,g_zscr,g_zfock_ao)
c     call ga_zgemm('N','N',n,n,n,z1,g_zinv_c,g_zdens_ao,z0,
c    1 g_ztemp1)

c     call ga_zgemm('N','N',n,n,n,z1,g_ztemp1,g_zinv_c,z0,
c    1 g_ztemp2)

c     call ga_print(g_ztemp2)
C--------------------------------------------------------------------------
      write(6,*)"here 0"

      if (.not.MA_Push_Get(MT_Dbl,n*n,'ADIPX',ladipx,iadipx))
     &   call errquit('rtsemiemp: cannot allocate adipx',0,
     1    MA_ERR)

       if (.not.MA_Push_Get(MT_Dbl,n*n,'ADIPY',ladipy,iadipy))
     &   call errquit('rtsemiemp: cannot allocate adipy',0,
     1    MA_ERR)

       if (.not.MA_Push_Get(MT_Dbl,n*n,'ADIPZ',ladipz,iadipz))
     &   call errquit('rtsemiemp: cannot allocate adipz',0,
     1    MA_ERR)

      if (got_orbspace) then
      call mo_proj(g_dipx,g_cgsa,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(iadipx),n)


      call mo_proj(g_dipy,g_cgsa,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(iadipy),n)

      call mo_proj(g_dipz,g_cgsa,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(iadipz),n)
      else

      call ga_trimat(g_dipx,Dbl_MB(iadipx),n)
      call ga_trimat(g_dipy,Dbl_MB(iadipy),n)
      call ga_trimat(g_dipz,Dbl_MB(iadipz),n)
      end if

      if (got_orbspace) then
        if (.not.MA_Push_Get(MT_Dbl,n*n,'BDIPX',lbdipx,ibdipx))
     &   call errquit('rtsemiemp: cannot allocate adipx',0,
     1    MA_ERR)

         if (.not.MA_Push_Get(MT_Dbl,n*n,'BDIPY',lbdipy,ibdipy))
     &   call errquit('rtsemiemp: cannot allocate adipy',0,
     1   MA_ERR)

         if (.not.MA_Push_Get(MT_Dbl,n*n,'BDIPZ',lbdipz,ibdipz))
     &   call errquit('rtsemiemp: cannot allocate adipz',0,
     1   MA_ERR)


      call mo_proj(g_dipx,g_cgsb,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(ibdipx),n)


      call mo_proj(g_dipy,g_cgsb,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(ibdipy),n)

      call mo_proj(g_dipz,g_cgsb,n,g_scr)
      call ga_trimat(g_scr,Dbl_MB(ibdipz),n)
       end if
       write(6,*)"here 3" 
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
C      g_zdens_ao)

      if (.not.rtdb_get(rtdb, "rt_semiemp:print_moocc",
     $     mt_log, 1, print_moocc))
     $ continue

       if (.not.rtdb_get(rtdb, "rt_semiemp:print_dipole",
     $     mt_log, 1, print_dipole)) call errquit (pname//
     $     "failed to read print_dipole from rtdb", 0, RTDB_ERR)

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
      if (lrestart) call  rt_semiemp_os_restart_load(it,n,t,g_zdensa,
     $ g_zdensb)
      ihomo = 0
      ilumo = 0

      call ga_add(z1,g_zdensa,z1,g_zdensb,g_ztemp1)

      if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $  call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)

      if (got_orbspace) call rt_activespace_os(rtdb,n)
      if (got_orbspace) then
       call convert_d2z (1d0, g_cgsa, 0d0, g_cgsa, g_zca)
       call convert_d2z (1d0, g_cgsb, 0d0, g_cgsb, g_zcb)
      end if

   50 continue

c     call ga_add(z1,g_zdensa,z1,g_zwa,g_ztemp1)
c     call ga_add(z1,g_zdensb,z1,g_zwb,g_ztemp2)
      if (inp_compare(.false.,field_type,'gaussian')) then
       fgauss = 0.0d0
       fgauss = rt_semiemp_gaussian_window(frequency,
     $ center,width,t)
c      write(6,*)"fgauss=",fgauss
      call ga_add(z1, g_zfocka, z1*fgauss, g_zscr, g_zfocka)
      call ga_add(z1, g_zfockb, z1*fgauss, g_zscr, g_zfockb)
      end if

      call rt_semiemp_avgs2(rtdb,g_zdensa,g_zdensb,g_densa,g_densb,
     1 g_scr,g_zscr,n,nel,t)

      if (print_moocc) then
      write(6,*)"||ALPHA OCCUPATIONS||"
      call convert_z2d(g_zdensa,g_densre,g_scr)
      call rt_semiemp_moocc(rtdb,g_densre,n,g_cgsa,t,1)
      write(6,*)"||BETA OCCUPATIONS||"
      call convert_z2d(g_zdensb,g_densre,g_scr)
      call rt_semiemp_moocc(rtdb,g_densre,n,g_cgsb,t,2)
      end if

      if (print_dipole) then 
      if (got_orbspace) then
      call ga_zero(g_zdensa_mo)
      call ga_zero(g_zdensb_mo)
c     call ga_print(g_ztemp1)
c     call ga_print(g_ztemp2)
      call mo_zproj(g_zdensa,g_zca,n,g_zdensa_mo)
      call mo_zproj(g_zdensb,g_zcb,n,g_zdensb_mo)
c     call ga_print(g_zdensa_mo)
c     call ga_print(g_zdensb_mo)
      call ga_add(z1,g_zdensa_mo,z1,g_zdensb_mo,g_ztemp1)
      call rt_semiemp_os_mo_dipole(rtdb,g_ztemp1,g_zdensa_mo,g_zdensb_mo
     1 ,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,
     1 Dbl_MB(iadipx),Dbl_MB(iadipy),Dbl_MB(iadipz),Dbl_MB(ibdipx),
     1 Dbl_MB(ibdipy),Dbl_MB(ibdipz),tag,t,g_s,g_ft,it,nrestarts,
     1 got_orbspace)
      else
 
      call rt_semiemp_os_dipole(rtdb,g_ztemp1,g_zdensa,g_zdensb,RAD,KAN,
     1 NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,n,NRR,ISW1,Z12,D12,NEL,E,
     1 Dbl_MB(iadipx),Dbl_MB(iadipy),Dbl_MB(iadipz),tag,t,g_s,g_ft,it,
     1 nrestarts,got_orbspace)
       end if 
       end if
      if (mod(it,nrestarts).eq.0) 
     $ call rt_semiemp_os_restart_save(it,n,t,g_zdensa,g_zdensb)


c     if (it.ge.1) then
     
      miter = 0

      call prop_chebyshev(g_zdensa,g_zfocka,GPARAMA,n,emina,g_zcheba,
     1 dt)

c     call mo_zproj(g_zfockb,g_zcb,n,g_zfockb)
c     call mo_zproj(g_zdensb,g_zcb,n,g_zdensb)

      call prop_chebyshev(g_zdensb,g_zfockb,GPARAMB,n,eminb,g_zcheba,
     1 dt)


      call ga_add(z1,g_zdensa,z1,g_zdensb,g_ztemp1)

      call semiemp_os_fock_build(g_ztemp1,g_zdensa,SUMELEL,QQ,LABI,
     1 LABJ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfocka,g_zfockb,2,GPARAMA,
     2 GPARAMB)



      t = t + dt
      it = it + 1      
      if (t.le.tmax) GO TO 50


C     == Cleanup GAs for time-dependent quantities and params ==
      if (got_orbspace) then
       if (.not.ma_pop_stack(lbdipz))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack bdipz',3,
     & MA_ERR)

      if (.not.ma_pop_stack(lbdipy))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack bdipz',3,
     & MA_ERR)

      if (.not.ma_pop_stack(lbdipx))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack bdipx',3,
     & MA_ERR)

       end if

      if (.not.ma_pop_stack(ladipz))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack adipz',3,
     & MA_ERR)

      if (.not.ma_pop_stack(ladipy))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack adipy',3,
     & MA_ERR)

      if (.not.ma_pop_stack(ladipx))
     & call errquit('rt_semiemp_os_cheb:cannot pop stack adipx',3,
     & MA_ERR)

      ok = .true.
      if (got_orbspace) then
      ok = ok .and. ga_destroy(g_zdensb_mo)
      ok = ok .and. ga_destroy(g_zdensa_mo)
      ok = ok .and. ga_destroy(g_zcb)
      ok = ok .and. ga_destroy(g_zca)
      end if
      ok = ok .and. ga_destroy(g_ztemp2)
      ok = ok .and. ga_destroy(g_ztemp1)
      ok = ok .and. ga_destroy(g_ztemp)
c     ok = ok .and. ga_destroy(g_zwb)
c     ok = ok .and. ga_destroy(g_zwa)
c     ok = ok .and. ga_destroy(g_zft)
      ok = ok .and. ga_destroy(g_zscr)
      ok = ok .and. ga_destroy(g_scr)
c     ok = ok .and. ga_destroy(g_zchebb)
      ok = ok .and. ga_destroy(g_zcheba)
      ok = ok .and. ga_destroy(g_zdensb)
      ok = ok .and. ga_destroy(g_zdensa)
      ok = ok .and. ga_destroy(g_zdens_ao)
      ok = ok .and. ga_destroy(g_zfockb)
      ok = ok .and. ga_destroy(g_zfocka)
c     ok = ok .and. ga_destroy(g_cb)
c     ok = ok .and. ga_destroy(g_ca)



      end subroutine
c $Id: rt_tddft_cs_prop.F 28240 2016-05-22 00:18:37Z klopata $
