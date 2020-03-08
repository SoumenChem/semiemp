      subroutine rtsemiemp(rtdb,SH,nb,ipol,FA,SUMELEL,QQ,LABI,LABJ,LABK
     1 ,LABL,LL,Q,NIN,NG,GAMMA,dipx,dipy,dipz,RAD,KAN,NA,CO,NU,NW,NP,
     2 ALPHA,H,BB,NEL,ISW1,Z12,D12,FT,RA,RB)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      integer, intent(in) :: rtdb !< Runtime database
C     == Variables ==
      type(rt_params_t) :: params

      character*80 theory
C     == Parameters ==
      character(*), parameter :: pname =
     $     "rt_semiemp: "

C     == External ==
      integer  ga_create_atom_blocked
      external ga_create_atom_blocked

      logical got_elh
C     == Variables ==
      integer i,j, istart,NEL,ISW1,NA
      double complex zval
      integer ocor
      integer me, ipol
      integer g_s, g_svecs, g_xpad,g_zXpad,g_c
      integer g_zX,g_zY,g_zs,g_densre,g_ft,g_temp,g_core
      integer g_dipx, g_dipy, g_dipz,g_scr,g_cgs
      integer isvals, lsvals,iovl,lovl
      integer g_densa, g_densb
      integer prop_method
      logical atmblk
      integer lseig,iseig,lsvecs,isvecs, inverse
      integer LABI(*), LABJ(*), LABK(*), LABL(*), LL, NIN(*), NG(*)
      double precision mult
      double precision QQ(*), SUMELEL(*), Q(*), GAMMA(*), FT(*)
C     (tmp)
      double precision SH(*),FA(*),RAD(*),Z12(*),D12(*),H(*),BB(*)
      double precision ALPHA(*),CO(*),RA(*),RB(*)
      integer NU(*),NW(*),NP(*),KAN(*)
      double precision dipx(*),dipy(*),dipz(*)
      double precision toll_s, D,pstrace, ga_trace_diag, pptrace
      double precision avgs2
      integer numdep,nb
      external dft_lindep
      integer dft_lindep, nlindmo 
      integer test
C     integer g_scr, g_s12, g_sn12
     

      me = ga_nodeid ()
      got_elh = .false.
C
C     Store rtdb in params struct, then determine the theory we are
C     using.
C
      params%rtdb = rtdb

      if (.not. ga_create(mt_dbl, nb, nb, 'AO ovl',
     &        0, 0, g_s))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'AO Dens',
     &        0, 0, g_densre))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'X DIPOLE',
     &        0, 0, g_dipx))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'Y DIPLOE',
     &        0, 0, g_dipy))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'Z DIPOLE',
     &        0, 0, g_dipz))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'FT',
     &        0, 0, g_ft))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'CORE',
     &        0, 0, g_core))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

            
      call ga_zero(g_s)
      call ga_zero(g_densre)
      call ga_zero(g_dipx)
      call ga_zero(g_dipy)
      call ga_zero(g_dipz)
      call ga_zero(g_ft)
      call ga_zero(g_core)
C
C Put density in ga
C

      if (.not. rtdb_cget(rtdb, "task:theory", 1, theory))
     $     call errquit(pname//"theory not specified",0, RTDB_ERR)

      if (theory .eq. "semiemp") then

         if (ipol.eq.0) then             !closed shell DFT
            if (me.eq.0) then
               write (luout, *) " "
               write (luout, *) "*************************"
               write (luout, *) "* Closed-shell RT-SEMIEMP *"
               write (luout, *) "*************************"
               write (luout, *) " "
               call util_flush (luout)
               call rt_semiemp_cs_init(params)
            endif
         elseif (ipol.eq.1) then         !open shell DFT
            if (me.eq.0) then
               write (luout, *) " "
               write (luout, *) "***********************"
               write (luout, *) "* Open-shell RT-SEMIEMP *"
               write (luout, *) "***********************"
               write (luout, *) " "
               call util_flush (luout)
               call rt_semiemp_os_init(params)
            endif
         endif
      else
         call errquit (pname//"Invalid theory: "//trim(theory), 0, 0)
      endif
 

            if (.not.rtdb_get(rtdb,'rt_semiemp:prop_method',
     $           mt_int,1,prop_method))
     $   prop_method = 0
c         write(6,*)"prop_method",prop_method

             if (.not.MA_Push_Get(MT_Dbl,nb*nb,'svecs',lsvecs,isvecs))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)


      if (.not. ga_create(mt_dbl, nb, nb, 'scr',
     &        0, 0, g_scr))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)


      if (.not. ga_create(mt_dbl, nb, nb, 'temp',
     &        0, 0, g_temp))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, 'ground state coefficient',
     &        0, 0, g_cgs))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      call ga_zero(g_scr)
      call ga_zero(g_temp)
      call ga_zero(g_cgs)

      call ga_sqmat(g_ft,FT,nb,1)
      call ga_sqmat(g_s,SH,nb,1)
c     call ga_sqmat(g_syo,SYO,nb,1)
c     call ga_sqmat(g_core,Q,nb,1)

c     test = ga_spd_invert(g_s)
      if (ipol.eq.0) then

      ocor = nel/2
      if (.not.rtdb_put(rtdb,'rtsemiemp:ocor',mt_int,1,ocor))
     1 call errquit('rtsemiemp: failed to write ocor',0, RTDB_ERR)

c     write(6,*)"test",test
c     call ga_print(g_s)
      if (.not. rtdb_get (rtdb, "rt_semiemp:got_elh",
     $    mt_log, 1, got_elh)) 
     1 continue

c     write(6,*)"got_elh",got_elh
      if (got_elh) then
        call eh_state(rtdb,FA,nb)

       call ga_sqmat(g_densre,FA,nb,1)
      end if
 
      if (.not.rtdb_get(rtdb,'semiemp:coefficient',mt_dbl,nb*nb,
     1 Dbl_MB(isvecs)))
     1call errquit('rtsemiemp: failed to write coefficient',0, RTDB_ERR)
      call ga_put(g_cgs,1,nb,1,nb,Dbl_MB(isvecs),nb)
      end if

      if (.not.got_elh) call ga_sqmat(g_densre,FA,nb,1)

      call ga_dgemm('N','N',nb,nb,nb,1.0d0,g_ft,g_densre,0.0d0,
     1 g_scr)


      call ga_dgemm('N','N',nb,nb,nb,1.0d0,g_scr,g_ft,0.0d0,
     1 g_temp)

      call ga_dgemm('N','N',nb,nb,nb,1.0d0,g_s,g_temp,0.0d0,
     1 g_scr)


C
      if (ga_nodeid().eq.0) write(6,*)"nel",nel
       pstrace=ga_trace_diag(g_scr)
      if (ga_nodeid().eq.0) write(6,*)"TRACE[PS]=", pstrace

      if (ipol.eq.0) then
      call ga_elem_multiply(g_temp,g_s,g_scr)
      call rt_semiemp_mullpop(rtdb,g_scr,0)
      end if

      call ga_sqmat(g_dipx,DIPX,nb,1)
      call ga_sqmat(g_dipy,DIPY,nb,1)
      call ga_sqmat(g_dipz,DIPZ,nb,1)
C     Store parameters and GA handles in params struct.
C
      params%nbf_ao = nb
      params%nbf_mo = nb
C
C     Super matrix size is just nbf_ao or nbf_mo (since this is vanilla
C     closed shell).
C
      params%ns_ao = params%nbf_ao
      params%ns_mo = params%nbf_mo
      params%tol_zero = 1d-8

      call rt_semiemp_print_initialization(params)

      call ga_zero(g_temp)

      if (ipol.eq.0.and.prop_method.eq.0) then             !closed shell DFT
      call rt_semiemp_cs_prop (rtdb,g_densre,SUMELEL,QQ,LABI,LABJ,LABK
     1,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,nb,RAD,KAN,NA,CO,NU,
     2 NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgs,g_core,
     3 g_temp)

       elseif (ipol.eq.0.and.prop_method.eq.2) then
      call rt_semiemp_cs_cheb(rtdb,g_densre,SUMELEL,QQ,LABI,LABJ,LABK
     1,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,nb,RAD,KAN,NA,CO,NU,
     2 NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgs,g_core,
     3 g_temp)



        elseif (ipol.eq.1.and.prop_method.eq.0) then    

      if (.not. ga_create(mt_dbl, nb, nb, 'Beta density',
     &        0, 0, g_densa))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
      if (.not. ga_create(mt_dbl, nb, nb, 'Alpha density',
     &        0, 0, g_densb))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      call ga_zero(g_densa)
      call ga_zero(g_densb)


      if (got_eh) eh_state_u(rtdb,RA,RB,nb)


      call ga_sqmat(g_densa,RA,nb,1)
      call ga_sqmat(g_densb,RB,nb,1)

       if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $       call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)

      call ga_dgemm('N','N',nb,nb,nb,1.0d0,g_densa,g_densb,0.0d0,
     1 g_scr)

       pptrace=ga_trace_diag(g_scr)
       avgs2 = ((mult-1)/2)**2+DBLE(nel)/2-pptrace

      if (ga_nodeid().eq.0) write(6,*)"<S2>", avgs2


      call rt_semiemp_os_prop (rtdb,g_densre,SUMELEL,QQ,LABI,LABJ,LABK
     1,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,nb,RAD,KAN,NA,CO,NU,
     2 NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgs,g_core,
     3 g_temp,g_densa,g_densb)

      if (.not. ga_destroy (g_densa))
     $     call errquit (pname//"Failed to destroy g_densa", 0, 0)
      if (.not. ga_destroy (g_densb))
     $     call errquit (pname//"Failed to destroy g_densb", 0, 0)

          elseif (ipol.eq.1.and.prop_method.eq.2) then

      if (.not. ga_create(mt_dbl, nb, nb, 'Beta density',
     &        0, 0, g_densa))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)
      if (.not. ga_create(mt_dbl, nb, nb, 'Alpha density',
     &        0, 0, g_densb))
     &        call errquit(pname//'Error creating ga',0,GA_ERR)

      call ga_zero(g_densa)
      call ga_zero(g_densb)

      call ga_sqmat(g_densa,RA,nb,1)
      call ga_sqmat(g_densb,RB,nb,1)

       if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $       call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)

      call ga_dgemm('N','N',nb,nb,nb,1.0d0,g_densa,g_densb,0.0d0,
     1 g_scr)

       write(6,*)"nel",nel
       write(6,*)"pstrace",pstrace

       pptrace=ga_trace_diag(g_scr)
       avgs2 = (((mult-1)/2)**2)+(DBLE(nel)/2)-pptrace
       

      if (ga_nodeid().eq.0) write(6,*)"<S2>", avgs2


      call rt_semiemp_os_cheb (rtdb,g_densre,SUMELEL,QQ,LABI,LABJ,LABK
     1,LABL,LL,NIN,Q,NG,GAMMA,g_dipx,g_dipy,g_dipz,nb,RAD,KAN,NA,CO,NU,
     2 NW,NP,ALPHA,H,BB,NEL,ISW1,Z12,D12,g_s,g_ft,g_cgs,g_core,
     3 g_temp,g_densa,g_densb)

      if (.not. ga_destroy (g_densa))
     $     call errquit (pname//"Failed to destroy g_densa", 0, 0)
      if (.not. ga_destroy (g_densb))
     $     call errquit (pname//"Failed to destroy g_densb", 0, 0)
            end if 

      if (.not. ga_destroy (g_cgs))
     $     call errquit (pname//"Failed to destroy g_temp", 0, 0)

      if (.not. ga_destroy (g_temp))
     $     call errquit (pname//"Failed to destroy g_temp", 0, 0)

      if (.not. ga_destroy (g_scr))
     $     call errquit (pname//"Failed to destroy g_scr", 0, 0)

      if (.not. ga_destroy (g_core))
     $     call errquit (pname//"Failed to destroy g_core", 0, 0)

      if (.not. ga_destroy (g_ft))
     $     call errquit (pname//"Failed to destroy g_ft", 0, 0)


      if (.not. ga_destroy (g_dipz))
     $     call errquit (pname//"Failed to destroy g_dipz", 0, 0)

      if (.not. ga_destroy (g_dipy))
     $     call errquit (pname//"Failed to destroy g_dipy", 0, 0)

      if (.not. ga_destroy (g_dipx))
     $     call errquit (pname//"Failed to destroy g_dipx", 0, 0)

      if (.not.ma_pop_stack(lsvecs))
     &   call errquit('rtsemeimp: cannot pop stack svecs',3, MA_ERR) 

      if (.not. ga_destroy (g_densre))
     $     call errquit (pname//"Failed to destroy g_densre", 0, 0)

      if (.not. ga_destroy (g_s))
     $     call errquit (pname//"Failed to destroy g_s", 0, 0)


      return
      END






