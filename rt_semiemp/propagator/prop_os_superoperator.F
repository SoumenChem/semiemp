      subroutine prop_os_superoperator(g_zdensa,g_zdensb,g_zfocka,
     1 g_zfockb,GPARAMA,GPARAMB,nbf,g_zwa,g_zwb,dt,SUMELEL,QQ,LABI,
     1 LABJ,LABK,LABL,LL,NIN,Q,NG,GAMMA,g_zfpa,g_zfpb,g_zcheba,g_zchebb)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
c#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      double precision GPARAMA, GPARAMB !(emax-emin)
      double precision emin
      double precision Aalpha, Balpha !scaling factor
      double precision dt !time step
      double precision Abeta, Bbeta !effective time step
      double precision  QQ(*), SUMELEL(*), Q(*), GAMMA(*)
      integer LABI(*), LABJ(*), LABK(*), LABL(*), NG(*), NIN(*)
      integer g_zdensa, g_zdensb
      integer g_zfocka, g_zfockb
      integer  g_ztempa, g_ztempb
      integer g_zta_n, g_ztb_n   !nth term in Chebyshev polynomial
      integer g_zta_n_1, g_ztb_n_1 !(n-1)th term in Chebyshev polynomial
      integer g_zta_n_2, g_ztb_n_2 !(n-2)th term in Chebyshev polynomial
      integer g_zwa, g_zwb !delta(P0(t))
      integer g_zfpa, g_zfpb ![F0,P0]
      integer i,itag
      integer g_zcheba, g_zchebb !new propagated density 
      integer g_zcheba_n, g_zchebb_n  !nth term in Chebyshev polynomial
      integer g_zcheba_n_1, g_zchebb_n_1 !(n-1)th term in Chebyshev polynomial
      integer IT,nbf,nrr,LL
      double precision TEST !convergence criteria
      double precision Abessfun, Bbessfun !bessel function
      double precision total_wall, total_cpu
      LOGICAL CONVERGED,status
      DATA TEST /1.0D-10/
      CONVERGED = .FALSE.

c     call ga_print(g_zt_n_1)
c     call ga_print(g_zfocka)
      nrr = nbf*(nbf+1)/2
      Aalpha = 1.5*GPARAMA*dt
      Abeta = 1.5*GPARAMA
      Balpha = 1.5*GPARAMB*dt
      Bbeta = 1.5*GPARAMB
      IT = 0
      call ga_zero(g_zcheba)
      call ga_zero(g_zchebb)
C-------------------------------------------------------------
C     1st term 
C------------------------------------------------------------
      call bessel_function(IT,Aalpha,Abessfun)
      call ga_add(z1,g_zcheba,z1*Abessfun,g_zwa,g_zcheba)
      call bessel_function(IT,Balpha,Bbessfun)
      call ga_add(z1,g_zchebb,z1*Bbessfun,g_zwb,g_zchebb)  
C-----------------------------------------------------------    
      IT = 1
      call bessel_function(IT,Aalpha,Abessfun)
      call bessel_function(IT,Balpha,Bbessfun)
      if (.not. ga_duplicate(g_zfocka, g_zcheba_n_1, "ACheb_n_1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_zchebb_n_1, "BCheb_n_1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_zta_n, "alpha zeta n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_ztb_n, "beta zeta n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_zta_n_1, "alpha zeta n-1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_ztb_n_1, "beta zeta n-1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_zta_n_2, "alpha zeta n-2"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_ztb_n_2, "beta zeta n-2"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      call ga_zero(g_zcheba_n_1)
      call ga_zero(g_zchebb_n_1)
      call ga_zero(g_zta_n)
      call ga_zero(g_ztb_n)
      call ga_zero(g_zta_n_1)
      call ga_zero(g_ztb_n_1)
      call ga_zero(g_zta_n_2)
      call ga_zero(g_ztb_n_2)
C-------------------------------------------------------------
C     2nd term
C------------------------------------------------------------
      call superopertor_os(g_zwa,g_zwb,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL
     1 ,NIN,nbf,Q,NG,GAMMA,GPARAMA,GPARAMB,g_zta_n_1,g_ztb_n_1,g_zfpa,
     1 g_zfpb,g_zdensa,g_zdensb,Abeta,Bbeta)

c     call ga_print(g_zcheba_n_1)
c     call ga_print(g_zta_n_1)

      call ga_add(z0,g_zcheba_n_1,Abessfun*z1,g_zta_n_1,g_zcheba_n_1)

      call ga_add(z1,g_zcheba,z1,g_zcheba_n_1,g_zcheba)

      call ga_add(z0,g_zchebb_n_1,Bbessfun*z1,g_ztb_n_1,g_zchebb_n_1)

      call ga_add(z1,g_zchebb,z1,g_zchebb_n_1,g_zchebb)
      
c     call ga_print(g_zcheba)
c     call ga_print(g_zchebb)
C---------------------------------------------------------------
      call ga_copy(g_zwa,g_zta_n_2)
      call ga_copy(g_zwb,g_ztb_n_2)

c      call ga_print(g_zta_n_1)
      if (.not. ga_duplicate(g_zfocka, g_zcheba_n, "alpha Cheb_n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_zchebb_n, "beta Cheb_n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      if (.not. ga_duplicate(g_zfocka, g_ztempa, "Temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfocka, g_ztempb, "Temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      call ga_zero(g_ztempa)
      call ga_zero(g_ztempb)
      call ga_zero(g_zcheba_n)
      call ga_zero(g_zchebb_n)
c     write(6,*)"here 1"
  10  CONTINUE
C-----------------------------------------------------------------
C     Iterations
C-----------------------------------------------------------------
       IT = IT + 1
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c      call ga_print(g_zfpb)
c      call ga_print(g_zfpa)

      call superopertor_os(g_zta_n_1,g_ztb_n_1,SUMELEL,QQ,LABI,LABJ,LABK
     1 ,LABL,LL,NIN,nbf,Q,NG,GAMMA,GPARAMA,GPARAMB,g_ztempa,g_ztempb,
     2 g_zfpa,g_zfpb,g_zdensa,g_zdensb,Abeta,Bbeta)
    
c     call ga_print(g_ztempa)
c     call ga_print(g_ztempb)

       call ga_add(2*z1,g_ztempa,z1,g_zta_n_2,g_zta_n)
       call ga_add(2*z1,g_ztempb,z1,g_ztb_n_2,g_ztb_n)
       call bessel_function(IT,Aalpha,Abessfun) 
       call bessel_function(IT,Balpha,Bbessfun)      

       call ga_add(z0,g_zcheba_n,z1*Abessfun,g_zta_n,g_zcheba_n)
       call ga_add(z0,g_zchebb_n,z1*Bbessfun,g_ztb_n,g_zchebb_n)

       call ga_add(z1,g_zcheba,z1,g_zcheba_n,g_zcheba)
       call ga_add(z1,g_zchebb,z1,g_zchebb_n,g_zchebb)

c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "chebyshev it cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

      if (IT.eq.2) GO TO 20

       if (( mat_is_same (g_zcheba_n_1,g_zcheba_n,test)).and.
     1 (mat_is_same (g_zchebb_n_1,g_zchebb_n,test))) then
        CONVERGED = .true.
       ELSE
        CONVERGED = .false.
       END IF
      
      IF (CONVERGED) GO TO 30

  20  CONTINUE
c     call ga_print(g_ztb_n_1)
c     call ga_print(g_zchebb_n)
      call ga_copy(g_ztb_n_1,g_ztb_n_2)
      call ga_copy(g_zta_n_1,g_zta_n_2)
      call ga_copy(g_ztb_n,g_ztb_n_1)
      call ga_copy(g_zta_n,g_zta_n_1)
      call ga_copy(g_zchebb_n,g_zchebb_n_1)
      call ga_copy(g_zcheba_n,g_zcheba_n_1)
c     write(6,*)"IT",IT
      GO TO 10

  30  CONTINUE
     
c     call ga_print(g_zcheb_n)
c     call ga_print(g_zcheb_n_1)
c     call ga_print(g_zcheb)
c     do i = 0 , nrr-1
c      write(6,*)"Tn-1",tre_n_1(i+1)
c     end do


      if (ga_nodeid().eq.0) write(6,*)"IT",IT
      status =.true.
      status=status.and.ga_destroy(g_ztempb)
      status=status.and.ga_destroy(g_ztempa)
      status=status.and.ga_destroy(g_zchebb_n)
      status=status.and.ga_destroy(g_zcheba_n)
      status=status.and.ga_destroy(g_ztb_n_2)
      status=status.and.ga_destroy(g_zta_n_2)
      status=status.and.ga_destroy(g_ztb_n_1)
      status=status.and.ga_destroy(g_zta_n_1)
      status=status.and.ga_destroy(g_ztb_n)
      status=status.and.ga_destroy(g_zta_n)
      status=status.and.ga_destroy(g_zchebb_n_1)
      status=status.and.ga_destroy(g_zcheba_n_1)
       return
       end


c------------------------------------------------------------------------------------------------------
c
c------------------------------------------------------------------------------------------------------

      subroutine superopertor_os(g_zwa,g_zwb,SUMELEL,QQ,LABI,LABJ,LABK,
     1 LABL,LL,NIN,n,Q,NG,GAMMA,GPARAMA,GPARAMB,g_zsupera,g_zsuperb,
     1 g_zfpa,g_zfpb,g_zdensa,g_zdensb,Abeta,Bbeta)
 
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
c#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"
 
       double precision QQ(*), GAMMA(*), SUMELEL(*), Q(*)
       integer NIN(*), NG(*), LABI(*), LABJ(*), LABK(*), LABL(*)
       integer LL, N
       integer g_zwa,g_zwb !(L^n)W
       integer g_zcomm ![[F,Pnew]-[F0,P0]]
       integer g_zfpa, g_zfpb ![F0,P0]
       integer g_zdensa, g_zdensb !P0
       integer g_zdens
       integer g_ztemp1, g_core 
       integer g_ztemp2
       integer g_zsupera, g_zsuperb !super operator
       integer g_zdensa_new,g_zdensb_new !P0+norm*g_zw*10^-5
       integer g_zwa_scale,g_zwb_scale !norm*g_zw*10^-5
       double precision GPARAMA, GPARAMB, emin  
       double precision Abeta, Bbeta, Anorm, Bnorm
       double precision alpha_scale, beta_scale
       integer g_zfocka, g_zfockb !new Fock matrix
       logical status
       double precision total_wall, total_cpu
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c     write(6,*)"in superoperator_os"

      call ga_zero(g_zsupera)
      call ga_zero(g_zsuperb)

      if (.not. ga_duplicate(g_zdensa, g_zdensa_new, "new density"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zdensb_new, "new density"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zfocka, "new alpha fock"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zfockb, "new beta fock"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_ztemp1, "temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_ztemp2, "temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zcomm, "comm"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zwa_scale, "scaled zwa"))
     $     call errquit ("couldnt duplicate g_zdens", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zwb_scale, "scaled zwb"))
     $     call errquit ("couldnt duplicate g_zdens", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdensa, g_zdens, "total density"))
     $     call errquit ("couldnt duplicate g_zdens", 0, GA_ERR)
      call ga_zero(g_zdensa_new)
      call ga_zero(g_zdensb_new)
      call ga_zero (g_zfocka)
      call ga_zero (g_zfockb)
      call ga_zero (g_ztemp1)
      call ga_zero (g_ztemp2)
      call ga_zero(g_zcomm)
      call ga_zero(g_zwa_scale)
      call ga_zero(g_zwb_scale)
      call ga_zero(g_zdens)


      call ga_norm_infinity(g_zwa,Anorm)
      call ga_copy(g_zwa,g_zwa_scale)

      call ga_norm_infinity(g_zwb,Bnorm)
      call ga_copy(g_zwb,g_zwb_scale)

      alpha_scale = 1.0D0
      beta_scale = 1.0D0

c     write(6,*)"Anorm,Bnorm",Anorm,Bnorm
c     write(6,*)"alpha_scale,beta_scale", alpha_scale,beta_scale

      call ga_scale(g_zwa_scale,alpha_scale)
      call ga_scale(g_zwb_scale,beta_scale)

c     call ga_print(g_zdensa)
c     call ga_print(g_zwa)
c     call ga_print(g_zwa_scale)
  
      call ga_add(z1,g_zdensa,z1,g_zwa_scale,g_zdensa_new)
      call ga_add(z1,g_zdensb,z1,g_zwb_scale,g_zdensb_new)
c     call ga_print(g_zdens)
c     call ga_print(g_zwa)
c----------------------------------------------------------------
C     Build the new Fock matrix
C----------------------------------------------------------------
      call ga_add(z1,g_zdensa_new,z1,g_zdensb_new,g_zdens)

c     call ga_print(g_zdens)
      call semiemp_os_fock_build(g_zdens,g_zdensa_new,SUMELEL,QQ,
     1 LABI,LABJ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfocka,g_zfockb,2,
     2 GPARAMA,GPARAMB)

c     call ga_print(g_zfocka)
c     call ga_print(g_zfockb)
c     call ga_print(g_zdensa_new)
c     call ga_print(g_zdensb_new)
C---------------------------------------------------------------
C--------------------------------------------------------------
C     Constructing super operator
C--------------------------------------------------------------
      call ga_zgemm('N','N',n,n,n,z1,g_zfocka,g_zdensa_new,z0,
     1 g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_zdensa_new,g_zfocka,z0,
     1 g_ztemp2)

c     call ga_print(g_ztemp1)
c     call ga_print(g_ztemp2)

      call ga_add(z1,g_ztemp1,-z1,g_ztemp2,g_zcomm)

c     call ga_print(g_zcomm)
           
      call ga_add(z1,g_zcomm,-z1,g_zfpa,g_zsupera)

      call ga_scale(g_zsupera,-zi/(Abeta*alpha_scale))

c     call ga_print(g_zsupera)
      call ga_zero(g_ztemp1)
      call ga_zero(g_ztemp2)
      call ga_zero(g_zcomm)

      call ga_zgemm('N','N',n,n,n,z1,g_zfockb,g_zdensb_new,z0,
     1 g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_zdensb_new,g_zfockb,z0,
     1 g_ztemp2)

c     call ga_print(g_ztemp1)
c     call ga_print(g_ztemp2)

      call ga_add(z1,g_ztemp1,-z1,g_ztemp2,g_zcomm)

c     call ga_print(g_zcomm)
      call ga_add(z1,g_zcomm,-z1,g_zfpb,g_zsuperb)
      call ga_scale(g_zsuperb,-zi/(Bbeta*beta_scale)) 

c     call ga_print(g_zsuperb)
c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "superoperator cpu time"!, total_cpu, "total time",
c    1 total_wall
c     end if

C--------------------------------------------------------------
      status = .true.
      status=status.and.ga_destroy(g_zdens)
      status=status.and.ga_destroy(g_zwb_scale)
      status=status.and.ga_destroy(g_zwa_scale)
      status=status.and.ga_destroy(g_zcomm)
      status=status.and.ga_destroy(g_ztemp2)
      status=status.and.ga_destroy(g_ztemp1)
      status=status.and.ga_destroy(g_zfockb)
      status=status.and.ga_destroy(g_zfocka)
      status=status.and.ga_destroy(g_zdensb_new)
      status=status.and.ga_destroy(g_zdensa_new)
      
      return 
      end 
