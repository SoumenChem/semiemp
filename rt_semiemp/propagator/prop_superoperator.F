      subroutine prop_superoperator(g_zdens,g_zfock,GPARAM,nbf,emin,g_zw
     1,dt,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,NIN,Q,NG,GAMMA,g_zfp,
     2 g_zcheb,g_core)

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

      double precision GPARAM !(emax-emin)
      double precision emin
      double precision alpha !scaling factor
      double precision dt !time step
      double precision beta !effective time step
      double precision  QQ(*), SUMELEL(*), Q(*), GAMMA(*)
      integer LABI(*), LABJ(*), LABK(*), LABL(*), NG(*), NIN(*)
      integer g_zdens, g_core
      integer g_zfock
      integer  g_ztemp
      integer g_zt_n   !nth term in Chebyshev polynomial
      integer g_zt_n_1 !(n-1)th term in Chebyshev polynomial
      integer g_zt_n_2 !(n-2)th term in Chebyshev polynomial
      integer g_zw !delta(P0(t))
      integer g_zfp ![F0,P0]
      integer i,itag
      integer g_zcheb !new propagated density 
      integer g_zcheb_n  !nth term in Chebyshev polynomial
      integer g_zcheb_n_1 !(n-1)th term in Chebyshev polynomial
      integer IT,nbf,nrr,LL
      integer me
      double precision TEST !convergence criteria
      double precision bessfun !bessel function
      double precision total_wall, total_cpu
      LOGICAL CONVERGED,status
      DATA TEST /1.0D-10/
      CONVERGED = .FALSE.

      me = ga_nodeid()
c     call ga_print(g_zfock)
c     if (me.eq.0) write(6,*)"in prop_superoperator"
      nrr = nbf*(nbf+1)/2
      alpha = 1.5*GPARAM*dt
      beta = 1.5*GPARAM
      IT = 0
      call ga_zero(g_zcheb)
c     call ga_print(g_zw)
C-------------------------------------------------------------
C     1st term 
C------------------------------------------------------------
      call bessel_function(IT,alpha,bessfun)
      call ga_add(z1,g_zcheb,z1*bessfun,g_zw,g_zcheb)  
C-----------------------------------------------------------    
      IT = 1
      call bessel_function(IT,alpha,bessfun)
      if (.not. ga_duplicate(g_zfock, g_zcheb_n_1, "Cheb_n_1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n, "zeta n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n_1, "zeta n-1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n_2, "zeta n-2"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      call ga_zero(g_zcheb_n_1)
      call ga_zero(g_zt_n)
      call ga_zero(g_zt_n_1)
      call ga_zero(g_zt_n_2)
C-------------------------------------------------------------
C     2nd term
C------------------------------------------------------------
c     call ga_print(g_zw)
      call superopertor(g_zw,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,NIN,nbf,Q
     1,NG,GAMMA,GPARAM,g_zt_n_1,g_zfp,g_zdens,beta,g_core)
c     call ga_print(g_zt_n_1)
      call ga_add(z0,g_zcheb_n_1,bessfun*z1,g_zt_n_1,g_zcheb_n_1)

      call ga_add(z1,g_zcheb,z1,g_zcheb_n_1,g_zcheb)
C---------------------------------------------------------------
      call ga_copy(g_zw,g_zt_n_2)

      if (.not. ga_duplicate(g_zfock, g_zcheb_n, "Cheb_n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      if (.not. ga_duplicate(g_zfock, g_ztemp, "Temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
       call ga_zero(g_zcheb_n)
       call ga_zero(g_ztemp)

  10  CONTINUE
C-----------------------------------------------------------------
C     Iterations
C-----------------------------------------------------------------
       IT = IT + 1
c      write(6,*)"IT",IT
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
      call superopertor(g_zt_n_1,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,NIN,
     1 nbf,Q,NG,GAMMA,GPARAM,g_ztemp,g_zfp,g_zdens,beta,g_core)

       call ga_add(2*z1,g_ztemp,z1,g_zt_n_2,g_zt_n)
       call bessel_function(IT,alpha,bessfun)       

       call ga_add(z0,g_zcheb_n,z1*bessfun,g_zt_n,g_zcheb_n)

       call ga_add(z1,g_zcheb,z1,g_zcheb_n,g_zcheb)
c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "chebyshev it cpu time"!, total_cpu, "total time",
c    1 total_wall
c     end if

      if (IT.eq.2) GO TO 20

       if ( mat_is_same (g_zcheb_n_1,g_zcheb_n,test)) then
        CONVERGED = .true.
       ELSE
        CONVERGED = .false.
       END IF
      
      IF (CONVERGED) GO TO 30

  20  CONTINUE
      call ga_copy(g_zt_n_1,g_zt_n_2)
      call ga_copy(g_zt_n,g_zt_n_1)
      call ga_copy(g_zcheb_n,g_zcheb_n_1)
      GO TO 10

  30  CONTINUE
     
c     call ga_print(g_zcheb_n)
c     call ga_print(g_zcheb_n_1)
c     call ga_print(g_zcheb)
c     do i = 0 , nrr-1
c      write(6,*)"Tn-1",tre_n_1(i+1)
c     end do
 
      if (me.eq.0) write(luout,*)"Chebyshev Iteration",IT
      status=.true.
      status=status.and.ga_destroy(g_ztemp)
      status=status.and.ga_destroy(g_zcheb_n)
      status=status.and.ga_destroy(g_zt_n_2)
      status=status.and.ga_destroy(g_zt_n_1)
      status=status.and.ga_destroy(g_zt_n)
      status=status.and.ga_destroy(g_zcheb_n_1)
      if (.not.status) call errquit ("prop_superoperator:destroy failed"
     1 , 0, GA_ERR)
       return
       end


c------------------------------------------------------------------------------------------------------
c
c------------------------------------------------------------------------------------------------------

      subroutine superopertor(g_zw,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,
     1 NIN,n,Q,NG,GAMMA,GPARAM,g_zsuper,g_zfp,g_zdens,beta,g_core)
 
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
       integer g_zw !(L^n)W
       integer g_zcomm ![[F,Pnew]-[F0,P0]]
       integer g_zfp ![F0,P0]
       integer g_zdens !P0
       integer g_ztemp1, g_core 
       integer g_ztemp2
       integer g_zsuper !super operator
       integer g_zdens_new !P0+norm*g_zw*10^-5
       integer g_zw_scale !norm*g_zw*10^-5
       double precision GPARAM, emin, beta, norm, scale, beta_scale
       integer g_zfock !new Fock matrix
       logical status
       double precision total_wall, total_cpu
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c     write(6,*)"i am in superoperator" 
      call ga_zero(g_zsuper)
     
      if (.not. ga_duplicate(g_zdens, g_zdens_new, "new density"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_zfock, "new density"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_ztemp1, "temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_ztemp2, "temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_zcomm, "comm"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_zw_scale, "scaled zw"))
     $     call errquit ("couldnt duplicate g_zdens", 0, GA_ERR)
       
      call ga_zero(g_zdens_new)
      call ga_zero(g_zfock)
      call ga_zero(g_ztemp1)
      call ga_zero(g_ztemp2)       
      call ga_zero(g_zcomm)
      call ga_zero(g_zw_scale)


      call ga_norm_infinity(g_zw,norm)
c     write(6,*)"norm",norm
c     call ga_dgop(374, norm, 1, '+')
      scale = norm*1D-5
      call ga_copy(g_zw,g_zw_scale)
      call ga_scale(g_zw_scale,scale)

      call ga_add(z1,g_zdens,z1,g_zw_scale,g_zdens_new)
c     call ga_print(g_zdens_new)
c     call ga_print(g_zw_scale)
c----------------------------------------------------------------
C     Build the new Fock matrix
C----------------------------------------------------------------
      call semiemp_cs_fock_build(g_zdens_new,SUMELEL,QQ,LABI,LABJ
     1 ,LABK,LABL,LL,NIN,n,Q,NG,GAMMA,g_zfock,2,GPARAM,emin,g_core)
C---------------------------------------------------------------
C--------------------------------------------------------------
C     Constructing super operator
C--------------------------------------------------------------
      call ga_zgemm('N','N',n,n,n,z1,g_zfock,g_zdens_new,z0,
     1 g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_zdens_new,g_zfock,z0,
     1 g_ztemp2)

      call ga_add(z1,g_ztemp1,-z1,g_ztemp2,g_zcomm)
      
      call ga_add(z1,g_zcomm,-z1,g_zfp,g_zsuper)

      beta_scale = beta*scale
      call ga_scale(g_zsuper,-zi/(beta_scale))

c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "superoperator cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

C--------------------------------------------------------------
      status=.true.
      status=status.and.ga_destroy(g_zw_scale)
      status=status.and.ga_destroy(g_zcomm)
      status=status.and.ga_destroy(g_ztemp2)
      status=status.and.ga_destroy(g_ztemp1)
      status=status.and.ga_destroy(g_zfock)
      status=status.and.ga_destroy(g_zdens_new)
      if (.not.status) call errquit ("prop_superoperator:destroy failed"
     1 , 0, GA_ERR)
      
      return 
      end 
