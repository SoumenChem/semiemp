      subroutine prop_chebyshev(g_zdens,g_zfock,GPARAM,nbf,emin,g_zcheb
     1 ,dt)

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

      double precision  GPARAM !(emax-emin)
      double precision scale
      double precision emin
      double precision alpha !scaling factor
      double precision dt !time step
      integer g_zdens
      integer g_zfock !complex Fock matrix
      integer g_tmp1, g_tmp2
      integer g_zt_n !nth term in Chebyshev polynomial
      integer g_zt_n_1 !(n-1)th term in Chebyshev polynomial
      integer g_zt_n_2 !(n-2)th term in Chebyshev polynomial
      integer limat, iimat
      integer g_I !Identity matrix
      integer g_zcheb !Total Chebyshev propagator
      integer g_zscr
      integer g_zw !Hermitian conjugate of g_cheb
      logical status
      integer i,itag
      integer nbf,nrr
      double precision param

         nrr = nbf*(nbf+1)/2
c     if (ga_nodeid().eq.0)  write(6,*)"in prop_chebyshev"

      if (.not. ga_create(mt_dbl, nbf, nbf, "Re[Dens]", 0, 0, g_I))
     $     call errquit ("couldnt create g_I", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n, "Tn"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n_1, "Tn-1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zt_n_2, "Tn-2"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_zscr, "scr"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      if (.not. ga_duplicate(g_zfock, g_zw, "Cheby"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      call ga_zero(g_I)
      call ga_zero(g_zt_n)
      call ga_zero(g_zt_n_1)
      call ga_zero(g_zt_n_2)
      call ga_zero(g_zcheb)
      call ga_zero(g_zscr)
      call ga_zero(g_zw)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'Imat',limat,iimat))
     &   call errquit('prop_chebyshev: cannot allocate imat',0,
     1    MA_ERR)



      call get_imat(Dbl_MB(iimat),nbf)

      call mat_sym(Dbl_MB(iimat),nbf,1)

C---------------------------------------------------------------------
C     Create T0
C---------------------------------------------------------------------

      call ga_put(g_I,1,nbf,1,nbf,Dbl_MB(iimat),nbf)

      call convert_d2z (1d0, g_I, 0d0, g_I, g_zt_n_2)

      scale = 1/(1.5*GPARAM)
  
c     write(6,*)"scale",scale
      call ga_scale(g_zfock,scale)

c     call ga_print(g_zfock)      

c     call convert_d2z (1d0, g_I, 0d0, g_I, g_zt_n_2)
C---------------------------------------------------------------------
C     Create T1
C-------------------------------------------------------------------- 
      call ga_copy(g_zfock,g_zt_n_1)
C-------------------------------------------------------------------
      alpha = 1.5*GPARAM *dt      

      call chebyshev_propagator(g_zt_n_1,g_zt_n_2,g_zt_n,g_zfock,g_zcheb
     1 ,nbf,alpha)

c       call ga_print(g_zcheb)
 
         if (.not. mat_is_unitary (g_zcheb, 1d-10))
     &   call errquit('prop_chebyshev: g_zcheb not unitary',0,
     1    MA_ERR)

      call mat_hermconj (g_zcheb,g_zw)
      
c     call ga_print(g_zdens)
      call ga_zgemm('N','N',nbf,nbf,nbf,z1,g_zdens,g_zw,z0,
     1 g_zscr)

c      call ga_print(g_zscr)
       call ga_zgemm('N','N',nbf,nbf,nbf,z1,g_zcheb,g_zscr,z0,
     1 g_zdens)

c     call ga_print(g_zdens)
      if (.not.ma_pop_stack(limat))
     &   call errquit('prop_chebyshev: cannot pop stack fock',3,
     1 MA_ERR)
      status = .true.
      status=status.and.ga_destroy(g_zw)
      status=status.and.ga_destroy(g_zscr) 
      status=status.and.ga_destroy(g_zt_n_2)
      status=status.and.ga_destroy(g_zt_n_1)
      status=status.and.ga_destroy(g_zt_n)
      status=status.and.ga_destroy(g_I)
      if (.not.status) call errquit ("prop_chebyshev:destroy failed", 0
     1 , GA_ERR)

      return
      end 

C-----------------------------------------------------------------------------------
C
C-----------------------------------------------------------------------------------
      subroutine chebyshev_propagator(g_zt_n_1,g_zt_n_2,g_zt_n,g_zfock,
     1 g_zcheb,nbf,alpha)

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
      integer g_zfock !complex Fock matrix
      integer g_ztemp 
      integer g_zt_n !nth term in Chebyshev polynomial
      integer g_zt_n_1 !(n-1)th term in Chebyshev polynomial
      integer g_zt_n_2 !(n-2)th term in Chebyshev polynomial
      integer g_zcheb !Total Chebyshev operator
      integer g_zcheb_n !nth term in chebyshev propagator
      integer g_zcheb_n_1 !(n-1)th term in chebyshev propagator
      integer IT,nbf,nrr,i
      double precision TEST !convergence criteria
      double precision alpha !scaling factor
      double precision bessfun !bessel function
c     double complex zi
      logical, external ::mat_is_same_semiemp

      LOGICAL CONVERGED,status
      DATA TEST /1.0D-10/
c     DATA TEST/1.0D-23/
c     DATA TEST /1.0D-30/
      CONVERGED = .FALSE.

c     call ga_print(g_zt_n_2)
c     call ga_print(g_zcheb)
      nrr = nbf*(nbf+1)/2
      IT = 0
C-------------------------------------------------------------
C     1st term
C------------------------------------------------------------
      Call bessel_function(IT,alpha,bessfun)
      call ga_add(z0,g_zcheb,z1*bessfun,g_zt_n_2,g_zcheb)
C-----------------------------------------------------------
C     2nd term
C-----------------------------------------------------------      
      IT = 1
      call bessel_function(IT,alpha,bessfun)
      if (.not. ga_duplicate(g_zfock, g_zcheb_n_1, "Cheb_n_1"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
      call ga_zero(g_zcheb_n_1)

      call ga_add(z0,g_zcheb_n_1,-zi*bessfun,g_zt_n_1,g_zcheb_n_1)

      call ga_add(z1,g_zcheb,z1,g_zcheb_n_1,g_zcheb)

      if (.not. ga_duplicate(g_zfock, g_zcheb_n, "Cheb_n"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)

      if (.not. ga_duplicate(g_zfock, g_ztemp, "Temp"))
     $     call errquit ("couldnt duplicate g_zfock", 0, GA_ERR)
       call ga_zero(g_zcheb_n)
       call ga_zero(g_ztemp)

  10  CONTINUE
C---------------------------------------------------------------
C    Iteration
C--------------------------------------------------------------
       IT = IT + 1
c      write(6,*)"IT",IT

       call ga_zgemm('N','N',nbf,nbf,nbf,2*z1,g_zfock,g_zt_n_1,z0,
     1 g_ztemp)

       call ga_add(z1,g_ztemp,-z1,g_zt_n_2,g_zt_n)

       call bessel_function(IT,alpha,bessfun)       

c      call ga_print(g_zcheb_n)
       call ga_add(z0,g_zcheb_n,((-zi)**IT)*bessfun,g_zt_n,g_zcheb_n)
       
c      call ga_print(g_zcheb_n)
       call ga_add(z1,g_zcheb,z1,g_zcheb_n,g_zcheb)

c     call ga_print(g_zcheb)
      if (IT.eq.2) GO TO 20

       if ( mat_is_same_semiemp (g_zcheb_n_1,g_zcheb_n,test)) then
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


      if (ga_nodeid().eq.0)  write(6,*)"Chebyshev Iteration",IT
      status = .true.
      status=status.and.ga_destroy(g_ztemp)
      status=status.and.ga_destroy(g_zcheb_n)
      status=status.and.ga_destroy(g_zcheb_n_1)
      if (.not.status) call errquit ("prop_chebyshev:destroy failed", 0
     1 , GA_ERR)
       return
       end
C------------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------------
      subroutine get_imat(imat,n)

      double precision imat(*)
      integer n,i,j,L

       do j = 1,n*n
        imat(j) =0.0d0
       end do
      L = 0
      do i = 1,n
       do j = 1,i
        L = L+1
        if (i.eq.j) then
        imat(L) = 1.0d0
        end if
       end do
      end do

      return
      end subroutine
C-----------------------------------------------------------------------------------------------
C
C----------------------------------------------------------------------------------------------
      subroutine bessel_function(IT,alpha,bessfun)
      implicit none        
      double precision FACT(35), P, SUM, TEST, PREV, TERM, bessfun,diff
      double precision alpha
      integer it, m

      data FACT/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
     X          ,3628800.0,39916800.0,479001600.0,6227020800.0,
     X          87178291200.0,1307674368000.0,20922789888000.0,
     X          355687428096000.0,6402373705728000.0,
     X          1.21645100408832000E+017,2.43290200817664000E+018,
     X          5.10909421717094400E+019,1.12400072777760768E+021,
     X          2.58520167388849782E+022,6.20448401733239410E+023,
     X          1.55112100433309861E+025,4.03291461126605650E+026,
     X          1.08888694504183519E+028,3.04888344611713837E+029,
     X          8.84176199373970077E+030,2.65252859812191032E+032,
     X          8.22283865417792243E+033, 2.63130836933693518E+035,
     X          8.68331761881188594E+036,2.95232799039604120E+038/
      DATA TEST /1.0D-10/

      IF (IT.eq.0) then
       P = 1.0d0
      ELSE
       P = 2.0d0
      END IF
      SUM = 0.0d0
      PREV = 0.0d0
      TERM = 0.0d0
      bessfun = 0.0d0
      m = 0
   10 CONTINUE
      TERM = P*((-1)**m)/(FACT(m+1)*FACT(m+IT+1))*(alpha/2)**(2*m+IT)
c     write(6,*)"TERM",TERM
      SUM = SUM +TERM
      diff = abs(TERM) - abs(PREV)
      IF (abs(diff).le.TEST)  GO TO 20
      PREV = TERM
      m = m + 1
      GO TO 10

   20 CONTINUE
c     write(6,*) "SUM",SUM
      bessfun = SUM
      RETURN
      END


 
