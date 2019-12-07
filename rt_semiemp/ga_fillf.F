      subroutine ga_fillf(g_dens,g_fock,nbf,LABI,LABJ,LABK,LABL,QQ,
     1 g_core,NINTGL,itag,NIN)
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

      integer LABI(*), LABJ(*), LABK(*), LABL(*), NIN(*) 
      integer NINTGL, LL, I, J, K, L, ITAG, IN, KN, nbf
      double precision QQ(*), V, alpha, beta
      integer g_dens,g_fock,g_core, me
      double precision dterm, dterm1, fterm, fterm1

      me = ga_nodeid()

      call ga_zero(g_fock)
      alpha = 1.0d0
      beta = 1.0d0
      call ga_sync()
      if (me.eq.0) then
c     write(6,*)"NINTGL",NINTGL
      LL = 0
   20 CONTINUE
C     BEGIN TWO _ ELECTRON PART
      LL = LL + 1
C     write(6,*)"LL,NINTGL",LL,NINTGL

      IF (LL.LE.NINTGL) THEN
      GO TO 30
      ELSE
      GO TO 170
      ENDIF
c     IF (NINTGL.LE.0) GO TO 170
c     I = 0
c     J = 0
c     K = 0
c     L = 0
c     LL = 1
   30 CONTINUE
      I = LABI(LL)
      J = LABJ(LL)
      K = LABK(LL)
      L = LABL(LL)
      V = QQ(LL)
      dterm = 0.0d0
      fterm = 0.0d0
      dterm1 = 0.0d0
      fterm1 = 0.0d0

      IF ((I-J).LT.0)   GO TO 90
      IF ((I-J).EQ.0)   GO TO 40
      GO TO 180


   40 IF ((K-L).LT.0)  GO TO 80
      IF ((K-L).EQ.0)  GO TO 50
      GOTO 180

   50 IF ((I-K).LT.0)  GO TO 70
      IF ((I-K).EQ.0) GO TO 60
      GOTO 180

c     (I,I/I,I)
   60 call ga_get(g_dens,i,i,i,i,dterm,1)

      call ga_get(g_fock,i,i,i,i,fterm,1)

      fterm = fterm + (dterm-0.5*dterm)*V

      call ga_put(g_fock,i,i,i,i,fterm,1)
      GO TO 160
C     (I,I/K,K) 
   70 call ga_get(g_dens,i,i,i,i,dterm,1)
      call ga_get(g_dens,k,k,k,k,dterm1,1)

      call ga_get(g_fock,i,i,i,i,fterm,1)
      call ga_get(g_fock,k,k,k,k,fterm1,1) 
     
      fterm = fterm + dterm1*V
      fterm1 = fterm1 + dterm*V

      call ga_put(g_fock,i,i,i,i,fterm,1)
      call ga_put(g_fock,k,k,k,k,fterm1,1)

      call ga_get(g_dens,k,k,i,i,dterm,1)
      call ga_get(g_fock,k,k,i,i,fterm,1)  

      fterm = fterm - 0.5*dterm*V

      call ga_put(g_fock,k,k,i,i,fterm,1)
      GO TO 160
C     (I,I/K,L)
   80 call ga_get(g_dens,i,i,i,i,dterm,1)
      call ga_get(g_dens,l,l,k,k,dterm1,1)

      call ga_get(g_fock,i,i,i,i,fterm,1)
      call ga_get(g_fock,l,l,k,k,fterm1,1)

      if (ITAG.eq.1)
     .fterm = fterm + 2.*dterm1*V
      fterm1 = fterm1 + dterm*V

      call ga_put(g_fock,i,i,i,i,fterm,1)
      call ga_put(g_fock,l,l,k,k,fterm1,1)

      call ga_get(g_dens,k,k,i,i,dterm,1)
      call ga_get(g_dens,l,l,i,i,dterm1,1)

      call ga_get(g_fock,k,k,i,i,fterm,1)
      call ga_get(g_fock,l,l,i,i,fterm1,1)
      fterm = fterm - 0.5*dterm1*V
      fterm1 = fterm1 - 0.5*dterm*V

      IF (I.NE.K) THEN

      call ga_put(g_fock,k,k,i,i,fterm,1)
      call ga_put(g_fock,l,l,i,i,fterm1,1)
     
      GO TO 160
      END IF
      IF (ITAG.eq.1)
     .fterm = fterm - 0.5*dterm1*V
      IF (ITAG.eq.2)
     .fterm = fterm + 0.5*dterm1*V

      call ga_put(g_fock,k,k,i,i,fterm,1)
      call ga_put(g_fock,l,l,i,i,fterm1,1)
      GO TO 160
   90 IF ((K-L).LT.0)  GO TO 110
      IF ((K-L).EQ.0)  GO TO 100
      GOTO 180
C     (I,J,/K,K)
  100 call ga_get(g_dens,j,j,i,i,dterm,1) 
      call ga_get(g_dens,k,k,k,k,dterm1,1)

      call ga_get(g_fock,j,j,i,i,fterm,1)
      call ga_get(g_fock,k,k,k,k,fterm1,1)
 
      fterm = fterm + dterm1*V
      IF (ITAG.eq.1)
     .fterm1 =fterm1 + 2.*dterm*V

      call ga_put(g_fock,j,j,i,i,fterm,1)
      call ga_put(g_fock,k,k,k,k,fterm1,1)
     
      call ga_get(g_dens,k,k,i,i,dterm,1)
      call ga_get(g_dens,k,k,j,j,dterm1,1)

      call ga_get(g_fock,k,k,i,i,fterm,1)
      call ga_get(g_fock,k,k,j,j,fterm1,1)

      IF (K.LT.J) call ga_get(g_dens,j,j,k,k,dterm1,1)
      IF (K.LT.J) call ga_get(g_fock,j,j,k,k,fterm1,1)

      fterm = fterm - 0.5*dterm1*V
      fterm1 = fterm1 - 0.5*dterm*V
      IF (J.NE.K) THEN
       call ga_put(g_fock,k,k,i,i,fterm,1)
       IF (K.LT.J) THEN
        call ga_put(g_fock,j,j,k,k,fterm1,1)
       ELSE
        call ga_put(g_fock,k,k,j,j,fterm1,1)
       END IF
      GO TO 160
      END IF
      IF (ITAG.eq.1)
     .fterm1 = fterm1 - 0.5*dterm*V
      IF (ITAG.eq.2)
     .fterm1 = fterm1 + 0.5*dterm*V

      call ga_put(g_fock,k,k,i,i,fterm,1)

      IF (K.LT.J) THEN
      call ga_put(g_fock,j,j,k,k,fterm1,1)
      ELSE
      call ga_put(g_fock,k,k,j,j,fterm1,1)
      END IF
      GO TO 160
C     (I,J/K,L)
  110 call ga_get(g_dens,j,j,i,i,dterm,1)
      call ga_get(g_dens,l,l,k,k,dterm1,1)

      call ga_get(g_fock,j,j,i,i,fterm,1)
      call ga_get(g_fock,l,l,k,k,fterm1,1)
      IF (ITAG.eq.1)
     .fterm = fterm + 2.*dterm1*V
      call ga_put(g_fock,j,j,i,i,fterm,1)
      IN = NIN(J) + I
      KN = NIN(L) + K
      IF (IN.EQ.KN) GO TO 120
      IF (ITAG.eq.1)
     .fterm1 = fterm1 + 2.*dterm*V
      call ga_put(g_fock,l,l,k,k,fterm1,1)

  120 call ga_get(g_dens,k,k,i,i,dterm,1)
      call ga_get(g_dens,l,l,j,j,dterm1,1)

      call ga_get(g_fock,k,k,i,i,fterm,1)
      call ga_get(g_fock,l,l,j,j,fterm1,1)
 
      IF (J.GT.L) call ga_get(g_dens,j,j,l,l,dterm1,1)
      IF (J.GT.L) call ga_get(g_fock,j,j,l,l,fterm1,1)

      fterm = fterm - 0.5*dterm1*V
      fterm1 = fterm1 - 0.5*dterm*V
      call ga_put(g_fock,k,k,i,i,fterm,1)
      IF (J.GT.L) THEN
      call ga_put(g_fock,j,j,l,l,fterm1,1)
      ELSE
      call ga_put(g_fock,l,l,j,j,fterm1,1)
      END IF

      IF (I.NE.K) GO TO 130
      IF (J.EQ.L) GO TO 140
      IF (ITAG.eq.1)
     .fterm = fterm - 0.5*dterm1*V
      IF (ITAG.eq.2)
     .fterm = fterm + 0.5*dterm1*V
      call ga_put(g_fock,k,k,i,i,fterm,1)
      GO TO 140
  130 IF (J.NE.L) GO TO 140
      fterm1 = fterm1 - 0.5*dterm*V
      IF (J.GT.L) THEN
      call ga_put(g_fock,j,j,l,l,fterm1,1)
      ELSE
      call ga_put(g_fock,l,l,j,j,fterm1,1)
      END IF

  140 call ga_get(g_dens,l,l,i,i,dterm,1)
      call ga_get(g_dens,j,j,k,k,dterm1,1)

      call ga_get(g_fock,l,l,i,i,fterm,1)
      call ga_get(g_fock,j,j,k,k,fterm1,1) 

      IF (K.GT.J) call ga_get(g_dens,k,k,j,j,dterm1,1)
      IF (K.GT.J) call ga_get(g_fock,k,k,j,j,fterm1,1)

      fterm = fterm - 0.5*dterm1*V

      call ga_put(g_fock,l,l,i,i,fterm,1)

      IN = NIN(L) + I
      KN = NIN(J) + K
      IF (IN.EQ.KN) GO TO 150
      fterm1 = fterm1 - 0.5*dterm*V

      IF (K.GT.J) THEN
      call ga_put(g_fock,k,k,j,j,fterm1,1)
      ELSE
      call ga_put(g_fock,j,j,k,k,fterm1,1)
      END IF

  150 IF (J.NE.K) GO TO 160
      IF (ITAG.eq.1)
     .fterm1 = fterm1 - 0.5*dterm*V
      IF (ITAG.eq.2)
     .fterm1 = fterm1 + 0.5*dterm*V
      call ga_put(g_fock,j,j,k,k,fterm1,1)
  160 CONTINUE

      GO TO 20

  170 CONTINUE
  175 CONTINUE

      call ga_sync()
      call ga_print(g_fock)
c     write(6,*)"HERE"
      do i = 2, nbf
       do j = 1, i-1
        call ga_get(g_fock,i,i,j,j,fterm,1)
        if (itag.eq.1)  call ga_put(g_fock,j,j,i,i,fterm,1)
        if (itag.eq.2)  call ga_put(g_fock,j,j,i,i,-fterm,1)
       end do
      end do
      end if
c     call ga_print(g_fock)
c     call ga_print(g_core)

      IF (ITAG.eq.2) GO TO 190
      call ga_add(alpha,g_fock,beta,g_core,g_fock)
  189 CONTINUE
c     WRITE (6,*)"BB", (BB(M),M=1,NRR)
      GO TO 190
  180 WRITE (6,9000) I, J, K, L
      CALL ZEXIT(3)
CNR
  190 CONTINUE
CNR
      call ga_print(g_fock)
      RETURN

 9000 FORMAT (' **** ERROR IN INTEGRAL LIST',4I5,'****')
      end subroutine



