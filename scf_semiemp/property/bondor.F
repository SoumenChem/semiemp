      SUBROUTINE BONDOR(P,NRR,NIN,NU,ISW1,NA,NB,ISWANN)
C
C     THIS SUBROUTINE CALC A BOND INDEX A-LA PARR ZERNER
C     THIS IS RELATED TO THE WYBERG INDEX?
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N) 
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
#include "inp.fh"
c     INCLUDE '../include/rhf.cmn'
C****************************************************************
C***      COMMON/RHF/
C***     .  AR(21),BR(21),FOP(8),NOP,NDN,MIM(14),ITAPE
C****************************************************************
c     INCLUDE '../include/switch.cmn'
C****************************************************************
C***      COMMON/SWITCH/
C***     .  ISWANN,III,IPRINT,IPUN,ISAVE,ISPEC,ISW1,ISW2,ISW3,ISW4
C***     .  JSW(10),IBETA,IDD1,IDD2
C****************************************************************
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     .  ATOM(100),NIN(MXB),
C***     .  NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     .  NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     .  Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D(MXA)
C****************************************************************
C
c     INCLUDE '../include/debug.cmn'
      DIMENSION    P(*),  NIN(*), NU(*)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0, IZERO = 0)

c     write(6,*)"i am in bondor"!,(P(I),I=1,NRR)
c     IF(ISW1.EQ.0.OR.ISW1.EQ.7.OR.ISW1.EQ.8.OR.ISW1.EQ.9) THEN

             if (.not.MA_Push_Get(MT_Dbl,nrr,'P2',lP2,iP2))
     &   call errquit('scfit: cannot allocate W',0,
     1    MA_ERR)
c     ELSE
c       IF(JSW(16).NE.1) CALL FINDP(P,NRR,ISWANN,NDN,1)
c     ENDIF
      IF (ISWANN .EQ. IZERO) THEN
       DO I = 1, NRR
c      write(6,*)"P",2*P(I)
        P(I) = 2*P(I)
       END DO
      END IF
      IF (ISW1.EQ.7) GO TO 50
      NAR = (NA*NA+NA)/2
      DO 10 I = 1, NAR
        Dbl_MB(iP2+I-1) = ZERO 
   10 CONTINUE
      DO 30 I = 1, NB
         NI = NU(I)
         DO 20 J = 1, I
            NJ = NU(J)
            IJ = NIN(I) + J
            NINJ = NIN(NI) + NJ
            Dbl_MB(iP2+NINJ-1) = Dbl_MB(iP2+NINJ-1) + P(IJ)**2
   20    CONTINUE
   30 CONTINUE
      IA = 0
      DO 40 I = 1, NA
         IA = IA + I
         Dbl_MB(iP2+IA-1) = ZERO 
   40 CONTINUE
c     IF (ISWANN .EQ. IZERO) THEN
c      DO I = 1, NRR
c      write(6,*)"P",2*P(I)
c       P(I) = P(I)/2
c      END DO
c     END IF

      if (ga_nodeid().eq.0) WRITE (6,9000)
c     if (ga_nodeid().eq.0)  write(6,*)"P2",(Dbl_MB(iP2+I),I=1,NRR)
 9000 FORMAT (11X,'ATOMIC BOND INDEX'//)

c     CALL PRNT(P2,NA,NAR)
   50 CONTINUE
c     IF (USEUNIT0) WRITE (0,*) 'BONDOR out'
      if (.not.ma_pop_stack(lP2))
     &   call errquit('scfit: cannot pop stack Q',3, MA_ERR)
      RETURN

      END
