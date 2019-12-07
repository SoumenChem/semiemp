      SUBROUTINE DIFOVL(LH,N1,L1,AMU,N2,L2,BMU,R,S,SP,SD,SF,DSS,DSP,DSD,
     .                  DSF,FACT,NIN,PP)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
C
C     EVALUATES OVERLAP AND ITS DERIVATIVE BETWEEN TWO SLATER TYPE ORBIT
C
C     N1,N2, ARE PRINC. Q.NO.-S, L1,L2, ARE SECONDARY Q.N.-S,, AMU AND
C     BMU ARE THE EXPONENTIAL CONSTANTS. R IS THE SEPARATION IN ANGSTROM
C     S IS THE SIGMA, SP THE PI, SD THE DELTA AND SF THE F COMPOENTS OF
C     THE OVERLAP, TO BE PUT TOGETHER USING LH AND SUB. GEOM.
C
C     FACT(I+1)= FACTORIAL I.
C
      PARAMETER (BOHRS    = 5.29177249D-1)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)

      DIMENSION    FACT(*), SS(4), DS(4), NIN(*)
      DIMENSION    A(35),B(35),DD(35)
C      REAL PP/0.0D0/,TT/0.0D0/
c     INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON/EV/
C***     .  A(35),B(35),DD(35)
C****************************************************************
c     INCLUDE '../include/pptt.cmn'
C****************************************************************
C***      COMMON/PPTT/PP
C****************************************************************
c     write(6,*)"i am in difovl"
      IF (PP.EQ.-10.D0) THEN
         TT = 0.D0
      END IF
      RR = R/BOHRS
      P = (AMU+BMU)*RR/TWO
      T = (AMU-BMU)/ (AMU+BMU)
      IF (R.LE.0.001) GO TO 30
      IF (ABS(T).LE.0.0001) T = ZERO
      IF (P.EQ.PP .AND. T.EQ.TT) GO TO 10
      PP = P
      TT = T
C      NM1=N1+N2+1
      CALL AUC(P,T,A,B)
   10 CONTINUE
      SS(2) = ZERO
      SS(3) = ZERO
      SS(4) = ZERO
      DS(2) = ZERO
      DS(3) = ZERO
      DS(4) = ZERO
      LMIN = L1
      IF (L2.LT.L1) LMIN = L2
      LMIN = LMIN + 1
      DO 20 II = 1, LMIN
         I = II - 1
         CALL MOLPD1(N1,N2,L1,L2,I,I,AMU,BMU,A,B,FACT,RR,SSSS,DSSSS,NIN)
c     write(6,*)"SSSS,DSSSS",SSSS,DSSSS
         SS(II) = SSSS
         DS(II) = DSSSS
   20 CONTINUE
      S = SS(1)
      SP = SS(2)
      SD = SS(3)
      SF = SS(4)
C
C  GET DERIVATE OF OVERLAP
      DSS = DS(1)
      DSP = DS(2)
      DSD = DS(3)
      DSF = DS(4)
      GO TO 40
C     ONE CENTER OVERLAP.
   30 S = ZERO
      SP = ZERO
      SD = ZERO
      SF = ZERO
      DSS = ZERO
      DSP = ZERO
      DSD = ZERO
      DSF = ZERO
      IF (L1.NE.L2) GO TO 40
      NN = 2*N1 + 1
      MM = 2*N2 + 1
      F = (ONE-T)**MM* (ONE+T)**NN
      S = FACT(N1+N2+1)*SQRT(F/ (FACT(NN)*FACT(MM)))
      IF (L1.GT.0) SP = (-ONE)** (L1+1)*S
      IF (L1.GT.1) SD = (-ONE)** (L1+2)*S
      IF (L1.GT.2) SF = (-ONE)** (L1+3)*S
      S = (-ONE)**L1*S
   40 RETURN

      END
