      SUBROUTINE DIFPEN(FACT,AMU,N,BMU,M,RR,DPEN,IP)

CGMP          IP equal to I in the penet.f

C
C     DERIVATIVES OF INTEGRALS X(A)*X(A)/R(B)  WHERE X(A) IS A.O.
C     IP IS TYPE
C     OF PENET. INT. IP=1 FOR INT. X(A)*Y(A)/R(B) WHERE X(A) IS ORB ON CE
C     NTER A, PRINC. QUANTUM NO. N AND EXPONENTIAL CONSTANT AMU. IP=2 IS
C     TYPE  INT.X(A)*X(B)/R(A), WHERE BMU AND M ARE ASSOCIATED WITH X(B)
C     RR IS SEPARATION BETWEEN ATOMIC CENTERS A AND B IN ANGSTROMS, FACT
C     (N+1) IS FACTORIAL N, AND DPEN IS THE SIGMA COMPONENT OF THE ANSWER
C
C  MCZ UG 1980
C  REVISED BY JDHD    FEB 3 85
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
C     CONST.CMN
      PARAMETER (BOHRS    = 5.29177249D-1,
     .           EVOLT    = 2.72113961D1,
     .           ENUBAR   = 2.1947463067D5,
     .           DEBYE    = 2.541748047D0)

c     INCLUDE '../include/binc.cmn'
C****************************************************************
C***      COMMON/BINC/BINCOE(465)
C****************************************************************
c     INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON/EV/
      DIMENSION A(35),B(35),DD(35)
C****************************************************************
      DIMENSION    FACT(*),BINCOE(465)

c     write(6,*)"i am in difpen"
C
C     FORM BINOMIAL COEFFECIENT
C
C     ( N )
C     ( M )  INDEX = N*(N+1)/2 + M+1
C
      BINCOE(1) = ONE
      IND = 2
      DO 70 I = 1, 29
       I1 = I + 1
       DO 60 J = 1, I1
        JJ = I1 - J + 1
        IF(J.LT.26) GO TO 40
        BIN = 1.D0
        IF(J.GE.JJ) THEN
         DO 20 K = J, I
          BIN = BIN*DFLOAT(K)
   20    CONTINUE
         BINCOE(IND) = BIN/FACT(JJ)
        ELSE
         DO 30 K = JJ, I
          BIN = BIN*DFLOAT(K)
   30    CONTINUE
         BINCOE(IND) = BIN/FACT(J)
        END IF
        GO TO 50

   40   BINCOE(IND) = FACT(I1)/FACT(J)/FACT(JJ)
   50   IND = IND + 1
   60  CONTINUE
   70 CONTINUE

      GO TO (100,300), IP

  100 T = ONE   
      R = RR/BOHRS
      Z = (AMU+BMU)/TWO
      IST = N + M
      P = Z*R
C     THIS SWITCH IS REQUIRED BECAUSE AUX HAS A 40.0 CUTOFF
      IF (P.GE.40.0D0) THEN
         NN = 2*N + 1
         MM = 2*M + 1
         T = ABS((AMU-BMU)/(AMU+BMU))
         F = (ONE-T)**MM* (ONE+T)**NN
         S = FACT(N+M +1)*SQRT(F/ (FACT(NN)*FACT(MM)))
         DPEN = -ONE*S/ (R*R)
C        THESE ARE ALL ONE CENTER AND MUST BE IN PHASE
C        IF (LA.GT.0) SP = (-ONE)** (LA+1)*S
C        IF (LA.GT.1) SD = (-ONE)** (LA+2)*S
C        IF (LA.GT.2) SF = (-ONE)** (LA+3)*S
C        S = (-ONE)**LA*S
C        IF (LA.GT.0) DPENPI = DPEN
C        IF (LA.GT.1) DPEND = DPEN
C        IF (LA.GT.2) DPEND = DPEN          
         RETURN

      END IF
      CALL AUC(P,T,A,B)
C
C  DERIVATIVE OF A AND B FUNCTIONS GIVEN BY
C  DA(I) = -A(I+1)   DB(I) = -T*B(I+1) WITH  T=1.
C
      IBINC = (IST-1)*IST/2
      SUM = 0.D0
      DSUM = 0.D0
      DO 200 I = 1, IST
         IBINC = IBINC + 1
         SUM = SUM + A(IST-I+1)*B(I)*BINCOE(IBINC)
         DSUM = DSUM - (A(IST-I+1)*B(I+1)+A(IST-I+2)*B(I))*BINCOE(IBINC)
  200 CONTINUE
      SCAL = (AMU**N)* (BMU**M)/FACT(2*M+1)*
     .       SQRT(AMU*BMU*FACT(2*M+1)/FACT(2*N+1))* (R**IST)
      PEN = SUM*SCAL
      DPEN = PEN*DFLOAT(IST)/R + DSUM*SCAL*Z
  300 RETURN

      END






