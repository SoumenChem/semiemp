      SUBROUTINE PENET(FACT,AMU,N,LA,BMU,M,LB,RR,PEN,PENPI,PEND,PENDF,I
     1 ,NIN)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     EVALUATES PENETRATION INTEGRALS BETWEEN  S.T.O.-S.      I IS TYPE
C     OF PENET. INT. I=1 FOR INT. X(A)*Y(A)/R(B) WHERE X(A) IS ORB ON 
C     CENTER A, PRINC. QUANTUM NO. N AND EXPONENTIAL CONSTANT AMU. I=2 IS
C     TYPE  INT.X(A)*X(B)/R(A), WHERE BMU AND M ARE ASSOCIATED WITH X(B)
C     RR IS SEPARATION BETWEEN ATOMIC CENTERS A AND B IN ANGSTROMS, FACT
C     (N+1) IS FACTORIAL N, AND PEN IS THE SIGMA COMPONENT OF THE ANSWER
C     PENPI THE PI COMPONENT, PEND THE DELTA COMPONENT IN HARTREES.
C     LA, LB = 0,1,2,3
C              S,P,D,F
C     M. C. ZERNER   UF.  MODIFIED 1988.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
C
c     INCLUDE '../include/binc.cmn'
C****************************************************************
C***      COMMON/BINC/
C***     .  BINCOE(465)
C****************************************************************
c     INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON/EV/
C***     .  A(35),B(35),DD(25)
C****************************************************************
C
      DIMENSION A(35),B(35),DD(25)
      DIMENSION    FACT(*), BINCOE(465),NIN(*)
C
c     write(6,*)"I am in penet",AMU,N,LA,BMU,M,LB,RR,NIN(1)
C
C     FORM BINOMIAL COEFFECIENT
C
C     ( N )
C     ( M )  INDEX = N*(N+1)/2 + M+1
C
      BINCOE(1) = ONE
      IND = 2
      DO 7 L = 1, 29
       I1 = L + 1
       DO 6 J = 1, I1
        JJ = I1 - J + 1
        IF(J.LT.26) GO TO 4
        BIN = 1.D0
        IF(J.GE.JJ) THEN
         DO 2 K = J, L
          BIN = BIN*DFLOAT(K)
   2    CONTINUE
         BINCOE(IND) = BIN/FACT(JJ)
        ELSE
         DO 3 K = JJ, L
          BIN = BIN*DFLOAT(K)
   3    CONTINUE
         BINCOE(IND) = BIN/FACT(J)
        END IF
        GO TO 5

   4   BINCOE(IND) = FACT(I1)/FACT(J)/FACT(JJ)
   5   IND = IND + 1
   6   CONTINUE
   7  CONTINUE

      N2 = N + M
      FN = FLOAT(N2)
      R = RR!/BOHRS
      GO TO (10,190), I

   10 CONTINUE
      PENPI = ZERO
      PEND = ZERO
C     (X(A)Y(A)/1/R(B))
      Z = (AMU+BMU)/2.0D0
      AR = Z*R
      IF (AR.GT.40.0D0) THEN
         IF (LA.EQ.LB) THEN
            PEN = 1.0D0/R
C
            NN = 2*N + 1
            MM = 2*M + 1
            T = ABS((AMU-BMU)/(AMU+BMU))
            F = (ONE-T)**MM* (ONE+T)**NN
            S = FACT(N+M +1)*SQRT(F/ (FACT(NN)*FACT(MM)))
C           THESE ARE ALL ONE CENTER AND MUST BE IN PHASE
C           IF (LA.GT.0) SP = (-ONE)** (LA+1)*S
C           IF (LA.GT.1) SD = (-ONE)** (LA+2)*S
C           IF (LA.GT.2) SF = (-ONE)** (LA+3)*S
C           S = (-ONE)**LA*S                               
            PEN = S*PEN
            IF (LA.GT.0) PENPI = PEN
            IF (LA.GT.1) PEND = PEN
            IF (LA.GT.2) PEND = PEN
         ELSE
            PEN = ZERO
         END IF
         GO TO 200

      END IF
      IF (R.LT.1.E-6) GO TO 180
      CALL AUC(AR,ONE,A,B)
      X = ZERO
      IF (LA.GE.LB) THEN
         IJ = (LA* (LA+1))/2 + LB + 1
      ELSE
         IJ = (LB* (LB+1))/2 + LA + 1
      END IF
      GO TO (20,40,60,80,100,120,130,150,160,
     .       170), IJ
C     S-S
   20 CONTINUE
      IBINC = ((N2-1)*N2)/2
      DO 30 J = 1, N2
         IBINC = IBINC + 1
         X = X + B(N2-J+1)*A(J)*BINCOE(IBINC)
   30 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)
      PEN = PEN*X/FACFAC
c     write(6,*)"PEN,PENPI,PEND,PENDF",PEN,PENPI,PEND,PENDF
      GO TO 200
C     S-P
   40 CONTINUE
      IN2 = N2 - 1
      IBINC = ((IN2-1)*IN2)/2
      DO 50 J = 1, IN2
         IBINC = IBINC + 1
         X = X + (B(N2-J)*A(J)+B(N2-J+1)*A(J+1))*BINCOE(IBINC)
   50 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)*SQRT(3.0D0)
      PEN = PEN*X/FACFAC
      GO TO 200

   60 CONTINUE
C     P-P  X IS THE SIGMA PART, Y IS THE PI PART.
      IN2 = N2 - 2
      IBINC = ((IN2-1)*IN2)/2
      Y = ZERO
      DO 70 J = 1, IN2
         IBINC = IBINC + 1
         X = X + (B(N2-J-1)*A(J)+2.0D0*B(N2-J)*A(J+1)+B(N2-J+1)*A(J+2))*
     .       BINCOE(IBINC)
         Y = Y + (B(J)-B(J+2))* (A(N2-J+1)-A(N2-1-J))*BINCOE(IBINC)
   70 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)*3.0D0
      PENPI = PEN*Y/ (2.0D0*FACFAC)
      PEN = PEN*X/FACFAC
      GO TO 200

   80 CONTINUE
C     S-D
      IN2 = N2 - 2
      IBINC = ((IN2-1)*IN2)/2
      DO 90 J = 1, IN2
         IBINC = IBINC + 1
         X = X + (A(J)* (3.0D0*B(N2-J-1)-B(N2-J+1))+
     .       A(J+2)* (3.0D0*B(N2-J+1)-B(N2-J-1))+4.0D0*A(J+1)*B(N2-J))*
     .       BINCOE(IBINC)
   90 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)*SQRT(3.0D0)
      PEN = PEN*X/FACFAC
      GO TO 200

  100 CONTINUE
C     P-D  X=SIGMA PART,  Y= PI PART
      IN2 = N2 - 3
      IBINC = ((IN2-1)*IN2)/2
      Y = ZERO
      DO 110 J = 1, IN2
         IBINC = IBINC + 1
         X = X + (3.0D0*A(J)*B(IN2-J+1)+
     .       7.0D0* (A(J+2)*B(IN2-J+3)+A(J+1)*B(IN2-J+2))-
     .       A(J+2)*B(IN2-J+1)-A(J)*B(IN2-J+3)+3.0D0*A(J+3)*B(IN2-J+4)-
     .       A(J+3)*B(IN2-J+2)-A(J+1)*B(IN2-J+4))*BINCOE(IBINC)
         Y = Y + ((A(J+2)-A(J))* (B(IN2-J+1)-B(IN2-J+3))+
     .       (A(J+3)-A(J+1))* (B(IN2-J+2)-B(IN2-J+4)))*BINCOE(IBINC)
  110 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)
C     NORMALIZERS  1.936.. = SQRT(15/4), 6.708...= SQRT(45)
      PENPI = 6.708203932*PEN*Y/ (2.0D0*FACFAC)
      PEN = 1.936491673*PEN*X/FACFAC
      GO TO 200

  120 CONTINUE
C     D - D   X = SIGMA, Y= PI, ZZ= DELTA PARTS.
  130 CONTINUE
      IN2 = N2 - 4
      IBINC = ((IN2-1)*IN2)/2
      Y = ZERO
      ZZ = ZERO
      DO 140 J = 1, IN2
         IBINC = IBINC + 1
         X = X + (A(J+4)* (B(IN2-J+1)-6.0D0*B(IN2-J+3)+9.0D0*B(IN2-J+
     .       5))+A(J)* (B(IN2-J+5)+9.0D0*B(IN2-J+1)-6.0D0*B(IN2-J+3))+
     .       6.0D0*A(J+2)* (6.0D0*B(IN2-J+3)-B(IN2-J+5)-B(IN2-J+1))+
     .       8.0D0*A(J+1)* (3.0D0*B(IN2-J+2)-B(IN2-J+4))+
     .       8.0D0*A(J+3)* (3.0D0*B(IN2-J+4)-B(IN2-J+2)))*BINCOE(IBINC)
         Y = Y + ((A(J+4)-A(J+2))* (B(IN2-J+3)-B(IN2-J+5))+
     .       2.0D0* (A(J+3)-A(J+1))* (B(IN2-J+2)-B(IN2-J+4))+
     .       (A(J+2)-A(J))* (B(IN2-J+1)-B(IN2-J+3)))*BINCOE(IBINC)
         ZZ = ZZ + ((A(J+4)-2.0D0*A(J+2)+A(J))*
     .        (B(IN2-J+1)-2.0D0*B(IN2-J+3)+B(IN2-J+5)))*BINCOE(IBINC)
  140 CONTINUE
      FACFAC = SQRT(FACT(N+N+1))*SQRT(FACT(M+M+1))
      PEN = (R**N2)* (AMU**N)* (BMU**M)*SQRT(AMU*BMU)
C     NORMALIZERS  3.75=15/4,  1.25 = 5/4
      PEND = 3.75D0*PEN*ZZ/ (2.0D0*FACFAC)
      PENPI = 15.0D0*PEN*Y/ (2.0D0*FACFAC)
      PEN = 1.25D0*PEN*X/FACFAC
      GO TO 200

  150 CONTINUE
  160 CONTINUE
C     YOU HAVE THESE FORMULA'S Z. AS FASCINATING AS THEY ARE
C     DON'T WORK THEM OUT AGAIN.
  170 CONTINUE
      WRITE (6,9000)

 9000 FORMAT (' ****ERROR IN PENET, INTEGRAL NOT AVAILABLE****')

      CALL ZEXIT(3)
  180 IF (LA.EQ.LB) THEN
         IF (AMU.NE.BMU) GO TO 190
         PEN = AMU*2.0D0/FN
         IF (LA.GT.0) PENPI = PEN
         IF (LA.GT.1) PEND = PEN
      ELSE
         PEN = ZERO
      END IF
      GO TO 200

  190 CONTINUE
      CALL OVLAP(LH,N-1,LA,AMU,M,LB,BMU,RR,SP,PP,PD,PF,FACT,-TEN,NIN)
      FACTOR = (TWO*AMU)/SQRT(FLOAT(2*N* (2*N-1)))
      PEN = FACTOR*SP
      PENPI = FACTOR*PP
      PEND = FACTOR*PD
      PENDF = FACTOR*PF
c     write(6,*)"PEN,PENPI,PEND,PENDF",PEN,PENPI,PEND,PENDF
C     AR=AMU*R
C      BR=BMU*R
C      M2=2*M+1
C     N2=2*N+1
C      P=(AMU+BMU)*R/TWO
C      T=(AMU-BMU)/(AMU+BMU)
C     CALL AUX(P,T,A,B)
C     T=FACT(M2)*FACT(N2)
C     P=SQRT (AR**N2*BR**M2/T)/R
C     M2=M+N
C      T=ZERO
C     DO 25 J=1,M2
C     N2=M2-J+1
C  25 T=T+DD(J)*B(J)*A(N2)
C     PEN=P*T
  200 CONTINUE
      RETURN

      END
