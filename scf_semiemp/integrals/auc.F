      SUBROUTINE AUC(PP,TT,A,B)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     THIS SUBROUTINE CALCULATES A AND B FNS. FOR MOLECULAR INTEGRALS.
C     B-FNS. OF ARGUMENT PP*TT AND WITH INDICES OF 0 TO IX, A FNS. OF
C     ARGUMENT PP AND INDICES FROM 0 TO IX.
C **  REFERENCES:
C **  (A) MULLIKAN/RIEKE/ORLOFF/ORLOFF (1949) JCP 17, 1248;
C **      PP (=P OR =RHO) = 1/2*[MU(A)+MU(B)]*R/A(H)
C **      TT (=T OR =TAU) = [MU(A)-MU(B)]/[MU(A)+MU(B)]
C **  (B) KOTANI/AMEMIYA/SIMOSE (1940) PROC. PHYS. MAT. SOC. JAPAN 22, 1
C **
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      DIMENSION    A(35), B(35)
      DATA         UPPLIM/1.0D+36/
      parameter   (ZERO=0.0D0, ONE = 1.0D0, TWO = 2.0D0)
C
C
      IX = 32
      RHO2 = PP*TT
      IR = ABS(RHO2+RHO2)
      IF (IR.GT.170) GO TO 250
      IS = MIN0(IR+1,15)
C
C     FIRST THE MORE DIFFICULT B-FN.
C
      IF (ABS(RHO2).LT.0.001) GO TO 210
      IF (PP.GT.40.0D0) GO TO 250
   10 D = EXP(RHO2)
      H = EXP(-RHO2)
      R = D - H
C
C     IF THE VALUE OF RHO IS TOO SMALL THE SINH MUST BE OBTAINED BY
C     SUMMING THE INFINITE SERIES RATHER THAN BY ADDITION OF TWO
C     EXPONENTIALS
C
c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in aux"
      IF (ABS(R)-0.1) 20, 50, 50
   20 RA = ZERO
      T = RHO2
      DO 30 I = 2, 25
         IF (ABS(T).LT.1.D-18) GO TO 40
         T = T*RHO2*RHO2/DFLOAT((I+I-1)* (I+I-2))
         RA = RA + T
   30 CONTINUE
   40 CONTINUE
      R = (RA+RHO2)*TWO
   50 B(1) = R/RHO2
C
C     AS MANY SUCCESSIVE B FUNCTIONS ARE GENERATED FROM B(0) BY THE
C     RECURSION FORMULA AS ACCURACY WILL PERMIT.
C
      IXS = IX
      DO 200 I = 2, IXS, IS
         IF (IR) 60, 100, 60
   60    IL = IS - 1
         IF ((I+IL).GT.IX) IL = IX - I + 1
         DO 90 J = 1, IL
            K = I + J - 1
C             IF((-1)**K) 29,29,30
            IF (MOD(K,2)) 80, 80, 70
   70       B(K) = (R+DFLOAT(K-1)*B(K-1))/RHO2
            GO TO 90

   80       B(K) = - (D+H-DFLOAT(K-1)*B(K-1))/RHO2
   90    CONTINUE
  100    IN = I + IS - 1
         IF (IN-IX) 110, 110, 230
C
C     AFTER THE RECURRENCE FORMULA HAS BEEN APPLIED AN APPROPRIATE NO.
C     OF TIMES THE NEXT B FUNCTION IS OBTAINED BY SUMMING THE INFINITE
C     SERIES.
C
C   39 IF((-1)**IN) 44,44,42
  110    IF (MOD(IN,2)) 120, 120, 160
  120    TR = RHO2
  130    B(IN) = - (TR+TR)/DFLOAT(IN+1)
         DO 150 J = 1, 500
            TR = TR*RHO2*RHO2/DFLOAT((J+J)* (J+J+1))
C      NOTE ACCURACY CRITERION
C      IF(ABS(TR/B(IN))-0.000011 ) 51,51,43
            IF (ABS(TR)-1.D-7*ABS(B(IN))) 200, 200, 140
  140       B(IN) = B(IN) - (TR+TR)/DFLOAT(IN+1+J+J)
  150    CONTINUE
  160    TR = ONE
  170    B(IN) = (TR+TR)/DFLOAT(IN)
         DO 190 J = 1, 500
            TR = TR*RHO2*RHO2/DFLOAT((J+J)* (J+J-1))
C      NOTE ACCURACY CRITERION
C      IF(ABS(TR/B(IN))-0.00001  ) 51,51,46
            IF (ABS(TR)-1.D-7*ABS(B(IN))) 200, 200, 180
  180       B(IN) = B(IN) + (TR+TR)/DFLOAT(IN+J+J)
  190    CONTINUE
  200 CONTINUE
      GO TO 230
C
C     IF THE ARGUMENT OF THE B-FN IS ZERO A SEPARATE FORMULA IS USED.
C
  210 CONTINUE
      JY = IX/2
      DO 220 I = 1, JY
         B(I+I) = ZERO
         B(I+I-1) = TWO/DFLOAT(I+I-1)
  220 CONTINUE
  230 CONTINUE
c     write(6,*)"B(IJ)",(B(IJ),IJ=1,35)

C
C     NOW THE A-FUNCTIONS
C
      C = EXP(-PP)
      A(1) = C/PP
      IXM = IX - 3
      DO 240 I = 2, IXM
         IF (A(I-1).GE.UPPLIM) THEN
            A(I) = A(I-1)
            GO TO 240

         END IF
         A(I) = (DFLOAT(I-1)*A(I-1)+C)/PP
  240 CONTINUE
c     write(6,*)"A(I)",(A(I),I=1,35)
      GO TO 270

  250 CONTINUE
      DO 260 I = 1, IX
         A(I) = 0.0D0
         B(I) = 0.0D0
  260 CONTINUE
  270 RETURN

      END

