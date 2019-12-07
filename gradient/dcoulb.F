      FUNCTION DCOULB(M,Z2,Z1,R,FACT,C,NIN)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C     SUBROUTINE TAKES THE DERIVATIVE OF THE COULOMB INTEGRAL WRT R.
C     M=1,2,3,4, ETC FOR (1S,1S]1S,1S), (1S,1S]2S,2S), (2S,2S]2S,2S) ETC
C     ELECTRON LABEL IS (1,1]2,2)
C     Z2 IS EXPONENT GOING WITH LOWER N VALUE
C     R IS IN ANGSTROM UNITS, ANSWER IN A.U.
C     C IS THE COULOMB INTEGRAL, CC THE DERIVATIVE=DCOULB
C     M.C.ZERNER   UNIV. OF GUELPH
C
      DIMENSION    FACT(*), AFAC(10),NIN(*)
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

c     INCLUDE '../include/pptt.cmn'
C****************************************************************
C***      COMMON/PPTT/PP
C****************************************************************
      DATA         AFAC/0.35355339059327D+00, 0.13975424859374D+00,
     .             0.82679728470768D-01, 0.57210228577644D-01,
     .             0.43066295528486D-01, 0.34175618247351D-01,
     .             0.28119734619391D-01, 0.23755299416522D-01,
     .             0.20475351001628D-01, 0.17929397691743D-01/
CGP      DATA         IONE, IZERO/1, 0/

c     write(6,*)"i am in dcoulb"
      RAD = R/BOHRS
      CC = 0.0D0
      S = Z1 + Z2
      D = Z1 - Z2
      ZZ1 = Z1*Z1
      ZZ2 = Z2*Z2
      T = D/S
      RO = S*RAD*0.5D0
      ST = 1.D-05
      IF (M.GT.3) GO TO 160
      GO TO (10,60,110,160,160,160,160,160,160,
     .       160), M

   10 CONTINUE
      IF (ABS(T).GT.ST) GO TO 20
      Y = 0.0D0
      GO TO 30

   20 Y = (ZZ1+ZZ2)/ (ZZ1-ZZ2)
   30 IF (ABS(RO).GT.ST) GO TO 40
      T2 = T*T
      C = (1.0D0-T2)* (5.0D0-T2)*S*0.0625D0
      GO TO 270

   40 IF (ABS(Y).GT.ST) GO TO 50
      RO2 = RO*RO
      RO3 = RO2*RO
      X = 1.0D0 + 1.375D0*RO + 0.75D0*RO2 + 0.16666666D0*RO3
      EX1 = EXP(-RO-RO)
      C = (1.0D0-X*EX1)/RAD
      CC = (-C+S*EX1* (X- (0.6875D0+0.75D0*RO+0.25D0*RO2)))/RAD
      GO TO 270

   50 R1 = Z1*RAD
      R2 = Z2*RAD
      X1 = 1.0D0 - Y
      X2 = 1.0D0 + Y
      X3 = Y*0.25D0
      X4 = 0.5D0 - X3
      X3 = 0.5D0 + X3
      EX1 = EXP(-R1-R1)
      EX2 = EXP(-R2-R2)
      X1 = X1*X1* (X3+0.25D0*R1)*EX1
      X2 = X2*X2* (X4+0.25D0*R2)*EX2
      C = (1.0D0-X1-X2)/RAD
      XD1 = Z1* ((1.D0-Y)**2*EX1/4.D0-2.D0*X1)
      XD2 = Z2* ((1.D0+Y)**2*EX2/4.D0-2.D0*X2)
      CC = (-C-XD1-XD2)/RAD
      GO TO 270

   60 CONTINUE
      IF (ABS(T).GT.ST) GO TO 70
      Y = 0.0D0
      GO TO 80

   70 Y = (ZZ1+ZZ2)/ (ZZ2-ZZ1)
   80 IF (ABS(RO).GT.ST) GO TO 90
      T2 = T*T
      T3 = T2*T
      T4 = T2*T2
      C = (1.0D0-T2)* (14.0D0-7.0D0*T-T2+3.0D0*T3-T4)*S*0.015625D0
      CC = 0.0D0
      GO TO 270

   90 IF (ABS(Y).GT.ST) GO TO 100
      RO2 = RO*RO
      RO3 = RO2*RO
      RO4 = RO2*RO2
      RO5 = RO3*RO2
      X = 1.0D0 + 1.5625D0*RO + 1.125D0*RO2 + 0.47916667D0*RO3 +
     .    0.125D0*RO4 + 0.01666667D0*RO5
      EX1 = EXP(-RO-RO)
      C = (1.0D0-X*EX1)/RAD
      CC = (-C+S*EX1* (X- (0.78125D0+1.125D0*RO+0.71875D0*RO2+
     .     0.25D0*RO3+0.0416667D0*RO4)))/RAD
      GO TO 270

  100 R1 = Z2*RAD
      R2 = Z1*RAD
      Y2 = Y*Y
      Y3 = Y2*Y
      R22 = R2*R2
      R23 = R22*R2
      X11 = 1.0D0 - Y
      EX1 = EXP(-R1-R1)
      X22 = 1.0D0 + Y
      EX2 = EXP(-R2-R2)
      X1 = X11**3* ((1.0D0-5.0D0*Y-4.0D0*Y2)*0.0625D0-Y*R1*0.125D0)*EX1
      X2 = X22*X22* ((15.0D0-22.0D0*Y+15.0D0*Y2-4.0D0*Y3)*0.0625D0+
     .     (3.0D0-3.0D0*Y+Y2)*R2*0.375D0+ (2.0D0-Y)*R22*0.25D0+
     .     R23*8.3333333D-02)*EX2
      DX1 = - (2.0D0*X1+EX1*X11**3*Y/8.0D0)*Z2
      DX2 = -2.0D0*Z1*X2 + Z1*X22**2*EX2*
     .      ((3.0D0-3.0D0*Y+Y2)*0.375D0+0.5D0* (2.0D0-Y)*R2+R22/4.0D0)
      C = (1.0D0-X1-X2)/RAD
      CC = (-C-DX1-DX2)/RAD
      GO TO 270

  110 CONTINUE
      IF (ABS(T).GT.ST) GO TO 120
      Y = 0.0D0
      GO TO 130

  120 Y = (ZZ1+ZZ2)/ (ZZ1-ZZ2)
  130 IF (ABS(RO).GT.ST) GO TO 140
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      C = (1.0D0-T2)* (93.0D0-47.0D0*T2+23.0D0*T4-5.0D0*T6)*S/512.0D0
      CC = 0.0D0
      GO TO 160

  140 IF (ABS(Y).GT.ST) GO TO 150
      RO2 = RO*RO
      RO3 = RO*RO2
      RO4 = RO2*RO2
      RO5 = RO3*RO2
      RO6 = RO3*RO3
      RO7 = RO3*RO4
      X = 1.0D0 + 1.63671875D0*RO + 1.2734375D0*RO2 + 0.61979167D0*RO3 +
     .    0.208333333D0*RO4 + 0.05D0*RO5 + 8.3333333D-03*RO6 +
     .    7.93650793D-04*RO7
      EX1 = EXP(-RO-RO)
      C = (1.0D0-X*EX1)/RAD
      DX1 = S* (1.63671875D0+2.546875D0*RO+1.859375D0*RO2+
     .      0.8333333D0*RO3+0.25D0*RO4+0.05D0*RO5+5.555556D-03*RO6)/
     .      2.0D0
      CC = (-C- (DX1-S*X)*EX1)/RAD
      GO TO 270

  150 Y2 = Y*Y
      Y3 = Y2*Y
      Y4 = Y2*Y2
      R1 = Z1*RAD
      R2 = Z2*RAD
      R12 = R1*R1
      R13 = R12*R1
      R22 = R2*R2
      R23 = R22*R2
      EX1 = EXP(-R1-R1)
      EX2 = EXP(-R2-R2)
      X1 = (8.0D0-Y-27.0D0*Y2-30.0D0*Y3-10.0D0*Y4)/16.0D0 +
     .     (11.0D0-19.0D0*Y-44.0D0*Y2-20.0D0*Y3)*R1/32.0D0 +
     .     (1.0D0-5.0D0*Y-4.0D0*Y2)*R12/16.0D0 - Y*R13/24.0D0
      X1 = X1*EX1
      X2 = (8.0D0+Y-27.0D0*Y2+30.0D0*Y3-10.0D0*Y4)/16.0D0 +
     .     (11.0D0+19.0D0*Y-44.0D0*Y2+20.0D0*Y3)*R2/32.0D0 +
     .     (1.0D0+5.0D0*Y-4.0D0*Y2)*R22/16.0D0 + Y*R23/24.0D0
      X2 = X2*EX2
      C = (1.0D0- (1.0D0-Y)**3*X1- (1.0D0+Y)**3*X2)/RAD
      DX1 = -2.D0*Z1*X1 + EX1* ((11.D0-19.D0*Y-44.D0*Y2-20.D0*Y3)/32.D0+
     .       (1.0D0-5.0D0*Y-4.0D0*Y2)*R1/8.D0-Y*R12/8.D0)*Z1
      DX2 = -2.D0*Z2*X2 + EX2* ((11.D0+19.D0*Y-44.D0*Y2+20.D0*Y3)/32.D0+
     .       (1.0D0+5.0D0*Y-4.0D0*Y2)*R2/8.0D0+Y*R22/8.0D0)*Z2
      CC = (-C- (1.0D0-Y)**3*DX1- (1.0D0+Y)**3*DX2)/RAD
      GO TO 270

  160 IF (M.LE.15) GO TO (170,180,180,190,190,190,200,
     .                    200,200,200,210,210,210,210,
     .                    210), M
      MS = M - 15
      GO TO (220,220,220,220,220,220,230,230,230,
     .       230,230,230,230), MS

  170 ML = 0
      NG = 1
      GO TO 240

  180 ML = 1
      NG = 2
      GO TO 240

  190 ML = 3
      NG = 3
      GO TO 240

  200 ML = 6
      NG = 4
      GO TO 240

  210 ML = 10
      NG = 5
      GO TO 240

  220 ML = 15
      NG = 6
      GO TO 240

  230 ML = 21
      NG = 7
  240 NN = 2* (M-ML)
      FNN = DFLOAT(NN)
      B2 = 2.0D0*Z1
      A2 = 2.0D0*Z2
      MM = 2*NG
      IF (RAD.LT.1.0D-5) GO TO 260
      NGG = MM - 1
      CALL DIFPEN(FACT,Z1,NG,Z1,NG,R,DPEN,1)
C
C  INITIALISE PP TO -10.DO FOR DIFOVL
C
      PP = -10.D0
      X = 0.0D0
      DO 250 I = 1, NN
         NL = NN - I
C
C  GET DERIV OF OVERLAP
C  SS? ARE THE OVERLAPS,  S? ARE THE DERIVS
C
         CALL DIFOVL(IONE,NGG,IZERO,B2,NL,IZERO,A2,R,SS,SSP,SSD,SSF,S,
     .               SP,SD,SF,FACT,NIN,PP)
         FI = DFLOAT(I)/FNN
         AAA = SQRT(FACT(2*NL+1))*FI/ (FACT(NN-I+1)* (2.0D0)** (NN-I))
         X = X + S*AAA
  250 CONTINUE
      AAA = SQRT(Z1**3/Z2)*AFAC(NG)
      CC = DPEN - AAA*X
      GO TO 270

  260 CC = 0.0D0
  270 DCOULB = CC
      RETURN

      END
