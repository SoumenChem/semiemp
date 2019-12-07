      logical FUNCTION COULB(N,AMU,BMU,R,FACT,INT,NIN)
C
C     THIS FUNCTION CALCULATES COULOMB INTEGRALS.
C     ONE AND TWO CENTER INTEGRALS  ONLY. CONVENTION IS (NS,NS/MS,MS) =
C     (NS/MS) WITH N LESS THAN OR EQ. TO M.  AMU IS EXP. FOR NS, BMU IS
C     EXP. FOR MS, RR IS DISTANCE IN ANGSTROMS. IF N = 1 (1S/1S), N=2
C     (1S/2S), 3 (2S/2S), 4 (1S/3S), 5 (2S/3S), 6 (3S/3S),7 ((1S/4S), 8
C     (2S/4S), 9 (3S/4S), 10 (4S/4S).
C
C     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      double precision INT,RR,R
      DIMENSION    A(35), B(35)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
C
      PARAMETER (BOHRS    = 5.29177249D-1)
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON/EV/
C***     .  A(35),B(35),DD(35)
C****************************************************************
c     INCLUDE '../include/pptt.cmn'
C****************************************************************
C***      COMMON/PPTT/
C***     .  PP
C****************************************************************
C
      Dimension    FACT(*), AFAC(10), NIN(*)
C
      DATA         AFAC/0.35355339059327D+00, 0.13975424859374D+00,
     .             0.82679728470768D-01, 0.57210228577644D-01,
     .             0.43066295528486D-01, 0.34175618247351D-01,
     .             0.28119734619391D-01, 0.23755299416522D-01,
     .             0.20475351001628D-01, 0.17929397691743D-01/
C      DATA         IZERO/0/IONE/1/
C
        COULB=.false.
        INT=0.0
        D = 0.0d0
c     write(6,*)"i am in coulb",N,AMU,BMU,RR
      RR = R!*BOHRS
      IF (ABS(R).LT.1.D-7) R = ZERO
      TT = (AMU-BMU)/ (AMU+BMU)
      T2 = TT**2
      IF (T2.LE.1.0D-6) THEN
         TT = ZERO
         T2 = ZERO
         AMU = BMU
      END IF
   10 CA = R*AMU
      CB = R*BMU
      C = 0.5D0* (AMU+BMU)*R
      IF (AMU-BMU) 20, 30, 20
   20 D = (AMU**2+BMU**2)/ (AMU**2-BMU**2)
c     write(6,*)" i am in 1"
   30 E = 0.5D0* (AMU+BMU)
      IF (N.GE.10) GO TO 290
      GO TO (40,90,120,170,200,230,260,270,280,
     .       290), N
c     write(6,*)" i am in 2"
   40 IF (R) 50, 80, 50
   50 IF (AMU-BMU) 60, 70, 60
   60 COUL = (E/C)* (ONE- ((ONE-D)**2* ((TWO+D)+CA)*EXP(-CA-CA)+ (ONE+
     .       D)**2* ((TWO-D)+CB)*EXP(-CB-CB))/FOUR)
c     write(6,*)" i am in 3"
      GO TO 410

   70 COUL = (E/C)* (ONE- (ONE+11.0D0*C/EIGHT+0.75D0*C**2+C**3/SIX)*
     .       EXP(-C-C))
c     write(6,*)" i am in 4"
      GO TO 410

   80 COUL = (ONE-T2)* (FIVE-T2)*E/EIGHT
c     write(6,*)" i am in 5"
      GO TO 410

   90 IF (R) 110, 100, 110
  100 COUL = (ONE-T2)*E* (14.0D0-SEVEN*TT-T2+THREE*T2*TT-T2*T2)/32.0D0
c     write(6,*)" i am in 6"
      GO TO 410

  110 IF (T2.LE.0.0001) GO TO 260
      COUL = (E/C)* (ONE- ((ONE-D)**3)*
     .       ((ONE-FIVE*D-FOUR*D**2)/16.D0-0.125D0*D*CA)*EXP(-CA-CA)-
     .       (ONE+D)**2* ((15.D0-22.0D0*D+15.0D0*D**2-FOUR*D**3)/16.0D0+
     .       THREE* (THREE-THREE*D+D**2)*CB/EIGHT+0.25D0* (TWO-D)*CB**2+
     .       CB**3/12.0D0)*EXP(-CB-CB))
c     write(6,*)" i am in 7"
      GO TO 410

  120 IF (R) 140, 130, 140
  130 COUL = (ONE-T2)* (93.0D0-47.0D0*T2+23.0D0*T2*T2-FIVE*T2**3)*E/
     .       256.0D0
c     write(6,*)" i am in 8"
      GO TO 410

  140 IF (AMU-BMU) 160, 150, 160
  150 COUL = (E/C)* (ONE- (ONE+419.0D0*C/256.0D0+163.0D0*C**2/128.0D0+
     .       119.0D0*C**3/192.0D0+FIVE*C**4/24.0D0+C**5/20.0D0+
     .       C**6/120.0D0+C**7/1260.0D0)*EXP(-C-C))
c     write(6,*)" i am in 9"
      GO TO 410

  160 COUL = (E/C)* (ONE- (ONE-D)**3* ((EIGHT-D-27.0D0*D**2-30.0D0*D**3-
     .       TEN*D**4)/16.0D0+ (11.0D0-19.0D0*D-44.0D0*D**2-
     .       20.0D0*D**3)*CA/32.0D0+ (ONE-FIVE*D-
     .       FOUR*D**2)*CA**2/16.0D0-D*CA**3/24.0D0)*EXP(-CA-CA)-
     .       ((ONE+D)**3)* ((EIGHT+D-27.0D0*D**2+30.0D0*D**3-
     .       TEN*D**4)/16.0D0+ (11.0D0+19.0D0*D-44.0D0*D**2+
     .       20.0D0*D**3)*CB/32.0D0+ (ONE+FIVE*D-
     .       FOUR*D**2)*CB**2/16.0D0+D*CB**3/24.0D0)*EXP(-CB-CB))
c     write(6,*)" i am in 10",COUL
      GO TO 410

  170 IF (R) 190, 180, 190
  180 AA = (BMU/ (AMU+BMU))**6
      COUL = BMU* (ONE-AA)/THREE - AMU*AA*BMU/ (AMU+BMU)
c     write(6,*)" i am in 11"
      GO TO 410

  190 CALL AUXC(A,B,CA)
      Z = 720.0D0* (A(2)*B(1)-A(1)*B(2))
      P = TWO*C
      TT = -TT
      CALL AUX(P,TT,A,B)
      X = CB
      Y = (X**5)* (A(7)*B(1)+FOUR*A(6)*B(2)+FIVE*A(5)*B(3)-
     .    FIVE*A(3)*B(5)-FOUR*A(2)*B(6)-A(1)*B(7)) +
     .    TEN* (X**4)* (A(6)*B(1)+THREE*A(5)*B(2)+TWO*A(4)*B(3)-
     .    TWO*A(3)*B(4)-THREE*A(2)*B(5)-A(1)*B(6)) +
     .    60.0D0* (X**3)* (A(5)*B(1)+TWO*A(4)*B(2)-TWO*A(2)*B(4)-
     .    A(1)*B(5)) + 240.0D0* (X**2)* (A(4)*B(1)+A(3)*B(2)-A(2)*B(3)-
     .    A(1)*B(4)) + 600.0D0*X* (A(3)*B(1)-A(1)*B(3)) +
     .    720.0D0* (A(2)*B(1)-A(1)*B(2))
      COUL = (Z-Y)* (R**2)* (AMU**3)/1440.0D0
      GO TO 410

  200 IF (R) 220, 210, 220
  210 AA = BMU/ (AMU+BMU)
      AA6 = AA**6
      ABU = AMU/BMU
      COUL = (ONE-AA6)*BMU/THREE - 0.50D0*AA6*AA*AMU*
     .       (28.0D0*ABU*ABU*AA*AA/THREE+SEVEN*ABU*AA+THREE)
c     write(6,*)" i am in 12"
      GO TO 410

  220 CALL AUXC(A,B,CA)
      Z = 720.D0* (A(4)*B(1)-THREE*A(3)*B(2)+THREE*A(2)*B(3)-A(1)*B(4))
      P = TWO*C
      TT = -TT
      CALL AUX(P,TT,A,B)
      X = CB
      Y = - (X**5)* (A(9)*B(1)+TWO*A(8)*B(2)-TWO*A(7)*B(3)-
     .    SIX*A(6)*B(4)+SIX*A(4)*B(6)+TWO*A(3)*B(7)-TWO*A(2)*B(8)-
     .    A(1)*B(9)) - TEN* (A(8)*B(1)+A(7)*B(2)-THREE*A(6)*B(3)-
     .    THREE*A(5)*B(4)+THREE*A(4)*B(5)+THREE*A(3)*B(6)-A(2)*B(7)-
     .    A(1)*B(8))*X**4 - 60.0D0* (A(7)*B(1)-THREE*A(5)*B(3)+
     .    THREE*A(3)*B(5)-A(1)*B(7))*X**3
      YY = -240.0D0* (A(6)*B(1)-A(5)*B(2)-TWO*A(4)*B(3)+TWO*A(3)*B(4)+
     .     A(2)*B(5)-A(1)*B(6))*X**2 - 600.0D0*
     .     (A(5)*B(1)-TWO*A(4)*B(2)+TWO*A(2)*B(4)-A(1)*B(5))*X -
     .     720.0D0* (A(4)*B(1)-THREE*A(3)*B(2)+THREE*A(2)*B(3)-
     .     A(1)*B(4))
      Y = Y + YY
      COUL = (Z+Y)* (R**4)* (AMU**5)/17280.0D0
c     write(6,*)" i am in 13"
      GO TO 410

  230 IF (R) 250, 240, 250
  240 T = BMU + AMU
      COUL = ((BMU**7)/ (THREE*T**11))*
     .       (T**11/ (BMU**6)-42.0D0*AMU**5-42.0D0*T*AMU**4-
     .       28.0D0*AMU**3* (T**2)-14.0D0* (AMU**2)* (T**3)-
     .       FIVE*AMU* (T**4)-T**5)
c     write(6,*)" i am in 14"
      GO TO 410

  250 CALL AUXC(A,B,CB)
      Z = A(6)*B(1) - FIVE*A(5)*B(2) + TEN*A(4)*B(3) - TEN*A(3)*B(4) +
     .    FIVE*A(2)*B(5) - A(1)*B(6)
      P = TWO*C
      CALL AUX(P,TT,A,B)
      X = CA
      U = X**3
      Y = (X**5)* (A(11)*B(1)-FIVE*A(9)*B(3)+TEN*A(7)*B(5)-
     .    TEN*A(5)*B(7)+FIVE*A(3)*B(9)-A(1)*B(11)) +
     .    TEN* (X**4)* (A(10)*B(1)-A(9)*B(2)-FOUR*A(8)*B(3)+
     .    FOUR*A(7)*B(4)+SIX*A(6)*B(5)-SIX*A(5)*B(6)-FOUR*A(4)*B(7)+
     .    FOUR*A(3)*B(8)+A(2)*B(9)-A(1)*B(10))
      YY = 60.0D0*U* (A(9)*B(1)-TWO*A(8)*B(2)-TWO*A(7)*B(3)+
     .     SIX*A(6)*B(4)-SIX*A(4)*B(6)+TWO*A(3)*B(7)+TWO*A(2)*B(8)-
     .     A(1)*B(9)) + 240.0D0* (X**2)* (A(8)*B(1)-THREE*A(7)*B(2)+
     .     A(6)*B(3)+FIVE*A(5)*B(4)-FIVE*A(4)*B(5)-A(3)*B(6)+
     .     THREE*A(2)*B(7)-A(1)*B(8)) + 600.0D0*X*
     .     (A(7)*B(1)-FOUR*A(6)*B(2)+FIVE*A(5)*B(3)-FIVE*A(3)*B(5)+
     .     FOUR*A(2)*B(6)-A(1)*B(7)) + 720.0D0*
     .     (A(6)*B(1)-FIVE*A(5)*B(2)+TEN*A(4)*B(3)-TEN*A(3)*B(4)+
     .     FIVE*A(2)*B(5)-A(1)*B(6))
      Y = Y + YY
      COUL = ((BMU**7)* (R**6)/ (518400.0D0))* (720.0D0*Z-Y)
c     write(6,*)" i am in 15"
      GO TO 410

  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
  290 CONTINUE
c     write(6,*)"i am in 16"
      NS = N - 15
      IF (N.LE.15) GO TO (300,310,310,320,320,320,330,
     .                    330,330,330,340,340,340,340,
     .                    340), N
      IF (N.GT.15) GO TO (350,350,350,350,350,350,360,
     .                    360,360,360,360,360,360), NS
  300 ML = 0
      NG = 1
      GO TO 370

  310 ML = 1
      NG = 2
      GO TO 370

  320 ML = 3
      NG = 3
      GO TO 370

  330 ML = 6
      NG = 4
      GO TO 370

  340 ML = 10
      NG = 5
      GO TO 370

  350 ML = 15
      NG = 6
      GO TO 370

  360 ML = 21
      NG = 7
  370 NN = 2* (N-ML)
      FNN = DFLOAT(NN)
      A2 = TWO*AMU
      B2 = TWO*BMU
      MM = 2*NG
      IF (R.LT.1.D-5) GO TO 390
      NGG = 2*NG - 1
      CALL PENET(FACT,BMU,NG,IZERO,BMU,NG,IZERO,RR,PEN,PENPI,PEND,
     x           PENDF,IONE,NIN)
C
C  SET PP TO -10.D0 SO AS TO INITIALISE OVLAP
C
c     write(6,*)"i am in 17"
      PP = -10.D0
      X = ZERO
      DO 380 I = 1, NN
         NL = NN - I
         CALL OVLAP(IONE,NGG,IZERO,B2,NL,IZERO,A2,RR,S,SP,SD,SF,FACT,
     1   PP,NIN)
         FI = DFLOAT(I)/FNN
         AAA = SQRT(FACT(2*NL+1))*FI/ (FACT(NN-I+1)* (TWO)** (NN-I))
         X = X + S*AAA
c     write(6,*)"after ovlap"
  380 CONTINUE
      AAA = SQRT(BMU**3/AMU)*AFAC(NG)
      COUL = PEN - X*AAA
      GO TO 410
C
C     GENERAL ONE CENTER INTEGRAL OF FORM (NS,NS/MS,MS), WHERE N=NN/2,
C     M=MM/2, AND AMU IS ASSOCIATED WITH N, BMU WITH M.
  390 A2 = TWO*AMU
      B2 = TWO*BMU
      P = A2** (NN+1)*B2** (MM+1)/FACT(MM+1)
      COUL = (FACT(MM)/A2** (NN+1))* (ONE/ (B2**MM)-ONE/ ((A2+B2)**MM))
      FNN = DFLOAT(NN)
      AA = A2* (A2+B2)** (MM+NN)
      DO 400 I = 2, NN
         II = I - 1
         T = DFLOAT(II)/ (FNN*FACT(NN-II+1))
         AA = AA*A2/ (A2+B2)
         COUL = COUL - T*FACT(MM+NN-II)/AA
  400 CONTINUE
      COUL = P*COUL
  410 INT = COUL
      END
c-------------------------------------------------------------------------------------
c   subroutine AUX
c------------------------------------------------------------------------------------
      SUBROUTINE AUX(PP,TT,A,B)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
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
      DIMENSION    A(35), B(35)
      DATA         UPPLIM/1.0D+36/
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
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
c     write(6,*)"i am in aux"
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
      GO TO 270

  250 CONTINUE
      DO 260 I = 1, IX
         A(I) = 0.0D0
         B(I) = 0.0D0
  260 CONTINUE
  270 RETURN

      END
c----------------------------------------------------------------------------------------
c  AUXC subroutine
c---------------------------------------------------------------------------------------
      SUBROUTINE AUXC(A,B,T)
C
C     THIS FUNCTION CALCULATES A(T) AND B(-T) FOR THE COULOMB FUNCTION.
C     AUXC0002
C
c      INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      DIMENSION    A(35), B(35)
C
      X = EXP(T)
      Y = EXP(-T)
      A(1) = Y/T
      B(1) = (X-Y)/T
C
C  ROUTINE ONLY CALLED FROM COULB
C  WHICH ONLY USES UPTO N=6
C
c     write(6,*)"i am in auxc"
      DO 10 N = 2, 6
         G = DFLOAT(N-1)
         A(N) = (Y+G*A(N-1))/T
         B(N) = - ((-1.D0)** (N-1)*Y-X+G*B(N-1))/T
   10 CONTINUE
      RETURN

      END
C----------------------------------------------------------------------------------------------
C COUL SUBROUTINE
C----------------------------------------------------------------------------------------------

                                       
