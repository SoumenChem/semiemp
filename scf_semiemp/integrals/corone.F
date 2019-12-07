      SUBROUTINE CORONE(SH,RAD,F,ALP,FACT,NB,NU,NA,NW,KAN,NP,NIN,ISW1)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      CHARACTER*10  ATOMK, ATOMKX, ATOMM, ATOMMX
C
C     THIS SUBROUTINE CALC. THE OVERLAP MATRIX ELEMENTS BETWEEN CORE A.O
C
C****************************************************************
C***      COMMON/PPTT/
C***     .  PP
C****************************************************************
C
      DIMENSION    D(8), E(8), SH(*), RAD(NA,NA), F(*), NU(*), NW(*),
     .             KAN(*), NP(*), ALP(*), FACT(*)
      DIMENSION    KFLAG(3), MFLAG(3)
      PARAMETER (ZERO =  0.0D0, ONE  =  1.0D0)
      DATA         MONE/-1.D0/


C
C
C  INTRODUCE ATOMK AND ATOMM TO TEST IF TREATING THE SAME SHELL
C
C    KFLAG(1)  --> K ATOM NO
C    KFLAG(2)  --> K PQN
C    KFLAG(3)  --> K LQN
C    MFLAG(1)  --> M ATOM NO
C    MFLAG(2)  --> M PQN
C    MFLAG(3)  --> M LQN
C  INITIALISE ATOMKX AND ATOMMX TO ZERO
      EDEL = 0.0d0
      EPI = 0.0d0
      MM = 0
      M = 0
      SS = 0.0d0
      SP = 0.0d0
      SD = 0.0d0
      SF = 0.0d0
      ATOMKX = '  0  0  0'
      ATOMMX = '  0  0  0'
      ATOMK  = '  0  0  0'
      ATOMM  = '  0  0  0'
c     write(6,*)"i am in corone"
C
C  SET PP TO -10.D0 SO AS TO INITIALISE OVLAP
      PP = -10.D0
      NA1 = NA*NA
      NA2 = NA1 + NA1
      IJ = 0
      NBT = NB
      DO 330 I = 1, NBT
         K = NW(I)
         KA = NU(I)
         KFLAG(1) = KA
         KX = KAN(KA)
         IF (K-1) 10, 20, 20
   10    LI = NP(I)
         NL = 0
         GO TO 70

   20    IF (K-3) 30, 30, 40
   30    LI = NP(I) + 3
         NL = 1
         GO TO 70

   40    IF (K-8) 50, 50, 60
   50    LI = 8
         NL = 2
         GO TO 70

   60    LI = 15
         NL = 3
   70    KFLAG(3) = NL
         NT1 = NP(I)
         KFLAG(2) = NT1
C
C      DO AN INTERNAL WRITE TO CONVERT THE SAVED FLAG VALUES TO
C      CHARACTER STRINGS.
C
c        WRITE (ATOMK,'(3I3)') KFLAG
         NLP1 = NL + 1
         DO 320 J = 1, I
            IJ = IJ + 1
            IF (I.EQ.J) THEN
               SH(IJ) = ONE
               GO TO 320

            END IF
            M = NW(J)
            IF (ISW1.EQ.7) M = 3
            MA = NU(J)
            MFLAG(1) = MA
            MX = KAN(MA)
            IF (MA-KA) 90, 80, 90
   80       CONTINUE
            SH(IJ) = 0.0D0
            GO TO 310

   90       CONTINUE
            IF (K-M) 100, 100, 110
  100       LA = M
            LH = K
            K1M = KA + (MA-1)*NA
            GO TO 120

  110       LA = K
            LH = M
            K1M = MA + (KA-1)*NA
  120       U1 = F(K1M)
            U2 = F(K1M+NA1)
            U3 = F(K1M+NA2)
            CALL transmat(U1,U2,U3,LH,D)
            CALL transmat(U1,U2,U3,LA,E)
            IF (M-1) 130, 140, 140
  130       LJ = NP(J)
            ML = 0
            GO TO 190

  140       IF (M-3) 150, 150, 160
  150       LJ = NP(J) + 3
            ML = 1
            GO TO 190

  160       IF (M-8) 170, 170, 180
  170       LJ = 8
            ML = 2
            GO TO 190

  180       LJ = 15
            ML = 3
  190       MFLAG(3) = ML
            NT2 = NP(J)
            MFLAG(2) = NT2
            MLP1 = ML + 1
            WRITE (ATOMM,'(3I3)') MFLAG
C
C  GET PHASE FOR E TRANSFORM
C     PHASE (-1)**(L-M) ON LA.GE.THE OTHER L VALUE
C
            IF (LA.EQ.M) MM = MLP1
            IF (LA.EQ.K) MM = NLP1
            GO TO (230,210,220,200), MM
C
C  F-ORBS
C
  200       EPHI = ONE
            EDEL = MONE
C
C P-ORBS
C
  210       EPI = ONE
            ESIG = MONE
            GO TO 240
C
C  D-ORBS
C
  220       EDEL = ONE
            EPI = MONE
C
C  S-ORBS
C
  230       ESIG = ONE
  240       CONTINUE
            KAMA = MA + (KA* (KA-1))/2
            R = RAD(KA,MA)
            IF (ATOMK.EQ.ATOMKX .AND. ATOMM.EQ.ATOMMX) GO TO 270
            ATOMKX = ATOMK
            ATOMMX = ATOMM
            SS = 0.0D0
            SP = 0.0D0
            SD = 0.0D0
            SF = 0.0D0
C      DO 38 ISS=1,2
C      AMU=ALP(I,ISS+1)
C      IF(AMU.LT.1.0D-4) GO TO 38
C      ACU=ALP(I,ISS+4)
C      DO 32 JSS=1,2
C      BMU=ALP(J,JSS+1)
C      IF(BMU.LT.1.0D-4) GO TO 32
C      BCU=ALP(J,JSS+4)
C
C     ASSUME SINGLE ZETA FOR CORE
C
            ACU = ONE
            BCU = ONE
            AMU = ALP(I)
            BMU = ALP(J)
            CALL OVLAP(LH,NT1,NL,AMU,NT2,ML,BMU,R,SSS,SSP,SSD,SSF,FACT,
     .      PP,NIN)
C           WRITE(6,2222) I, J, NT1, NT2, NL,ML,AMU,BMU,R,FACT(3),SSS
C2222       FORMAT(' CORONE', 6I4, 4F10.6,F12.8)
            ABC= 0.0d0
            ABC = ACU*BCU
            SS = SS + SSS*ABC
            SP = SP + SSP*ABC
            SD = SD + SSD*ABC
            SF = SF + SSF*ABC
  250       CONTINUE
  260       CONTINUE
C
C
C  FORM SIGMA OVERLAP
C
  270       SH(IJ) = E(1)*D(1)*SS
C           WRITE(6,*) ' CORONE', IJ, SH(IJ)
            MINM = MIN0(NL,ML) + 1
C
C  CHECK HOW MANY TERMS IN OVERLAP
C
            GO TO (310,300,290,280), MINM
C
C  PHI OVER
C
  280       OV = (E(6)*D(6)+E(7)*D(7))*SF*EPHI
            SH(IJ) = SH(IJ) + OV
C
C  DEL OVER
C
  290       OV = (E(4)*D(4)+E(5)*D(5))*SD*EDEL
            SH(IJ) = SH(IJ) + OV
C
C  PI OVER
C
  300       OV = (E(2)*D(2)+E(3)*D(3))*SP*EPI
            SH(IJ) = SH(IJ) + OV
C
  310       CONTINUE
  320    CONTINUE
  330 CONTINUE
      RETURN

      END
