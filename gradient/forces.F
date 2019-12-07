      SUBROUTINE DFORCES(NA,DG,RAD,FA,H,F,FALPH,FZ,DPEN,ISWANN,NIN,FACT,
     1 NTYP,KAN,NP,NU,NW,NG,CO,Z12,ETS,ETP,ETD,ETF,NB,ALPHA,GRAD,ISPEC)

c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     LOGICAL DRYONE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
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
C     COMPLETES THE CALCULATION OF THE FORCES.
C     ON ENTRY FZ CONTAINS THE ONE ELECTRON PART OF THE FORCES.
C     F IS THE TOTAL DENSITY,FALPH IS THE ALPHA DENSITY
C     THIS IS SET IN DELEN WHICH CALLS FORCES.
C
c     INCLUDE '../include/cf.cmn'
c     INCLUDE '../include/dipgrn.cmn'
c     INCLUDE '../include/rfield.cmn'
c     INCLUDE '../include/rhf.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
C
      LOGICAL DRYONE

      DIMENSION DG(*),RAD(NA,NA),FA(*),H(*),F(*),FZ(*),FALPH(*)
      DIMENSION DPEN(NB,NB), ALPHA(NB,6),GRAD(3,NA)
      DIMENSION FACT(*), NIN(*), NTYP(*), KAN(*), NP(*), NU(*), NW(*)
      DIMENSION NG(*),CO(3,NA),Z12(*),ETS(*),ETP(*),ETD(*),ETF(*)
      DATA      THRESH/1.0D-13/
C
C     NOW EVALUATE THE FORCES
C
C     F IS THE TOTAL DENSITY,FALPH IS THE ALPHA DENSITY
c     write(6,*)"i am in forces"
      NA1 = NA
      NA2 = 2*NA
      NA11 = NA*NA
      NA12 = 2*NA11
      NRR = (NB*(NB+1))/2
      NDN = 0
      KA = 0
      IPRINT = 0
      XX = 0.0d0
C     TRUST NO_ONE  GET THE TOTAL DENSITY
c     CALL FINDP(F,NRR,ISWANN,NDN,IONE)
C     GET FALPH FROM UNIT FOR THE UHF OPEN SHELL CASES
      IF (ISWANN.EQ.0) GO TO 10
      IF (ISWANN.EQ.3) GO TO 10
c     CALL FINDP(FALPH,NRR,ISWANN,NDN,ITHREE)
   10 CONTINUE
C
      JD = 0
      NTA = 0
      DO 250 J = 1, NB
         JA = NU(J)
         JD = JD + J
         KD = 0
         LJ = NW(J)
         X = ZERO
         Y = ZERO
         Z = ZERO
         XXX = ZERO
         YYY = ZERO
         ZZZ = ZERO
         NTB = 0
         NGJ = NG(J)
C***********
         NPQN = NP(J)
C***********
         IF ((NTYP(JA).GT.3) .AND. LJ.EQ.0) NTA = NTA + 1
         DO 150 K = 1, NB
            KA = NU(K)
            KD = KD + K
            LK = NW(K)
            NGK = NG(K)
            IF ((NTYP(KA).GT.3) .AND. LK.EQ.0) NTB = NTB + 1
            IF (KA-JA) 20, 150, 30
   20       KAKA = (JA* (JA-1))/2 + KA
            IJ = JD - J + K
            KC = (NGJ* (NGJ-1))/2 + NGK
            GO TO 40

   30       KAKA = (KA* (KA-1))/2 + JA
            IJ = KD - K + J
            KC = (NGK* (NGK-1))/2 + NGJ
   40       R = RAD(KA,JA)!/BOHRS
c           write(6,*)"R in force",R
            ITJ = JA + (KA-1)*NA
            FAC1 = 0.5D0*F(IJ)**2
            IF (ISWANN.NE.0 .AND. ISWANN.NE.3) FAC1 = FALPH(IJ)**2 +
     .          (F(IJ)-FALPH(IJ))**2
            FAC = -FAC1*DG(KC) + F(KD)*F(JD)*DG(KC)
            CT = FA(ITJ)
            ST = ONE - CT*CT
            IF (THRESH-ST) 60, 50, 50
   50       ST = ZERO
            GO TO 70

   60       ST = SQRT(ST)
   70       CF = FA(ITJ+NA11)
            SF = FA(ITJ+NA12)
            XXX = XXX + FAC*ST*CF
            YYY = YYY + FAC*ST*SF
            ZZZ = ZZZ + FAC*CT
            IF (LK.NE.0) GO TO 120
C           NEW TREATMENT OF PENETRATION INTEGRALS  OCT 6 83
C           KP LOCATES THE "P" SYMMETRY INTEGRAL, KC THE "S".
C           IF(NGJ.LT.NGK)KP=NGJ+((NGK+1)*NGK)/2
C           IF(NGJ.GT.NGK)KP=KC+1
            IF (NTYP(KA).LT.4) THEN
               IF (Z12(KA).LE.2.002) THEN
                  KP = KC
                  XX = -F(JD)*DPEN(NGJ,NGK)*Z12(KA)
               ELSE
                  XX = -F(JD)* (2.D0*DPEN(NGJ,NGK)+
     .                 (Z12(KA)-2.D0)*DPEN(NGJ,NGK+1))
               END IF
            END IF
C
C           TREAT D(N-2)S(2) CONFIG
C
            IF (NTYP(KA).NE.4) GO TO 80
C           IF(NGJ.LT.NGK) KDD=NGJ+((NGK+2)*(NGK+1))/2
C           IF(NGJ.GT.NGK) KDD=KC+2
            XX = -F(JD)* (ETS(NTB)*DPEN(NGJ,NGK)+
     .           ETD(NTB)*DPEN(NGJ,NGK+2))
            GO TO 90
C
C           TREAT S(2)D(1)F(N-3) CONFIG
   80       IF (NTYP(KA).NE.5) GO TO 90
C           IF(NGJ.LT.NG(K)) THEN
C              KDD=NGJ+((NG(K)+2)*(NG(K)+1))/2
C              KF=NGJ+(NG(K)+3)*(NG(K)+2)/2
C           ELSE
C              KDD=KC+2
C              KF=KC+3
C           END IF
            IF (KAN(KA).LT.89) THEN
               XX = -F(JD)* (ETS(NTB)*DPEN(NGJ,NGK)+
     .              ETD(NTB)*DPEN(NGJ,NGK+2)+ETF(NTB)*DPEN(NGJ,NGK+3))
            ELSE
               XX = -F(JD)* (ETS(NTB)*DPEN(NGJ,NGK)+
     .              ETP(NTB)*DPEN(NGJ,NGK+1)+ETD(NTB)*DPEN(NGJ,NGK+2)+
     .              ETF(NTB)*DPEN(NGJ,NGK+3))
            END IF
   90       CONTINUE
            XXX = XXX + XX*ST*CF
            YYY = YYY + XX*ST*SF
            ZZZ = ZZZ + XX*CT
CRDB        IF(LJ.NE.0) GO TO 210
            IF (LJ.NE.0) GO TO 150
            IF (ISPEC.EQ.0 .OR. ISPEC.EQ.3) GO TO 100
            XX = Z12(JA)*Z12(KA)*DG(KC)
            GO TO 110

  100       CONTINUE
            IF (R.LT.1.D-4) GO TO 110
            XX = -Z12(JA)*Z12(KA)/ (R*R)
  110       CONTINUE
            XXX = XXX + XX*ST*CF
            YYY = YYY + XX*ST*SF
            ZZZ = ZZZ + XX*CT
  120       CONTINUE
            IF (LJ.NE.0) GO TO 150
C           NEW TREATMENT OF PENETRATION INTEGRAL
C           IF(NGJ.LT.NGK)KP=KC+1
C           IF(NGJ.GT.NGK)KP=NGK+((NGJ+1)*NGJ)/2
            IF (NTYP(JA).LT.4) THEN
               IF (Z12(JA).LE.2.002) THEN
C                 KP=KC
                  XX = F(KD)*DPEN(NGK,NGJ)*Z12(JA)
               ELSE
                  XX = F(KD)* (2.D0*DPEN(NGK,NGJ)+
     .                 (Z12(JA)-2.D0)*DPEN(NGK,NGJ+1))
               END IF
            END IF
C
C           TREAT S(2)D(N-2) CONFIG
            IF (NTYP(JA).NE.4) GO TO 130
C           IF(NGJ.LT.NG(K)) KDD=KC+2
C           IF(NGJ.GT.NG(K)) KDD=NG(K)+((NGJ+2)*(NGJ+1))/2
            XX = F(KD)* (ETS(NTA)*DPEN(NGK,NGJ)+
     .           ETD(NTA)*DPEN(NGK,NGJ+2))
            GO TO 140
C
C           TREAT S(2)D(1)F(N-3) CONFIG
  130       IF (NTYP(JA).NE.5) GO TO 140
            IF (NGJ.LT.NG(K)) THEN
               KDD = KC + 2
               KF = KC + 3
            ELSE
               KDD = NG(K) + ((NGJ+2)* (NGJ+1))/2
               KF = NG(K) + (NGJ+3)* (NGJ+2)/2
            END IF
            IF (KAN(JA).LT.89) THEN
               XX = F(KD)* (ETS(NTA)*DPEN(NGK,NGJ)+
     .              ETD(NTA)*DPEN(NGK,NGJ+2)+ETF(NTA)*DPEN(NGK,NGJ+3))
            ELSE
               XX = F(KD)* (ETS(NTA)*DPEN(NGK,NGJ)+
     .              ETP(NTA)*DPEN(NGK,NGJ+1)+ETD(NTA)*DPEN(NGK,NGJ+2)+
     .              ETF(NTA)*DPEN(NGK,NGJ+3))
            END IF
  140       CONTINUE
            Z = Z + XX*CT
            Y = Y + XX*ST*SF
            X = X + XX*ST*CF
  150    CONTINUE
C***********************************************************************
C
C  ADDS IN POINT CHARGES
C
         ZZN = 0.0D0
         YYN = 0.0D0
         XXN = 0.0D0
         DO 230 N1 = 1, NA
            IF (KAN(N1).NE.0) GO TO 230
            IF (N1-JA) 160, 160, 170
  160       KAKA = (JA* (JA-1))/2 + N1
            GO TO 180

  170       KAKA = (N1* (N1-1))/2 + JA
  180       RR = RAD(KA,JA)*BOHRS
            IF (RR.LT.1.D-4) GO TO 230
            DD = ZERO
            DO 210 JJ = 2, 3
               AMU = ALPHA(J,JJ)
               IF (AMU.LT.1.D-4) GO TO 210
               CMU = ALPHA(J,JJ+3)
               DO 200 JJJ = 2, 3
                  BMU = ALPHA(J,JJJ)
                  IF (BMU.LT.1.D-4) GO TO 200
                  DMU = ALPHA(J,JJJ+3)
  190             CALL DIFPEN(FACT,AMU,NPQN,BMU,NPQN,RR,DPENA,IONE)
                  DD = DD + DPENA*CMU*DMU
  200          CONTINUE
  210       CONTINUE
            ITJ = JA + (N1-1)*NA
            DDD = -DD*F(JD)*Z12(N1)
            IF (LJ.NE.0) GO TO 220
            RR = RR!/BOHRS
            DDD = DDD - Z12(JA)*Z12(N1)/ (RR*RR)
  220       CT = FA(ITJ)
            ST = 1.D0 - CT*CT
            IF (ST.LT.THRESH) ST = 0.D0
            ST = SQRT(ST)
            CF = FA(ITJ+NA11)
            SF = FA(ITJ+NA12)
            XXN = XXN + DDD*ST*CF
            YYN = YYN + DDD*ST*SF
            ZZN = ZZN + DDD*CT
C
C  FORCES ON PT CHARGES
            FZ(N1) = FZ(N1) + DDD*ST*CF
            FZ(N1+NA1) = FZ(N1+NA1) + DDD*ST*SF
            FZ(N1+NA2) = FZ(N1+NA2) + DDD*CT
  230    CONTINUE
C***********************************************************************
         IF (LJ.NE.0) GO TO 240
         XXX = XXX - X
         YYY = YYY - Y
         ZZZ = ZZZ - Z
  240    CONTINUE
c        write(6,*)"FZ",FZ(JA),FZ(JA+NA1),FZ(JA+NA2)
c        write(6,*)"XXX,YYY,ZZZ",XXX,YYY,ZZZ
c        write(6,*)"XXN,YYN,ZZN",XXN,YYN,ZZN
         FZ(JA) = FZ(JA) - XXX - XXN
         GRAD(1,JA) = FZ(JA)
         FZ(JA+NA1) = FZ(JA+NA1) - YYY - YYN
         GRAD(2,JA) = FZ(JA+NA1)
         FZ(JA+NA2) = FZ(JA+NA2) - ZZZ - ZZN
         GRAD(3,JA) = FZ(JA+NA2)
  250 CONTINUE
C
      IF (ISWANN.NE.3) GO TO 410
C
C     ADD THE OPEN SHELL CONTRIBUTIONS
C
      DO 400 NIJ = 2, NDN
         DO 260 I = 1, NRR
            F(I) = ZERO
            FA(I) = ZERO
  260    CONTINUE
C        THIS READS THE DENSITY OF THE VARIOUS OPEN SHELL OPERATORS
         NIJM1 = NIJ - 1
         DO 310 NIK = 2, NDN
            NIKM1 = NIK - 1
            IF (NIJ-NIK) 270, 280, 280
  270       IZJZ = NIN(NIKM1) + NIJM1
            GO TO 290

  280       IZJZ = NIN(NIJM1) + NIKM1
  290       CONTINUE
c           CALL FINDP(H,NRR,ISWANN,NIK,ITHREE)
C           F WILL BE THE COULOMB DENSITY, FA THE EXCHANGE DENSITY
            DO 300 I = 1, NRR
               F(I) = F(I) + H(I)*AR(IZJZ)
               FA(I) = FA(I) + H(I)*BR(IZJZ)
  300       CONTINUE
  310    CONTINUE
c        CALL FINDP(H,NRR,ISWANN,NIJ,ITHREE)
         JD = 0
         DO 390 J = 1, NB
            JA = NU(J)
            JD = JD + J
            KD = 0
            LJ = NW(J)
            XXX = ZERO
            YYY = ZERO
            ZZZ = ZERO
            NGJ = NG(J)
            DO 380 K = 1, NB
               KA = NU(K)
               KD = KD + K
               LK = NW(K)
               NGK = NG(K)
               IF (KA-JA) 320, 380, 330
  320          IJ = JD - J + K
               KAJA = NIN(JA) + KA
               KC = (NGJ* (NGJ-1))/2 + NGK
               GO TO 340

  330          IJ = KD - K + J
               KAJA = NIN(KA) + JA
               KC = (NGK* (NGK-1))/2 + NGJ
  340          CONTINUE
               FAC = -0.5D0* (H(KD)*F(JD)+H(JD)*F(KD)-H(IJ)*FA(IJ))*
     .               DG(KC)
               CT = (CO(3,KA)-CO(3,JA))/RAD(KA,JA)
               ST = ONE - CT*CT
               IF (ST-THRESH) 350, 360, 360
  350          ST = ZERO
               CF = ONE
               SF = ZERO
               GO TO 370

  360          ST = SQRT(ST)
               CF = (CO(1,KA)-CO(1,JA))/ (ST*RAD(KA,JA))
               SF = (CO(2,KA)-CO(2,JA))/ (ST*RAD(KA,JA))
  370          CONTINUE
               XXX = XXX + FAC*ST*CF
               YYY = YYY + FAC*ST*SF
               ZZZ = ZZZ + FAC*CT
  380       CONTINUE
c           write(6,*)"FZ",FZ(JA),FZ(JA+NA1),FZ(JA+NA2)
c           write(6,*)"XXX,YYY,ZZZ",XXX,YYY,ZZZ
            FZ(JA) = FZ(JA) - XXX
            GRAD(1,JA) = FZ(JA)
            FZ(JA+NA1) = FZ(JA+NA1) - YYY
            GRAD(2,JA) = FZ(JA+NA1)
            FZ(JA+NA2) = FZ(JA+NA2) - ZZZ
            GRAD(3,JA) = FZ(JA+NA2)
  390    CONTINUE
  400 CONTINUE
  410 CONTINUE
C
C      IF(NFZ.EQ.0) GO TO 411
C      DO 403 I=1,NFZ
C      J=NFA(I)+(NFD(I)-1)*NA
C  403 FZ(J)=ZERO
C      WRITE(6,2006) (NFA(I),I=1,NFZ)
C 2006 FORMAT(1H0,10X,'FORCE CONSTRAINTS ON ATOMS',8I4//)
C  411 CONTINUE
C
C     FOR SELF CONSISTENT REACTION FIELD THE DIPOLE CONTRIBUTION MUST
C     BE ADDED
C
c     IF (JSW(6).GT.0) THEN
c        FCSCRF = GAMME(1)*BOHRS*2.0D0
c        DO 420 I = 1, NA
c           FZ(I) = FZ(I) + (Z12(I)-D12(I))*DIPGR(1)*FCSCRF
c           FZ(I+NA1) = FZ(I+NA1) + (Z12(I)-D12(I))*DIPGR(2)*FCSCRF
c           FZ(I+NA2) = FZ(I+NA2) + (Z12(I)-D12(I))*DIPGR(3)*FCSCRF
c 420    CONTINUE
c     END IF
C
C     IF(INTERNALS(1).GE.0) THEN
C       DO 777 I=1,NA
C       IF(CQV(I).LE.0.D0) THEN
C          FZ(I)=0.D0
C          FZ(I+NA1)=0.D0
C          FZ(I+NA2)=0.D0
C       ENDIF
C777    CONTINUE
C     ENDIF
C
      DRYONE=.FALSE.
c     IF(JSW(12).EQ.1) CALL DIFCOR(RAD,H,FALPH,FZ,FA,WA,DRYONE)

      IF (IPRINT.GE.0.and.ga_nodeid().eq.0) THEN
         WRITE (6,9000) ISW3

 9000    FORMAT (' GEOMETRY',I5/10X,'ATOM:ENERGY GRADIENT(A.U.)'/20X,
     .          'X',10X,'Y',10X,'Z')

         DO 440 I = 1, NA
            WRITE (6,9010) I, FZ(I), FZ(I+NA1), FZ(I+NA2)
  440    CONTINUE

 9010    FORMAT (10X,I5,3F10.6)

      END IF
      RETURN

      END
