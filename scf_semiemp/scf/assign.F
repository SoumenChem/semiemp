      SUBROUTINE ASSIGN(BB,NC,NROOT,JJ,IT,NASS,MIM,FOP,ISWANN,
     1 ISW3,IPUN,NB,EIG,NEL,MULT,W,NCORE)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      double precision mult,nmult,FW,FOP(*)
      INTEGER      RHFCLO, RHFOPN, UHFANN, UHFNOA
      INTEGER      AOCHAR, AUFBAU, REARRG, SYMASN
      PARAMETER    (RHFCLO=0,UHFNOA=1,UHFANN=2,RHFOPN=3)
      PARAMETER    (SYMASN=-2,REARRG=-1,AUFBAU=0,AOCHAR=1)
      LOGICAL      IBGIV, DRYRUN
      PARAMETER (MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000,MXB=1000000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   ( ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
       Dimension
     x W(*),WW(NB),WWW(NB),TE(NB),TE1(NB),TE2(NB),XT(NB),
     x KK(NB),KK1(NB),KK2(NB),KKF(NB),KKF1(NB),KKF2(NB),
     x NASS(7),KKSYM(8,5)
      DIMENSION    NN(23), NIG(NB),OCCNUM(NB),MIM(*)
      DIMENSION    BB(*),EIG(*)
      DATA         IQ, NAST/0, 0/
      MAXNAS = NB
      NELL = NEL
c     write(6,*)"i am in assign",NELL
c     write(6,*)"MIM",(MIM(I),I=1,10)
c     write(6,*)"FOP",(FOP(I),I=1,10)
C     SEE THE COMMON BLOCK assel.cmn
      ICHK = 0
      IXS = 0
      NN(1) = 1
      NAS = NASS(1)
      FW = 0.0D0
      IF (MIM(JJ).NE.0) FW = FOP(JJ)/DFLOAT(MIM(JJ))
c     write(6,*)"i am here 1",FW
      NELL = NC + MIM(JJ)
c     write(6,*)"in assign",NELL,NROOT
      NE1 = NELL - 14
      IF (NE1.LT.NC+1) NE1 = NC + 1
      NE2 = NE1 + 21
      IF (NE2.GE.NC+NROOT) NE2 = NC + NROOT
C
      IF (NAS) 10, 80, 200
C   ********************************************************************
C   *                                                                  *
C   *                                                                  *
C   * NAS.LT.0                                                         *
C   * REARRANGE VECTORS ACCORDING TO NEXT CARD                         *
C   * (IF UHF, NEXT TWO CARDS)                                         *
C   * TYPICAL INPUT FOR NAS=-1                                         *
C   *   4 7 5 6 8                                                      *
C   * RESULTS IN ORBITAL SEVEN MOVING TO FIFTH POSITION                *
C   * AFTER REARRANGMENT, ORBITALS ARE FILLED BY AUFBAU                *
C   * THE FIRST ORBITAL GIVEN (I.E. 4 ABOVE) IS THE REFERENCE, THEN A  *
C   * MAX OF MAXNAS LABELS CAN BE GIVEN.                               *
C   * THIS PROCEDURE IS USED ONLY ON THE FIRST CYCLE!!!                *
C   *                                                                  *
   10 IF (NAS.EQ.REARRG) THEN
c     write(6,*)"i am here"
C
C REARRANGE STARTING VECTORS
C
         IF (IT.NE.0) GO TO 80
         IF (ISWANN.EQ.RHFOPN .AND. JJ.GT.1) GO TO 80
         IF (ISWANN.EQ.RHFOPN) GO TO 30
         IF (JJ.NE.2) GO TO 30
C
C MOVE ALPHA UHF ORDER INTO KK
C
         DO 20 I = 1, MAXNAS
            KK(I) = KK1(I)
   20    CONTINUE
C
   30    CONTINUE
         NBOT = KK(1)
         NN(1) = KK(1)
         DO 40 I = 2, MAXNAS
            NN(I) = NBOT - 1 + I
            IF (KK(I).EQ.0) GO TO 50
   40    CONTINUE
   50    I = I - 1
         NTOP = KK(I)
         DO 70 I = NBOT, NTOP
            J = I - NBOT + 1
            IF (KK(J).EQ.NN(J)) GO TO 70
            IS = NB* (I-1)
            JS = NB* (KK(J)-1)
            KJ = KK(J) - NBOT + 1
            NN(KJ) = NN(J)
            DO 60 K = 1, NB
               X = BB(IS+K)
               BB(IS+K) = BB(JS+K)
               BB(JS+K) = X
   60       CONTINUE
   70    CONTINUE
         NDIF = NTOP - NBOT + 1
c        IF (ISW3.EQ.0 .OR. IXS.EQ.1) THEN
c           WRITE (6,9060) (KK(I),I=1,NDIF)
c        END IF
         NN(1) = 1
         GO TO 80

      ELSE
C
C IF NAS LESS THAN ZERO AND NOT EQUAL TO REARRG, THEN
C NAS EQUALS SYMASN.  SO,
C
C ASSIGN ORBITAL OCCUPANCIES BY SYMMETRY OF M.O.'S
C
c        CALL SYMASS(BB,W,NC,NROOT,JJ,IT,IXS)
         GO TO 610
C
      END IF
C **********************************************************************
C *                                                                    *
C *                  *****  *    *  *****   ****   *****  *   *        *
C * NAS.EQ.0         *   *  *    *  *       *   *  *   *  *   *        *
C *                  *****  *    *  *****   *****  *****  *   *        *
C *                  *   *  *    *  *       *   *  *   *  *   *        *
C *                  *   *  ******  *       ****   *   *  *****        *
   80 CONTINUE
      DO 90 I = 1, NB
         W(I) = 0.0D0
   90 CONTINUE
C
C     FILL M.O.'S IN ORDER
C
      DO 100 I = NC + 1, NC + MIM(JJ)
c      write(6,*)"FW",FW
         W(I) = FW
  100 CONTINUE
C
C
  110 DO 120 I = 1, MAXNAS
         WW(I) = 0.0D0
  120 CONTINUE
c     write(6,*)"i am here 3",NE1,NE2
C
C     CHECK FOR DEGENERACIES ,WHICH GIVE FRACTIONAL OCCUPANCIES
C
      IF (ISWANN.EQ.1) GO TO 480
      IF (ISWANN.EQ.2) GO TO 480
      IF (IPUN.EQ.3) GO TO 480
      IF (IPUN.EQ.5) GO TO 480
      IF (IPUN.EQ.8) GO TO 480
      IF (IPUN.EQ.10) GO TO 480
      IF (NE1.EQ.NE2) GO TO 480
C
C *** BETA OR CLOSED SHELL
C
  130 CONTINUE
      DO 140 I = 1, 23
         NN(I) = 1
  140 CONTINUE
      L = 1
c     write(6,*)"i am here 4",(EIG(I),I=1,10)
      DO 170 I = NE1, NE2
         WW(L) = W(I) + WW(L)
         IF (I+1.GT.NB) GO TO 160
         IF (ABS(EIG(I)-EIG(I+1))-1.0D-4) 150, 150, 160
  150    CONTINUE
         IF (I.EQ.NE2) GO TO 160
         NN(L) = NN(L) + 1
C        IF(I.EQ.NE2) WW(L)=WW(L)/DFLOAT(NN(L))
         GO TO 170

  160    WW(L) = WW(L)/REAL(NN(L))
         L = L + 1
  170 CONTINUE
C     W(NE1)=WW(1)
C     NE3=NE1+1
C     L=1
C     DO 25 I=NE3,NE2
C     IF(ABS(EIG(I)-EIG(I-1)).GT.1.0D-3) L=L+1
C     IF(W(I).NE.WW(L).AND.ISWANN.EQ.3) GO TO 4000
C     25 W(I)=WW(L)
C     JDHD
C     CHECK FOR DEGENERACIES
c     write(6,*)"NN",NN
      L = 1
      LTEST = 0
      DO 180 I = NE1, NE2
c     write(6,*)"WW(L),I",WW(L),I
         W(I) = WW(L)
         IF (ABS(EIG(I)-EIG(I+1)).GT.1.D-4) L = L + 1
         IF (W(I+1).NE.WW(L) .AND. ISWANN.EQ.3) LTEST = 1
  180 CONTINUE

      IF (LTEST.EQ.1) WRITE (6,9040) JJ
C
      IF (ISWANN.EQ.3) GO TO 480
C
C    CHECK MULTIPLICITY
C
      S = 0.0D0
      DO 190 I = 1, NB
c        write(6,*)"W(I)",W(I)
         IF (W(I).GT.1.0D0) S = S + (W(I)-1.0D0)
  190 CONTINUE
c     write(6,*)"S,NEL",S,NEL
      NMULT = (DFLOAT(NEL)-2.0D0*S) + 1.0D0
c     write(6,*)"here ??",(NMULT-MULT)
              IF ((NMULT-MULT).EQ.0) GO TO 480
              if (ga_nodeid().eq.0) WRITE (6,9110)
              ICHK = 1
              GO TO 480
C **********************************************************************
C *                                                                    *
C *                                                                    *
C *   SPECIAL M.O.'S ARE DESIGNATED WITH SPECIFIC OCCUPANCY.           *
C *   TYPICAL INPUT FOR NAS=3, NASS(2)=3 , NASS(3)=2                   *
C *     1.0        22  0.7         1 -1.0         5 -1.0               *
C *     1.0        23  0.7                                             *
C *     1.0        24  0.6                                             *
C *                                                                    *
C *                                                                    *
  200  IF (IQ.EQ.1) GO TO 220
C     THE DIMENSIONS ARE SET IN assel.cmn
              IF (NAS.GT.100) GO TO 590
              IQ = 1
              WRITE (6,9000)

 9000 FORMAT ('  ASSIGN TABLE........................................')

              WRITE (6,9010) NASS
              DO 210 J = 1, NAS, 18
                 NASTOP = J + 18
                 WRITE (6,9020) (XT(I),I=J,NASTOP)
                 WRITE (6,9010) (KK(I),I=J,NASTOP)
                 WRITE (6,9010) (KKF(I),I=J,NASTOP)
                 WRITE (6,9020) (TE(I),I=J,NASTOP)
                 WRITE (6,9010) (KK1(I),I=J,NASTOP)
                 WRITE (6,9010) (KKF1(I),I=J,NASTOP)
                 WRITE (6,9020) (TE1(I),I=J,NASTOP)
                 WRITE (6,9010) (KK2(I),I=J,NASTOP)
                 WRITE (6,9010) (KKF2(I),I=J,NASTOP)
                 WRITE (6,9020) (TE2(I),I=J,NASTOP)
                 WRITE (6,9130)
  210 CONTINUE

 9010 FORMAT (4X,19I4)
 9020 FORMAT (4X,19F4.1)

  220 CONTINUE
C
C     RESET OCCUPATION NUMBERS AND OCCUPIED FLAGS (NIG(I))
C
              DO 250 I = 1, NB
C        IF(ISWANN.NE.3) GO TO 43
C        IF(JJ.EQ.1)     GO TO 43
C        IF(IPUN.EQ.5)   GO TO 44
C        IF(IPUN.EQ.6)   GO TO 44
C        IF(IPUN.EQ.8)   GO TO 44
C        IF(IPUN.EQ.10)  GO TO 44
  230    NIG(I) = 0
  240    W(I) = 0.0D0
  250 CONTINUE
C
              FNNN = 0.0D0
              IF (JJ.GT.1) GO TO 270
C
C      BETA OR CLOSED SHELL
C
  260 NAST = NASS(1) - NASS(3)
      IF (ISWANN.EQ.3) NAST = NASS(3)
        NBEG = 1
        NTOP = NB
        IF (NAST.EQ.0) GO TO 430
        GO TO 290
C
C      ALPHA OR OPEN SHELL
C
  270 NBEG = NAST + 1
      NAST = NASS(1)
      IF (ISWANN.NE.3) GO TO 280
        NAST = NBEG - 1 + NASS(JJ+2)
C     NTOP=NTOP-MIM(JJ-1)
C     IF(IT.NE.0) GO TO 57
C     IF(IPUN.EQ.5) GO TO 557
C     IF(IPUN.EQ.6) GO TO 557
C     IF(IPUN.EQ.8) GO TO 557
C     IF(IPUN.EQ.10) GO TO 557
C     GO TO 57
C     557 CONTINUE
C     NTOP=NB
  280 IF (NBEG.GT.NAST) GO TO 430
  290 CONTINUE
C
C
C
        NRRR = NROOT
       IF (ISWANN.EQ.3 .AND. JJ.NE.1) NRRR = NROOT + NC
              DO 420 J = NBEG, NAST
                 TEST = TE(J)
                 TEST1 = TE1(J)
                 TEST2 = TE2(J)
                 K = KK(J)
                 K1 = KK1(J)
                 K2 = KK2(J)
                 KF = KKF(J)
                 KF1 = KKF1(J)
                 KF2 = KKF2(J)
                 X = XT(J)
C        THE NO. OF ELECTRONS IS EITHER 0 OR THE AVERAGE NO.
C        OF ELECTRONS FOR A GIVEN OPERATOR - UNLESS EHT.
                 IF (ABS(X-FW).LT.0.000002) GO TO 300
                 IF (X.EQ.0.0D0) GO TO 300
c                IF (ISW1.EQ.0) GO TO 300
c        IF(INDEX(CONTRLO,'SCF_IT').EQ.0.AND.IT.NE.0) GO TO 300
                 WRITE (6,9050) X, JJ, FW
                 X = FW
  300    CONTINUE
C
C
         DO 410 I = 1, MAXNAS
C
C           CHECK ORBITALS ABOVE NELL FIRST
C
            JB = NELL + I
            IF (JB.GT.NRRR) GO TO 390
  310       IF (NIG(JB).GE.1) GO TO 390
C
C           FIRST TEST ===========================
C
C           LOOK FOR ORBITAL COEFICIENT
C
            IB = (JB-1)*NB
            IF (KF.GT.0) GO TO 320
            XX = BB(IB+K)**2 - TEST
            IF (TEST.LT.-0.5D0) XX = 1.001D0 - XX
            IF (XX.LT.0.0D0) GO TO 390
            GO TO 330
C
C           LOOK FOR ORBITAL PAIR
C
  320       T1 = BB(IB+K)*BB(IB+KF)/TEST
            IF (T1.LT.1.0D0) GO TO 390
C
C           SECOND TEST==========================
C
C           LOOK FOR ORBITAL COEFICIENT
C
C
  330       IF (K1.EQ.0) GO TO 370
            IF (KF1.GT.0) GO TO 340
            XX = BB(IB+K1)**2 - TEST1
            IF (TEST1.LT.-0.5D0) XX = 1.001D0 - XX
            IF (XX.LT.0.0D0) GO TO 390
            GO TO 350
C
C           LOOK FOR ORBITAL PAIR
C
  340       T1 = BB(IB+K1)*BB(IB+KF1)/TEST1
            IF (T1.LT.1.0D0) GO TO 390
C
C           THIRD TEST ==========================
C
C           LOOK FOR ORBITAL COEFFICIENTS
C
  350       IF (K2.EQ.0) GO TO 370
            IF (KF2.GT.0) GO TO 360
            XX = BB(IB+K2)**2 - TEST2
            IF (TEST2.LT.-0.5D0) XX = 1.001D0 - XX
            IF (XX.LT.0.0D0) GO TO 390
            GO TO 370
C
C           LOOK FOR ORBITAL PAIR
C
  360       T1 = BB(IB+K2)*BB(IB+KF2)/TEST2
            IF (T1.LT.1.0D0) GO TO 390
  370       CONTINUE
C
C           COUNT UP ASSIGNED ELECTRONS
C
  380       W(JB) = X
            FNNN = FNNN + X
            NIG(JB) = 2
            GO TO 420

  390       IF (JB.LE.NELL) GO TO 400
C
C           NOW LOOK AT ORBITALS BELOW NELL
C
            JB = NELL - I + 1
            IF (JB.LE.NC) GO TO 400
            GO TO 310

  400       CONTINUE
  410    CONTINUE
  420 CONTINUE
C
C
C
C
  430 CONTINUE
      RNE3 = 1.0D-8
      RNE4 = FOP(JJ) - FNNN
      IF (ABS(RNE4).LT.1.0D-2) GO TO 460
C
C     ASSIGN REMAINING ELECTRONS BY AUFBAU
C
      DO 450 I = NC + 1, NB
         IF (NIG(I).NE.0) GO TO 450
C
         FW = 0.0D0
         IF (MIM(JJ).NE.0) FW = FOP(JJ)/DFLOAT(MIM(JJ))
         W(I) = FW
         NIG(I) = 2
         RNE3 = RNE3 + FW
         IF (RNE3-RNE4) 450, 460, 440
  440    IF (ABS(RNE3-RNE4).LT.1.0D-2) GO TO 460
C
C        ASSIGN LAST ELECTRON
C
         W(I) = RNE3 - RNE4
         NIG(I) = 2
         GO TO 460

  450 CONTINUE
C
C     REMOVE FLAGS IF NO ELECTRON WAS ASSIGNED
C
  460 CONTINUE
      DO 470 I = NC + 1, NB
         IF (NIG(I).NE.2) GO TO 470
         IF (W(I).LT.1.0D-2) NIG(I) = 1
         NIG(I) = NIG(I) - 1
  470 CONTINUE
C
C     NOW WRITE OUT ASSIGNED ELECTRONS
C
  480 CONTINUE
c     write(6,*)"i am here 5",ICHK,ISWANN,NAS,JJ
      IF (ICHK.EQ.1) GO TO 490
      IF (ISWANN.EQ.0 .AND. NAS.EQ.0) GO TO 600
      IF (ISWANN.EQ.3) GO TO 530
      IF (JJ.GT.1) GO TO 510
C
C     CLOSED SHELL RHF AND BETA UHF
C
  490 CONTINUE
      IF ((ISW3.EQ.0 .OR. IXS.EQ.1).and.ga_nodeid().eq.0) THEN
c       IF(IT.EQ.0.OR.INDEX(CONTRLO,'SCF_IT').NE.0) THEN
        IF (IT.EQ.0) THEN
         WRITE (6,9070) NCORE
c        WRITE (6,9070) NE1 - 1
         DO 500 J = NE1, NE2, 12
            NE22 = J + 11
            IF (NE22.GT.NE2) NE22 = NE2
            WRITE (6,9080) (I,I=J,NE22)
            IF (ISWANN.EQ.0)  WRITE (6,9120) (W(I),I=J,NE22)
            IF (ISWANN.NE.0)  WRITE (6,9090) (W(I),I=J,NE22)
            WRITE (6,9130)
  500    CONTINUE
         ENDIF
        ENDIF
c     END IF
      GO TO 610
C
C     ALPHA UHF
C
  510 CONTINUE
      IF ((ISW3.EQ.0 .OR. IXS.EQ.1.AND.ga_nodeid().eq.0)) THEN
c       IF(IT.EQ.0.OR.INDEX(CONTRLO,'SCF_IT').NE.0) THEN
        IF (IT.EQ.0) THEN
c        WRITE (6,9070) NCORE
c        WRITE (6,9070) NE1 - 1
         DO 520 J = NE1, NE2, 12
            NE22 = J + 11
            IF (NE22.GT.NE2) NE22 = NE2
            IF (ga_nodeid().eq.0) WRITE (6,9080) (I,I=J,NE22)
            IF (ga_nodeid().eq.0) WRITE (6,9100) (W(I),I=J,NE22)
            WRITE (6,9130)
  520    CONTINUE
         ENDIF
        ENDIF
c     END IF
      GO TO 610
C
C     OPEN SHELL RHF---CLOSED SHELL OPERATOR
C
  530 CONTINUE
      IF (JJ.GT.1) GO TO 560
  540 CONTINUE
      IF (ISW3.EQ.0 .OR. IXS.EQ.1) THEN
c       IF(IT.EQ.0.OR.INDEX(CONTRLO,'SCF_IT').NE.0) THEN
c        WRITE (6,9070) NE1 - 1
c        DO 550 J = NE1, NE2, 12
c           NE22 = J + 11
c           IF (NE22.GT.NE2) NE22 = NE2
c           WRITE (6,9080) (I,I=J,NE22)
c           WRITE (6,9120) (W(I),I=J,NE22)
c           WRITE (6,9130)
c 550    CONTINUE
c       ENDIF
      END IF
      GO TO 610
C
C     OPEN SHELL RHF---OPEN SHELL OPERATOR
C
  560 CONTINUE
      IF (ISW3.EQ.0 .OR. IXS.EQ.1) THEN
c       IF(IT.EQ.0.OR.INDEX(CONTRLO,'SCF_IT').NE.0) THEN
c         JJ1 = JJ - 1
c         WRITE (6,9030) JJ1

 9030     FORMAT (//,'  OPEN SHELL OCCUPANCY',I5)

          DO 570 J = NE1, NE2, 12
            NE22 = J + 11
            IF (NE22.GT.NE2) NE22 = NE2
            WRITE (6,9080) (I,I=J,NE22)
            WRITE (6,9120) (W(I),I=J,NE22)
            WRITE (6,9130)
  570     CONTINUE
c       END IF
      END IF
cmgc
      DO 111 MGC=NE1,NE2
 111  OCCNUM(MGC)=W(MGC)
cmgc
      GO TO 610
C
C     ERROR MESSAGES*******
C

  590 WRITE (6,9140) MAXNAS,NAS
      CALL ZEXIT(3)

 9040 FORMAT(' *** WARNING - SPLIT DEGENERACY IN OPEN SHELL',I5,' ***')
 9050 FORMAT (' *** INVALID ASSIGNMENT OF ',F5.2,' ELECTRONS FOR ',
     .       ' OPERATOR ',I2,3X,'...',F5.2,' ELECTRONS ASSIGNED')
 9060 FORMAT ('  STARTING VECTORS HAVE BEEN REARRANGED--NEW ORDER IS',
     .       //,20I4)
 9070 FORMAT (//' THE FIRST',I5,' M.O.''S ARE FULLY OCCUPIED,THEN:'/)
 9080 FORMAT (' M.O.:',12X,12 (I5))
 9090 FORMAT (' BETA  OCCUP:',5X,12F5.1)
 9100 FORMAT (' ALPHA OCCUP:',5X,12F5.1)
 9110 FORMAT (/,
     .     ' *** THE HOMO AND LUMO ARE DEGENERATE, SO THIS CYCLE PRO'
     .     ,'CEEDS'/' AS AN AVERAGED EXCITED SINGLET-TRIPLET STATE ***'/
     .     )
 9120 FORMAT (' ','OCCUPANCY:  ',5X,12F5.1)
 9130 FORMAT ()
 9140 FORMAT (' ***** MAX. OF',I4,' ELECTRONS CAN BE ASSIGNED, LE.',I5)

  600 IF (IT.EQ.0 .AND. ICHK.EQ.0) THEN
         IF ((ISW3.EQ.0 .OR. IXS.EQ.1).and.ga_nodeid().eq.0) THEN
            WRITE (6,9150)
         END IF
      END IF

 9150 FORMAT (/,' RESTRICTED HARTREE-FOCK CALC.- ',/,
     .  ' ALL OCC. ORBITALS CONTAIN 2 ELECTRONS AND ARE FILLED IN ORDER'
     .       )
 9200 FORMAT (/,' UNRESTRICTED HARTREE-FOCK CALC.- ',/
     .       )
  610 CONTINUE
      RETURN

      END

