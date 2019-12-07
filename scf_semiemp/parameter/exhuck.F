      SUBROUTINE EXHUCK(H,F,EIS,EIP,EID,EIF,RAD,NA,WCOR,NRR,LQNO,NB,NIN,
     1 KAN,NW,NU,BETA,BETB,ISW1,NP,ALPHA,BKA,BKB)
C     H ENTERS THIS SUBROUTINE AS THE WEIGHTED OVERLAP, AND
C     RETURNS AS THE STARTING MATRIX.  F RETURNS AS THE OVERLAP
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/datain.cmn'
c     INCLUDE '../include/factor.cmn'

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      Parameter ( EVOLT = 2.72113961D1, MXA = 1200, MPE = 105) 
      parameter (IBETA = 0, MXB = 4000,ONE=1.0D0, ISW3 = 0.0d0)
C
c     INCLUDE '../include/switch.cmn'
C****************************************************************
C***      COMMON/SWITCH/
C***     .   ISWANN,III,IPRINT,IPUN,ISAVE,ISPEC,ISW1,ISW2,ISW3,ISW4,
C***     .   JSW(10),IBETA,IDD1,IDD2
C****************************************************************
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     .   ATOM(100),NIN(MXB),
C***     .   NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     .   NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     .   Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
C
      DIMENSION    BETA(MPE,10),BETB(MPE,4)
      DIMENSION    F(*),  H(*), IEX(NB), WCOR(*), LQNO(*), NIN(*)
      DIMENSION    KAN(NB),NW(NB),NU(NB),RAD(NA,NA),ALPHA(NB,6),NP(NB)
C     DIMENSION    BETA(MPE,10), BETB(MPE,4)
      DIMENSION    EIS(MPE), EIP(MPE), EID(MPE), EIF(MPE), FIP(4,NA)
      DIMENSION    BKA(NB), BKB(NB)
C
c     write(6,*)"i am in exhuck",(F(IJ),IJ=1,nb*(nb+1)/2)
      HTEMP = 0.0d0
      X = 0.0d0 
      do i = 1, NA
       FIP(1,I) = 0.0d0
       FIP(2,I) = 0.0d0
       FIP(3,I) = 0.0d0
       FIP(4,I) = 0.0d0
      end do
      NEXTR = 0
      HUKFAC = 1.00D0
      JBETA = 0
      DO I = 1, NRR
      F(I) = H(I)
      END DO

C     WRITE(6,*) ' EXHUCK  NEXTR =',(NEXTR(I),I=1,NB)
      IF (IBETA.EQ.2 .OR. IBETA.EQ.3) JBETA = 1
C
      IF (ISW3.EQ.0) THEN
         IF (IBETA.EQ.0.AND.ga_nodeid().eq.0) WRITE (6,9000)

 9000    FORMAT (/4X,'>>>>> SUM FORM FOR  BETA <<<<< '//)

         IF (IBETA.EQ.1) WRITE (6,9010)

 9010    FORMAT (/4X,'>>>>> PRODUCT FORM FOR BETA <<<<< '//)

         IF (IBETA.EQ.2 .OR. IBETA.EQ.3) WRITE (6,9020)

 9020    FORMAT (/4X,'>>>>> BETA = A + B/R     <<<<< '//)

      END IF
C
C     write(6,*) "i am in 1"
      DO 10 I = 1, NB
         IEX(I) = 0
   10 CONTINUE

      IF (IBETA.EQ.5) WRITE(6,9025) 
 9025 FORMAT(/4X,'>>>> BETA = A + K*(Z1-Z2)**2 <<<<'//)
C     write(6,*) "i am in 2"
      DO 90 I = 1, NB
         II = (I* (I+1))/2
         K = NW(I)
         KU = NU(I)
         KA = KAN(I)
         KL = LQNO(I)
c      write(6,*)"K,KU,KA,KL",K,KU,KA,KL
C     INNER SHELL ORBITALS
c        IF (NEXTR(I).LT.0) THEN
c           CALL PARAMC(HHHH,I,NP(I),KU,KA,K)
c           H(II) = HHHH
c           WCOR(I) = HHHH
c           GO TO 70

c        END IF
C       write(6,*)"i am in 3"
C
C        IF (KA.EQ.0) GO TO 60
c        IF (NEXTR(I).EQ.1) THEN
c        write(6,*)" i am in 4"
C        POLARIZATION FUNCTIONS
C        THE ADDED 20. UNCOUPLES THESE FUNCTIONS FOR THE FIRST SCF CYCLE.
c           WCOR(I) = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN) 
c           H(II) = WCOR(I) + 20.0D0
c           GO TO 70

c        END IF
c        IF (NEXTR(I).EQ.2) THEN
c        write(6,*)"i am in 5"
c           HTEMP = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN)
c        END IF
C     NOTE HOW THIS TREATS 3D A.O.-S.
C
         IF (K.EQ.0) GO TO 20
         IF (K.LE.3) GO TO 30
         GO TO 40

   20    IF (ABS(FIP(1,KU)).LT.0.0001) THEN
C        write(6,*)"i am in 3.4",KA,EIS(KA),EVOLT
            H(II) = EIS(KA)/EVOLT
c          write(6,*)"H(II) in exhuck",H(II)
         ELSE
            WRITE (6,9030) I, KU, FIP(1,KU)

 9030       FORMAT (' #### IP FOR ORBITAL',I5,' ATOM',I5,' RESET TO ',
     .             F10.6,' ####')

            H(II) = FIP(1,KU)/EVOLT
         END IF
C        write(6,*)" i am in 3.5",H(II),IEX(I)
         WCOR(I) = H(II)

         IF (IEX(I).NE.0 .AND. EIS(KA).LT.-99.0D0) THEN
            WCOR(I) = HTEMP
            H(II) = HTEMP + 20.0D0
         END IF
         IF (ISW1.NE.7) GO TO 70
         H(II) = EIP(KA)/EVOLT
         WCOR(I) = H(II)
         GO TO 70

   30    IF (ABS(FIP(2,KU)).LT.0.0001) THEN
            H(II) = EIP(KA)/EVOLT
         ELSE
            WRITE (6,9030) I, KU, FIP(2,KU)
            H(II) = FIP(2,KU)/EVOLT
         END IF
         WCOR(I) = H(II)

         IF (IEX(I).NE.0 .AND. EIP(KA).LT.-99.0D0) THEN
c           HTEMP = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN)
            H(II) = HTEMP + 20.0D0
            WCOR(I) = HTEMP
         END IF
         GO TO 70

   40    CONTINUE
         IF (K.GT.8) GO TO 50
         IF (ABS(FIP(3,KU)).LT.0.0001) THEN
            H(II) = EID(KA)/EVOLT
c           write(6,*)"H(II) in exhuck",H(II)
         ELSE
            WRITE (6,9030) I, KU, FIP(3,KU)
            H(II) = FIP(3,KU)/EVOLT
         END IF
         WCOR(I) = H(II)

         IF (IEX(I).NE.0 .AND. EID(KA).LT.-99.0D0) THEN
c           HTEMP = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN)
            H(II) = HTEMP + 20.0D0
            WCOR(I) = HTEMP
         END IF
         GO TO 70

   50    CONTINUE
         IF (ABS(FIP(4,KU)).LT.0.0001) THEN
            H(II) = EIF(KA)/EVOLT
         ELSE
            WRITE (6,9030) I, KU, FIP(4,KU)
            H(II) = FIP(4,KU)/EVOLT
         END IF
         WCOR(I) = H(II) 

         IF (KA.LT.52) THEN
c           HTEMP = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN)
            H(II) = HTEMP + 20.0D0
            WCOR(I) = HTEMP
         END IF
         GO TO 70

   60    CONTINUE
C     UNCOUPLE RYDBERG OR OTHER EXTRA ORBITALS
c        H(II) = UPRIME(KA,NP(I),ALPHA(I,1),K,Z12(KU),POT,FKIN)
         WCOR(I) = H(II) 
   70    CONTINUE
         LL = I - 1
         IF (LL.EQ.0) GO TO 90
         ID = II - I
         DO 80 J = 1, LL
            IJ = ID + J
            IF (ISW1.EQ.0) THEN
C       EXTENDED HUCKEL PART
               II = NIN(I) + I
               JJ = NIN(J) + J
               FACEX = 1.0D0
c              IF (NEXTR(I).NE.0 .OR. NEXTR(J).NE.0) FACEX = 0.300D0
c              BBA = HUKFAC*FACEX* (H(II)+H(JJ))/2.0D0
c              X = BBA
            ELSEIF (ISW1.GT.50.AND.ISW1.LT.54) THEN
               MA = KAN(J)
               ML = LQNO(J)
c           write(6,*)"need to get HAM1",MA,ML
               CALL HAM1(KA,KL,MA,ML,X,HSP,HSD,HSF,ISW1)
            ELSE IF (IBETA.LE.3) THEN
C           write(6,*)"i am here"
               CALL HOFF(I,J,X,BBA,BBB,RAD,H,IEX,JBETA,IONE,NW,NU,KAN,
     1               BETA,BETB,NIN,NA,NP,ISW1,ALPHA,NB,BKA,BKB)
C           write(6,*)"I,J,X,BBA,BBB,RAD,H",I,J,X,BBA,BBB
            ELSE IF (IBETA.EQ.5) THEN
c           write(6,*)"need to get HOFFP"
c              CALL HOFFP(I,J,X,BBA,BBB,RAD,IEX,JBETA,ITWO,NW,NU,KAN)
            ELSE IF (IBETA.EQ.6) THEN
               X = ONE
            ELSE
               WRITE (6,*) ' ****ERROR IN IBETA, SUBROUTINE EXHUCK ****'
               CALL ZEXIT(3)
            END IF
c          write(6,*)"H(IJ)",H(IJ),"X",X
            H(IJ) = X*H(IJ)
c          write(6,*)"H(IJ)",H(IJ)
   80    CONTINUE
   90 CONTINUE
c     WRITE(6,*) ' IN EXHUCK'
C     CALL PRNT(H,NB,NRR)
c     do i = 1, nb*(nb+1)/2
c     write(6,*)"H(I) in exhuck",H(I)
c     end do
      RETURN

      END
c-------------------------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------------------------
       SUBROUTINE HOFF(IO,JO,X,BBA,BBB,RAD,H,IEX,JBETA,ITIME,NW,NU,KAN
     1  ,BETA,BETB,NIN,NA,NP,ISW1,ALPHA,NB,BKA,BKB)
C
C     OFF DIAGONAL ONE ELECTRON MATRIX "BETA" BETWEEN AO'S
C     IO AND JO, JO . LE. IO.
C
C     ITIME = 1, FIRST TIME:  ITIME = 2, LOOK IT UP.
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      parameter (EVOLT = 2.72113961D1)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105, IBETA = 0)
      Dimension FKAPPA(7,7),  BETA(MPE,10),BETB(MPE,4),BETC(MPE,4)
      DIMENSION    KAN(*),NW(*),NU(*),NP(*)
      DIMENSION    H(*), RAD(NA,NA), IEX(NB), NIN(*)
      DIMENSION    BKA(NB), BKB(NB),ALPHA(NB,6)
      DIMENSION    FNP(8)
      DATA         FNP/0.500D0, 0.50D0, 0.50D0, 0.50D0, 0.500D0,
     .             0.500D0, 0.500D0, 0.500D0/
c     SAVE         BKA, BKB
C
c     write(6,*)"BKA",BKA
c     write(6,*)"BKB",BKAB
      do i =1, 7
       do j = 1, 7
         FKAPPA(i,j)=0.0d0
       end do 
      end do 
c     write(6,*)"i am in hoff",IO,JO
      IF (JO.GT.IO) THEN
         WRITE (6,*) ' ***** ERROR IN SUBROUTINE HOFF ****'
         CALL ZEXIT(3)
      END IF
      HUKFAC = 1.00D0
C
      I = IO
      J = JO
      IFLAG = 1
      XXX = FKAPPA(NP(I),NP(J))*(ALPHA(I,1)-ALPHA(J,1))**2/EVOLT
C     GO TO (10,70), ITIME
C
   10 CONTINUE
      K = NW(I)
      KU = NU(I)
      KA = KAN(I)
      BBB = ZERO
C     INNER SHELL ORBITALS
c     IF (NEXTR(I).LT.0) THEN
c        BKA(I) = BETP(KA,DK(KU))
c        BKB(I) = BETQ(KA,DK(KU))
c        IEX(I) = 1
c     END IF
C
c     IF (KA.EQ.0) GO TO 60
c     IF (NEXTR(I).EQ.1) THEN
C        POLARIZATION FUNCTIONS
c        BKA(I) = BETP(KA,DK(KU))
c        BKB(I) = BETQ(KA,DK(KU))
c        IEX(I) = 1
c        GO TO 70
c     END IF
C     NOTE HOW THIS TREATS 3D A.O.-S.
C
      IF (K.EQ.0) GO TO 20
      IF (K.LE.3) GO TO 30
      GO TO 40

   20 CONTINUE
      IF (ISW1.EQ.7) THEN
         BKA(I) = BETA(KA,2)
         BKB(I) = BETB(KA,2)
      ELSE
         IF (IBETA.EQ.6) THEN
            BKA(I) = EVOLT
            BKB(I) = EVOLT
         ELSE
            BKA(I) = BETA(KA,1)
            BKB(I) = BETB(KA,1)
c        write(6,*)"BKA(I),BKA(J)",BKA(I),BKA(J)
         END IF
      END IF
      GO TO 70

   30 CONTINUE
      IF (IBETA.EQ.6) THEN
         BKA(I) = EVOLT
         BKB(I) = EVOLT
      ELSE
         BKA(I) = BETA(KA,2)
         BKB(I) = BETB(KA,2)
      END IF
      GO TO 70

   40 CONTINUE
      IF (K.GT.8) GO TO 50
      IF (IBETA.EQ.6) THEN
         BKA(I) = EVOLT
         BKB(I) = EVOLT
      ELSE
         BKA(I) = BETA(KA,3)
         BKB(I) = BETB(KA,3)
      END IF
      GO TO 70

   50 CONTINUE
      IF (IBETA.EQ.6) THEN
         BKA(I) = EVOLT
         BKB(I) = EVOLT
      ELSE
         BKA(I) = BETA(KA,4)
         BKB(I) = BETB(KA,4)
      END IF
      GO TO 70

c  60 CONTINUE
C     UNCOUPLE RYDBERG OR OTHER EXTRA ORBITALS
c     BKA(I) = BETP(KA,DK(KU))
c     BKB(I) = BETQ(KA,DK(KU))
c     IEX(I) = 1
   70 CONTINUE

C     FIRST TIME THROUGH IO = 2, JO = 1, and need to set AO 1 PARAMETERS
      IF (I.EQ.2) THEN
         I = 1
         GO TO 10

      END IF
      I = IO
      IF (I.EQ.2) THEN
         K = NW(I)
         KU = NU(I)
         KA = KAN(I)
      END IF
      M = NW(J)
      MU = NU(J)
      MA = KAN(J)
      IF (ISW1.NE.0) GO TO 80
C
C  EXTENDED HUCKEL PART
      II = NIN(I) + I
      JJ = NIN(J) + J
      FACEX = ONE
      IF (IEX(I).EQ.1 .OR. IEX(J).EQ.1) THEN
C       FACEX=0.300D0
         FACEX = (FNP(NP(I))+FNP(NP(J)))/TWO
      END IF
      BBA = HUKFAC*FACEX* (H(II)+H(JJ))/TWO
      X = BBA
      GO TO 100
C
   80 CONTINUE
   90 CONTINUE
      FACEX = 1.0D0
C     THE NEXT LINE CAN BE USED TO UNCOUPLE EXTRA ORBITALS
      IF (IEX(I).EQ.1 .OR. IEX(J).EQ.1) THEN
C       FACEX=0.300D0
         FACEX = (FNP(NP(I))+FNP(NP(J)))/TWO
      END IF
      IF (JBETA.EQ.0) THEN
         IF (IBETA.EQ.0 .OR. IBETA.EQ.2 .OR. IBETA.EQ.6) THEN
            BBA = 0.5D0* (BKA(I)+BKA(J))*FACEX/EVOLT + XXX*FACEX
            X = BBA
c         write(6,*)"X1",X,"BKA(I)",BKA(I),"BKA(J)",BKA(J),FACEX,XXX
         ELSE IF (IBETA.EQ.1 .OR. IBETA.EQ.3) THEN
            BBA = -SQRT(BKA(I)*BKA(J))*FACEX/EVOLT + XXX*FACEX
            X = BBA
         ELSE
            WRITE (6,*) ' ***** THERE IS AN ERROR IN INPUT ON RESONANCE'
     .        , '  INTEGRAL OPTION, IBETA = ', IBETA, '*****'
            CALL ZEXIT(3)
         END IF
      ELSE IF (JBETA.EQ.1) THEN
C       WRITE(6,*) ' HOFF,I,J,KU,MU,BKA(I),BKB(I)',I,J,KU,MU,
C    .              BKA(I),BKB(I)
c        KUMU = NIN(KU) + MU
         R = RAD(KU,MU)
C        write(6,*)"RAD",R
         IF (IBETA.EQ.0 .OR. IBETA.EQ.2) THEN
C        write(6,*)"i am here 1"
            BBA = 0.5D0* (BKA(I)+BKA(J))*FACEX/EVOLT + XXX*FACEX
            IF (KU.EQ.MU) THEN
               BBB = 0.0D0
            ELSE
               BBB = 0.5D0* (BKB(I)+BKB(J))*FACEX/ (EVOLT*R)
            END IF
            X = BBA + BBB
         ELSE IF (IBETA.EQ.1 .OR. IBETA.EQ.3) THEN
            BBA = -SQRT(BKA(I)*BKA(J))*FACEX/EVOLT + XXX*FACEX
            IF (KU.EQ.MU) THEN
               BBB = 0.0D0
            ELSE
               BBB = SQRT(BKB(I)*BKB(J))*FACEX/ (EVOLT*R)
            END IF
            X = BBA + BBB
         ELSE
            WRITE (6,*) ' ***** THERE IS AN ERROR IN INPUT ON RESONANCE'
     .        , '  INTEGRAL OPTION, IBETA = ', IBETA, '*****'
            CALL ZEXIT(3)
         END IF
      ELSE
         WRITE (6,*) ' ***** THERE IS AN ERROR IN INPUT ON RESONANCE',
     .     '  INTEGRAL OPTION, JBETA = ', JBETA, '*****'
         CALL ZEXIT(3)
      END IF
C     IF(KA.EQ.0.OR.MA.EQ.0) X=-SQRT(BKA(I)*BKA(J))/EVOLT
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
c     write(6,*)"X",X
      RETURN

      END
C-------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------
      SUBROUTINE HAM1(KA,KL,MA,ML,HSS,HSP,HSD,HSF,ISW1)
C     THIS SUBROUTINE FOR AM1 and PM3.
C     ISW1=51,52,53 FOR MNDO,AM1,PM3,RESPECTIVELY
C     R IS IN ANGSTROMS
C
C     KA IS ATOMIC NUMBER OF FIRST ATOM
C     MA IS ATOMIC NUMBER OF SECOND ATOM,
C     KL,ML = 0,1 2 = s p d f
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/switch.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (EVOLT    = 2.72113961D1)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)

      DIMENSION BMNDO(86,3),BAM1(86,3),BPM3(86,3)
      DATA BMNDO/
     . -6.989064,  0.000000, -1.350040, -4.017096, -8.252054,
     .-18.985044,-20.495758,-32.688082,-48.290466,  0.000000,
     .  0.000000, -2.586000, -2.670284, -9.086804, -6.791600,
     .-10.761670,-14.262320,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -1.000000,
     .  0.000000, -4.516479,  0.000000,  0.000000, -8.917107,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -3.235147,
     .  0.000000,  0.000000, -7.414451,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -0.404525,
     .  0.000000, -8.042387,  0.000000,  0.000000,  0.000000,
     .  0.000000,
     .  0.000000,  0.000000, -1.350040, -4.017096, -8.252054,
     . -7.934122,-20.495758,-32.688082,-36.508540,  0.000000,
     .  0.000000, -2.842000, -2.670284, -1.075827, -6.791600,
     .-10.108433,-14.262320,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -2.000000,
     .  0.000000, -1.755517,  0.000000,  0.000000, -9.943740,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -4.290416,
     .  0.000000,  0.000000, -6.196781,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -6.206683,
     .  0.000000, -3.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  86*0.D0/
      DATA BAM1 /
     . -6.173787,  0.000000,  0.000000,  0.000000, -9.599114,
     .-15.715783,-20.299110,-29.272773,-69.590277,  0.000000,
     .  0.000000,  0.000000, -3.866822, -3.784852, -6.353764,
     . -3.920566,-24.594670,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -1.997429,
     .  0.000000, -4.356607,  0.000000,  0.000000,-19.399880,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000, -9.861000,  0.000000,  0.000000, -2.921157,
     . -2.540030,  0.000000,  0.000000,  0.000000, -3.235000,
     .  0.000000,  0.000000, -8.443327,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -4.927437,
     .  0.000000,  0.000000,  0.000000,  0.000000, -0.908657,
     .   6*0.0D0,
     .  0.000000,  0.000000,  0.000000,  0.000000, -6.273757,
     . -7.719283,-18.238666,-29.272773,-27.922360,  0.000000,
     .  0.000000,  0.000000, -2.317146, -1.968123, -6.590709,
     . -7.905278,-14.637216,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -4.758119,
     .  0.000000, -0.991091,  0.000000,  0.000000, -8.957195,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000, -4.940000,  0.000000,  0.000000, -1.182660,
     . -1.155416,  0.000000,  0.000000,  0.000000, -2.577890,
     .  0.000000,  0.000000, -6.323405,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -0.482089,
     .  0.000000,  0.000000,  0.000000,  0.000000, -4.909384,
     .  6*0.0D0,   86*0.0D0/
      DATA BPM3/
     . -5.626512,  0.000000, -0.550000, -3.962053,  0.000000,
     .-11.910015,-14.062521,-45.202651,-48.405939,  0.000000,
     .  0.000000,  0.000000, -0.594301, -2.862145,-12.615879,
     . -8.827465,-27.528560,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -0.715578,
     . -4.945618, -5.325002, -8.232165, -6.157822,-31.171342,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000, -8.581944, -2.993319, -2.785802,
     .-14.794217, -2.665146,-14.494234,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -3.101365,
     . -1.084495, -6.126024,  0.000000,  0.000000,  0.000000,
     .  0.000000,
     . -5.626512,  0.000000, -1.500000, -2.780684,  0.000000,
     . -9.802755,-20.043848,-24.752515,-27.744660,  0.000000,
     .  0.000000,  0.000000, -0.956550, -3.933148, -4.160040,
     . -8.091415,-11.593922,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -6.351864,
     . -0.407053, -2.250157, -5.017386, -5.493039, -6.814013,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000, -0.601034, -1.828908, -2.005999,
     . -2.817948, -3.895430, -5.894703,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,
     .  0.000000,  0.000000,  0.000000,  0.000000, -3.464031,
     . -7.946799, -1.395430,  0.000000,  0.000000,  0.000000,
     .  0.000000,  86*0.D0/
C
c     write(6,*)"i am in ham1"
      HSS=ZERO
      HSP=ZERO
      HSD=ZERO
      HSF=ZERO
      I=KL+1
      J=ML+1
      IF(ISW1.EQ.51)THEN
         HSS=(BMNDO(KA,I)+BMNDO(MA,J))*HALF
      ELSEIF(ISW1.EQ.52)THEN
         HSS=(BAM1(KA,I)+BAM1(MA,J))*HALF
      ELSEIF(ISW1.EQ.53)THEN
         HSS=(BPM3(KA,I)+BPM3(MA,J))*HALF
      ENDIF
      HSS=HSS/EVOLT
c     write(6,*)"HERE"
      IF(KL.GE.1.AND.ML.GE.1)HSP=HSS
      IF(KL.GE.2.AND.ML.GE.2)HSD=HSS
C
      RETURN
      END






