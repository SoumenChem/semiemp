      SUBROUTINE NUCREP(NA,RAD,KAN,Z12,REP,PO,AALP,ISW1,NB,NG
     1  ,GAMMA,ISPEC)
C  NUCLEAR REPULSION FOR AM1, MNDO, PM3 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION   Z12(*), RAD(NA,NA), KAN(*),PO(9,105),AALP(105)
      DIMENSION NB(*), NG(*), GAMMA(*)
      PARAMETER (BOHRS    = 5.29177249D-1, ZERO = 0.0d0, IZERO = 0)

c     WRITE(6,*)"i am in nucrep"
      IQMMM = IZERO
      IPRINT = IZERO
      REP = ZERO
c     ADDLIN = .FALSE.
       IF(ISW1.GT.50.AND.ISW1.LT.54) THEN
C      MNDO, AM1, PM3
         DO 5 I = 1, NA

C     Linear scaling - calculate global nuclear repulsion
C     Add it if one of the atoms belongs to the fragment 
c         IF (IFRAG(I).EQ.NOWF.AND.ITGL.EQ.IONE) THEN
c          ADDLIN = .TRUE.
c         ELSE
c          ADDLIN = .FALSE.
c         END IF
          II = I + 1
          IF (II .GT. NA) GO TO  6
          DO  3 J=II,NA
           IJ = I + (J * (J - 1)) / 2
           IF ((KAN(I) .EQ .0) .AND. (KAN(J) .EQ. 0)) GO TO  3
           IF ((Z12(I) .EQ. 0) .OR. (Z12(J) .EQ. 0)) GO TO  3
           IF (ABS(RAD(I,J)) .LT. 1.0D-2) THEN
             WRITE (6,9000) I,J,RAD(I,J)
             PRINT*,KAN(I),KAN(J)
             CALL ZEXIT(3)
           ENDIF
           NI = KAN(I)
           NJ = KAN(J)
           CALL CCREP(KAN(I),KAN(J),Z12(I),Z12(J),RAD(I,J),ZZ,PO(9,NI),
     1          PO(9,NJ),AALP(NI),AALP(NJ),ISW1)
c         write(6,*)"REP,ZZ",REP,ZZ
           REP = REP + ZZ
C     Linear scaling
c          IF (ADDLIN) REPGL = REPGL + ZZ
    3   CONTINUE
    5   CONTINUE
    6   CONTINUE
       IF (IQMMM.EQ.IZERO.AND.IPRINT.NE.-2) WRITE (6,9010) REP
       RETURN
       ENDIF
c       DO I = 1, NA
c       write(6,*)"NB",NB(I)
c       END DO
       IF (ISPEC.EQ.1) THEN
C     SPECTROSCOPIC CALCULATION USES GAMMA
c      write(6,*)"ISPEC",ISPEC
       
       ID = 1
       L = 0
       M = 0
       DO 50, I=2,NA
        L = L+NB(I-1)
c      write(6,*)"NFIRST",L+1
        IORB = L+1
        NGD  = (NG(IORB) * (NG(IORB) - 1)) / 2
        DO 30 J = 1, I-1
        IF (J.EQ.1) THEN 
         M = 1
        ELSE
         M = M + NB(J-1)
        END IF
         JORB = M
C     Linear scaling - see above
c        IF (IFRAG(J).EQ.NOWF.AND.ITGL.EQ.IONE) THEN
c         ADDLIN = .TRUE.
c        ELSE
c         ADDLIN = .FALSE.
c        END IF

         IJ   = ID + J
         IF ((KAN(I) .EQ. 0) .AND. (KAN(J) .EQ. 0)) GO TO 30
         IF (ABS(RAD(I,J)) .LT. 1.0D-2) THEN
          WRITE (6,9000) I,J,RAD(I,J)
          CALL ZEXIT(3)
         END IF
         NGA = NGD + NG(JORB)
         IF(KAN(I) .EQ. 0 .OR. KAN(J) .EQ. 0) THEN
          REP = REP + BOHRS* Z12(I) * Z12(J) / RAD(I,J)
C     Linear scaling
c         IF (ADDLIN) REPGL = REPGL + BOHRS*Z12(I)*Z12(J)/RAD(IJ)
         ELSE
c         write(6,*)"I,J,GAMMA",Z12(I),Z12(J),GAMMA(NGA)
          REP = REP + Z12(I) * Z12(J) * GAMMA(NGA)
c         write(6,*)"REP 2",REP
C     Linear scaling
c         IF (ADDLIN) write(6,*)"i am in addlin"
c        REPGL = REPGL + Z12(I)*Z12(J)*GAMMA(NGA)
         ENDIF
   30   CONTINUE
        ID = ID + I
   50  CONTINUE
       END IF
c      IF (IQMMM.EQ.IZERO.AND.IPRINT.NE.-2) WRITE (6,9010) REP
 9000          FORMAT (' ERROR IN GEOMETRY, ATOM ',I3,' AND',I3,
     .         ' AT',F10.6,'A')
 9010 FORMAT (//' CORE REPULSION ENERGY =',F14.6//)
       RETURN
       END
C--------------------------------------------------------------------
C
C--------------------------------------------------------------------
      SUBROUTINE CCREP(NI,NJ,ZI,ZJ,RIJ,ENUCLR,POI,POJ,AALPI,AALPJ,ISW1)
C     From voityuk@theochem.tu-muenchen.de Thu Jul 31 17:11:13 1997
C     CORE-CORE REPULSIONS FOR MNDO,AM1,PM3   
C     NI AND NJ ARE ATOMIC NUMBERS.  ZI AND ZJ ARE CORE CHARGES
C      R IS IN ANGSTROMS
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/am1.cmn'
     
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (ZERO =0.0d0, ONE = 1.0D0, BOHRS    = 5.29177249D-1) 
      LOGICAL MOD
c     write(6,*)"i am in ccrep"
      RAU    = RIJ
      R = RIJ*BOHRS
c     write(6,*)"POJ,POK",POI,POJ
c     write(6,*)"AALPI,AALPJ",AALPI,AALPJ
c     GAB    = ONE/SQRT(RAU**2+(PO(9,NI)+PO(9,NJ))**2)
      GAB    = ONE/SQRT(RAU**2+(POI+POJ)**2)
      ENUC   = ZI*ZJ *GAB
      ADD    = ZERO
c
      ENI    = EXP(-AALPI*R)
      IF(NI.EQ.NJ) THEN
         ENJ = ENI
      ELSE
         ENJ = EXP(-AALPJ*R)
      ENDIF
      SCALE  = ONE+ENI+ENJ
C
      MOD=(NI.EQ.7 .OR. NI .EQ.8) .AND. NJ.EQ.1 .OR.
     .    (NJ.EQ.7 .OR. NJ .EQ.8) .AND. NI.EQ.1
      IF(MOD) THEN
         IF(NJ.EQ.1) SCALE=SCALE+(R-ONE)*ENI
         IF(NI.EQ.1) SCALE=SCALE+(R-ONE)*ENJ
      ENDIF
C     SPECIAL CORE-CORE REPULSION TERMS FOR AM1 AND PM3.
      IF(ISW1.EQ.52 .OR. ISW1.EQ.53) THEN
        CALL GFUNC (ISW1,NI,NJ,ZI,ZJ,R,ADD)
C       WRITE(6,*) ' CCREP ',ISW1,NI,NJ,ZI,ZJ,ADD
      ENDIF
C
      ENUCLR = ENUC*SCALE+ADD
C
      RETURN
      END
C------------------------------------------------------------------------------
C
C-----------------------------------------------------------------------------
      SUBROUTINE GFUNC(ISW1,NI,NJ,ZI,ZJ,RAU,ADD)
C     *
C     ADDITIONAL TERMS TO CORE-CORE REPULSION IN AM1 AND PM3
C     *
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/am1.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (ZERO =0.0d0, ONE = 1.0D0,  EVOLT    = 2.72113961D1,
     1 BOHRS    = 5.29177249D-1)
      DIMENSION BORON1(3,4),BORON2(3,4),BORON3(3,4)
      DIMENSION IMP(105), GK(105,4), GL(105,4), GM(105,4)

      DATA BORON1/  0.182613D0,  0.118587D0, -0.073280D0,
     1              0.412253D0, -0.149917D0,  0.000000D0,
     2              0.261751D0,  0.050275D0,  0.000000D0,
     3              0.359244D0,  0.074729D0,  0.000000D0/
      DATA BORON2/  6.0D0,  6.0D0,  5.0D0,
     1             10.0D0,  6.0D0,  0.0D0,
     2              8.0D0,  5.0D0,  0.0D0,
     3              9.0D0,  9.0D0,  0.0D0/
      DATA BORON3/  0.727592D0,  1.466639D0,  1.570975D0,
     1              0.832586D0,  1.186220D0,  0.000000D0,
     2              1.063995D0,  1.936492D0,  0.000000D0,
     3              0.819351D0,  1.574414D0,  0.000000D0/
      DATA CUTOFF/25.0D0/
    
c     R = RAU/BOHRS
      R = RAU
      VIJ    = ZI*ZJ/R
      ADD    = ZERO
      CALL GFPRM(ISW1,GK,GL,GM,IMP)

C      AM1 parameters
c      write(6,*)"i am in gfunc"
      IF(NI.EQ.5.OR.NJ.EQ.5) THEN
C ***  ATOM PAIRS INVOLVING BORON.
         NK  = NI+NJ-5
         NL  = 1
         IF(NK.EQ.1) NL=2
         IF(NK.EQ.6) NL=3
         IF(NK.EQ.9.OR.NK.EQ.17.OR.NK.EQ.35.OR.NK.EQ.53) NL=4
C
         II=3
         DO  I=1,II
         GK(5,I)=BORON1(I,NL)
         GL(5,I)=BORON2(I,NL)
         GM(5,I)=BORON3(I,NL)
         ENDDO
C
      ENDIF
C
c     write(6,*)"IMP(NI)",IMP(NI)
c     write(6,*)"IMP(NJ)",IMP(NJ)
      DO  IG=1,IMP(NI)
      XX     = GL(NI,IG)*(R-GM(NI,IG))**2
c     write(6,*)"XX",XX
      IF(XX.LT.CUTOFF) ADD = ADD+VIJ*GK(NI,IG)*EXP(-XX)
c     write(6,*)"ADD",ADD
      ENDDO
      IF(NI.EQ.NJ) THEN
         ADD = ADD+ADD
      ELSE
         DO  IG=1,IMP(NJ)
         XX  = GL(NJ,IG)*(R-GM(NJ,IG))**2
c     write(6,*)"XX",XX
         IF(XX.LT.CUTOFF) ADD = ADD+VIJ*GK(NJ,IG)*EXP(-XX)
c     write(6,*)"ADD",ADD
        ENDDO
      ENDIF
      ADD=ADD/EVOLT
c     write(6,*)"ADD",ADD
      RETURN
      END











