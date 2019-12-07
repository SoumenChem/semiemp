      SUBROUTINE SYMORB(C,EGN,NU,NW,IIMO,NORB,JTEN,IFLAG)
C-----------------------------------------------------------------------
C     THIS SUB. ANALYSES THE SYM. OF VECTORS BY PROJECTION.
C     (NOTE THAT THIS SUB CAN ALSO BE USED TO SYM. BLOCK ANY
C     MATRIX)
C     JTEN = 3 FOR CORE ORBITALS
C-----------------------------------------------------------------------
      INCLUDE '../include/param.cmn'
      INCLUDE '../include/const.cmn'
      INCLUDE '../include/mark1.cmn'
      INCLUDE '../include/switch.cmn'
      INCLUDE '../include/sym.cmn'
      DIMENSION    G2(2,2), G4(4,4), G8(8,8), JOPER(8,8), SGN2(16,7)
      DIMENSION    NU(*), NW(*), EGN(*), C(*), IWO(MXA), IIMO(MXB)
      DIMENSION    NFIRST(MXA)
      DATA         G2/1.0D0, 1.0D0, 1.0D0, -1.0D0/
      DATA         G4/6*1.0D0, -1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0,
     .             -1.0D0, 1.0D0, -1.0D0, -1.0D0, 1.0D0/
      DATA         G8/12*1.0D0, 4*-1.0D0, 1.0D0, 1.0D0, -1.0D0, -1.0D0,
     .             1.0D0, 1.0D0, -1.0D0, -1.0D0, 1.0D0, 1.0D0, 4*-1.0D0,
     .             1.0D0, 1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0,
     .             -1.0D0, -1.0D0, 1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0,
     .             -1.0D0, 1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0, -1.0D0,
     .             1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0,
     .             -1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0/
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
      DATA         SGN2/1.0D0, -1.0D0, -1.0D0, 1.0D0, 3*1.0D0, 2*-1.0D0,
     .             1.0D0, 2*-1.0D0, 2*1.0D0, 2*-1.0D0,

     .             3*1.0D0, -1.0D0, 3*1.0D0, 2*-1.0D0, -1.0D0, 2*1.0D0,
     .             2*-1.0D0, 2*1.0D0,

     .             1.0D0, 3*-1.0D0, 5*1.0D0, 7*-1.0D0,

     .             1.0D0, 1.0D0, -1.0D0, 3*1.0D0, -1.0D0, 1.0D0, -1.0D0,
CCC  X 2*1.0D0,2*-1.0D0,2*1.0D0,-1.0D0,
     .             2*1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0,

     .             1.0D0, -1.0D0, 4*1.0D0, -1.0D0, -1.0D0, 1.0D0,
CCC  X 1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0,
     .             1.0D0, -1.0D0, 1.0D0, 1.0D0, -1.0D0, -1.0D0, 1.0D0,

     .             1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0, 1.0D0, -1.0D0,
     .             1.0D0, -1.0D0,
CCC  X 2*-1.0D0,2*1.0D0,2*-1.0D0,1.0D0,
     .             2*-1.0D0, 1.0D0, -1.0D0, 1.0D0, -1.0D0, 1.0D0,

     .             1.0D0, 1.0D0, -1.0D0, -1.0D0, 1.0D0, 1.0D0, 2*-1.0D0,
     .             1.0D0,
CCC  X -1.0D0,1.0D0,-1.0D0,1.0D0,-1.0D0,1.0D0,-1.0D0/
     .             -1.0D0, 1.0D0, -1.0D0, -1.0D0, 1.0D0, 1.0D0, -1.0D0/
C
C-----------------------------------------------------------------------
C
C     SGN2 IS THE PHASE (SIGN) OF ORBITAL J UNDER OPERATOR I.
C     SGN2(J,I): I= OPERATOR = 1 FOR C2,=2 FOR SIGH, =3 FOR I,
C      =4 FOR SIG(XZ), =5 FOR SIG(YZ),=6 FOR C2(Y),=7 FOR C2(X).
C     J=2 A PX, =3 A PY,=4 A PZ, =5 A D(Z2),=6 A D(X2-Y2),
C      = 7 A D(XY), = 8 A D(XZ), = 9 A D(YZ), = 10 A F(Z3),
C      = 11 A F(XZ2), = 12 A F(YZ2), = 13 A F(XYZ), = 14 A F(Z(X2-Y2)),
C      = 15 A F(X(X2-3Y2), = 16 A F(Y(3X2-Y2)
C
C-----------------------------------------------------------------------
      DATA         JOPER/1, 7*0, 2, 7*0, 3, 7*0, 1, 4, 5, 5*0, 1, 2, 3,
     .             5*0, 1, 6, 7, 5*0, 1, 6, 7, 3, 2, 4, 5, 0, 8*0/
C-----------------------------------------------------------------------
C
C     JOPER(J,I) THE J  OPERATOR AS GIVEN UNDER SGN2 ABOVE FOR THE
C     I'TH GROUP
C
C-----------------------------------------------------------------------
C
C
C
      write(6,*)"i am in symorb"
      IFLAG = ITSYM
      MB = NORB 
C     WRITE(6,*) ' VECTORS IN SYMORB'
C     IF (JTEN.EQ.3) CALL PRNT(C,MB,MB*MB)
      IF (ITSYM.EQ.0 .OR. ISYM.LE.1) GO TO 380
C     CALCULATE NFIRST
      NOLD = 0
      DO 3 I = 1, MXA
      NFIRST(I) = 0
    3 CONTINUE
      DO 5 I = 1, MB   
      IF(NU(I).NE.NOLD) THEN
         NFIRST(NU(I)) = I
         NOLD = NU(I)
      ENDIF
    5 CONTINUE              
      NAT = NOLD
      DO 10 I = 1, MB
         IIMO(I) = 0
   10 CONTINUE
CMK
C      WRITE(6,*) ' IN SYMORB, ITSYM, ISYM = ', ITSYM, ISYM
C      CALL SYMMK(C, MB, ITSYM, ISYM, IIMO, EGN, IERR)
C      IF (IERR .EQ. 0) GOTO 500
CMK
      ISSS = ISYM
      IF ((ISYM/2)-2) 20, 30, 40
   20 CONTINUE
      CALL SGRPOR(C,EGN,NU,NW,NFIRST,JOPER,IWO,SGN2,G2,ISSS,NORB,NAT,
     x            JTEN,IFLAG)
      GO TO 50

   30 CONTINUE
      CALL SGRPOR(C,EGN,NU,NW,NFIRST,JOPER,IWO,SGN2,G4,ISSS,NORB,NAT,
     x            JTEN,IFLAG)
      GO TO 50

   40 CONTINUE
      CALL SGRPOR(C,EGN,NU,NW,NFIRST,JOPER,IWO,SGN2,G8,ISSS,NORB,NAT,
     x            JTEN,IFLAG)
   50 CONTINUE

      IF (IFLAG.LT.0) GO TO 370
      IF (ITSYM.EQ.0) GO TO 380
      LL = 0
      DO 100 I = 1, MB
         IF (IIMO(I).NE.0) LL = LL + 1
  100 CONTINUE
      DO 110 I = 1, NAT
         IWO(I) = 0
  110 CONTINUE
      GO TO (120,120,120,200,200,200,290), ITSYM
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
  120 IM = JOPER(1,ITSYM)
      DO 190 I = 1, NAT
         IF (IWO(I).NE.0) GO TO 190
         I2 = NCON(I,1)
         JT1 = NFIRST(I)
         IF(JT1.EQ.0) GO TO 190
         JT2 = NFIRST(I2)
         IWO(I) = 1
         IWO(I2) = 1
         NT = 0
  130    CONTINUE
         NT = NT + 1
         IF (I.EQ.NU(NT+JT1)) GO TO 130
         DO 180 IB = 1, NT
            IBB = NW(JT1+IB-1) + 1
            DO 170 J = 1, MB
               IF (IIMO(J).NE.0) GO TO 160
               IST = (J-1)*MB + IB - 1
               CALL SCLASS(C(IST+JT1),K)
               IF (K.NE.1) GO TO 160
               DO 140 IIREP = 1, 2
                  CCC = ABS(C(IST+JT1)+G2(IIREP,2)*SGN2(IBB,IM)*
     .                  C(IST+JT2))
                  IF (CCC.LT.ABS(C(IST+JT1))*1.998D+00) GO TO 140
                  IIMO(J) = IIREP
                  LL = LL + 1
                  GO TO 150

  140          CONTINUE
  150          CONTINUE
  160          IF (LL.GE.MB) GO TO 400
  170       CONTINUE
  180    CONTINUE
  190 CONTINUE
      GO TO 370
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
  200 IM = JOPER(1,ITSYM)
      IN = JOPER(2,ITSYM)
      IP = JOPER(3,ITSYM)
      DO 270 I = 1, NAT
         IF (IWO(I).NE.0) GO TO 270
         I2 = NCON(I,1)
         I3 = NCON(I,2)
         I4 = NCON(I,3)
         JT1 = NFIRST(I)
         IF(JT1.EQ.0) GO TO 270
         JT2 = NFIRST(I2)
         JT3 = NFIRST(I3)
         JT4 = NFIRST(I4)
         IWO(I) = 1
         IWO(I2) = 1
         IWO(I3) = 1
         IWO(I4) = 1
         NT = 0
  210    CONTINUE
         NT = NT + 1
         IF (I.EQ.NU(NT+JT1)) GO TO 210
         DO 260 IB = 1, NT
            IBB = NW(JT1+IB-1) + 1
            DO 250 J = 1, MB
               IF (IIMO(J).NE.0) GO TO 240
               IST = (J-1)*MB + IB - 1
               CALL SCLASS(C(IST+JT1),K)
               IF (K.NE.1) GO TO 240
               DO 220 IIREP = 1, 4
                  CCC = ABS(C(IST+JT1)+G4(IIREP,2)*C(IST+JT2)*
     .                  SGN2(IBB,IM)+G4(IIREP,3)*C(IST+JT3)*
     .                  SGN2(IBB,IN)+G4(IIREP,4)*C(IST+JT4)*
     .                  SGN2(IBB,IP))
                  IF (CCC.LT.ABS(C(IST+JT1))*3.98D+00) GO TO 220
                  IIMO(J) = IIREP
                  LL = LL + 1
                  GO TO 230

  220          CONTINUE
  230          CONTINUE
  240          IF (LL.GE.MB) GO TO 400
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
      GO TO 370
C-----------------------------------------------------------------------
C
C
C    D2H
C-----------------------------------------------------------------------
  290 IM = 1
      IN = 6
      IP = 7
      IQ = 3
      IR = 2
      IS = 4
      IT = 5
      DO 360 I = 1, NAT
         IF (IWO(I).NE.0) GO TO 360
         I2 = NCON(I,1)
         I3 = NCON(I,2)
         I4 = NCON(I,3)
         I5 = NCON(I,4)
         I6 = NCON(I,5)
         I7 = NCON(I,6)
         I8 = NCON(I,7)
         JT1 = NFIRST(I)
         IF(JT1.EQ.0) GO TO 360
         JT2 = NFIRST(I2)
         JT3 = NFIRST(I3)
         JT4 = NFIRST(I4)
         JT5 = NFIRST(I5)
         JT6 = NFIRST(I6)
         JT7 = NFIRST(I7)
         JT8 = NFIRST(I8)
         IWO(I) = 1
         IWO(I2) = 1
         IWO(I3) = 1
         IWO(I4) = 1
         IWO(I5) = 1
         IWO(I6) = 1
         IWO(I7) = 1
         IWO(I8) = 1
         NT = 0
  300    CONTINUE
         NT = NT + 1
         IF (I.EQ.NU(NT+JT1)) GO TO 300
         DO 350 IB = 1, NT
            IBB = NW(JT1+IB-1) + 1
            DO 340 J = 1, MB
               IF (IIMO(J).NE.0) GO TO 330
               IST = (J-1)*MB + IB - 1
               CALL SCLASS(C(IST+JT1),K)
               IF (K.NE.1) GO TO 330
               DO 310 IIREP = 1, 8
                  CCC = ABS(C(IST+JT1)+C(IST+JT2)*G8(IIREP,2)*
     .                  SGN2(IBB,IM)+C(IST+JT3)*G8(IIREP,3)*
     .                  SGN2(IBB,IN)+C(IST+JT4)*G8(IIREP,4)*
     .                  SGN2(IBB,IP)+C(IST+JT5)*G8(IIREP,5)*
     .                  SGN2(IBB,IQ)+C(IST+JT6)*G8(IIREP,6)*
     .                  SGN2(IBB,IR)+C(IST+JT7)*G8(IIREP,7)*
     .                  SGN2(IBB,IS)+C(IST+JT8)*G8(IIREP,8)*
     .                  SGN2(IBB,IT))
                  IF (CCC.LT.ABS(C(IST+JT1))*7.98D+00) GO TO 310
                  IIMO(J) = IIREP
                  LL = LL + 1
                  GO TO 320

  310          CONTINUE
  320          CONTINUE
  330          IF (LL.GE.MB) GO TO 400
  340       CONTINUE
  350    CONTINUE
  360 CONTINUE
  370 CONTINUE
      WRITE (6,9000) ISYM
      WRITE (6,9010) (IIMO(I),I=1,MB)
      WRITE (6,9020)
C     BECAUSE OF DEGENRACIES IN THE CORE ORBOTALS THE
C     SYMMETRY MIGHT NOT BE ASSIGNED CORRECTLY. DO NOT
C     RESET INPUT SYMMETRY!
      IFLAG = -1
      IF (JTEN.NE.3) ISYM = 1
  380 CONTINUE
      DO 390 I = 1, MB
         IIMO(I) = 1
  390 CONTINUE
  400 CONTINUE
      IF (ISWANN.EQ.0 .OR. ISWANN.EQ.3) GO TO 430
      GO TO (410,430,450), JTEN

  410 CONTINUE
      DO 420 J = 1, MB
         IMA(J) = IIMO(J)
  420 CONTINUE
  430 CONTINUE
      IF (JTEN.LT.3) THEN
         DO 440 I = 1, NORB
            IMO(I) = IIMO(I)
  440    CONTINUE
      END IF
  450 RETURN

 9000 FORMAT (
     .      '  ***** M.O.S CANNOT BE CLASSIFIED ACCORDING TO SYMMETRY',
     .       I5,' *****')
 9010 FORMAT (20I5)
 9020 FORMAT (' M.O.S ABOVED FLAGGED WITH ZERO ARE NOT CLASSIFIED')

      END
C***********************************************************************
      SUBROUTINE SYMMK(VEC1,NB,ITSYM,ISYM,IIMO,EIG,IERR)
C***********************************************************************
C
      INCLUDE '../include/param.cmn'
      INCLUDE '../include/units.cmn'
      INCLUDE '../include/const.cmn'
C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION    VEC1(*), VEC2(1,1), EIG(MXB)
C     DIMENSION    VEC1(*), VEC2(MXB,MXB), EIG(MXB)
      INTEGER      IST(16), IEN(16), IIMO(MXB), MAT(8,10)
      CHARACTER*6  BSYM
      CHARACTER*4  LAB(16), LABEL(MXB), BLANK
      CHARACTER*3  HH(10)
C
      DATA         HUNDR/100.0D0/
      DATA         HH/' C2', 'C1H', ' CI', 'C2V', 'C2H', ' D2', 'D2H',
     .             ' D2', 'D2H', 'C1H'/
      DATA         BLANK/'    '/
      DATA         MAT/1, 2, 6*0, 1, 2, 6*0, 1, 2, 6*0, 1, 2, 3, 4, 4*0,
     .             1, 4, 2, 3, 4*0, 1, 2, 3, 4, 4*0, 1, 2, 3, 4, 5, 6,
     .             8, 7, 1, 2, 3, 4, 4*0, 1, 2, 3, 4, 5, 6, 8, 7, 1, 2,
     .             6*0/
C
      IF(MXB.GT.1) THEN
       WRITE(6,*) ' ***** SYMMK ERROR *****FIX DIMENSION VEC2***'
       CALL ZEXIT(3)
      ENDIF
      IO_SYMORB=-1
C     This unit (98 in original code) has no connection in any otherpart
C     of Zindo. Please note that it is also not in units.cmn
      IF (IO_SYMORB.LE.0) THEN
         PRINT*,'IO_SYMORB IS NOT ASSIGNED IN SYMMK (SYMROB.F)'
         CALL ZEXIT (3)
      END IF
      REWIND (IO_SYMORB)
      READ (IO_SYMORB,END=60) BSYM, NIR
      READ (IO_SYMORB) (IST(JJJ),JJJ=1,NIR)
      READ (IO_SYMORB) (IEN(JJJ),JJJ=1,NIR)
      READ (IO_SYMORB) (LAB(JJJ),JJJ=1,NIR)
      READ (IO_SYMORB) ((VEC2(JJ1,JJ2),JJ1=1,NB),JJ2=1,NB)
C
      INDX = INDEX(BSYM,HH(ITSYM))
      IF (INDX.EQ.0) THEN
         WRITE (6,*) ' ***** WRONG SYMMETRY MATRIX ON UNIT 98 *****'
         GO TO 60

      END IF
      IF (NIR.NE.ISYM) THEN
         WRITE (6,*) ' ***** NUMBER OF IRREPS NOT CONSISTENT *****'
         GO TO 60

      END IF
C
      DO 40 I = 1, NB
         LABEL(I) = BLANK
         DO 30 L = 1, NIR
            SUM = ZERO
            DO 20 J = IST(L), IEN(L)
C
               X = ZERO
               KK = (I-1)*NB
               DO 10 K = 1, NB
                  KK = KK + 1
                  X = X + VEC1(KK)*VEC2(K,J)
C                  X  = X + VEC1(K, I) * VEC2(K, J)
   10          CONTINUE
               SUM = SUM + X*X
C
   20       CONTINUE
            IF (SUM.GT.0.9D0) THEN
               LABEL(I) = LAB(L)
               IIMO(I) = MAT(L,ITSYM)
            ELSE IF (SUM.GT.0.1) THEN
               WRITE (6,*) ' ERROR: SUM(', I, ') = ', SUM, ' LAB = ',
     .           LAB(L)
            END IF
   30    CONTINUE
   40 CONTINUE
C
      NLEFT = NB
      DO 50 J = 1, NB, 6
         NT = MIN0(6,NLEFT) + J - 1
         WRITE (6,9030) (JJJ,JJJ=J,NT)
         WRITE (6,9040) (LABEL(JJJ),JJJ=J,NT)
         EEIG = EIG(J)*EVOLT
         IF (ABS(EEIG).LT.1000.0D0) THEN
            WRITE (6,9050) (EIG(JJJ)*EVOLT,JJJ=J,NT)
         ELSE
            WRITE (6,9060) (EIG(JJJ)*EVOLT,JJJ=J,NT)
         END IF
         NLEFT = NLEFT - 6
   50 CONTINUE
C
cmk
      WRITE (6,9000)
      WRITE (6,9010) (I,IIMO(I),I=1,NB)

 9000 FORMAT ('      M.O. AND M.O. SYMMETRY'//)
 9010 FORMAT (10 (2X,2I3))
cmk
 9020 FORMAT (10 (I5,2X,A4))
 9030 FORMAT (' NO.',16X,6 (I5,5X))
 9040 FORMAT (' SYMMETRY',8X,6 (6X,A4))
 9050 FORMAT (' VALUE(EV)',7X,6F10.2,//)
 9060 FORMAT (' VALUE(EV)',7X,6F10.1,//)
C
      IERR = 0
      RETURN

   60 IERR = 1
      RETURN

      END
