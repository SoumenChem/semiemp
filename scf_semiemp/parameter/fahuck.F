C----------------------------------------------------------------------

      SUBROUTINE FAHUCK(NB,NRR,NW,NU,FA,H,S,IAPX,UCOR)
C-----------------------------------------------------------------------
C     INITIALIZE THE ONE-ELECTRON MATRIX TO THE EXTENDED HUCKEL
C     HAMILTONIAN.
C     H IS THE HUEKEL-LIKE STARTING MATRIX, FA THE ONE-EL MATRIX
C     S IS THE WEIGHTED OVERLAP
C-----------------------------------------------------------------------

c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/z.cmn'
c     INCLUDE '../include/datain.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/factor.cmn'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION    S(*),FA(*),H(*),NW(NB),NU(NB),UCOR(4,86)
C-----------------------------------------------------------------------
C
C     AB-INITIO DON'T MESS WITH ONE EL MATRIX
C     BUT SCALE THE OFF DIAGONAL OF THE STARTING MATRIX
c     write(6,*)"i am in fahuck" ,(S(I),I=1,nrr)
      IF(IAPX.EQ.9) THEN
         ID = 1
         DO 5 I = 2, NB
            DO 4 J = 1, I-1
               IJ = ID+J
C              THIS SCALING OF THE OFF DIAG. STARTING MATRIX IS NEEDED
C              TO SET THE PHASES FOR METHODS WITH OVERLAP
               H(IJ) = H(IJ)*1.89D0
    4       CONTINUE
            ID = ID + I
    5    CONTINUE
         RETURN
      ENDIF
      IF(IAPX.NE.8) THEN
         DO 10 L = 1, NRR
c      write(6,*)"H,FA",H(L),FA(L)
            FA(L) = H(L) + FA(L)
   10    CONTINUE
      ELSE
C        XINDO  DO NOT RESET THE OFF-DIAGONAL ELEMENTS OF ONE -EL.
C        BUT SCALE THE STARTING MATRIX
         LD = 0
         DO 12 L = 1, NB
            LD = LD + L
            FA(LD) = H(LD) + FA(LD)
   12    CONTINUE
C        ID = 1
C        DO 15 I = 2, NB
C           DO 14 J = 1, I-1
C              IJ = ID+J
C              THIS SCALING OF THE OFF DIAG. STARTING MATRIX IS NEEDED
C              TO SET THE PHASES FOR METHODS WITH OVERLAP
C              H(IJ) = H(IJ)*1.89D0
C  14       CONTINUE
C           ID = ID + I
C  15    CONTINUE
      ENDIF
C     
C     WORRY ABOUT THE DIAGONAL
C     
      ID = 0
      DO 20 I = 1,NB
         ID=ID+I
         FAC = 0.0d0
c        IF(NEXTR(I).EQ.1) THEN
c           FA(ID) = H(ID) -20.0D0
c           FAC = 20.0D0
c        ENDIF
         NATOM = NU(I)
         NORB = NW(I)
         IF(NORB.GT.8) THEN
            NL = 4
         ELSEIF(NORB.GT.3) THEN
            NL = 3
         ELSEIF(NORB.GT.0) THEN
            NL = 2
         ELSE
            NL = 1
         ENDIF
c        IF(ABS(UCOR(NL,NATOM)).GT.1.0D-7) THEN 
c           THE CORE INTEGRAL HAS BEEN INPUTTED IN A TAB
c           write(6,*)"FA,H,FAC1",FA(ID),H(ID),FAC
c           FA(ID) = FA(ID) - H(ID) + FAC
c        ELSEIF(IAPX.GE.50.AND.IAPX.LT.54) THEN
          IF(IAPX.GE.50.AND.IAPX.LT.54) THEN
C           AM1 HAS ITS OWN CORE INTEGRAL
c           write(6,*)"FA,H,FAC2",FA(ID),H(ID),FAC
            FA(ID) = FA(ID) - H(ID) + FAC
         ENDIF
   20 CONTINUE

C     UNCOUPLES FROM THE PARAMETERS FOR THE STARTING MATRIX
C     PARAMETERS FOR THE MODEL HAMILTONIAN AS BETB
C
C     THIS IS DONE IN UFORM ---

C     IF(ISW1.EQ.8) CALL PARAMX

C     ID = 1
C     DO 40 I = 2, NB
C     KI = KAN(NU(I))
C     IF(KI.EQ.0) GO TO 40
C     DO 30 J = 1, I-1
C     KJ = KAN(NU(J))
C     IF(KJ.EQ.0) GO TO 30
C     IJ = ID + J
C     X = (BETB(KI,LQNO(I)+1) + BETB(KJ,LQNO(J)+1))/(TWO*EVOLT)
C     FA(IJ) = FA(IJ) + S(IJ)*X
C  30 CONTINUE
C     ID = ID + I
C  40 CONTINUE
c     write(6,*)"H(I)",(H(I),I=1,nrr)
c     write(6,*)"FA(I)",(FA(I),I=1,nrr)
      RETURN
      END
