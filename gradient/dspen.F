      FUNCTION DSPEN(I,RR,CINTB,XX,NP,ALPHA,FACT,IONE,NB)
C     CALCULATES THE DERIVATIVE OF THE ELECTRON NUCLEAR ATTRACTION
C     I IS THE AO, AND THE NUCLEUS IS RR ANGSTROMS AWAY.
C     DIFPEN IS THE EXACT INTEGRAL, XX IS THE ZDO APPROXIMATION
C     AVERAGED IN.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
      DIMENSION NP(*), ALPHA(NB,6), FACT(*)
c     write(6,*)"i am in dspen"
      IDD1 = 0
      IF (IDD1.EQ.1) THEN
         YY = 0.0D0
         NK = NP(I)
         NL = NP(I)
         DO 20 JJ = 2, 3
            AMU = ALPHA(I,JJ)
            CMU = ALPHA(I,JJ+3)
            IF (AMU.LT.1.D-4) GO TO 20
            DO 10 JJJ = 2, 3
               BMU = ALPHA(I,JJJ)
               DMU = ALPHA(I,JJJ+3)
               IF (BMU.LT.1.D-4) GO TO 10
               CALL DIFPEN(FACT,AMU,NK,BMU,NL,RR,PEN,IONE)
               YY = YY + PEN*CMU*DMU
   10       CONTINUE
   20    CONTINUE
         DPDR = 0.3D0*YY + 0.7D0*XX
      ELSE
         DPDR = XX
      END IF
      DSPEN = DPDR
      RETURN

      END
