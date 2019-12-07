C-----------------------------------------------------------------------

      SUBROUTINE IMPROV(NB,NU,DK,H,ISW1)
C-----------------------------------------------------------------------
C     THIS ROUTINE IMPROVES THE STARTING HAMILTONIAN FOR THE SCF.
C-----------------------------------------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/geop.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
C---------------------------------------------------------------------
c     LOGICAL      EVEN
      DIMENSION    H(*),NU(*),DK(*)
C---------------------------------------------------------------------
c     EVEN(N) = MOD(N,2) .EQ. 0
C---------------------------------------------------------------------
c     write(6,*)"i am in improv"!,(H(I),I=1,40)
      LL = 0
      DO 10 L = 1, NB
         LL = LL + L
         KA = NU(L)
c        IF (NEXTR(L).EQ.1) THEN
c           H(LL) = H(LL) - 14.0
c        END IF
C       IF ((DK(KA) .NE. 0) .AND. (ISW1 .NE. 7)) THEN
c        IF ((ABS(DK(KA)-0.1D-13).GE.0.1D-5) .AND. (ISW1.NE.7)) THEN
c           H(LL) = H(LL)*Z12(KA)/DK(KA)
c        END IF
   10 CONTINUE

C-----------------------------------------------------------------------
C     FOR BETTER STARTING MATRIX DIVIDE ION. POT. BY 2
C-----------------------------------------------------------------------
      IF (ISW1.EQ.0.OR.ISW1.EQ.2.OR.ISW1.EQ.4) GO TO 5000
         LL = 0
         DO 20 L = 1, NB
            LL = LL + L
            H(LL) = H(LL)/2.0D0
   20    CONTINUE

 5000 RETURN

      END
