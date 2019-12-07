      SUBROUTINE SCONS(D12,D13,IT,ITMAX,NA,SCFTOL,EE,JSLIT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
      COMMON       /DIISX/DIIS, DIISCHK, DDFLAG
      DIMENSION    D12(*),D13(*),JSLIT(4)
      SAVE         PREVS
      LOGICAL      CONVERGED,DIIS,DIISCHK,DDFLAG

C     
C     This is the energy convergence criterion, close to
C     machine accuracy:
      DATA TEST /1.0D-8/
      IQMMM = 0
c     DATA TEST/1.0D-23/
c     DATA TEST /1.0D-30/
      CONVERGED = .FALSE.
c     write(6,*)"i am in scons" ,scftol
C    
      IF (IT .GT. 1) THEN
C     
C     Form the maximum change in density - we may need it for print
C     
       IMAX = IZERO
       WMAX = ZERO
       DO 20 I = 1, NA
        W = ABS(D12(I)-D13(I))
        IF (W.GT.WMAX) THEN
c        write(6,*)"i am here"
         IMAX = I
         WMAX = W
        END IF
   20  CONTINUE
C     
C     Check for number of iterations
C     
       IF (IT.GE.ITMAX) THEN
        if (ga_nodeid().eq.0) WRITE (6,9030)
c       write(6,*)"1st criteria"
        CONVERGED=.TRUE.
       END IF

C     
C     Test for energy
C     
        DIFF = PREVS - EE
c       write(6,*)"DIFF/EE",DIFF/EE
       IF (ABS(DIFF/EE).LE.TEST .AND. DIFF.GE.0.0D0) THEN
        if (ga_nodeid().eq.0) WRITE (6,9040)
       CONVERGED=.TRUE.
c        write(6,*)"2nd criteria"
C     The next ones can be turned off by setting a small SCFTOL:
        ELSE IF (ABS(SCFTOL).GT.1.0D-6) THEN
C     
C     Check commutator in DIIS
C     
c        IF (SMALL.LE.1.0D-10) THEN
c          WRITE (6,9060)
c         CONVERGED=.TRUE.
c        write(6,*)"3rd criteria"
c        ELSE
C     
C     Check density
C     
          IF (WMAX.LT.ABS(SCFTOL)) THEN
           if (ga_nodeid().eq.0) WRITE (6,9050)
           CONVERGED=.TRUE.
c         write(6,*)"4th criteria"
          END IF
         END IF
        END IF
c      END IF
C     
C     Tests done, how are we doing?
C     
      IF (CONVERGED) THEN
C     
C     Turn on sense light and turn off DIIS.
C     
C      THE CALCULATION HAS CONVERGED
       CALL SLITE(2,JSLIT)
C      PRINT THE RESULTS
       CALL SLITE(3,JSLIT)
       ELSE
C     
C     All tests show that there is no convergence yet.
C     Go for the next cycle.
C     
C     Record energy to be used at next iteration:
       PREVS = EE
      END IF
      IF (IQMMM .EQ. IZERO.AND.ga_nodeid().eq.0) THEN
       IF ((NA .GT. 30) .AND. (IT .GT. 1))  WRITE (6,9020) WMAX,IMAX
      ENDIF
   30 CONTINUE
 5000 RETURN

 9000   FORMAT (/,'  SCONS: D12 (NEW)',/,5 (2X,10F10.6,/))
 9010   FORMAT (/,'  SCONS: D13 (OLD)',/,5 (2X,10F10.6,/))
 9020 FORMAT (/,' MAXIMUM SCF DENSITY DIFFERENCE=',F10.6,
     .   ' FOR ATOM',I5)
 9030 FORMAT (//' SENSE LIGHT 2 IS ON (TOO MANY ITERATIONS)'//)
 9040 FORMAT (//' SENSE LIGHT 2 IS ON (ENERGY)'//)
 9050 FORMAT (//' SENSE LIGHT 2 IS ON (DENSITY)'//)
 9060 FORMAT (//' SENSE LIGHT 2 IS ON (UDIIS)'//)


      END



