      SUBROUTINE MOLPAB(N1,N2,L1,L2,M1,M2,SK1,SK2,A,B,FACT,R,VEST,NIN)
C
C     THIS SUBROUTINE CALCULATES TWO CENTERED OVERLAP INTEGRALS
C
c     INCLUDE '../include/param.cmn'
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      DIMENSION    A(35), B(35), FACT(35), NIN(*)
c
c      if (ga_nodeid().eq.0)
c    .write(6,*) "i am in molpab"
c     write(6,*)N1,N2,L1,L2,M1,M2,SK1,SK2
      VEST = 0.0D0
      IF (M1.NE.M2) GO TO 30
      STRAD = 0.0D0
      RHOA = R*SK1
      RHOB = R*SK2
      F1 = FACT(2*N1+1)
      F2 = FACT(2*N2+1)
      IFF = L1 - M1 + 1
      F3 = FACT(IFF)
      IFF = L2 - M2 + 1
      F4 = FACT(IFF)
      IFF = L1 + M1 + 1
      F5 = FACT(IFF)
      IFF = L2 + M2 + 1
      F6 = FACT(IFF)
C      TERM=2.0D0**(N1+N2-L1-L2)*SQRT(DFLOAT((2*L1+1)*(2*L2+1))*
C     1   (F3/F1/F5)*(F4/F2/F6)*(SK1/SK2)**(2*N1+1))
      Q1 = (SK1**N1/SQRT(F1))*R**N1
      Q2 = (SK2**N2/SQRT(F2))*R**N2
      TERM = SQRT(DFLOAT((2*L1+1)* (2*L2+1))* (F3/F5)* (F4/F6)*RHOA*
     .       RHOB)/2.D0** (L1+L2+1)
      TERM = Q1*TERM*Q2
      JEND = 1 + ((L1-M1)/2)
      KEND = 1 + ((L2-M2)/2)
      DO 20 J = 1, JEND
         JU = J - 1
         IFF = 2*L1 - 2*JU + 1
         F11 = FACT(IFF)
         IFF = L1 - M1 - 2*JU + 1
         F13 = FACT(IFF)
         F15 = FACT(JU+1)
         IFF = L1 - JU + 1
         F17 = FACT(IFF)
         DO 10 K = 1, KEND
            KU = K - 1
            IFF = 2*L2 - 2*KU + 1
            F12 = FACT(IFF)
            IFF = L2 - M2 - 2*KU + 1
            F14 = FACT(IFF)
            F16 = FACT(KU+1)
            IFF = L2 - KU + 1
            F18 = FACT(IFF)
c        write(6,*)"cfunct",M1,RHOA,RHOB
            CALL CFUNCT(N1-L1+2*JU,N2-L2+2*KU,L1-M1-2*JU,L2-M2-2*KU,M1,
     .                  RHOA,RHOB,VALUE,FACT,A,B,NIN)
c            write(6,*)"VALUE",VALUE
            STRAD = STRAD + VALUE* (F11/F13/F15/F17)* (F12/F14/F16/F18)*
     .               (-1)** (JU+KU)
   10    CONTINUE
   20 CONTINUE
      VEST = TERM*STRAD
   30 RETURN

      END
