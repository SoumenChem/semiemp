      SUBROUTINE OVLAP(LH,N1,L1,AMU,N2,L2,BMU,R,S,SP,SD,SF,FACT,PP,NIN)
C
C     EVALUATES OVERLAP BETWEEN TWO SLATER TYPE ORBITALS.
C
C     N1,  N2  are the principal quantum numbers
C     L1,  L2  are the secondary quantum numbers
C     AMU, BMU are the exponential constants
C     R  is the separation in angstroms
C
C     S, SP, SD and SF are the sigma, pi, delta and phi components
C     of the overlap, to be put together using LH and
C     subroutine GEOM.
C
C     FACT(I+1)= FACTORIAL I.
C
c      INCLUDE '../include/param.cmn'
C      INCLUDE '../include/const.cmn'
C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C
c      INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON / EV /  A(35), B(35), DD(35)
C****************************************************************
c      INCLUDE '../include/pptt.cmn'
C****************************************************************
C***      COMMON / PPTT / PP
C****************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (BOHRS    = 5.29177249D-1,
     .           EVOLT    = 2.72113961D1,
     .           ENUBAR   = 2.1947463067D5,
     .           DEBYE    = 2.541748047D0,
     .           EVTOCM   = 8065.541D0,
     .           ESUQD    = 1.345035248D0,
     .           AVOGAD   = 6.0221367D0,
     .           PI       = 3.14159265358979323846D0,
     .           AUTOKC   = 627.510D0,
     .           ECHARG   = 1.60217733D0,
     .           AUTOJ =  4.3597482D-18      )
      DIMENSION    FACT(35), SS(4)
      parameter   (ZERO=0.0D0, ONE = 1.0D0, TWO = 2.0D0)
      dimension A(35),B(35),NIN(*)

C
c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in ovlap",NIN(1)
      IF (PP.EQ.-10.D0) THEN
         TT = ZERO
      END IF
      RR = R
      P = (AMU+BMU)*RR/TWO
      T = (AMU-BMU)/ (AMU+BMU)
c     write(6,*)"AMU,BMU,RR",AMU,BMU,RR
c     write(6,*)"here 1",R
C
C     If the two STO's are on the same center,
C     evaluate one center overlap
C
      IF (R.LE.0.001) GO TO 30
C
c      if (ga_nodeid().eq.0)
c    .write(6,*) "i am here R" ,R
      IF (ABS(T).LE.0.0001) T = ZERO
C
C     This checking using PP only saves time;
C     AUX.F is a very time consuming routine.
C
C     IF (P.EQ.PP .AND. T.EQ.TT) GO TO 10
      PP = P
      TT = T
C
C     NM1 = N1 + N2 + 1
C
c     write(6,*)"P,T",P,T
      CALL AUC(P,T,A,B)
C
   10 CONTINUE
C
C     Initialize matrix SS(4)
C
      SS(2) = ZERO
      SS(3) = ZERO
      SS(4) = ZERO
C
      LMIN = L1
      IF (L2.LT.L1) LMIN = L2
      LMIN = LMIN + 1
C
      DO 20 II = 1, LMIN
         I = II - 1
C
C     Evaluate overlap in subroutine MOLPAB
C
c     write(6,*)N1,N2,L1,L2,I,I,AMU,BMU         
         CALL MOLPAB(N1,N2,L1,L2,I,I,AMU,BMU,A,B,FACT,RR,SSSS,NIN)
c        write(6,*) "SSSS",SSSS
         SS(II) = SSSS
   20 CONTINUE
C
C     Update matrix SS(4)
C
      S = SS(1)
      SP = SS(2)
      SD = SS(3)
      SF = SS(4)
      GO TO 40
C
C     One center overlap ( R = 0 )
C
   30 S = ZERO
      SP = ZERO
      SD = ZERO
      SF = ZERO
C
C     For two AO's of different symmetries (e.g. sigma, pi)
C     centered on the same atom, the overlap is S=0 --> RETURN
C
      IF (L1.NE.L2) GO TO 40
C
      NN = 2*N1 + 1
      MM = 2*N2 + 1
      F = (ONE-T)**MM* (ONE+T)**NN
      S = FACT(N1+N2+1)*SQRT(F/ (FACT(NN)*FACT(MM)))
      IF (L1.GT.0) SP = (-ONE)** (L1+1)*S
      IF (L1.GT.1) SD = (-ONE)** (L1+2)*S
      IF (L1.GT.2) SF = (-ONE)** (L1+3)*S
      S = (-ONE)**L1*S
C
   40 RETURN

      END

