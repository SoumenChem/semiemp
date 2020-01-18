      SUBROUTINE twoint_semiemp(nat,TB,NB,ALPHA,RAD,GAMMA,NW,NP,KAN,NU
     1 ,NG,ISW1,ISPEC,PO,NIN)
C
C     This subroutine organizes the Coulomb matrix GAMMA(I,J)
C

C     JSW(1) contains switch  ITRIP
C
      integer nat,ISPEC,ISW1
      double precision  GAMMA(*),INT,R
      double precision RAD(nat,nat)
      integer NB(nat)
      integer NP(nat)
      integer TB !total number of basis
      integer NU(TB) !type of atom array
      double precision ALPHA(TB,6)
      integer I,J,IJ,LA,KA,L,MA,KG,MN,MMM
      integer M
      double precision AMW,AMU,BMU
      integer NW(TB)
      integer NG(TB)
      integer KAN(TB),NIN(*)
      logical COULB,COUL, status
      double precision F0GAM
      double precision FACT(40), PO(*)
      data FACT/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
     X          ,3628800.0,39916800.0,479001600.0,6227020800.0,
     X          87178291200.0,1307674368000.0,20922789888000.0,
     X          3.55687428096000E+014,6.402373705728000E+015,
     X          1.21645100408832000E+017,2.43290200817664000E+018,
     X          5.10909421717094400E+019,1.12400072777760768E+021,
     X          2.58520167388849782E+022,6.20448401733239410E+023,
     X          1.55112100433309861E+025,4.03291461126605650E+026,
     X          1.08888694504183519E+028,3.04888344611713837E+029,
     X          8.84176199373970077E+030,2.65252859812191032E+032,
     X          8.22283865417792243E+033, 2.63130836933693518E+035,
     X          8.68331761881188594E+036,2.95232799039604120E+038,
     X        1.0333148D+040,3.7199333D+041,1.3763753D+043,5.2302262D+044
     X          ,2.0397882D+046/

c      write(6,*)"ISPEC",ISPEC
C
c      WRITE (6,9000)
c9000 FORMAT ('  >>> Theoretical integrals for Rydbergs <<< '//)
      MN = 0
      KG = 0
C
C     For ground state CNDO
C
      DO 120 I = 1, TB
C
C  NW(I) is the symmetry of AO I:  0    1    2    3    4  ...
C                                  s   px   py   pz   dz2 ...
c        write(6,*) "i am in twoint"
         LA = NW(I)
C
C        IF(LA.GT.1.AND.LA.NE.4.AND.LA.NE.9) GO TO 69
         IF (LA.LE.1) GO TO 10
         IF (LA.EQ.4) GO TO 10
         IF (LA.NE.9) GO TO 110
C
C        NP(I) is the principal quantum number of atomic orbital I
C
   10    K = NP(I)
c        write(6,*)"K",NP(I)
         KA = NU(I)
         L = KAN(I)
         KG = KG + 1
         AMW = ALPHA(I,1)
C
C
         DO 100 J = 1, I
            LB = NW(J)
C        IF(LB.GT.1.AND.LB.NE.4.AND.LB.NE.9) GO TO 68
            IF (LB.LE.1) GO TO 20
            IF (LB.EQ.4) GO TO 20
            IF (LB.NE.9) GO TO 100
   20       MA = NU(J)
            N = KAN(J)
            AMU = AMW
            BMU = ALPHA(J,1)
            M = NP(J)
c           write(6,*)"M",M
C           K1K2=(MA*(MA-1))/2+KA
c           K1K2 = KAKA + MA
            R = RAD(KA,MA)
            MN = MN + 1
C           AM1, PM3, MNDO
            IF(ISW1.GT.50.AND.ISW1.LT.54) THEN
c             write(6,*)"before f0gam"
              GAMMA(MN) =  F0GAM(N,L,R,PO)
c             write(6,*)"after F0GAM",GAMMA(MN)
              GO TO 90
            ENDIF
C
            IF (L.EQ.47 .AND. N.EQ.47) ISPEC = 0
            IF (ISPEC.LE.0 .OR. ISPEC.EQ.2) GO TO 50
c           IF (L.EQ.47 .OR. N.EQ.47) GO TO 50
c           IF (JSW(1).EQ.1) GO TO 40
c           IF (ISPEC.EQ.4) GO TO 40
c           IF (ISPEC.EQ.3) GO TO 30
C
C   For ISPEC = 1  use Mataga-Nishimoto gammas (spectroscopic theory)
C             = 5  use modified Warshel gammas
C
            IF (L.EQ.0 .OR. N.EQ.0) GO TO 50
c           IF (NEXTR(I).LT.0 .OR. NEXTR(J).LT.0) GO TO 50
C
C        IF(NEXTR(I).EQ.1.OR.NEXTR(J).EQ.1) GOTO 63
C        IF EITHER CENTRE N OR L IS A RYDBERG CENTRE OR IF EITHER IS A
C        CORE ORBITAL, CALC THE TWO ELECTRON INTEGRAL THEORETICALLY
C
c           write(6,*)"before coul",L,N
            status = COUL(KA,MA,L,N,R,LA,LB,AMU,BMU,K,M,
     x                  0,0,FACT,ISPEC,INT)   
c           write(6,*)"coul INT",INT
            GAMMA(MN) = INT   
            GO TO 90
C
C   For ISPEC = 3   use Ohno-Klopman  gammas (CNDO, INDO, PPP only!)
C
   30       IF (L.EQ.0 .OR. N.EQ.0) GO TO 50
c           IF (NEXTR(I).LT.0 .OR. NEXTR(J).LT.0) GO TO 50
c           IF (NEXTR(I).EQ.1 .OR. NEXTR(J).EQ.1) GO TO 50
C
C        IF EITHER CENTRE N OR L IS A RYDBERG CENTRE CALC THE TWO ELECTR
C        INTEGRAL THEORETICALLY
c          write(6,*)"before coulok"
c
c           GAMMA(MN) = COULOK(KA,MA,L,N,R,LA,LB,AMU,BMU,K,M,NEXTR(I),
c    x                         NEXTR(J),FACT)
c           GO TO 90
C
C  For JSW(1) = 1 or ISPEC = 4, use Pariser-Parr gammas
C
   40       CONTINUE
            IF (N.EQ.0 .OR. L.EQ.0) GO TO 50
c           IF (NEXTR(I).LT.0 .OR. NEXTR(J).LT.0) GO TO 50
c           IF (NEXTR(I).EQ.1 .OR. NEXTR(J).EQ.1) GO TO 50
c            GAMMA(MN) = COULP(L,N,R)
c           GO TO 90

   50       CONTINUE
c           write(6,*)"for Ag"
C
            IF (K-M) 80, 70, 60
   60       MMM = M + (K* (K-1))/2
c           write(6,*)"MMM,BMU,AMU,R,INT"
c           write(6,*)MMM,BMU,AMU,R,INT
            status = COULB(MMM,BMU,AMU,R,FACT,INT,NIN)
c           write(6,*)"coulb INT",INT
            GAMMA(MN)=INT
c           write(6,*)"GAMMA",GAMMA(MN)
            GO TO 90

   70       IF (AMU.LT.BMU) GO TO 60
   80       MMM = K + (M* (M-1))/2
c           write(6,*)"MMM,BMU,AMU,R,INT"
c           write(6,*)MMM,BMU,AMU,R,INT,NIN(1)
            status = COULB(MMM,AMU,BMU,R,FACT,INT,NIN)
c           write(6,*)"MN",MN
c           write(6,*)"coulb INT 1",INT
            GAMMA(MN)=INT
c           write(6,*)"GAMMA",GAMMA(MN)
   90       CONTINUE
  100    CONTINUE
  110     NG(I) = KG
  120 CONTINUE
  130 CONTINUE
c      DO 140 I = 1, NA
c        NTYP(I) = IABS(NTYP(I))
  140 CONTINUE
C
      RETURN

      END
C----------------------------------------------------------------------------
C
C----------------------------------------------------------------------------
      FUNCTION F0GAM (NI,NJ,R,PO)
C     SPERHICAL PART OF AM1 INTEGRALS
      double precision R, RAU, R2,EE, AEE
      double precision PO(9,105)
      double precision F0GAM
      parameter (ONE = 1.0D0)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/am1.cmn'
C     R ENTERS IN ANGS
C     *
C *** INITIALIZATION.
c     write(6,*)"i am in fogam"
      RAU = R!/BOHRS
      R2     = RAU*RAU
      AEE    = (PO(1,NI)+PO(1,NJ))**2
C ***
         EE     = ONE/SQRT(R2+AEE)
      F0GAM = EE
C ***
      write(6,*)"EE",EE
      RETURN
      END
