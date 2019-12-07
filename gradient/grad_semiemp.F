      subroutine grad_semiemp(NA,NB,ALPHA,NP,NW,KAN,NU,RAD,FA,NIN,F,H
     1 ,BETA,BETB,ISW1,BKA,BKB,NRR,NG,RA,RB,Z12,NTYP,ISWANN,CO,ETS,ETP,
     1 ETD,ETF,FZ,GRAD,ISPEC)
C
C     DIFOVL NOW CALCULATES THE OVERLAP AND ITS DERIVATIVE
C
C     SUBROUTINE DELEN ORGANIZES THE DERIVATIVES OF INTEGRALS.
C     SUBROUTINE FORCES EVALUATES THE FORCE AT EACH NUCLEI
C
C     ***** D2LEN IS AN ENTRY FOR QM/MM THAT BYBASSES THE 
C           DYNAMIC ALLOCATION.  IT IS CALLED FROM quantum.f
C           AND IS NEARLY IDENTICAL TO D1LEN.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "sym.fh"
#include "util.fh"
#include "dftps.fh"

C
C  SET UP MASK TO AVOID REPEATED INTEGRAL CALCULATIONS
C  SEE ONEINT FOR AN EXPLANATION
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
C     CONST.CMN
      PARAMETER (BOHRS    = 5.29177249D-1,
     .           EVOLT    = 2.72113961D1,
     .           ENUBAR   = 2.1947463067D5,
     .           DEBYE    = 2.541748047D0)
C
      INTEGER      KFLAG(3), MFLAG(3)
C     DIMENSION EIS(MPE),EIP(MPE),EID(MPE),EIF(MPE),BETA(MPE,10),
C    . BETB(MPE,4)
      DIMENSION     DPEN(NB,NB), WA(NRR), NG(*), NTYP(*),GRAD(*)
c     DIMENSION    DG(*), RAD(*), FA(*), H(*), F(*), FZ(*), FALPH(*)
      DIMENSION    NP(*), NW(*), KAN(*), NU(*), RAD(NA,NA), NIN(*), H(*)
      DIMENSION    FA(*), ALPHA(NB,6), FKAPPA(7,7), F(*), FACBLK(5,105)
      DIMENSION    DT(10), DP(10), ET(10), EP(10), D(10), E(10)
      DIMENSION    FSPDF(4,4),QQ(3*NA), FZ(3*NA),IEX(NB)
      DIMENSION    BETA(105,10),BETB(105,4),RA(*),RB(*),Z12(*)
      DIMENSION    FACTZ(10), IBB(5), BETEXP(12,105),FACT(35)
      DIMENSION    IZDSIG(136), IZDPI(120), IZDDEL(78), IZDPHI(28)
      DIMENSION    BKA(NB), BKB(NB), DG(NRR),CO(*)
      DIMENSION    ETS(*),ETP(*),ETD(*),ETF(*)
c     COMMON /NEWGEOM/ HZIJ,SIGG,SS2,PHIG,SF2,DELG,SD2,PIG,SP2,KA
C     DIMENSION    p1(*),pg(*),ph(*),veco(l1,nbond),vecv(l1,nvirt),
C    .             terms(nbond,2),xocc(*),lstocc(*)
C
C  THE ARRAYS DT,DP,ET,EP COULD BE EQUIVALENCED WITH QQ
      DATA         MONE/-1.0D0/
      DATA         IBB/0, 1, 3, 6, 10/
      DATA         IZDSIG/0, 0, 0, -1, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 2, 0, 0, -2, 2, 0, 0, 0, 0, 0, 0, 1, 0, 0,
     .             0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3,
     .             0, 0, 3, -3, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/
      DATA         IZDPI/0, 0, 0, 0, 4, 0, 0, -5, 0, 0, 0, -4, 0, 0, 0,
     .             -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -4, 5,
     .             4, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0,
     .             0, -1, 0, 0, 0, 0, 0, 0, 4, -5, -4, 0, 0, 0, 6, 0, 0,
     .             0, 7, 0, 0, 0, 0, 0, -7, 0, 0, 7, 0, 8, 0, 0, 0, 0,
     .             0, -8, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0/
      DATA         IZDDEL/0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 1, 0, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, -7, 0, 0, 0, 0, 0, -8, 0, 0,
     .             0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 8, 0, 0, 0, 0, -4,
     .             0, 0, 7, 0, 0, 0, 0, 0, -9, 0, 0, 0, 0, 0, 0, 9, 0,
     .             0, -10, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0/
      DATA         IZDPHI/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     .             0, 0, 0, 0, 9, 0, 0, 0, 0, 9, 0, 0, 0/
      DATA         FACTZ/1.0000000000D0, 1.7320508070D0, 2.4494897440D0,
     .             -1.0000000000D0, -1.7320508070D0, -2.4494897440D0,
     .             -1.5811388300D0, 1.5811388300D0, -1.2247448720D0,
     .             1.2247448720D0/
      data FACT/1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0
     X          ,3628800.0,39916800.0,479001600.0,6227020800.0,
     X          87178291200.0,1307674368000.0,20922789888000.0,
     X          355687428096000.0,6402373705728000.0,
     X          1.21645100408832000E+017,2.43290200817664000E+018,
     X          5.10909421717094400E+019,1.12400072777760768E+021,
     X          2.58520167388849782E+022,6.20448401733239410E+023,
     X          1.55112100433309861E+025,4.03291461126605650E+026,
     X          1.08888694504183519E+028,3.04888344611713837E+029,
     X          8.84176199373970077E+030,2.65252859812191032E+032,
     X          8.22283865417792243E+033, 2.63130836933693518E+035,
     X          8.68331761881188594E+036,2.95232799039604120E+038/

C     GET THE OVERLAP
C      REWIND 16
C      READ  (16) (H(I),I=1,NRR)
C
c     write(6,*)"i am in grad_semiemp"!,(F(I),I=1,21)
      FMBET = ZERO
      SWGT = ZERO
      SS2 = ZERO
      SS = ZERO 
      SP2 = ZERO
      SP =ZERO
      SIGG = ZERO
      SF2 =ZERO
      SF =ZERO
      SD2 = ZERO
      SD =ZERO
      PIG =ZERO
      PHIG = ZERO
      PHWGT = ZERO
      MM = IZERO
      IB1 = IZERO
      IB2 = IZERO
      PWGT = ZERO
      GGG =ZERO
      EPI = ZERO
      EDEL =ZERO      
      DELG = ZERO
      DWGT = ZERO
      DSWGT = zERO
      DSIGG = ZERO
      DPIG = ZERO
      DPWGT = ZERO
      DRAD = ZERO
      DPHIG = ZERO
      DPHWGT = ZERO
      DGGG = ZERO
      DDWGT = ZERO
      DDELG = ZERO
      BETAS = ZERO
      BETAP = ZERO
      BETAD = ZERO
      BETAF = ZERO
      S2T = ZERO
      C2P = ZERO
      C2T = ZERO




      JBETA = 0
      do i =1, 7
       do j = 1, 7
         FKAPPA(i,j)=0.0d0
       end do
      end do
      do i = 1, nb
        IEX(I) = 0
      end do
C     IF JBETA = 1 then form is a +b/R
      IF(IBETA.EQ.2 .OR. IBETA.EQ.3) JBETA = 1
      ITIME = 1
      IATMKX = 0
      IATMMX = 0
C 
C
C  INITIALISE PP TO -10.D0 TO SET UP DIFOVL
C
      HIJ = 0.0D0
      HALPIJ = 0.0d0
      HZIJ = 0.0D0
      PP = -10.D0
C
C
C     THIS PORTION OF THE PROGRAM IS CONCERNED WITH SETTING UP VARIOUS
C     GEOMETRICAL MATRICES
C
C     GG VALUE DEFINED IN ONEINT
      GG = FOUR
C
      DO 20 I = 1, 4
         DO 10 J = 1, 4
            FSPDF(I,J) = ONE
   10    CONTINUE
   20 CONTINUE
      if (ISPEC.eq.1) then
       FSPDF(1,1) = 1.0
       FSPDF(2,1) = 1.267
       FSPDF(2,2) = 0.585
       FSPDF(3,1) = 1.0
       FSPDF(3,2) = 1.0
       FSPDF(3,3) = 1.0
      end if

c     FSPDF(1,1) = FSSI1
c     FSPDF(2,1) = FPSI1
c     FSPDF(2,2) = FPPI1
c     FSPDF(3,1) = FDSI1
c     FSPDF(3,2) = FDPI1
c     FSPDF(3,3) = FDDE1
      NA1 = NA*NA
      NA2 = 2*NA1
   30 CONTINUE
C
C
c     CALL GEOMAT(FA,RAD)
C
C  ZERO THE FORCE ARRAYS
C  FORCES FROM THE OVERLAP AND RESONANCE INTEGRALS ARE
C  INITIALLY SET UP IN Q NUM
C
      NA11 = NA
      NA12 = NA + NA
      DO 40 I = 1, NA
         QQ(I) = ZERO
         QQ(I+NA11) = ZERO
         QQ(I+NA12) = ZERO
         FZ(I) = ZERO
         FZ(I+NA11) = ZERO
         FZ(I+NA12) = ZERO
   40 CONTINUE
C
C  H = DS/DR
C  FALPH= DS/DTHETA
C
C  TSIGN SETS THE APROPRIATE SIGN FOR DS/DTHETA
C  TSIGN POSITIVE WHEN GOING FROM ATOM 1 TO 2
C  TSIGN NEGATIVE FOR 2 TO 1
C  I.E. DS/DTHETA AN ANTISYMMETRIC ARRAY
C
C  FZ= DS/DPHI  ALL IN THE LOCAL COORDINATE FRAMEWORK
C
C
C     EVALUATES THE DERIVATIVE OF THE OVERLAP
C
C  AND INCLUDES THE RESONANCE INTEGRAL TERM
C
C     INITIALIZE THE PARAMETERS
C
      IJ = 0
C     GET THE TOTAL DENSITY AND CALL IT F
C     IF(JSW(19).GE.1)THEN
C       CALL gvbdnz(f,p1,pg,ph,veco,vecv,terms,xocc,lstocc,nocc1,nbond,
C    . nvirt,ngvb,l1)
C     ELSE

c       CALL FINDP(F,NRR,ISWANN,NDN,IONE)

C     ENDIF
C
C     DERIVATIVE OF OVERLAP TERMS
C
      DO 520 I = 1, NB
         K = NW(I)
         KA = NU(I)
         KFLAG(1) = KA
         IZ = KAN(I)
         NT1 = NP(I)
         IF(K-1) 50, 60, 60
   50    LI = NP(I)
         NL = 0
         GO TO 110

   60    IF(K-3) 70, 70, 80
   70    LI = NP(I) + 3
         NL = 1
         GO TO 110

   80    IF(K-8) 90, 90, 100
   90    LI = 8
         NL = 2
         GO TO 110

  100    LI = 15
         NL = 3
  110    KFLAG(3) = NL
C
  120    CONTINUE
         DO 510 J = 1, I
            IJ = IJ + 1
            M = NW(J)
            MA = NU(J)
            MFLAG(1) = MA
            NT2 = NP(J)
            JZ = KAN(J)
            IF(MA.EQ.KA) GO TO 510
            KAMA = NIN(KA) + MA
            R = RAD(KA,MA)*BOHRS
            RR = R/BOHRS
            IF(K-M) 130, 130, 140
C  ATOM KA TO MA
  130       LA = M
            LH = K
            K1M = KA + (MA-1)*NA
            TSIGN = -ONE
            GO TO 150
C  ATOM MA TO KA
  140       LA = K
            LH = M
            K1M = MA + (KA-1)*NA
            TSIGN = ONE
  150       CT = FA(K1M)
            COSP = FA(K1M+NA1)
            SINP = FA(K1M+NA2)
            CALL transmat(CT,COSP,SINP,LH,D)
            CALL transmat(CT,COSP,SINP,LA,E)
C  GET DIFFERENTIALS OF THE COORDINATE TRANSFORMATIONS
            CALL DGEOM(CT,COSP,SINP,LH,DT,DP)
            CALL DGEOM(CT,COSP,SINP,LA,ET,EP)
            IF(M-1) 160, 170, 170
  160       LJ = NP(J)
            ML = 0
            GO TO 220

  170       IF(M-3) 180, 180, 190
  180       LJ = NP(J) + 3
            ML = 1
            GO TO 220

  190       IF(M-8) 200, 200, 210
  200       LJ = 8
            ML = 2
            GO TO 220

  210       LJ = 15
            ML = 3
  220       MFLAG(3) = ML
C
C  GET BETA PAIR
C
            IF(IBETA.LT.4) THEN
               CALL HOFF(I,J,X,BBA,BBB,RAD,H,IEX,JBETA,ITIME,NW,NU,KAN,
     1               BETA,BETB,NIN,NA,NP,ISW1,ALPHA,NB,BKA,BKB)
C              X = B + C/R
               SIGG = X
               PIG = X
               DELG = X
               PHIG = X
               BBQ = ZERO
            ELSE IF(IBETA.EQ.5) THEN
c              CALL HOFFP(I,J,X,BBA,BBB,RAD,IEX,JBETA,
c    .                    ITWO)
               SIGG = X
               PIG = X
               DELG = X
               PHIG = X
               BBQ = ZERO
            ELSE IF(IBETA.EQ.6) THEN
               BBQ= FKAPPA(NT1,NT2)*(ALPHA(I,1)-ALPHA(J,1))**2
     x             /EVOLT
              IB1 = IBB(NL+1)
              IB2 = IBB(ML+1)
              BETAS = (BETA(IZ,IB1+1)+BETA(JZ,IB2+1))/ (TWO*EVOLT) +BBQ
              BETAP = (BETA(IZ,IB1+2)+BETA(JZ,IB2+2))/ (TWO*EVOLT) +BBQ
              BETAD = (BETA(IZ,IB1+3)+BETA(JZ,IB2+3))/ (TWO*EVOLT) +BBQ
              BETAF = (BETA(IZ,IB1+4)+BETA(JZ,IB2+4))/ (TWO*EVOLT) +BBQ
              SIGG = RR**FMBET*EXP(- (BETEXP(IB1+1,IZ)+BETEXP(IB2+1,
     .               JZ))*RR*RR)*BETAS
              PIG = RR**FMBET*EXP(-(BETEXP(IB1+2,IZ)+BETEXP(IB2+2,JZ))
     .              *RR*RR)*BETAP
              DELG = RR**FMBET*EXP(- (BETEXP(IB1+3,IZ)+BETEXP(IB2+3,
     .               JZ))*RR*RR)*BETAD
              PHIG = RR**FMBET*EXP(- (BETEXP(IB1+4,IZ)+BETEXP(IB2+4,
     .               JZ))*RR*RR)*BETAF
              BBB = ZERO
            ELSE
             WRITE(6,2266)  IBETA
 2266        FORMAT(' ***** ERROR IN OPTIONS, NO IBETA =',I5)
             CALL ZEXIT(3)
            END IF
            IF(LA.EQ.M) MM = ML + 1
            IF(LA.EQ.K) MM = NL + 1
            GO TO (260,240,250,230), MM
C  GET PHASES FOR AO TRANSFORMATIONS
C  F-ORBS
  230       EPHI = ONE
            EDEL = MONE
C  P-ORBS
  240       EPI = ONE
            ESIG = MONE
            GO TO 270
C  D-ORBS
  250       EDEL = ONE
            EPI = MONE
C  S-ORBS
C
  260       ESIG = ONE
  270       CONTINUE
            IF(R.LT.1.D-4) GO TO 510
            KFLAG(2) = NT1
            MFLAG(2) = NT2
            IATMK = KFLAG(1) + 10*KFLAG(2) + 100*KFLAG(3)
            IATMM = MFLAG(1) + 10*MFLAG(2) + 100*MFLAG(3)
            IF(IATMK.EQ.IATMKX .AND. IATMM.EQ.IATMMX) GO TO 300
            IATMKX = IATMK
            IATMMX = IATMM
            SS = ZERO
            SP = ZERO
            SD = ZERO
            SF = ZERO
            SS2 = ZERO
            SP2 = ZERO
            SD2 = ZERO
            SF2 = ZERO
            DO 290 I1 = 1, 2
               C1 = ALPHA(I,I1+4)
               AMU = ALPHA(I,I1+1)
               IF(AMU.LT.1.0D-4) GO TO 290
               DO 280 J1 = 1, 2
                  C2 = ALPHA(J,J1+4)
                  BMU = ALPHA(J,J1+1)
                  IF(BMU.LT.1.0D-4) GO TO 280
                  CALL DIFOVL(LH,NT1,NL,AMU,NT2,ML,BMU,R,TSS2,TSP2,TSD2,
     .                        TSF2,TSS,TSP,TSD,TSF,FACT,NIN,PP)
c      write(6,*)"TSS,TSP,TSD,TSF,TSS2,TSP2,TSD2,TSF2,C1,C2",TSS,TSP,TSD
c    1     ,TSF,TSS2,TSP2,TSD2,TSF2,C1,C2

C
C  GET LOCAL OVERLAP FOR ANGLE VARIATION AND D-ORBITALS
C
C      CALL OVLAP(LH,NT1,NL,AMU,NT2,ML,BMU,R,TSS2,TSP2,TSD2,TSF2,FACT)
                  C3 = C1*C2
                  SS = SS + TSS*C3
                  SP = SP + TSP*C3
                  SD = SD + TSD*C3
                  SF = SF + TSF*C3
                  SS2 = SS2 + TSS2*C3
                  SP2 = SP2 + TSP2*C3
                  SD2 = SD2 + TSD2*C3
                  SF2 = SF2 + TSF2*C3
C     SS,SP ETC. IS D(S)/D(R), SS2, SP2, ETC IS THE OVERLAP S.
  280          CONTINUE
  290       CONTINUE
C ----------------------------------------------------------------
C
c  INCLUDE DISTANCE DEPENDENT TERMS OF THE WEIGTING FACTORS
C       JAN 17 83
C
            GGG = EXP(- (ALPHA(I,1)+ALPHA(J,1))*R*GG)
            DGGG = - (ALPHA(I,1)+ALPHA(J,1))*GG*BOHRS*GGG
  300       CONTINUE
C
C
C
C
C -------------------------------------------------------------
C  ------------------------------------------------
            IF(IBETA.EQ.6) THEN
               DGGG = ZERO
               DRAD = ZERO
               DSIGG = (FMBET*RR** (FMBET-ONE)-
     .                 TWO*RR** (FMBET+ONE)* (BETEXP(IB1+1,
     .                 IZ)+BETEXP(IB2+1,JZ)))*EXP(- (BETEXP(IB1+1,
     .                 IZ)+BETEXP(IB2+1,JZ))*RR*RR)
               DSIGG = DSIGG * BETAS
               DPIG = (FMBET*RR** (FMBET-ONE)-
     .              TWO*RR**(FMBET+ONE)*(BETEXP(IB1+2,IZ)+BETEXP(IB2+2,
     .              JZ)))*EXP(- (BETEXP(IB1+2,IZ)+BETEXP(IB2+2,JZ))*
     .              RR*RR)*EPI
               DPIG = DPIG * BETAP
               DDELG = (FMBET*RR** (FMBET-ONE)-
     .                 TWO*RR** (FMBET+ONE)* (BETEXP(IB1+3,
     .                 IZ)+BETEXP(IB2+3,JZ)))*EXP(- (BETEXP(IB1+3,
     .                 IZ)+BETEXP(IB2+3,JZ))*RR*RR)*EDEL
               DDELG = DDELG*BETAD
               DPHIG = (FMBET*RR** (FMBET-ONE)-
     .                 TWO*RR** (FMBET+ONE)* (BETEXP(IB1+4,
     .                 IZ)+BETEXP(IB2+4,JZ)))*EXP(- (BETEXP(IB1+4,
     .                 IZ)+BETEXP(IB2+4,JZ))*RR*RR)*EPHI
               DPHIG = DPHIG*BETAF
            ELSE IF(IBETA.EQ.2.OR.IBETA.EQ.3) THEN
C              INCLUDE THE d(1/r)/dr term
C
C              BBB = B/R
               DRAD = -BBB/RR
               DSIGG = DRAD
               DPIG = DRAD
               DDELG = DRAD
               DPHIG = DRAD
            ELSE IF(IBETA.EQ.1.OR.IBETA.EQ.0) THEN
               DRAD = ZERO
               DSIGG = ZERO
               DPIG = ZERO
               DDELG = ZERO
               DPHIG = ZERO
C
            ENDIF
            IF(IBETA.LT.6) THEN
              IF(ABS(FACBLK(1,1)).LT.1.0D-8) THEN

                  IF(NL.EQ.ML) GO TO 310
                  IF(ISW1.EQ.0) GO TO 310
                  SWGT = ONE
                  PWGT = ONE
                  DWGT = ONE
                  PHWGT = ONE
                  DSWGT = ZERO
                  DPWGT = ZERO
                  DDWGT = ZERO
                  DPHWGT = ZERO
                  GO TO 320
  310             CONTINUE
                  NLP1 = NL + 1
                  MLP1 = ML + 1
                  FSWGT = (FSPDF(NLP1,1)+FSPDF(MLP1,1))/TWO
                  FPWGT=  (FSPDF(NLP1,2)+FSPDF(MLP1,2))/TWO
                  FDWGT = (FSPDF(NLP1,3)+FSPDF(MLP1,3))/TWO
                  FPHWGT = (FSPDF(NLP1,4)+FSPDF(MLP1,4))/TWO
                  SWGT=FSWGT+(ONE-FSWGT)*GGG
                  PWGT=FPWGT+(ONE-FPWGT)*GGG
                  DWGT=FDWGT+(ONE-FDWGT)*GGG
                  PHWGT=FPHWGT+(ONE-FPHWGT)*GGG
                  DSWGT = (ONE-FSWGT)*DGGG
                  DPWGT = (ONE-FPWGT)*DGGG
                  DDWGT = (ONE-FDWGT)*DGGG
                  DPHWGT = (ONE-FPHWGT)*DGGG
  320             CONTINUE

C  10             SIGMA FACTORS, THEN 10 PI, THEN 10 DELTA, THEN 10 PHI.
C ---------------------------------------------------------------
C                 IF(NLP1.GT.MLP1) THEN
C                    KLML = NIN(NLP1) + MLP1
C                 ELSE
C                    KLML = NIN(MLP1) + NLP1
C                 END IF
C                 SIGG = FACBLK(1,KLML)
C                 PIG = FACBLK(2,KLML)
C                 DELG = FACBLK(3,KLML)
C                 PHIG = FACBLK(4,KLML)
C                 DSIGG = DRAD
C                 DPIG = DRAD
C                 DDELG = DRAD
C                 DPHIG = DRAD
C
C --------END FIRST CHOICE---------------------------------------
C ---------------------------------------------------------------
C                 P-SIGMA, P-PI, D-DELTA, S-PSIG FOR FIRST 10 ATOMS.
C
C                 IZ=KAN(KA)
C                 JZ=KAN(MA)
C                 IF(IZ.GT.10.OR.JZ.GT.10) THEN
C                   WRITE(6,*) '  ***** ERROR IN FACBLK, SUB ONEINT****'
C                   CALL ZEXIT(3)
C                 ENDIF
C                 IF(NL.EQ.ML) GO TO 45
C                 IF(ISW1.EQ.0) GO TO 45
C                 SIGG=SQRT(FACBLK(4,IZ)*FACBLK(4,JZ))
C                 PIG=ONE
C                 DELG=ONE
C                 PHIG=ONE
C                 GO TO 455
C  45             CONTINUE
C --------SQRT ROOT FORMULAS ------------------------------------
C                 SIGG=SQRT(FACBLK(1,JZ)*FACBLK(1,IZ))
C                 PIG= SQRT(FACBLK(2,JZ)*FACBLK(2,IZ))
C                 DELG=SQRT(FACBLK(3,JZ)*FACBLK(3,IZ))
C                 PHIG=SQRT(FACBLK(4,JZ)*FACBLK(4,IZ))
C---------SUM FORMULAS ------------------------------------------
C                 SIGG=(FACBLK(1,IZ)+FACBLK(1,JZ))/TWO
C                 PIG= (FACBLK(2,IZ)+FACBLK(2,JZ))/TWO
C                 DELG=(FACBLK(3,IZ)+FACBLK(3,JZ))/TWO
C                 PHIG=(FACBLK(4,IZ)+FACBLK(4,JZ))/TWO
C --------END COMBINATION FORMULAS -----------------------------
C 455             CONTINUE
C                 WRITE(6,*) ' IZ,JZ,SIG,PIG', IZ,JZ,SIGG,PIG
               ENDIF
C ---------------------------------------------------------
C -------------------------------------------------------------
           ELSE
C            IBETA = 6
             SWGT = ONE
             PWGT = ONE
             DWGT = ONE
             PHWGT = ONE
             DSWGT = ZERO
             DPWGT = ZERO
             DDWGT = ZERO
             DPHWGT = ZERO
         ENDIF
C
C
C FORM SIGMA INTERACTION
C
            DSIGG = (DRAD*SWGT+SWGT*DSIGG + SIGG*DSWGT)*ESIG
c           write(6,*)"DSIGG",DSIGG
            SIGG = SIGG*SWGT*ESIG
c           write(6,*)"SIGG",SIGG
            SIJ = E(1)*D(1)*SS
c           write(6,*)"SIJ",SIJ
            SALPIJ = (ET(1)*D(1)+DT(1)*E(1))*SS2
c           write(6,*)"SALPIJ",SALPIJ
            SZIJ = (EP(1)*D(1)+DP(1)*E(1))*SS2
c           write(6,*)"SZIJ",SZIJ
            HIJ = E(1)*D(1)*SS*SIGG+E(1)*D(1)*SS2*DSIGG
c           write(6,*)"HIJ",HIJ
            HALPIJ = (ET(1)*D(1)+DT(1)*E(1))*SS2*SIGG
c           write(6,*)"HALPIJ",HALPIJ
            HZIJ = (EP(1)*D(1)+DP(1)*E(1))*SS2*SIGG
c           write(6,*)"HZIJ",HZIJ
C              HIJP = E(1)*D(1)*SS2*SIGG*DRAD
C
C  FIND MAX LZ VALUE USED IN INTEGRAL
C
            MINM = MIN0(NL,ML) + 1
            GO TO (360,350,340,330), MINM
C PHI INTERACTION
  330       CONTINUE
c           write(6,*)"here 1"
            DPHIG = (DRAD*PHWGT + PHWGT*DPHIG +  PHIG*DPHWGT)*EPHI
            PHIG = PHWGT*PHIG*EPHI
            HIJ = HIJ + (E(7)*D(7)+E(6)*D(6))*SF*PHIG +
     .            (E(7)*D(7)+E(6)*D(6))*SF2*DPHIG
            HALPIJ = HALPIJ + (ET(7)*D(7)+ET(6)*D(6)+DT(7)*E(7)+
     .               DT(6)*E(6))*SF2*PHIG
            HZIJ = HZIJ + (EP(7)*D(7)+EP(6)*D(6)+DP(7)*E(7)+DP(6)*E(6))*
     .             SF2*PHIG
            SIJ = SIJ + (E(7)*D(7)+E(6)*D(6))*SF
            SALPIJ = SALPIJ + (ET(7)*D(7)+ET(6)*D(6)+DT(7)*E(7)+
     .               DT(6)*E(6))*SF2
            SZIJ = SZIJ + (EP(7)*D(7)+EP(6)*D(6)+DP(7)*E(7)+DP(6)*E(6))*
     .             SF2
C              HIJP = HIJP + (E(7)*D(7)+E(6)*D(6))*SF2*PHIG*DRAD
C DELTA INTERACTION
  340       CONTINUE
c           write(6,*)"here 2"
            DDELG =( DRAD*DWGT + DWGT*DDELG +  DDWGT*DELG)*EDEL
            DELG = DWGT*DELG*EDEL
            HIJ = HIJ + (E(4)*D(4)+E(5)*D(5))*SD*DELG +
     .            (E(4)*D(4)+E(5)*D(5))*SD2*DDELG+DRAD
            HALPIJ = HALPIJ + (ET(4)*D(4)+ET(5)*D(5)+DT(4)*E(4)+
     .               DT(5)*E(5))*SD2*DELG
            HZIJ = HZIJ + (EP(4)*D(4)+EP(5)*D(5)+DP(4)*E(4)+DP(5)*E(5))*
     .             SD2*DELG
            SIJ = SIJ + (E(4)*D(4)+E(5)*D(5))*SD
            SALPIJ = SALPIJ + (ET(4)*D(4)+ET(5)*D(5)+DT(4)*E(4)+
     .               DT(5)*E(5))*SD2
            SZIJ = SZIJ + (EP(4)*D(4)+EP(5)*D(5)+DP(4)*E(4)+DP(5)*E(5))*
     .             SD2
C              HIJP = HIJP + (E(4)*D(4)+E(5)*D(5))*SD2*DELG*DRAD
C  PI INTERACTION
  350       CONTINUE
c           write(6,*)"here 3"
            DPIG = (DRAD*PWGT + PWGT*DPIG + PIG*DPWGT)*EPI
            PIG = PWGT*PIG*EPI
            HIJ = HIJ + (E(2)*D(2)+E(3)*D(3))*SP*PIG +
     .            (E(2)*D(2)+E(3)*D(3))*SP2*(DPIG)
            HALPIJ = HALPIJ + (ET(2)*D(2)+ET(3)*D(3)+DT(2)*E(2)+
     .               DT(3)*E(3))*SP2*PIG
            HZIJ = HZIJ + (EP(2)*D(2)+EP(3)*D(3)+DP(2)*E(2)+DP(3)*E(3))*
     .             SP2*PIG
            SIJ = SIJ + (E(2)*D(2)+E(3)*D(3))*SP
            SALPIJ = SALPIJ + (ET(2)*D(2)+ET(3)*D(3)+DT(2)*E(2)+
     .               DT(3)*E(3))*SP2
            SZIJ = SZIJ + (EP(2)*D(2)+EP(3)*D(3)+DP(2)*E(2)+DP(3)*E(3))*
     .             SP2
C              HIJP = HIJP + (E(2)*D(2)+E(3)*D(3))*SP2*PIG*DRAD
C
C  SCALE THE OVERLAP BY THE RESONANCE PARAMETER
C
  360       CONTINUE
c           write(6,*)"here 4",TSIGN,HALPIJ
C           IF(IBETA.NE.6) THEN
C              IF(JBETA.EQ.0) THEN
C                 HIJ = HIJ*BBA
C                 HALPIJ = HALPIJ*BBA
C                 HZIJ = HZIJ*BBA
C              ELSE
C                 HIJ = HIJ* (BBA+BBB) + HIJP
C                 HALPIJ = HALPIJ* (BBA+BBB)
C                 HZIJ = HZIJ* (BBA+BBB)
C              END IF
C           END IF
C  SET THIS PHASE
            HALPIJ = TSIGN*HALPIJ
C
C  FORM GRADIENT TERMS
C
c           write(6,*)"F(IJ),HIJ",F(IJ),HIJ,HALPIJ
            FAC = TWO*F(IJ)*HIJ
            FACTH = -TWO*F(IJ)*HALPIJ
            ST = ONE - CT*CT
            IF(ABS(ST).GT.1.D-7) THEN
               ST = SQRT(ST)
            ELSE
C
               YSUM = ZERO
C
               IND = (LA+1)*LA/2 + LH + 1
               LABEL = IZDSIG(IND)
                  BKC = ONE
               IF(LABEL) 370, 390, 380
C
  370          YSUM = YSUM + FACTZ(-LABEL)*COSP*SS2*ESIG*SIGG*BKC
               GO TO 390
C
  380          YSUM = YSUM + FACTZ(LABEL)*COSP*CT*SS2*ESIG*SIGG*BKC
C
  390          GO TO (490,460,430,400), MINM
C
C        PHI INTERACTION
C
C        IND=(LA-8)*(LA-7)/2+LH-8
  400          IND = (LA-5)* (LA-4)/2 + LH - 4
               LABEL = IZDPHI(IND)
                  BKC = ONE
               IF(LABEL) 410, 430, 420
C
  410          YSUM = YSUM + FACTZ(-LABEL)*COSP*SF2*EPHI*PHIG*BKC
               GO TO 430
C
  420          YSUM = YSUM + FACTZ(LABEL)*COSP*CT*SF2*EPHI*PHIG*BKC
C
C        DEL INTERACTION
C
  430          IND = (LA-3)* (LA-2)/2 + LH - 2
               LABEL = IZDDEL(IND)
                  BKC = ONE
               IF(LABEL) 440, 460, 450
C
  440          YSUM = YSUM + FACTZ(-LABEL)*COSP*SD2*EDEL*DELG*BKC
               GO TO 460
C
  450          YSUM = YSUM + FACTZ(LABEL)*COSP*CT*SD2*EDEL*DELG*BKC
C
C        PI INTERACTION
C
  460          IND = (LA-1)*LA/2 + LH
               LABEL = IZDPI(IND)
                  BKC = ONE
               IF(LABEL) 470, 490, 480
C
  470          YSUM = YSUM + FACTZ(-LABEL)*COSP*SP2*EPI*PIG*BKC
               GO TO 490
C
  480          YSUM = YSUM + FACTZ(LABEL)*COSP*CT*SP2*EPI*PIG*BKC
C
C        SCALE BY BETA AND THE DENSITY
C
  490          CONTINUE
C              IF(IBETA.NE.6) YSUM = YSUM*X
               IF(TSIGN.EQ.ONE) THEN
                  COSP = -COSP
                  CT = -CT
               END IF
               YY = -COSP*YSUM/RR*TWO*F(IJ)
               XX = -FACTH*CT*COSP/RR
               GO TO 500

            END IF
            IF(TSIGN.EQ.ONE) THEN
               CT = -CT
               SINP = -SINP
               COSP = -COSP
            END IF
            FACP = TWO*F(IJ)*HZIJ
c           write(6,*)"FAC,ST,COSP,FACTH,CT,COSP,FACP,SINP,ST,RR",FAC,ST
c    1 ,COSP,FACTH,CT,COSP,FACP,SINP,ST,RR
            XX = -FAC*ST*COSP - (FACTH*CT*COSP-FACP*SINP/ST)/RR
            YY = -FAC*ST*SINP - (FACTH*CT*SINP+FACP*COSP/ST)/RR
  500       ZZ = -FAC*CT + FACTH*ST/RR
c           write(6,*)"XX,YY,ZZ",XX,YY,ZZ
c           write(6,*)"QQ(KA),QQ(KA+NA11),QQ(KA+NA12)",QQ(KA),QQ(KA+NA11
c    1     ),QQ(KA+NA12)
            QQ(KA) = QQ(KA) + XX
            QQ(KA+NA11) = QQ(KA+NA11) + YY
            QQ(KA+NA12) = QQ(KA+NA12) + ZZ
            QQ(MA) = QQ(MA) - XX
            QQ(MA+NA11) = QQ(MA+NA11) - YY
            QQ(MA+NA12) = QQ(MA+NA12) - ZZ
  510    CONTINUE
  520 CONTINUE
C
C
C  TRANSFER FORCES TO THE FZ ARRAY
      DO 530 I = 1, NA
         FZ(I) = QQ(I)
         FZ(I+NA11) = QQ(I+NA11)
         FZ(I+NA12) = QQ(I+NA12)
c     write(6,*)"FZ",FZ(I),FZ(I+NA11),FZ(I+NA12)
  530 CONTINUE
C     EVALUATES THE DERIVATIVE OF THE COULOMB MATRIX
C
C  INITIALISE PP TO -10.DO FOR DIFOVL
C
      PP = -10.D0
      KG = 0
C      FOR G.S. CNDO
      MN = 0
      DO 640 I = 1, NB
         LA = NW(I)
C      IF(LA.GT.1.AND.LA.NE.4) GO TO 69
         IF(LA.LE.1) GO TO 540
         IF(LA.EQ.4) GO TO 540
         IF(LA.NE.9) GO TO 630
  540    K = NP(I)
         KA = NU(I)
         L = KAN(I)
         KG = KG + 1
         AMW = ALPHA(I,1)
         DO 620 J = 1, I
            LB = NW(J)
            IF(LB.LE.1) GO TO 550
            IF(LB.EQ.4) GO TO 550
            IF(LB.NE.9) GO TO 620
  550       MA = NU(J)
            MN = MN + 1
            IF(MA.EQ.KA) THEN
               DG(MN) = ZERO
               GO TO 610

            END IF
            N = KAN(J)
            AMU = AMW
            BMU = ALPHA(J,1)
            M = NP(J)
            K1K2 = KA* (KA-1)/2 + MA
            R = RAD(KA,MA)*BOHRS
c          write(6,*)"R",R
            IF(R.LT.1.D-4) THEN
               DG(MN) = ZERO
               GO TO 610

            END IF
            IF(ISPEC.EQ.0) GO TO 570
c           IF(JSW(1).EQ.1) GO TO 560
C           MATAGA NISHIMOTO WEISS, OHNO KLOPMANN OR WARSHELL
            IF(N.EQ.0 .OR. L.EQ.0) GO TO 570
            DG(MN) = DCOUL(L,N,R,LA,LB,AMU,BMU,K,M,0,0,
     .               FACT,ISPEC)
c           write(6,*)"DG",DG(MN)
            GO TO 610
C           PARISER PARR
  560       CONTINUE
            IF(N.EQ.0 .OR. L.EQ.0) GO TO 570
c           DG(MN) = DCOULP(N,L,R)
c           write(6,*)"DG",DG(MN)
            GO TO 610

  570       CONTINUE
            IF(K-M) 600, 590, 580
  580       MMM = M + (K* (K-1))/2
            DG(MN) = DCOULB(MMM,BMU,AMU,R,FACT,CINTB,NIN)
c           write(6,*)"DG",DG(MN)
            GO TO 610

  590       IF(AMU.LT.BMU) GO TO 580
  600       MMM = K + (M* (M-1))/2
            DG(MN) = DCOULB(MMM,AMU,BMU,R,FACT,CINTB,NIN)
c           write(6,*)"DG",DG(MN)
  610       CONTINUE
            DPEN(NG(I),NG(J)) = DG(MN)
            DPEN(NG(J),NG(I)) = DG(MN)
  620    CONTINUE
  630    NG(I) = KG
  640 CONTINUE
C     DEBUG
C      CALL PRNT (DG,KG,MN)
C      CALL PRNT (FALPH,NB,NRR)
C     DEBUG
  650 CONTINUE
C     CALCULATED DERIVATIVE OF THE ELECTRON NUCLEAR TERM.
C     DPEN(II,JJ), II = AO, JJ = NUCLEUS (OR AO ON NUCLEUS)
      IDD1 = 0
      IF(IDD1.EQ.1) THEN
         DO 690 I = 1, NB
            LA = NW(I)
            IF(LA.LE.1) GO TO 660
            IF(LA.EQ.4) GO TO 660
            IF(LA.NE.9) GO TO 690
  660       NPI = NP(I)
            KA = NU(I)
            IZ = KAN(KA)
            NGI = NG(I)
            AMW = ALPHA(I,1)
            DO 680 J = 1, NB
               LB = NW(J)
               IF(LB.LE.1) GO TO 670
               IF(LB.EQ.4) GO TO 670
               IF(LB.NE.9) GO TO 680
  670          NPJ = NP(J)
               NGJ = NG(J)
               KB = NU(J)
               JZ = KAN(KB)
               IF(KB.EQ.KA) THEN
                  DPEN(NGI,NGJ) = ZERO
                  GO TO 680

               END IF
               AMU = AMW
               IF(KA.GT.KB) THEN
                  K1K2 = NIN(KA) + KB
               ELSE
                  K1K2 = NIN(KB) + KA
               END IF
               R = RAD(KA,KB)*BOHRS
               IF(R.LT.1.D-4) THEN
                  DPEN(NGI,NGJ) = ZERO
                  GO TO 680

               END IF
               XX = DPEN(NGI,NGJ)
               DPEN(NGI,NGJ) = DSPEN(I,R,CINTB,XX,NP,ALPHA,FACT,IONE,NB)
c            write(6,*)"DPEN",DPEN(NGI,NGJ)
  680       CONTINUE
  690    CONTINUE
      END IF
CEND DEBUG
c     IF(JSW(19).GE.1)THEN
c       CALL GVBFOR(DG,RAD,FA,H,F,FALPH,FZ,DPEN,p1,pg,ph,veco,vecv,
c    . terms,nbond,nvirt,ngvb,l1)
c       ELSE
c      write(6,*)"i am here"
       CALL DFORCES(NA,DG,RAD,FA,H,F,RA,FZ,DPEN,ISWANN,NIN,FACT,NTYP,KAN
     1 ,NP,NU,NW,NG,CO,Z12,ETS,ETP,ETD,ETF,NB,ALPHA,GRAD,ISPEC)
c     ENDIF

c      write(6,*)"in grad_semiemp"
      return
      end

















