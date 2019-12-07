      subroutine oneint_semiemp(rtdb,NA,NB,ALPHA,NP,NW,KAN,NU,RAD,F,FA,H
     1 ,SH,NIN,ISW1,inttyp,FT)

C    This subroutine construct one electron part of the hamiltonian.
c    It will construct overlap matrix and use that to get weigted
c    overlap matrix.    
C     H is the overlap, SH is the weihted overlap
C     FA is the kinetic energy

c
c
c
      implicit double precision (E,S,T,P,B,D,O)
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
      dimension D(8),E(8)
      integer rtdb
      integer NA
      integer NB
      integer LA,LH,NA1,NA2,KA,KB,IZA
      integer NT1,KL,KLP1,MA,MB,IZB,NT2,ML
      integer MLP1,MM
      integer NU(NB)
      double precision AMU,ACU,BMU,BCU,ABC,RR,R
      double precision U1,U2,U3
      integer  NP(NB)
      integer  NW(NB)
      integer inttyp
      double precision KAN(NB)
      integer i,j,ij,ISS,JSS
      double precision ALPHA(NB,6)
      double precision F(*)
      double precision RAD(NA,NA)
      double precision FA(*) 
      double precision SH(*)
      double precision H(*)
      integer NIN(*)
      logical status
      double precision PP,GG
      double precision FACBLK(5,105)
      integer IDD1
      double precision infa(6)      
      CHARACTER*9  ATOMK, ATOMKX, ATOMM, ATOMMX
      integer    LANGL(16)
      double precision     FSPDF(4,4)
      integer    KFLAG(3), MFLAG(3), IBB(5)
      parameter (BOHRS = 5.29177249D-1)
      DATA         LANGL/0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3/
      DATA FACBLK /525 * 0.0D0/
C     THE ARRAY IBB ALLOWS FOR A TRIANGULAR MATRIX OF FACTORS BETWEEN
C     DIFFERENT ANGULAR MOMENTA.
      parameter (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DATA         IBB/0, 1, 3, 6, 10/
      double precision FACT(35)
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

        MM = 0   
        EPI = ZERO
        EDEL = ZERO
C
c      if (ga_nodeid().eq.0)
c    . write(6,*) "constructing oneint"
C     INTRODUCE ATOMK AND ATOMM TO TEST IF TREATING THE SAME SHELL
C
C     KFLAG(1)  --> K atom number;
C     KFLAG(2)  --> K principal quantum number;
C     KFLAG(3)  --> K angular   quantum number;
C     MFLAG(1)  --> M atom number;
C     MFLAG(2)  --> M principal quantum number;
C     MFLAG(3)  --> M angular quantum number.
C
C     Initialize ATOMKX and ATOMMX to zero
C
      ATOMKX = '  0  0  0'
      ATOMMX = '  0  0  0'
C     FSPDF(I,J) is the weighting factor.
C     I = ANGULAR QUANTUM NUMBER + 1,
C
c     if (.not.rtdb_get(rtdb,'semiemp:infa',mt_dbl,6,infa)) then
c     call errquit('oneint_semiemp: failed to read infa',0, RTDB_ERR)
c      endif
c     write(6,*) "i am in oneint", infa
C   ISW1  =   0      1       2      3      4      5      6     7
C            EHT   CNDO/1 CNDO/2  INDO/1 INDO/2 NDDO/1 NDDO/2 PPP
C         =  51     52      53
C           MNDO    AM1    PM3
      DO 20 i = 1, 4
C
C     Loop over the orbital types:
C          J = 1,2,3,4 for sigma, pi, delta, phi.
C
         DO 10 j = 1, 4
C
C     Standard theory: input
C
            FSPDF(i,j) = ONE  
   10    continue
   20 continue
      if (inttyp.eq.1) then
       FSPDF(1,1) = 1.0
       FSPDF(2,1) = 1.267
       FSPDF(2,2) = 0.585
       FSPDF(3,1) = 1.0
       FSPDF(3,2) = 1.0
       FSPDF(3,3) = 1.0
      end if
      IF(ISW1.NE.8) THEN
        XFAC = ONE
      ELSE
        XFAC = ONE
      ENDIF
c    calculate the one eletron overlap matrix
C     GG is just a number big enough to make
C     exp(-GG) the smallest of  all f's in...
C
C     h(i,j) = f * [( beta(A) + beta(B)) / 2 ] * (delta(i,j))
C

      GG = 4.0D0
C      GET BETA   THIS IS NOW DONE IN PARAM
C     IF (IBETA.EQ.6) CALL PARAMZ
C
C     Set PP (in COMMON PPTT) to -10.D0 so as to initialize OVLAP
C
      PP = -10.D0
c
      NA1 = NA*NA
      NA2 = NA1 + NA1
      IJ = 0
C     Loop over the atomic basis orbitals I
c       write(6,*) "number of basis function", NB
      DO 210 I = 1, NB
C
C       NU(I) is the atom to which orbital I belongs
         KA = NU(I)
C
C       NW(I) is the symmetry of orbital I:  0    1    2    ...
C                                            s   px   py    ...
C
         KB = NW(I)
C
C       KAN(KA) is the atomic number of center KA
C
CSG
         IZA = KAN(I)
C
C       For PPP wavefunction (ISW1 = 7), use only pz orbitals
C
c         IF (ISW1.EQ.7) KB = 3

C       NP(I) is the principal quantum number of atomic orbital I
C
         NT1 = NP(I)
         KL = LANGL(KB+1)
         KLP1 = KL + 1
         KFLAG(1) = KA
         KFLAG(2) = NT1
         KFLAG(3) = KL
C
C       Do an internal write to convert the saved flag values
C       to character strings
          IDD1= 0
c        WRITE (ATOMK,'(3I3)') KFLAG
         DO 200 J = 1, I
            IJ = IJ + 1
            IF (I.NE.J) GO TO 60
C       Diagonal elements of overlap, special (weighted) overlap
C       and kinetic energy
C
            H(IJ) = 1.0
            SH(IJ) = 1.0
            TS     = 0.0
C
C       For IDD1=1, calculate also nuclear-electronic attraction
C
 
            IF (IDD1.LT.9) GO TO 200
          write(6,*)"IDD1",IDD1
C       Do double zeta
C
            DO 40 ISS = 1, 2
               R = ZERO
               AMU = ALPHA(I,ISS+1)
C
C       If exponent is very small, consider it zero -->
C       exp(0) = 1 --> do single zeta
C
               IF (AMU.LT.1.0D-4) GO TO 40
               ACU = ALPHA(I,ISS+4)  
C
C
               DO 30 JSS = 1, 2
                  BMU = ALPHA(I,JSS+1)
C
C       If exponent is very small, consider it zero -->
C       exp(0) = 1 --> do single zeta
C
                  IF (BMU.LT.1.0D-4) GO TO 30
                  BCU = ALPHA(I,JSS+4)
c
   30          CONTINUE
   40       CONTINUE
C
C       ------ End of double zeta loop ---------
C
C       Update matrix FA
C
   50       FA(IJ) = ABS(TS) + FA(IJ)
c           write(6,*) "FA(IJ)=", FA(IJ)
            TTS = 0.0
            GO TO 200
C
   60       CONTINUE
C           TWO CENTER INTEGRALS
C       NU(J) is the atom to which atomic orbital J belongs
C
            MA = NU(J)
C
C       NW(J) is the symmetry of atomic orbital J:
C
C                        0    1    2    3    4   ...
C                        s   px   py   pz   dz2  ...
                                                   
C
            MB = NW(J)
C
C       KAN(MA) is the atomic number of center MA
CSG
            IZB = KAN(J)
C
C       For PPP wavefunction (ISW1 = 7), use only pz orbitals
C
c           IF (ISw1.EQ.7) MB = 3
C
            NT2 = NP(J)
            ML = LANGL(MB+1)
            MLP1 = ML + 1
            MFLAG(1) = MA
C
C       NT2 is the principal quantum number
C
            MFLAG(2) = NT2
            MFLAG(3) = ML
            WRITE (ATOMM,'(3I3)') MFLAG
C
C       One center overlaps (one center kinetic energy is zero)
C
            IF (MA.EQ.KA) THEN
               H(IJ) = 0.0
               SH(IJ) = 0.0
               GO TO 200
             END IF
C
            IF (KB.LE.MB) THEN
               LA = MB
               LH = KB
               KORM = KA + (MA-1)*NA
            ELSE
               LA = KB
               LH = MB
               KORM = MA + (KA-1)*NA
            END IF
C
C     Directional cosines previously calcd. in subroutine GEOMAT
C
            U1 = F(KORM)
            U2 = F(KORM+NA1)
            U3 = F(KORM+NA2)
C
C     The U's are the Eulerian factors
c           write(6,*)"before geom",U1,U2,U3,LH,D
c           write(6,*)"before geom",U1,U2,U3,LA,E
            CALL transmat(U1,U2,U3,LH,D)
            CALL transmat(U1,U2,U3,LA,E)
c           write(6,*)"E",E(1)
C.......................................................
C
C     Get phase for E transform
C
C     PHASE (-1)**(L-M) ON LA.GE.THE OTHER L VALUE
C
            IF (LA.EQ.MB) MM = MLP1
            IF (LA.EQ.KB) MM = KLP1
            GO TO (100,80,90,70), MM
C
C     f orbitals
C
   70       EPHI = ONE
            EDEL = -ONE
C
C     p orbitals
C
   80       EPI = ONE
            ESIG = -ONE
            GO TO 110
C
C     d orbitals
C
   90       EDEL = ONE
            EPI = -ONE
C
C     s orbitals
C
  100       ESIG = ONE
  110       CONTINUE
C
C.......................................................
C
            KAMA = MA + (KA* (KA-1))/2
            R = RAD(KA,MA)
C           IF (ATOMK.EQ.ATOMKX .AND. ATOMM.EQ.ATOMMX) GO TO 150
            ATOMKX = ATOMK
            ATOMMX = ATOMM
C
C     Initialize overlap
C
            SS = ZERO
            SP = ZERO
            SD = ZERO
            SF = ZERO
C
C     Initialize kinetic energy
C
            TS = ZERO
            TP = ZERO
            TD = ZERO
            TF = ZERO
C 
C     Initialize the nuclear - electronic intgral
C
            PS = ZERO
            PPI= ZERO
            PD = ZERO
            PF = ZERO

            PS1 = ZERO
            PPI1= ZERO
            PD1 = ZERO
            PF1 = ZERO
C     Do double zeta for each atom. When exponent is zero, quit.
C
            DO 130 ISS = 1, 2
               AMU = ALPHA(I,ISS+1)
               IF (AMU.LT.1.0D-4) GO TO 130
               ACU = ALPHA(I,ISS+4)
C
               DO 120 JSS = 1, 2
                  BMU = ALPHA(J,JSS+1)
                  IF (BMU.LT.1.0D-4) GO TO 120
                  BCU = ALPHA(J,JSS+4)
C
C     Evaluate product of MO coefficients
C
                  ABC = ACU*BCU
C
C     Evaluate overlap between STO's
C
                  CALL OVLAP(LH,NT1,KL,AMU,NT2,ML,BMU,R,SSS,SSP,SSD,SSF,
     .                       FACT,PP,NIN)
c                 WRITE(6,*) ' ONEINT AFTER OVLAP',i,j,R,SSS,SSP,SSD,SSF
c                 write(6,*)"SS,SSS,ABC",SS,SSS,ABC 
C

C
C     Update overlap
C
                  SS = SS + SSS*ABC
                  SP = SP + SSP*ABC
                  SD = SD + SSD*ABC
                  SF = SF + SSF*ABC
  120          CONTINUE
  130       CONTINUE
C
             RR = R/BOHRS
             GGG = EXP(- (ALPHA(I,1)+ALPHA(J,1))*RR*GG)
C
C
               BETAS = ONE
               BETAP = ONE
               BETAD = ONE
               BETAF = ONE


C
               IF (ABS(FACBLK(1,1)).LT.1.0D-10) THEN
c                 write(6,*) "in ABS(FACBLK)"
c              write(6,*)"KL",KL,"ML",ML
                  IF (KL.EQ.ML) GO TO 140
                  IF (ISW1.EQ.0) GO TO 140
C
                  SIGG = ONE
                  PIG = ONE
                  DELG = ONE
                  PHIG = ONE
                  GO TO 150

  140             CONTINUE
C
C
C      Standard theory: sigma = ( f + f' )/2
C
                  SIGG = (FSPDF(KLP1,1)+FSPDF(MLP1,1))/TWO
                  PIG = (FSPDF(KLP1,2)+FSPDF(MLP1,2))/TWO
                  DELG = (FSPDF(KLP1,3)+FSPDF(MLP1,3))/TWO
                  PHIG = (FSPDF(KLP1,4)+FSPDF(MLP1,4))/TWO
C
                  SIGG = SIGG + (ONE-SIGG)*GGG
                  PIG = PIG + (ONE-PIG)*GGG
                  DELG = DELG + (ONE-DELG)*GGG
                  PHIG = PHIG + (ONE-PHIG)*GGG
               ELSE
C
C         10 SIGMA FACTORS, THEN 10 PI, THEN 10 DELTA, THEN 10 PHI.
C
C             s-s   s-p  p-p  s-d  p-d  d-d  s-f   p-f  d-f  f-f
C -----------------------------------------------------------------
C  FACTR(1) = 1.26 1.26 1.26 1.26 1.26 1.26  1.26  1.26  1.26 1.26
C             0.56 0.56 0.56 0.56 0.56 0.56  0.56  0.56  0.56 0.56
C             1.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0  1.0
C             1.0  1.0  1.0  1.0  1.0  0.87  0.75  0.60  1.0  1.0
C -----------------------------------------------------------------
                  IF (KLP1.GT.MLP1) THEN
                     KLML = NIN(KLP1) + MLP1
                  ELSE
                     KLML = NIN(MLP1) + KLP1
                  END IF
C
                  SIGG = FACBLK(1,KLML)
                  PIG = FACBLK(2,KLML)
                  DELG = FACBLK(3,KLML)
                  PHIG = FACBLK(4,KLML)
C
C ------END FIRST CHOICE -----------------------------------------
C       P-SIGMA, P-PI, D-DELTA, S-PSIG FOR FIRST 10 ATOMS
C -----------------------------------------------------------------
C
C  ATM NO.     H   He   Li    Be    B    C     N     O     F
C  FACTR(1) = 1.0  1.00 1.00 1.0  1.0  1.10  1.25  1.50  1.0  1.0
C             1.0  1.0  1.00 1.0  1.0  0.585 0.72  0.90  1.0  1.0
C             1.0  1.0  1.0  1.0  1.0  1.0   1.0   1.0   1.0  1.0
C             1.0  1.0  1.0  1.0  1.0  0.87  0.75  0.60  1.0  1.0
C -----------------------------------------------------------------
C.......................................................
C     The overlap is SH and the special overlap is H.
C
C
C     The overlap is SH and the special overlap is H.
C
C
C     Form sigma overlap and kinetic energy
         ENDIF
  150       SH(IJ) = E(1)*D(1)*SS*ESIG
            H(IJ) = SH(IJ)*SIGG*BETAS
c           write(6,*)"SH(IJ),SIGG,BETAS,IJ",SH(IJ),SIGG,BETAS,IJ

c           write(6,*)"E,D,SS,ESIG",E(1),D(1),SS,ESIG
c           write(6,*)"H,SH,SIGG,BETAS",H(IJ),SH(IJ),SIGG,BETAS
            EK = ZERO
C
            FA(IJ) = FA(IJ) + XFAC*EK
c            WRITE(6,*) ' ONEINT', FA(IJ)
C            WRITE(6,*) ' ONEINT', ESIG, TS, I, J
C
C    Check how many terms in overlap
C
            MINKM = MIN0(KL,ML) + 1
c           write(6,*)"MINKM",MINKM
            GO TO (190,180,170,160), MINKM
C
C   Form phi overlap and kinetic energy
C
C
c 160       write(6,*) "in 160"
  160       OV = (E(6)*D(6)+E(7)*D(7))*SF*EPHI
c           IF (IDD1.EQ.9.OR.ISW1.EQ.8) THEN
c              EK = (E(6)*D(6)+E(7)*D(7))*TF*EPHI
c              PN(IJ) = PN(IJ) + (E(6)*D(6)+E(7)*D(7))*(PF*Z12(NU(I)) +
c    x                            PF1*Z12(NU(J)))*EPHI
c           ENDIF

            SH(IJ) = SH(IJ) + OV
            FA(IJ) = FA(IJ) + XFAC*EK
            H(IJ) = H(IJ) + OV*PHIG*BETAF
c           write(6,*)"H(IJ) here 1", H(IJ)
C
C  Form delta overlap and kinetic energy
C
c 170      write(6,*)"in 170"
  170       OV = (E(4)*D(4)+E(5)*D(5))*SD*EDEL
c           IF (IDD1.EQ.9.OR.ISW1.EQ.8)  THEN
c              EK = (E(4)*D(4)+E(5)*D(5))*TD*EDEL
c              PN(IJ) = PN(IJ) + (E(4)*D(4)+E(5)*D(5))*(PD*Z12(NU(I)) +
c    x                            PD1*Z12(NU(J)))*EDEL
c           ENDIF
            SH(IJ) = SH(IJ) + OV
            FA(IJ) = FA(IJ) + XFAC*EK
            H(IJ) = H(IJ) + OV*DELG*BETAD
c           write(6,*)"H(IJ) here 2",H(IJ)
C
C  Form pi overlap and kinetic energy
C
c 180       write(6,*)"in 180"
  180       OV = (E(2)*D(2)+E(3)*D(3))*SP*EPI
c           IF (IDD1.EQ.9.OR.ISW1.EQ.8)  THEN
c              EK = (E(2)*D(2)+E(3)*D(3))*TP*EPI
c              PN(IJ) = PN(IJ)+(E(2)*D(2)+E(3)*D(3))*(PPI*Z12(NU(I))+
c    x                          PPI1*Z12(NU(J)))*EPI
c           ENDIF
            SH(IJ) = SH(IJ) + OV
            FA(IJ) = FA(IJ) + XFAC*EK
            H(IJ) = H(IJ) + OV*PIG*BETAP
c           write(6,*)"H(IJ) here 3",H(IJ)
C.......................................................
  190       CONTINUE
  200    CONTINUE
  210 CONTINUE
c     write(6,*)"i am here",(SH(I),I=1,nb*nb)
c     do ij = 1, nb*(NB+1)/2
c       write(6,*)"H(IJ) in oneint",H(IJ)
c     end do 
      RETURN
      END




