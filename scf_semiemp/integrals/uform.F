      SUBROUTINE UFORM(NA,NB,NW,F,H,WCOR,GAMMA,KAN,LQNO,NU,RAD,ISW1,NG,
     1 ALPHA,COFSQ2,NIN,NTYP,DK,Z12,ISW2,NP,SCG1SP,SCF2PP,
     1 SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,SCG3SF,SCG2PF,SCG4PF,
     2 SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,SCF2FF,SCF4FF,SCF6FF,
     3 R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,R2PFFF,R2PPPF,R2SDFF,
     4 R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF,SUMELEL,ETS,ETP,ETD,ETF
     5 ,ISPEC) 
C
C     SUBROUTINE INITIALLY FORMS THE ONE-CENTRE CORE ELEMENTS AND THEN
C     THE COMPLETE ONE ELECTRON HAMILTONIAN MATRIX
C
C     DOES NOT ASSUME THE S AND P ORBITAL FUNCTIONS HAVE THE SAME RADIAL
C     PART.  APRIL 1976.
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/units.cmn'
      DOUBLE PRECISION L, M, N, NVE, NSE, NPE, NDE, NF, NFE
      PARAMETER (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, 
     +  THREE =3.0D0,FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, 
     +  SEVEN = 7.0D0,EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0,
     +   MXA=1200,MXB=4000,MPE=105) 
C
C****************************************************************
      Dimension SCG1SP(*),SCF2PP(*),SCG2SD(*),SCG1PD(*),
     . SCF2PD(*),SCG3PD(*),SCF2DD(*),SCF4DD(*),SCG3SF(*),
     . SCG2PF(*),SCG4PF(*),SCG1DF(*),SCG3DF(*),SCG5DF(*),
     . SCF2PF(*),SCF2DF(*),SCF4DF(*),SCF2FF(*),SCF4FF(*),
     . SCF6FF(*),
C***C D type R-integrals
     .   R1SPPD(*),R2SDDD(*),R2SDPP(*),
C***C F type R-integrals
     .   R1PDDF(*),R1SPDF(*),R2PFDD(*),R2PFFF(*),
     .   R2PPPF(*),R2SDFF(*),R2SDPF(*),R3PDDF(*),
     .   R3SFDF(*),R3SFPD(*),R4PFDD(*),R4PFFF(*)
      Dimension BETB(MPE,4),ALPHA(NB,6)
      Dimension ETS(NA), ETP(NA), ETD(NA), ETF(NA),DK(*),KAN(*)
      Dimension SUMNUCEL(NB), SUMELEL(*),NP(*),Z12(*),NG(NB)
      PARAMETER (EVOLT    = 2.72113961D1)
c     INCLUDE '../include/cf.cmn'
C****************************************************************
C***      COMMON/FACTOR/
C***     .  FPSI1,FDSI1,FPPI1,FDPI1,FDDE1,FSSI1,FACBLK(5,100)
C****************************************************************
C
      DIMENSION COFSQ2(105,5), NIN(*),NTYP(*),NU(NB),NW(NB)
      DIMENSION  WCOR(4000),FKO(8),F(*),H(*),GAMMA(*),RAD(NA,NA)
CGP      DATA         IZERO/0/, IONE/1/
         DATA         FKO/ 0.62500000, 0.36328125, 0.25813802,
     x         0.20090485, 0.16476059, 0.13980329, 0.12150772,
     x         0.10750626/
         
       XX = 0.0d0
       KA = 0
       IPRINT = 0
       IDD1 = 0
       IQMM = 0
       PEN = 0
       KL = 0
C
c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in uform 1"!,(H(I),I=1,40)
      NN = NG(NB)
      ISCM = (NN* (NN+1))/2
      LCOUNT = 0
      GO TO (10,20,10,20,10,20,10,10), ISW1
C     .............................................................
C
C     For ISW1 = 1 (CNDO/1) or ISW1 = 3 (INDO/1), evaluate
C     one center core elements from IP's only (and not IP's + EA's)
C
   10 CONTINUE
      CND = ONE
      GO TO 30

   20 CND = 0.5D0
   30 CONTINUE
C
      II = 0
      KB = 0
C
C     Loop over the basis functions
C
      DO 190 I = 1, NB
         II = II + I
C
C        NU(I) is the atomic center to which basis function I belongs
C
         MA = NU(I)
C
C        NW(I) is the type of ao basis function I is.
C
         MO = NW(I)
C
C        KAN(MA) is the atomic number of atom MA
C
         KA = KAN(I)
         IF(MO.GT.8) THEN
C           F ORBITAL
            JK = 4
         ELSEIF(MO.GT.3) THEN
C           D ORBITAL
            JK = 3
         ELSEIF(MO.GT.0) THEN
C           P ORBITAL
            JK = 2
         ELSE
C           S ORBITAL
            JK = 1
         ENDIF
C        IF(ABS(UCOR(JK,MA)).GT.1.0D-7) THEN
C           THE CORE INTEGRAL HAS BEEN INPUTTED IN A TAB
C           H(II) = H(II) + UCOR(JK,MA)
C           WCOR(I) = UCOR(JK,MA)
C           WRITE(6,2040) MA,I,JK-1, UCOR(JK,MA)
C2040       FORMAT (' CORE INTEGRAL INPUT FOR ATOM',I5,' ORBITAL =',I5,
C    x          ' L =',I5, 'IS ',F10.6)
C           GO TO 190
C        ENDIF
C        Branch for Rydberg or other extra orbitals
C
C        SUBTRACT OF THE IP or (IP-EA)/2 BEFORE CONSTRUCTING WCOR
C        WCOR IS THE CORE INTEGRAL
         KI = NG(I)
C        THE DIAGONAL IS CALCULATED
         IF (KA.EQ.0) GO TO 190
C
C        A polarization function and U is ab-initio
c        IF (NEXTR(I).EQ.1) GO TO 190

C        A CORE ORBITAL AND U IS A LOOK UP.
c        IF (NEXTR(I).LT.0) GO TO 190

C        U IS SEMI EMPIRICAL
c        write(6,*)"H(II),WCOR",H(II),WCOR(I)
         H(II) = H(II) - WCOR(I)
c        write(6,*)"H(II),WCOR",H(II),WCOR(I)
c       write(6,*)"i am in uform 3",(H(K),K=1,NB)
C
C        NVE = number of valence electrons
C
         CALL VALOC(I,DK,NVE,MA,L,M,N,NF,KAN,ISW1,ISW2,COFSQ2,NTYP,
     1   Z12,NU)
c        write(6,*)"L,M,N,NF",L,M,N,NF
C
         IF (NTYP(MA).GT.3 .AND. MO.EQ.0) THEN
            LCOUNT = LCOUNT + 1
            ETS(LCOUNT) = L
            ETP(LCOUNT) = M
            ETD(LCOUNT) = N
            ETF(LCOUNT) = NF
         END IF
C
Crdb
C
C        This section sets the Slater-Condon factors for atom MA;
C        NTYP is the type of basis set for atom MA;
C
C     NTYP(MA) = 1       basis set type: s
C              = 2                       sp
C              = 3                       ns, np, nd
C              = 4                       ns, np, (n-1)d
C              = 5                       ns, np, nd, nf
C
C    F0(sp), F0(sd), F0(pd) are calcd. in subroutine TWOINT
C
C
         IF (KB.EQ.0) GO TO 40
C        IF (KA.EQ.KB) GO TO 50
   40    CONTINUE
         KB = KA
C
C        Slater-Condon factors for sp basis:
C
C        ( s px | s px ) = G1/3, etc.
C
c        write(6,*)"MA",MA
c        write(6,*)"SCG1SP",SCG1SP(MA)
c        write(6,*)"SCF2PP",SCF2PP(MA)
         G1SP = SCG1SP(MA)*THREE
         F2PP = SCF2PP(MA)*25.0D0
C
         IF (NTYP(MA).GE.3) THEN
C
C           Slater-Condon factors for spd basis: G2(sd), G1(pd) etc.
C
            G2SD = SCG2SD(MA)*FIVE
            G1PD = SCG1PD(MA)*15.0D0
            F2PD = SCF2PD(MA)*35.0D0
            G3PD = SCG3PD(MA)*245.0D0
            F2DD = SCF2DD(MA)*49.0D0
            F4DD = SCF4DD(MA)*441.0D0
c        write(6,*)"G2SD,G1PD,F2PD,G3PD,F2DD,F4DD",
c    1   G2SD,G1PD,F2PD,G3PD,F2DD,F4DD
         ELSE
C
C           The remaining factors are obviously zero for a sp basis
C
            G2SD = ZERO
            G1PD = ZERO
            F2PD = ZERO
            G3PD = ZERO
            F2DD = ZERO
            F4DD = ZERO
         END IF
C
         IF (NTYP(MA).GE.5) THEN
C
C
C           Slater-Condon factors for a spdf basis
C
            G3SF = SCG3SF(MA)*SEVEN
            G2PF = SCG2PF(MA)*175.0D0
            G4PF = SCG4PF(MA)*189.0D0
            G1DF = SCG1DF(MA)*35.0D0
            G3DF = SCG3DF(MA)*315.0D0
            G5DF = SCG5DF(MA)*15246.0D0
            F2PF = SCF2PF(MA)*75.0D0
            F2DF = SCF2DF(MA)*105.0D0
            F4DF = SCF4DF(MA)*693.0D0
            F2FF = SCF2FF(MA)*225.0D0
            F4FF = SCF4FF(MA)*1089.0D0
            F6FF = SCF6FF(MA)*736164.0D0
         ELSE
C
C           The remaining factors are obviously zero for a
C           sp or spd basis
C
            G3SF = ZERO
            G2PF = ZERO
            G4PF = ZERO
            G1DF = ZERO
            G3DF = ZERO
            G5DF = ZERO
            F2PF = ZERO
            F2DF = ZERO
            F4DF = ZERO
            F2FF = ZERO
            F4FF = ZERO
            F6FF = ZERO
         END IF
   50    CONTINUE
C
C        Determine the number of s (=L), p (=M) and d (=N)
C        valence electrons the atom has
C
         IF (ISW1.NE.7) GO TO 60
C
C        PPP section: use only p functions, therefore L = N = NF = 0.
C
         L = ZERO
         M = Z12(MA)
         N = ZERO
         NF = ZERO
         ISS = (KI* (KI+1))/2
         ISP = ISS + KI
         IF (ISP.GT.ISCM) ISP = 1
 
C     Update PPP Hamiltonian matrix
C
         WCOR(I) = WCOR(I)  - (M-CND)*GAMMA(ISS)
         WCOR(I) = WCOR(I) + (M-CND)*TWO*F2PP/25.0D0 + L*G1SP/SIX +
     .           N* (G1PD/15.0D0+THREE*G3PD/70.0D0) +
     .           NF* (THREE*G2PF/70.0D0+TWO*G4PF/63.0D0)
         H(II) = H(II) + WCOR(I)
c        write(6,*)"H(II) in 1",H(II)
         GO TO 170

   60    CONTINUE
C        RYDBERG ORBITALS HAVE CORE INTEGRAL CALCULATED IN EXHUCK.
c        IF (NEXTR(I).EQ.1) GO TO 190
C
C        CORE ORBITALS ENERGY FROM EXHUCK, AND ARE ASSUMED OCC.
c        IF (NEXTR(I).LT.0) GO TO 190
C
         IF (NW(I).GT.8) GO TO 140
         KIP = (NW(I)+2)/3 + 1
         GO TO (70,80,120,120), KIP
   70    CONTINUE
C
C        Evaluate one-center, one-electron U(ss) matrix element
C
         ISS = (KI* (KI+1))/2
         ISP = ISS + KI
         IF (ISP.GT.ISCM) ISP = 1
         ISD = ((KI+1)* (KI+2))/2 + KI
         IF (ISD.GT.ISCM) ISD = 1
         ISF = NIN(KI+3) + KI
         IF (ISF.GT.ISCM) ISF = 1
C
C        Update Hamiltonian
C
c      if (ga_nodeid().eq.0)
c    .   write(6,*)"WCOR(I) before",WCOR(I)
         WCOR(I) = WCOR(I) - (L-CND)*GAMMA(ISS) - M*GAMMA(ISP) -
     .                N*GAMMA(ISD) - NF*GAMMA(ISF)
c        write(6,*)"WCOR(I) after",WCOR(I)
         IF (ISW1.LE.2) GO TO 75
C
C        Update Hamiltonian
C
c        write(6,*)"WCOR",WCOR(I)
         WCOR(I) = WCOR(I) + M*G1SP/SIX + N*G2SD/TEN + NF*G3SF/14.0D0
   75    H(II) = H(II) + WCOR(I)
c        write(6,*)"H(II) in 2",H(II)
         GO TO 170
C
C        Evaluate one-center, one-electron U(pp) matrix elements
C
   80    CONTINUE
         IPP = NIN(KI) + KI
         ISP = IPP - 1
         IPD = IPP + KI
         IPF = IPD + 1 + KI
         IF (IPD.GT.ISCM) IPD = 1
         IF (IPF.GT.ISCM) IPF = 1
C
C        Number of valence electrons in the ground state config.
C
         NSE = L
         NPE = M
         NDE = N
         NFE = NF
C
         IF (NTYP(MA).GT.3) GO TO 90
         IF (M.GT.ZERO) GO TO 110
         NPE = ONE
         NSE = L - ONE
         GO TO 110
C
C
   90    CONTINUE
         IF (KA.EQ.20) GO TO 100
         NSE = L - ONE
         NPE = ONE
         IF (NTYP(NA).EQ.5) THEN
C
C           For Lanthanides, use ground state config. s(2)f(n-2)d(0)
C
            IF (KA.GE.58 .AND. KA.LE.71) THEN
               NDE = ZERO
               NFE = NVE - TWO
            END IF
C
C           For Actinides, use ground state config. s(2)p(1)f(n-3)d(0)
C
            IF (KA.GE.89 .AND. KA.LE.103) THEN
               NSE = TWO
               NPE = ONE
               NDE = ZERO
               NFE = NVE - THREE
            END IF
         ELSE
            NDE = N
         END IF
         GO TO 110
C
C        Calcium is an exception: experimental ground state is s(2),
C        but one must use here s(1)p(1)d(0) configuration because the
C        experimental ionization potential is obtained for p electron.
C
  100    CONTINUE
         NSE = ONE
         NPE = ONE
         NDE = ZERO
C
C        For M > 0, there is at least one valence "p" orbital
C
  110    CONTINUE
         WCOR(I) = WCOR(I) - NSE*GAMMA(ISP) - (NPE-CND)*GAMMA(IPP) -
     .                NDE*GAMMA(IPD) - NFE*GAMMA(IPF)
C
C        For INDO or NDDO wavefunction ( ISW1 > 2 )
C
         IF (ISW1.LE.2) GO TO 115
C
C        Update Hamiltonian matrix
C
c        write(6,*)"NPE,CND,F2PP,NSE,G1SP,NDE,G1PD,G3PD,NFE,G2PF,G4PF"
c        write(6,*)NPE,CND,F2PP,NSE,G1SP,NDE,G1PD,G3PD,NFE,G2PF,G4PF
         WCOR(I) =  WCOR(I) + (NPE-CND)*TWO*F2PP/25.0D0 + NSE*G1SP/SIX +
     .           NDE* (G1PD/15.0D0+THREE*G3PD/70.0D0) +
     .           NFE* (THREE*G2PF/70.0D0+TWO*G4PF/63.0D0)
c       write(6,*)"WCOR?, H(II)",WCOR(I),H(II)
  115   H(II) = H(II) + WCOR(I)
c        write(6,*)"H(II) in 3",H(II)
         GO TO 170
C
C        Evaluate one-center, one-electron U(dd) matrix element
C
  120    CONTINUE
         NSE = L
         NPE = M
         NDE = N
         NFE = NF
         IF (NTYP(MA).GE.4) GO TO 130
         NDE = ONE
         NPE = M - ONE
  130    CONTINUE
C
C        For first row polarization functions
C
         IF (KA.LE.10) THEN
C
C        If calculated, then skip
C
            EID = ALPHA(I,1)* (ALPHA(I,1)/TWO-Z12(MA)/DFLOAT(NP(I)))
            IF (H(II).LT.0.9*EID) GO TO 170
C
C           Otherwise, estimate from atomic spectra
C
            NDE = CND
         END IF
C
         IDD = NIN(KI) + KI
         ISD = IDD - 2
         IPD = IDD - 1
         IFD = IDD + KI
         IF (NTYP(MA).EQ.5) THEN
C
C           For Lanthanides (58<KA<71) use config. d(1)f(n-3)
C
            IF (KA.GE.58 .AND. KA.LE.71) THEN
               NDE = 1.D0
               NFE = NVE - 3.D0
C
C             For Actinium (KA=89), use configuration s(2)p(0)d(1)f(n-3)
C
            ELSE IF (KA.EQ.89) THEN
               NSE = TWO
               NPE = ZERO
               NDE = ONE
               NFE = NVE - THREE
C
C     For all other actinides, use configuration-averaged  ground
C     states in which  COFSQ2(KA,I) is the contribution (in terms
C     of the square of the coefficient) of the I-th configuration
C     for atomic center KA.
C
            ELSE IF (KA.GE.90 .AND. KA.LE.103) THEN
               COFN = COFSQ2(KA,1) + COFSQ2(KA,3) + COFSQ2(KA,4)
               NSE = (TWO*COFSQ2(KA,1)+TWO*COFSQ2(KA,3)+COFSQ2(KA,4))/
     .               COFN
               NPE = ZERO
               NDE = (COFSQ2(KA,1)+TWO*COFSQ2(KA,3)+THREE*COFSQ2(KA,4))/
     .               COFN
               NFE = NVE - NSE - NPE - NDE
            END IF
         END IF
C
         IF (IFD.GT.ISCM) IFD = 1
C
C        Update CNDO Hamiltonian matrix
C
         WCOR(I) =  WCOR(I)  - NSE*GAMMA(ISD) - NPE*GAMMA(IPD) -
     .           (NDE-CND)*GAMMA(IDD) - NFE*GAMMA(IFD)
         IF (ISW1.LE.2) GO TO 135
C
C
C        Update INDO Hamiltonian matrix
C
         WCOR(I) = WCOR(I) + (NDE-CND)*TWO* (F2DD+F4DD)/63.0D0 +
     .           NSE*G2SD/TEN + NPE* (G1PD/15.0D0+THREE*G3PD/70.0D0) +
     .           NFE* (THREE*G1DF/70.0D0+TWO*G3DF/105.0D0+
     .           FIVE*G5DF/231.0D0)
c        write(6,*)"WCOR(I)??,H(II)",WCOR(I),H(II)
  135    H(II) = H(II) + WCOR(I)
c        write(6,*)"H(II) in 4",H(II)

C
C        DECOUPLE D ORBITALS FOR FIRST ROW STARTING MATRIX.
C
         GO TO 170
C
C        Evaluate one-center, one-electron U(ff) matrix elements
C
  140    CONTINUE
         NSE = L
         NPE = M
         NDE = N
         NFE = NF
         IF (NTYP(MA).EQ.5) GO TO 160
         IF (N.GE.0.99999) GO TO 150
         NPE = M - ONE
         NFE = ONE
         GO TO 160

  150    NDE = N - ONE
         NFE = ONE
C
C        For spdf basis set (NTYP=5)
C
  160    CONTINUE
C
C        For Lanthanides (58<KA<71)
C
         IF ((KA.GE.58) .AND. (KA.LE.71)) THEN
            NDE = COFSQ2(KA,1)
            NFE = NVE - L - M - NDE
C           ETS(LCOUNT) = NSE
C           ETP(LCOUNT) = NPE
C           ETD(LCOUNT) = NDE
C           ETF(LCOUNT) = NFE
C
C           For actinium (KA=89), use configuration s(2)p(0)d(0)f(n-2)
C
         ELSE IF (KA.EQ.89) THEN
            NSE = TWO
            NPE = ZERO
            NDE = ZERO
            NFE = NVE - TWO
C           ETS(LCOUNT) = NSE
C           ETP(LCOUNT) = NPE
C           ETD(LCOUNT) = NDE
C           ETF(LCOUNT) = NFE
C
C           For thorium (KA=90), use configuration-averaged
C           ground state with s(2)
C
         ELSE IF (KA.EQ.90) THEN
            NSE = TWO
            COFN = COFSQ2(KA,1) + COFSQ2(KA,2) + COFSQ2(KA,5)
            npE = COFSQ2(KA,2)/COFN
            NDE = COFSQ2(KA,1)/COFN
            NFE = NVE - NSE - NPE - NDE
C
C           ETS(LCOUNT) = NSE
C           ETP(LCOUNT) = NPE
C           ETD(LCOUNT) = NDE
C           ETF(LCOUNT) = NFE
C
C           For all other actinides (91<KA<103)
C
         ELSE IF ((KA.GE.91) .AND. (KA.LE.103)) THEN
            NSE = TWO* (COFSQ2(KA,1)+COFSQ2(KA,2)+COFSQ2(KA,3)+
     .            COFSQ2(KA,5)) + COFSQ2(KA,4)
            NPE = COFSQ2(KA,2)
            NDE = COFSQ2(KA,1) + TWO*COFSQ2(KA,3) + THREE*COFSQ2(KA,4)
C
C           Evaluate number of electrons in the f valence shell
C
            NFE = NVE - NSE - NPE - NDE
C
C           Error: NFE cannot be negative for actinides
C
            IF (NFE.LT.ZERO) THEN
               WRITE (6,*) ' ***** ERROR IN UFORM FOR ACTINIDE *****'
               CALL ZEXIT(3)
            END IF
C
C           ETS(LCOUNT) = NSE
C           ETP(LCOUNT) = NPE
C           ETD(LCOUNT) = NDE
C           ETF(LCOUNT) = NFE
         ELSE
            NFE = ZERO
C           ETS(LCOUNT) = NSE
C           ETP(LCOUNT) = NPE
C           ETD(LCOUNT) = NDE
C           ETF(LCOUNT) = NFE
         END IF
         IFF = NIN(KI) + KI
         ISF = IFF - 3
         IPF = IFF - 2
         IDF = IFF - 1
C
C        Update EHT or CNDO Hamiltonian
C
         WCOR(I) = WCOR(I)  - NSE*GAMMA(ISF) - NPE*GAMMA(IPF) -
     .              NDE*GAMMA(IDF) - (NFE-CND)*GAMMA(IFF)
         IF (ISW1.LE.2) GO TO 165
C
C        Update INDO Hamiltonian matrix diagonal element
C
         WCOR(I) = WCOR(I) + NSE*G3SF/14.0D0 +
     .           NPE* (THREE*G2PF/70.0D0+TWO*G4PF/63.0D0) +
     .           NDE* (THREE*G1DF/70.0D0+TWO*G3DF/105.0D0+
     .           FIVE*G5DF/231.0D0) + (NFE-CND)*
     .           (FOUR*F2FF/195.0D0+TWO*F4FF/143.0D0+
     .           100.0D0*F6FF/5577.0D0)
C
C        For CNDO/1 or CNDO/2 wavefunction, simply return
C
c        write(6,*)"WCOR(I)???,H(II)",WCOR(I),H(II)
  165    H(II) = H(II) + WCOR(I)
c        write(6,*)"H(II) in 5",H(II)
  170    CONTINUE
C
C        CORRECTION IF VALENCE CHARGE NOT NEUTRAL ATOM
C        TRANSITION STATE THEORY FOR CORE IONIZATION
C        GOSCINSKI ET AL
C
C        Z12 = number of valence electrons in neutral atom
C           ( e.g. 4 for carbon; 5 for nitrogen etc. )
C
         QQQQ = DK(MA) - Z12(MA)
C
C        For PPP wavefunction, set variable QQQQ equal to zero
C
         IF (ISW1.NE.7) GO TO 180
         QQQQ = ZERO
  180    CONTINUE
C
C        Update diagonal element of the Hamiltonian matrix
C
c       write(6,*)"DFLOAT",NP(I),H(II),QQQQ,ALPHA(I,1)
         H(II) = H(II) + QQQQ*ALPHA(I,1)/DFLOAT(NP(I))
c        write(6,*)"H(II) in 6",H(II)
         WCOR(I) = WCOR(I) + QQQQ*ALPHA(I,1)/DFLOAT(NP(I))
c         write(6,*)"WCOR(I)????,H(II)",WCOR(I),H(II)
  190 CONTINUE
C
C     ...............................
C
      II = 1
C
C     Loop over atomic centers
C
c     write(6,*)"NA",NA
      DO 220 I = 1, NA
C
C     NTYP(I) is the type of basis set for atom I
C
         NN = NTYP(I)
         IF (NN.LT.3) GO TO 210
C
C     For spd (NTYP=3) or spd(n-1) (NTYP=4) basis sets
C
         IS = (II* (II+1))/2
         IP = IS + II + 1
         ID = IS + 4*II + 10
         IF (NN.EQ.4) NN = 3
         IF (NN.EQ.5) GO TO 200
C
C        For specific purposes, print the core integral values
C
         IF (IPRINT.GT.0) WRITE (6,9000) I, H(IS), H(IP), H(ID)

 9000    FORMAT (' CORE INTEGRALS FOR ATOM     ',I5,'  ARE',4F10.6//)

         GO TO 210
C
C        For spdf basis sets (NTYP=5)
C
  200    NN = 4
         IFF = IS + 9*II + 45
         IF (IPRINT.GT.0) WRITE (6,9000) I, H(IS), H(IP), H(ID), H(IFF)
C
C     For s (NTYP=1) or sp (NTYP=2) basis sets
C
  210    II = II + NN*NN
  220 CONTINUE
C
      IF(ISW1.EQ.5.OR.ISW1.EQ.6) GO TO 390
C
C *** THIS PART DONE IN GFORM1 FOR NDDO
C
C      NOW THE COMPLETE ONE ELECTRON HAMILTONIAN
C      ADDING THE NUCLEAR - ELECTRONIC ATTRACTION.
C
C      H(II) = U(II) - SUM(DK(B) * SPEN(AB)
C
C      NOTE HOW THIS TREATS THE TWO CENTRE TERMS WHERE D ELECTRONS ARE
C      INVOLVED
C
      II = 0
C
C     Loop over the basis functions
C
      DO 380 I = 1, NB
         II = II + I
c        write(6,*)"II",II
         MA = NU(I)
         KI = NG(I)
         SUMNUCEL(I) = ZERO
         SUMELEL(I)  = ZERO
         KIKI = NIN(KI) + KI
         GAMA = GAMMA(KIKI)
         Z = ALPHA(I,1)
         AMU = ALPHA(I,1)
         IF (I.EQ.1) GO TO 240
         K = I - 1
         KIM = NG(K)
  240    CONTINUE
C
C +++++++++++++++++++++++++++++++++++++++
C
         DO 310 J = 1, NB
            MB = NU(J)
C
            IF (MA-MB) 250, 300, 260
  250       MAMB = NIN(MB) + MA
            GO TO 270

  260       MAMB = NIN(MA) + MB
  270       RR = RAD(MA,MB)
C
C     IF (MB.NE.MBP) THEN
C       ZLEFT=Z12(MB)
C       MBP=MB
C
c           IF (NEXTR(J).EQ.0) THEN
c        write(6,*)"J,NVE,MB",J,NVE,MB
         CALL VALOC(J,DK,NVE,MB,L,M,N,NF,KAN,ISW1,ISW2,COFSQ2,NTYP,
     1   Z12,NU)
c        write(6,*)"L,M,N,NF",L,M,N,NF
               KK = NW(J)
c              write(6,*)"KK",KK
               ZVAL = NF/SEVEN
               IF (KK.LT.9) ZVAL = N/FIVE
               IF (KK.LT.4) ZVAL = M/THREE
               IF (KK.LT.1) ZVAL = L
c           ELSE
c              IF (NEXTR(J).GT.0) THEN
c                 ZVAL = ZERO
c              ELSE
C
C                 This atomic orbital is a core atomic orbital
C
c                 ZVAL = TWO
c              END IF
c           END IF
C
            KJ = NG(J)
            KJKJ = NIN(KJ) + KJ
            GAMB = GAMMA(KJKJ)
            IF (KI.GT.KJ) KA = KJ + (KI* (KI-1))/2
            IF (KI.LT.KJ) KA = KI + (KJ* (KJ-1))/2
            XX = ZERO
C
C           For IDD1=1, evaluate nuclear electronic attraction
C
            IF (IDD1.EQ.1) THEN
C
C              Loop over the first double zeta function
C
               DO 290 JJ = 2, 3
                  IF(ISW1.EQ.8) THEN
C                   THIS GENERATES THE EFFECTIVE EXPONENT FROM THE
C                   COULOMB INTEGRAL
                    AMU = GAMA /FKO(NP(I))
                    CMU = ONE
                    IF(JJ.EQ.3) AMU = ZERO
                  ELSE
                    AMU = ALPHA(I,JJ)
                    CMU = ALPHA(I,JJ+3)
                  ENDIF
                  IF (AMU.LT.1.D-4) GO TO 290
C
C       Loop over the second double zeta function
C
                  DO 280 JJJ = 2, 3
                  IF(ISW1.EQ.8) THEN
                    BMU = GAMB /FKO(NP(J))
                    DMU = ONE
                    IF(JJJ.EQ.3) AMU = ZERO
                  ELSE
                     BMU = ALPHA(I,JJJ)
                     DMU = ALPHA(I,JJJ+3)
                  ENDIF
C
C       If exponent is very small, set it equal to zero -->
C         exp(0) = 1 ----> do single zeta
C
                     IF (BMU.LT.1.D-4) GO TO 280
C
C       Evaluate penetration integral through subroutine PENET
C
c                    CALL PENET(FACT,AMU,NP(I),IZERO,BMU,NP(I),IZERO,RR,
c    .                          PEN,PENPI,PEND,PENF,IONE)
                     XX = XX + PEN*CMU*DMU
  280             CONTINUE
  290          CONTINUE
            END IF
C
C           Update diagonal element of the Hamiltonian matrix
C
            KL=1
c           write(6,*)"ZVAL",ZVAL
            ZZZ = SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC,KL
     X       )*ZVAL
c           write(6,*)"ZZZ,H(II)",ZZZ,H(II)
            SUMNUCEL(I) = SUMNUCEL(I) - ZZZ
            H(II) = H(II) - ZZZ
c           write(6,*)"H(II) in 8",H(II)
C
  300       CONTINUE
  310    CONTINUE
C
C        This part is for expanded basis set
C
         DO 370 J = 1, NA
            IF (KAN(J).NE.0) GO TO 370
            IF (ABS(Z12(J)).LT.0.0001) GO TO 370
            IF (MA-J) 320, 370, 330
  320       IJ = NIN(J) + MA
            GO TO 340

  330       IJ = NIN(MA) + J
  340       RR = RAD(MA,J)
C
C.......................................................
C
C     DOUBLE ZETA FOR PENETRATION INTEGRALS
C
            SS = 0.D0
C
C     Loop over the first double zeta function
C
            DO 360 JJ = 2, 3
               AMU = ALPHA(I,JJ)
               CMU = ALPHA(I,JJ+3)
C
C     If exponent is very small, consider it zero --->
C           exp(0) = 1 ----> do single zeta
C
               IF (AMU.LT.1.D-4) GO TO 360
C
C     Loop over the second double zeta function
C
               DO 350 JJJ = 2, 3
                  BMU = ALPHA(I,JJJ)
                  DMU = ALPHA(I,JJJ+3)
C
C     If exponent is very small, consider it zero ---->
C           exp(0) = 1 ----> do single zeta
C
                  IF (BMU.LT.1.D-4) GO TO 350
C
C     Evaluate penetration integrals through subroutine PENET
C
c                 CALL PENET(FACT,AMU,NP(I),IZERO,BMU,NP(I),IZERO,RR,
c    .                       PEN,PENPI,PEND,PENF,IONE)
                  SS = SS + PEN*CMU*DMU
  350          CONTINUE
  360       CONTINUE
C
C           Update diagonal element of the Hamiltonian matrix
C
c           write(6,*)"Z12,SS,H(II)",Z12(J),SS,H(II)
            H(II) = H(II) - Z12(J)*SS
c           write(6,*)"H(II) in 9",H(II)
            SUMNUCEL(I) = SUMNUCEL(I) - Z12(J)*SS
C
C.....................................................................*
C
  370    CONTINUE
  380 CONTINUE
  390 CONTINUE
C     WRITE(6,*) ' UFORM', NVE,L, M, N
C     WRITE(6,*) ' UFORM', DK(1), Z12(1), D12(1)
C
C     IF XINDO , CONSIDER THIS!!
      IF(ISW1.EQ.8) THEN
C       GET THE WEIGHTED OVERLAP IN F
C       REWIND 16
c       CALL ZREWIND (IO_SCF_MATRICES, 'UFORM')
C       READ (16) (F(I),I=1,NRR)
c       CALL ZRREAD (IO_SCF_MATRICES, F, 1, NRR, 'UFORM')
C
C       GET BETB
        CALL PARAMD(BETB)
C
        RR = ZERO

        DO 440 I = 2, NB
        ID = NIN(I)
        KA = NU(I)
        NGI = NIN(NG(I))
        ZI = ONE
        KII = NGI+NG(I)
        GAMA = GAMMA(KII)
        KL=2
        POTI = SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC,KL)
        XXX = BETB(KAN(NU(I)),LQNO(I)+1)/ EVOLT
C
        DO 420 J = 1, I-1
        IJ = ID + J
        KB = NU(J)
c       POT = -SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC)
c    X  *Z12(KA)
        IF(KA.EQ.KB) GO TO 420
        KAKB = NIN(KA) + KB
C       RBOHRS = RAD(KAKB)/BOHRS

C       WRITE(6,*) ' UFORM I, J, SUMNUCEL',I,J,SUMNUCEL(I),SUMNUCEL(J)
C       SUMNUCEL(I) DOES NOT CONTAIN THE ONE-CENTER TERM
C       write(6,*)"H(IJ)",H(IJ)
C       write(6,*)"SUMNUCEL(I,J),F(IJ)",SUMNUCEL(I),SUMNUCEL(J),F(IJ)
        H(IJ) = H(IJ) +
     .           (SUMNUCEL(I) + SUMNUCEL(J))*F(IJ)/TWO
C         FOR XINDO1  ADD  BACK THE TWO CENTER TERMS THAT ARE PART OF
C             WCOR
C         write(6,*)"i am here ?"
          RRR = RAD(KA,KB)
          NGJ = NIN(NG(J))
          GAMA = GAMMA(NGI+NG(I))
          GAMB = GAMMA(NGJ+NG(J))
          KIJ = NGI + NG(J)
          KJJ = NIN(NG(J)) + NG(J)
          KL=3
          SUMFIX=SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC,KL)
     X    *Z12(KB)+SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC,
     X    KL)*Z12(KA)
C         WRITE(6,*) ' UFORM', WCOR(I), WCOR(J), F(IJ)
C         write(6,*)"H(IJ),SUMFIX,F(IJ)",H(IJ),SUMFIX,F(IJ) 
          H(IJ) = H(IJ)+ SUMFIX * F(IJ)/TWO
          ZJ = ONE
          GAMB = GAMMA(KJJ)
          KL=4
          POTJ =SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,ISPEC,KL)
          XX = XXX + BETB(KAN(NU(J)),LQNO(J)+1)/EVOLT
c         write(6,*)"H(IJ),XX,F(IJ)",H(IJ),XX,F(IJ)
          H(IJ) = H(IJ) + XX*F(IJ)/TWO
c         write(6,*)"H(IJ)",H(IJ)
c         H(IJ) = H(IJ) - PN(IJ)

  420   CONTINUE
  440   CONTINUE
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1  , 'UFORM')
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1  , 'UFORM')
      ENDIF
C
C     WRITE(6,*) ' UFORM ETS '
C     IF(LCOUNT.GT.0) THEN
C     WRITE (6,*) (ETS(I),I=1,LCOUNT)
C     WRITE (6,*) (ETP(I),I=1,LCOUNT)
C     WRITE (6,*) (ETD(I),I=1,LCOUNT)
C     WRITE (6,*) (ETF(I),I=1,LCOUNT)
C     ENDIF
c        write(6,*)"LCOUNT",LCOUNT
c     write(6,*)"leaving uform"
c     write(6,*)"i am here 8",(H(I),I=1,40)
      RETURN

      END


C
C  +++++++++++++++++++++++++  VALOC  ++++++++++++++++++++++++++++++++++
C
      SUBROUTINE VALOC(I,DK,NVE,MA,L,M,N,NF,KAN,ISW1,ISW2,COFSQ2,NTYP,
     1 Z12,NU)
C
C     Determines the number of  s(=L), p(=M) and d(=N)
C     valence electrons the atom has.
C
C      INCLUDE '../include/param.cmn'
C      INCLUDE '../include/units.cmn'
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/units.cmn'
      PARAMETER (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0,
     +  THREE =3.0D0,FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0,
     +  SEVEN = 7.0D0,EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0,
     +   MXA=1200,MXB=4000)
      DOUBLE PRECISION L, M, N, NVE, NF
C****************************************************************
C      INCLUDE '../include/switch.cmn'
      Dimension  DK(*),NTYP(*),KAN(*),Z12(*),COFSQ2(105,5)
      Dimension NU(*)

c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in valoc"!,KAN(I)
C
C     IF(NEXTR(I).NE.0) THEN
C       WRITE(6,2000) I
C2000   FORMAT ( ' *****ERROR IN VALOC*****, ORBITAL',I5,
C    . '  NOT A VALENCE TYPE')
C       CALL ZEXIT(3)
C     ENDIF
      KA = KAN(I)
C     write(6,*)"KA",KA
C
C     Are the core orbitals included ???
C
C     IF (JSW(11).GT.1) THEN
C
C       CORRECT THE NUCLEAR NUMBER FOR THE NO. OF VALENCE ELECTRONS
C
C        NWW = NW(I)
C  10    NPS = NP(I-NWW) - 1
C        NVE = DK(MA) - TWO* (NPS)**2
C        IF (NPS.GT.1) GO TO 10
C     ELSE
         NVE = DK(MA)
C     END IF
C
      L = ZERO
      M = ZERO
      N = ZERO
      NF = ZERO
C
C     For hydrogen (KA=1), there is only one s orbital
C
      IF (KA.EQ.1) THEN
         L = ONE
C
C     For 1st and 2nd row elements (including Potassium)
C
      ELSE IF (KA.LT.20) THEN
         L = MIN(TWO,NVE)
         M = MAX(NVE-TWO,ZERO)
C
C     For 1st row transition elements (including Copper)
C
C
C     For 1st row transition elements (including Copper)
C
      ELSE IF (KA.LT.30) THEN
         IF (NTYP(MA).EQ.1) THEN
            IF(Z12(MA).GT.TWO) THEN
              WRITE(6,*) ' ****ERROR IN BASIS SET DEFN*****'
              CALL ZEXIT(3)
            ENDIF
            L = NVE
            M = ZERO
            N = ZERO
         ELSEIF (NTYP(MA).EQ.2) THEN
            L = MIN(TWO,NVE    )
            M = MAX(ZERO,NVE    -TWO)
            N = ZERO
         ELSEIF (ISW2.EQ.4) THEN
            M = ZERO
            L = COFSQ2(KA,1) + TWO*COFSQ2(KA,2)
            N = NVE - L
         ELSE
            M = ZERO
            L = ONE
            IF (ISW2.EQ.2) L = TWO
            IF (ISW2.EQ.3) L = ZERO
            N = NVE - L
         END IF
C
C     For zinc (KA=30), ground state config. is s(2)p(0)d(10)
C
      ELSE IF (KA.EQ.30) THEN
         IF (NTYP(MA).EQ.4) THEN
            L = TWO
            M = ZERO
            N = TEN
         ELSE
            L = TWO
            M = ZERO
         END IF
C
C     For 2nd row transition elements (31<KA<38), use config.
C         s(2)p(n-2)       or       s(n)p(0)
C
      ELSE IF (KA.LE.38) THEN
         L = MIN(TWO,NVE)
         M = MAX(NVE-TWO,ZERO)
C
C
C     For 2nd row transition elements (39<KA<47),
C     use configuration-averaged ground states
C
      ELSE IF (KA.LE.47) THEN
         IF (NTYP(MA).EQ.1) THEN
            IF(Z12(MA).GT.TWO) THEN
              WRITE(6,*) ' ****ERROR IN BASIS SET DEFN*****'
              CALL ZEXIT(3)
            ENDIF
            L = NVE
            M = ZERO
            N = ZERO
         ELSEIF (NTYP(MA).EQ.2) THEN
            L = MIN(TWO,NVE    )
            M = MAX(ZERO,NVE    -TWO)
            N = ZERO
         ELSE
         L = COFSQ2(KA,1) + TWO*COFSQ2(KA,2) + COFSQ2(KA,4)
         M = COFSQ2(KA,4)
         N = (NVE-ONE)*COFSQ2(KA,1) + (NVE-TWO)*COFSQ2(KA,2) +
     .       NVE*COFSQ2(KA,3) + (NVE-ONE)*COFSQ2(KA,4)
C
C     Cadmium (KA=48) is an exception: use p(0) configuration
C
        ENDIF
      ELSE IF (KA.EQ.48) THEN
         IF (ISW2.EQ.4) THEN
            M = ZERO
            L = COFSQ2(KA,1) + TWO*COFSQ2(KA,2)
            N = NVE - L
         ELSE
            M = ZERO
            L = ONE
            IF (ISW2.EQ.2) L = TWO
            IF (ISW2.EQ.3) L = ZERO
            N = NVE - L
         END IF
C
C     Remaining elements (until Barium)
C
      ELSE IF (KA.LE.56) THEN
         L = MIN(TWO,NVE)
         M = MAX(NVE-TWO,ZERO)
C
C     For Lanthanides ( 57 < atomic number < 71 ),
C     the s(2)p(0)d(1)f(n-3) configuration dominates
C
      ELSE IF (KA.GE.57 .AND. KA.LE.71) THEN
         L = TWO
         M = ZERO
         N = ONE
         NF = NVE - THREE
C
C     Remaining elements (until Radium)
C
      ELSE IF (KA.LE.80) THEN
C     third transition metal series
         IF (NTYP(MA).EQ.1) THEN
            IF(Z12(MA).GT.TWO) THEN
              WRITE(6,*) ' ****ERROR IN BASIS SET DEFN*****'
              CALL ZEXIT(3)
            ENDIF
            L = NVE
            M = ZERO
            N = ZERO
         ELSEIF (NTYP(MA).EQ.2) THEN
            L = MIN(TWO,NVE    )
            M = MAX(ZERO,NVE    -TWO)
            N = ZERO
         ELSE
         L = (TWO*COFSQ2(KA,2)+COFSQ2(KA,1)+COFSQ2(KA,4))/(COFSQ2(KA,2)+
     X        COFSQ2(KA,1)+COFSQ2(KA,3)+COFSQ2(KA,4))
         M = COFSQ2(KA,4)+COFSQ2(KA,5)
         N = NVE - L - M
         ENDIF
      ELSE IF (KA.LE.88) THEN
         L = MIN(TWO,NVE)
         M = MAX(NVE-TWO,ZERO)
      ELSE
C
C     Use configuration-averaged ground states
C
         L = TWO* (COFSQ2(KA,1)+COFSQ2(KA,2)+COFSQ2(KA,3)+
     .       COFSQ2(KA,5)) + COFSQ2(KA,4)
         M = COFSQ2(KA,2)
         N = COFSQ2(KA,1) + TWO*COFSQ2(KA,3) + THREE*COFSQ2(KA,4)
         NF = NVE - L - M - N
      END IF
C
      RETURN

      END
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
      FUNCTION SPEN(RR,GAMA,GAMB,Z,GAMMA,KA,XX,II,JJ,NW,ISW1,
     1       ISPEC,KL)
C
C     FUNCTION EVALUATES THE NUCLEAR ELECTRONIC REPULSION
C     KA PULLS OUT GAMMA(KA)
C
      implicit none
C
C
      integer  la, kd,ld
      integer IPEC,ISPEC,KL
      Double precision   GAMMA(*),SPEN
      integer NW(*)
      integer ISW1,ISW2,KA,IDD1,II,JJ
      double precision W,S,ayw,rr,r,GAMA,GAMB,z,XX,ONE 
      PARAMETER (IDD1= 0,ONE=1.0D0)

C
      ayw = 0.0d0
      SPEN=0.0d0
      r = 0.0d0
      IPEC = ISPEC + 1
c     write(6,*)"i am in spen" ,ISPEC,IPEC,KL
cywh
C     if (isw1.eq.8.OR.JSW(16).EQ.1) go to 50
      if (isw1.eq.8) go to 50
cywh
      IF (IPEC.EQ.0.OR.IDD1.EQ.1) THEN
C     THIS HAS BEEN CALCULATED THEORETICALLY FROM PENET
        SPEN = XX
      ELSE
        IF(IPEC.LE.0) IPEC = 1
        GO TO (10,20,10,30,30), IPEC

   10   CONTINUE
c       write(6,*)"KA",KA
c       write(6,*)"GAMMA(KA)",GAMMA(KA)
        S = GAMMA(KA)
        GO TO 40

   20   CONTINUE
   30   CONTINUE
C       SPECTROSCOPIC INTEGRALS
C       B=1.0D0
C       C=2.224D0
C       R=RR/BOHRS
C       W=1.0D0+B*EXP(-C*(GAMA+GAMB)*R/2.0D0)
        W = 1.0D0
c       write(6,*)"GAMMA(KA)",GAMMA(KA)
        S = GAMMA(KA)*W
   40   CONTINUE
        IF (IDD1.EQ.1) THEN
C        SPEN = 0.3*XX + 0.7D0*S
         SPEN = S
        ELSE
         SPEN = S
        END IF
      ENDIF
      go to 1000
   50 continue
c
cyw   or mimic the way in coul.f to use PENA, but now I am not gonna totally
c     do that.
c
c
c     nw(i): symmetry of ao I
c
         la = nw(ii)
         if (LA .LE. 0) then
            kd = 1
         elseif (LA .LE. 3) then
            kd = 2
         elseif (LA .LE. 8) then
            kd = 3
         else
            kd = 4
         endif
c        write(6,*) "kd = ", kd
         la = nw(jj)
         if (LA .LE. 0) then
            ld = 1
         elseif (LA .LE. 3) then
            ld = 2
         elseif (LA .LE. 8) then
            ld = 3
         else
            ld = 4
         endif
c        write(6,*) "ld = ", ld
c
c kan: atomic number, ayw: (uu|1/R) = (ayw**2 + r**2)**(-0.5)
c
C        ayw = ( pena(kd,kan(nu(ii))) +  pena(ld,kan(nu(jj))))/TWO
c        WRITE(6,*) ' SPEN ', GAMA, GAMB, GAMMA(KA)
c        IF (ABS( pena(KD,NU(II))).LT.1.0D-6)
c    x      pena(KD,NU(II)) = GAMA*EVOLT
c        IF(ABS( pena(LD,NU(JJ))).LT.1.0D-6)
c    x      pena(LD,NU(JJ)) = GAMB*EVOLT
c        ayw = ( pena(kd,(nu(ii))) +  pena(ld,(nu(jj))))/TWO
c        r = rr/bohrs
c        ayw = evolt/ayw
C        write(6,*) "ayw= ", ayw
         IF(ISPEC.EQ.3) THEN
            spen = ONE/sqrt(ayw**2 + r**2)
         ELSEIF(ISPEC.EQ.1) THEN
            spen = ONE/(ayw + r)
C        ELSEIF(ISPEC.EQ.0) THEN
C           CALL PENET(FACT,AMU,N,LA,BMU,M,LB,RR,PEN,PENPI,PEND,PENDF,
C    x                IONE)
C           SPEN = PS
         ENDIF
   60 continue
c        write(6,*) 'spen= ', spen
 1000 continue
c        write(6,*) 'spen= ', spen

      END





