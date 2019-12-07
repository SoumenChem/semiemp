      logical FUNCTION COUL(M1,M2,K,L,R,LA,LB,AMU,BMU,NP1,NP2,NEXTA,
     x             NEXTB,FACT,ISPP,TAG)
C
C     SPECTROSCOPIC TWO ELECTRON INTEGRALS FROM
C     ISPP = 1 WEISS MATAGA NISHIMOTO METHOD
C     ISPP = 3 FROM  OHNO KLOPMAN METHOD
C     ISPP = 5 FROM MODIFIED WARSHELL METHOD
C
C      K AND L ARE THE ATOMIC NUMBERS, R THE DISTANCE IN ANGSTROMS
C      LA AND LB THE ORBITAL TYPE (S, PX, ETC)
C
C       WSS CONTAINS GAMMA SS
C       WPP CONTAINS GAMMA PP
C       WSP CONTAINS GAMMA SP
C       WDD CONTAINS METAL GAMMA DD
C       WSD CONTAINS METAL GAMMA SD,PD
C       WSF CONTAINS METAL GAMMA SF,PF
C       WDF CONTAINS METAL GAMMA DF
C       WFF CONTAINS METAL GAMMA FF
C
C       WDDP CONTAINS THE GAMMA DD FOR POLARIZATION FUNCTIONS
C       WSDP CONTAINS THE GAMMA SD FOR POLARIZATION FUNCTIONS
C       WSPP CONTAINS THE GAMMA SP FOR POLARIZATION FUNCTIONS
C       WPPP CONTAINS THE GAMMA PP FOR POLARIZATION FUNCTIONS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
C
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/datain.cmn'
c     INCLUDE '../include/const.cmn'
C******************************************************************
      PARAMETER (MXB = 700, MXA = 100, MPE = 105)
C******************************************************************
C MXB IS THE BASIS SET SIZE.  MXB IS THE MAX. NUMBER OF ATOMS ALLOWED
C********************************************************************
c     INCLUDE '../include/cofexp.cmn'
C********************************************************************
      COMMON/COFEXP/
     .              AS(MPE,3),AP(MPE,3),AD(MPE,3),AF(MPE,3),
     .              CS(MPE,3),CP(MPE,3),CD(MPE,3),CF(MPE,3)
C********************************************************************
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)

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
      DIMENSION FACT(*)
      DIMENSION WSS(MPE), WSP(MPE), WPP(MPE), WSD(MPE), WDD(MPE),
     +          WSF(MPE), WDF(MPE), WFF(MPE)
      DIMENSION WSPP(MPE), WPPP(MPE), WSDP(MPE), WDDP(MPE),
     +          WSFP(MPE), WPFP(MPE), WDFP(MPE), WFFP(MPE)
C
      DIMENSION NPPP(MPE)
C     DIMENSION WPF(MPE),WPDP(MPE),WPD(MPE)
C     ----------------------------------------------------------
      DATA WSS /
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   4.57D0,  6.78D0,
C         B        C        N        O        F        Ne
     .   8.68D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
C    .   3.31D0,  4.79D0,    Original smoothed values.
     .   4.31D0,  5.25D0,
C         Al       Si       P        S        Cl       Ar
     .   6.21D0,  7.57D0,  8.86D0, 10.09D0, 11.25D0,  0.00D0,
C         K        Ca
C    .   3.18D0,  3.25D0,   original smoothed values
     .   3.84D0,  3.87D0,
C         Sc       Ti       V        Cr       Mn
     .   3.89D0,  4.50D0,  5.07D0,  5.60D0,  6.09D0,
C         Fe       Co       Ni       Cu       Zn
     .   6.54D0,  6.96D0,  7.34D0,  7.68D0,  7.98D0,
C         Ga       Ge       As       Se       Br       Kr
     .   5.70D0,  0.00D0,  0.00D0,  8.35D0,  9.08D0,  0.00D0,
C         Rb       Sr
C    .   1.50D0,  2.05D0,   original smoothed values
     .   3.70D0,  3.14D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ      WALTER THIEL Zr = 6.061
CMCZ .   2.57D0,  3.06D0,  3.52D0,  3.98D0,  4.40D0,
     .   2.57D0,  3.06D0,  3.52D0,  5.57D0,  4.40D0,
C         Ru       Rh       Pd       Ag       Cd
     .   4.81D0,  5.19D0,  5.54D0,  0.00D0,  6.19D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  5.76D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   3.13D0,  4.28D0,
C         La
     .   5.420D0,
C**  .   5.76D0,  5.98D0,  6.00D0,  6.14D0,  6.22D0,  6.30D0,  6.38D0,
C**  .   6.47D0,  6.56D0,  6.64D0,  6.75D0,  6.81D0,  6.89D0,  6.98D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.504D0, 5.585D0, 5.663D0, 5.740D0, 5.815D0, 5.889D0, 5.969D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.047D0, 6.126D0, 6.205D0, 6.284D0, 6.363D0, 6.444D0, 6.524D0,
C                  Hf       Ta       W        Re
     x            5.12D0,  5.30D0,  5.48D0,  5.88D0,
C         Os       Ir       Pt       Au       Hg
     x   5.88D0,  6.08D0, 6.29D0, 6.51D0, 0.0D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     x  25*0.00D0/
C     ----------------------------------------------------------
      DATA WSP /
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   4.57D0,  6.78D0,
C         B        C        N        O        F        Ne
     .   8.68D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
C    .   3.31D0,  4.79D0,    Original smoothed values.
     .   4.31D0,  5.25D0,
C         Al       Si       P        S        Cl       Ar
     .   6.21D0,  7.57D0,  8.86D0, 10.09D0, 11.25D0,  0.00D0,
C         K        Ca
C    .   3.18D0,  3.25D0,   original smoothed values
     .   3.84D0,  3.87D0,
C         Sc       Ti       V        Cr       Mn
     .   3.89D0,  4.50D0,  5.07D0,  5.60D0,  6.09D0,
C         Fe       Co       Ni       Cu       Zn
     .   6.54D0,  6.96D0,  7.34D0,  7.68D0,  7.98D0,
C         Ga       Ge       As       Se       Br       Kr
     .   5.70D0,  0.00D0,  0.00D0,  8.35D0,  9.08D0,  0.00D0,
C         Rb       Sr
C    .   1.50D0,  2.05D0,   original smoothed values
     .   3.70D0,  3.14D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ .   2.57D0,  3.06D0,  3.52D0,  3.98D0,  4.40D0,
     .   2.57D0,  3.06D0,  3.52D0,  5.57D0,  4.40D0,
C         Ru       Rh       Pd       Ag       Cd
     .   4.81D0,  5.19D0,  5.54D0,  0.00D0,  6.19D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  5.76D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   3.13D0,  4.28D0,
C         La
     .   6.08D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.504D0, 5.585D0, 5.663D0, 5.740D0, 5.815D0, 5.889D0, 5.969D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.047D0, 6.126D0, 6.205D0, 6.284D0, 6.363D0, 6.444D0, 6.524D0,
C                  Hf       Ta       W        Re
     x            5.12D0,  5.30D0,  5.48D0,  5.88D0,
C         Os       Ir       Pt       Au       Hg
     x   5.88D0,  6.08D0, 6.29D0, 6.51D0, 0.0D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .  25*0.00D0/
C     ----------------------------------------------------------
      DATA WSPP /
C         H       He
     .   7.45D0,  0.00D0,
CWSP .  12.85D0,  0.00D0,
C         Li      Be
     .   4.57D0,  6.78D0,
C         B        C        N        O        F        Ne
     .   8.68D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
C    .   3.31D0,  4.79D0,    Original smoothed values.
     .   4.31D0,  5.25D0,
C         Al       Si       P        S        Cl       Ar
     .   6.21D0,  7.57D0,  8.86D0, 10.09D0, 11.25D0,  0.00D0,
C         K        Ca
C    .   3.18D0,  3.25D0,   original smoothed values
     .   3.84D0,  3.87D0,
C         Sc       Ti       V        Cr       Mn
     .   3.89D0,  4.50D0,  5.07D0,  5.60D0,  6.09D0,
C         Fe       Co       Ni       Cu       Zn
     .   6.54D0,  6.96D0,  7.34D0,  7.68D0,  7.98D0,
C         Ga       Ge       As       Se       Br       Kr
     .   5.70D0,  0.00D0,  0.00D0,  8.35D0,  9.08D0,  0.00D0,
C         Rb       Sr
C    .   1.50D0,  2.05D0,   original smoothed values
     .   3.70D0,  3.14D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ .   2.57D0,  3.06D0,  3.52D0,  5.57D0,  4.40D0,
     .   2.57D0,  3.06D0,  3.52D0,  3.98D0,  4.40D0,
C         Ru       Rh       Pd       Ag       Cd
     .   4.81D0,  5.19D0,  5.54D0,  0.00D0,  0.00D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  5.76D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   3.13D0,  4.28D0,
C         La
     .   6.08D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.504D0, 5.585D0, 5.663D0, 5.740D0, 5.815D0, 5.889D0, 5.969D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.047D0, 6.126D0, 6.205D0, 6.284D0, 6.363D0, 6.444D0, 6.524D0,
C                  Hf       Ta       W        Re
     x            5.12D0,  5.30D0,  5.48D0,  5.88D0,
C         Os       Ir       Pt       Au       Hg
     x   5.88D0,  6.08D0,  6.29D0,  6.51D0,  0.0D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .  25*0.00D0/
C     ----------------------------------------------------------
      DATA WPP/
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   4.57D0,  6.78D0,
C         B        C        N        O        F        Ne
     .   8.68D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
C    .   3.31D0,  4.79D0,    Original smoothed values.
     .   4.31D0,  5.25D0,
C         Al       Si       P        S        Cl       Ar
     .   6.21D0,  7.57D0,  8.86D0, 10.09D0, 11.25D0,  0.00D0,
C         K        Ca
C    .   3.18D0,  3.25D0,   original smoothed values
     .   3.84D0,  3.87D0,
C         Sc       Ti       V        Cr       Mn
     .   3.89D0,  4.50D0,  5.07D0,  5.60D0,  6.09D0,
C         Fe       Co       Ni       Cu       Zn
     .   6.54D0,  6.96D0,  7.34D0,  7.68D0,  7.98D0,
C         Ga       Ge       As       Se       Br       Kr
     .   5.70D0,  0.00D0,  0.00D0,  8.35D0,  9.08D0,  0.00D0,
C         Rb       Sr
C    .   1.50D0,  2.05D0,   original smoothed values
     .   3.70D0,  3.14D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ .   2.57D0,  3.06D0,  3.52D0,  3.98D0,  4.40D0,
     .   2.57D0,  3.06D0,  3.52D0,  5.57D0,  4.40D0,
C         Ru       Rh       Pd       Ag       Cd
     .   4.81D0,  5.19D0,  5.54D0,  0.00D0,  6.19D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  5.76D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   3.13D0,  4.28D0,
C         La
     .   6.08D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.504D0, 5.585D0, 5.663D0, 5.740D0, 5.815D0, 5.889D0, 5.969D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.047D0, 6.126D0, 6.205D0, 6.284D0, 6.363D0, 6.444D0, 6.524D0,
C                  Hf       Ta       W        Re
     x            5.12D0,  5.30D0,  5.48D0,  5.88D0,
C         Os       Ir       Pt       Au       Hg
     x   5.88D0,  6.08D0,  6.29D0,  6.51D0,  0.0D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .  25*0.00D0/
C     ----------------------------------------------------------
      DATA WPPP /
C         H       He
     .   7.45D0,  0.00D0,
CWPP .  12.85D0,  0.00D0,
C         Li      Be
     .   4.57D0,  6.78D0,
C         B        C        N        O        F        Ne
     .   8.68D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
C    .   3.31D0,  4.79D0,    Original smoothed values.
     .   4.31D0,  5.25D0,
C         Al       Si       P        S        Cl       Ar
     .   6.21D0,  7.57D0,  8.86D0, 10.09D0, 11.25D0,  0.00D0,
C         K        Ca
C    .   3.18D0,  3.25D0,   original smoothed values
     .   3.84D0,  3.87D0,
C         Sc       Ti       V        Cr       Mn
     .   3.89D0,  4.50D0,  5.07D0,  5.60D0,  6.09D0,
C         Fe       Co       Ni       Cu       Zn
     .   6.54D0,  6.96D0,  7.34D0,  7.68D0,  7.98D0,
C         Ga       Ge       As       Se       Br       Kr
     .   0.00D0,  0.00D0,  0.00D0,  8.35D0,  9.08D0,  0.00D0,
C         Rb       Sr
C    .   1.50D0,  2.05D0,   original smoothed values
     .   3.70D0,  3.14D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ .   2.57D0,  3.06D0,  3.52D0,  3.98D0,  4.40D0,
     .   2.57D0,  3.06D0,  3.52D0,  5.57D0,  4.40D0,
C         Ru       Rh       Pd       Ag       Cd
     .   4.81D0,  5.19D0,  5.54D0,  0.00D0,  0.00D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  5.76D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   3.13D0,  4.28D0,
C         La
     .   6.08D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.504D0, 5.585D0, 5.663D0, 5.740D0, 5.815D0, 5.889D0, 5.969D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.047D0, 6.126D0, 6.205D0, 6.284D0, 6.363D0, 6.444D0, 6.524D0,
C                  Hf       Ta       W        Re
     x            5.12D0,  5.30D0,  5.48D0,  5.88D0,
C         Os       Ir       Pt       Au       Hg
     x   5.88D0,  6.08D0,  6.29D0,  6.51D0,  0.0D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .  25*0.00D0/
C     ----------------------------------------------------------
      DATA WDD/
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   0.00D0,  0.00D0,
C         B        C        N        O        F        Ne
     .   0.00D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
     .   1.67D0,  2.43D0,
C         Al       Si       P        S        Cl       Ar
     .   3.42D0,  4.63D0,  6.09D0,  7.77D0,  9.68D0,  0.00D0,
C         K        Ca
     .   5.03D0,  6.03D0,
C         Sc       Ti       V        Cr       Mn
     .   7.02D0,  7.98D0,  8.91D0,  9.81D0, 10.68D0,
C         Fe       Co       Ni       Cu       Zn
     .  11.52D0, 12.32D0, 13.10D0, 13.84D0, 14.55D0,
C         Ga       Ge       As       Se       Br       Kr
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Rb       Sr
     .   2.48D0,  3.38D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ     WALTER THIEL Zr = 11.44
CMCZ .   4.28D0,  5.22D0,  6.12D0,  7.01D0,  7.88D0,
     .   4.28D0,  5.22D0,  6.12D0,  9.25D0,  7.88D0,
C         Ru       Rh       Pd       Ag       Cd
     .   8.73D0,  9.55D0, 10.36D0, 11.15D0, 11.93D0,
C       IODINE IS A GUESS
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   0.00D0,  0.00D0,
C         La
     .   6.911D0,
C**  .   6.700D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   7.073D0, 7.211D0, 7.332D0, 7.439D0, 7.534D0, 7.619D0, 7.638D0,
C**  .   6.930D0, 7.067D0, 7.190D0, 7.30D0,  7.400D0, 7.496D0, 7.543D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   7.646D0, 7.643D0, 7.629D0, 7.604D0, 7.568D0, 7.520D0, 7.459D0,
C**  .   7.575D0, 7.605D0, 7.628D0, 7.650D0, 7.650D0, 7.650D0, 7.650D0,
C                  Hf       Ta       W        Re
     x            8.65D0,  9.00D0,  9.45D0,  10.00D0,
C         Os       Ir       Pt       Au       Hg
     x  10.65D0, 11.40D0, 12.25D0, 13.20D0, 0.00D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   25*0.00D0/
C     ----------------------------------------------------------
C     THE D BLOCK FOR POLARIZATION FUNCTIONS
      DATA WDDP/
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   0.00D0,  0.00D0,
C         B        C        N        O        F        Ne
CWDD .   0.00D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
     .   0.00D0, 7.893D0, 8.532D0, 9.235D0, 9.946D0,  0.00D0,
C         Na       Mg
     .   1.67D0,  2.43D0,
C         Al       Si       P        S        Cl       Ar
     .   3.42D0,  4.63D0,  6.09D0,  7.77D0,  9.68D0,  0.00D0,
C         K        Ca
     .   5.03D0,  6.03D0,
C         Sc       Ti       V        Cr       Mn
     .   7.02D0,  7.98D0,  8.91D0,  9.81D0, 10.68D0,
C         Fe       Co       Ni       Cu       Zn
     .  11.52D0, 12.32D0, 13.10D0, 13.84D0, 14.55D0,
C         Ga       Ge       As       Se       Br       Kr
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Rb       Sr
     .   2.48D0,  3.38D0,
C         Y        Zr       Nb       Mo       Tc
     .   4.28D0,  5.22D0,  6.12D0,  7.01D0,  7.88D0,
C         Ru       Rh       Pd       Ag       Cd
     .   8.73D0,  9.55D0, 10.36D0, 11.15D0, 11.93D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Cs       Ba
     .   0.00D0,  0.00D0,
C         La
     .   6.911D0,
C**  .   6.700D0,
C**  .   6.930D0, 7.067D0, 7.190D0, 7.30D0,  7.400D0, 7.496D0, 7.543D0,
C**  .   7.575D0, 7.605D0, 7.628D0, 7.650D0, 7.650D0, 7.650D0, 7.650D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   7.073D0, 7.211D0, 7.332D0, 7.439D0, 7.534D0, 7.619D0, 7.638D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   7.646D0, 7.643D0, 7.629D0, 7.604D0, 7.568D0, 7.520D0, 7.459D0,
C                  Hf       Ta       W        Re
     x            8.65D0,  9.00D0,  9.45D0,  10.00D0,
C         Os       Ir       Pt       Au       Hg
     x  10.65D0, 11.40D0, 12.25D0, 13.20D0, 0.00D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   25*0.00D0/
C     ----------------------------------------------------------
      DATA WSD/
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   0.00D0,  0.00D0,
C         B        C        N        O        F        Ne
     .   0.00D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
C         Na       Mg
     .   2.35D0,  3.41D0,
C         Al       Si       P        S        Cl       Ar
     .   4.61D0,  5.92D0,  7.35D0,  8.85D0, 10.43D0,  0.00D0,
C         K        Ca
     .   4.00D0,  4.00D0,
C         Sc       Ti       V        Cr       Mn
     .   4.71D0,  5.38D0,  6.01D0,  6.60D0,  7.16D0,
C         Fe       Co       Ni       Cu       Zn
     .   7.68D0,  8.16D0,  8.61D0,  9.01D0,  9.39D0,
C         Ga       Ge       As       Se       Br       Kr
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Rb       Sr
     .   1.56D0,  1.97D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ     WALTER THIEL Zr = 6.94
CMCZ .   2.41D0,  2.93D0,  3.45D0,  3.99D0,  4.54D0,
     .   2.41D0,  2.93D0,  3.45D0,  6.38D0,  4.54D0,
C         Ru       Rh       Pd       Ag       Cd
     .   5.12D0,  5.71D0,  6.33D0,  6.96D0,  7.60D0,
C        IODINE IS A MERE GUESS
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  8.15D0,  0.00D0,
C         Cs       Ba
     .   0.00D0,  0.00D0,
C         La
     .   5.516D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.609D0, 5.696D0, 5.777D0, 5.856D0, 5.932D0, 6.005D0, 6.074D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.140D0, 6.203D0, 6.264D0, 6.322D0, 6.377D0, 6.429D0, 6.478D0,
C                  Hf       Ta       W        Re
     x            5.936D0, 6.142D0, 6.348, 6.577D0,
C         Os       Ir       Pt       Au       Hg
     x   6.806D0, 7.033D0, 7.271D0, 7.5204D0, 0.00D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   25*0.00D0/
C     ----------------------------------------------------------
      DATA WSDP/
C         H       He
     .  12.85D0,  0.00D0,
C         Li      Be
     .   0.00D0,  0.00D0,
C         B        C        N        O        F        Ne
CWSD .   0.00D0, 11.11D0, 12.01D0, 13.00D0, 14.00D0,  0.00D0,
     .   0.00D0, 7.893D0, 8.532D0, 9.235D0, 9.946D0,  0.00D0,
C         Na       Mg
     .   2.35D0,  3.41D0,
C         Al       Si       P        S        Cl       Ar
     .   4.61D0,  5.92D0,  7.35D0,  8.85D0, 10.43D0,  0.00D0,
C         K        Ca
     .   4.00D0,  4.00D0,
C         Sc       Ti       V        Cr       Mn
     .   4.71D0,  5.38D0,  6.01D0,  6.60D0,  7.16D0,
C         Fe       Co       Ni       Cu       Zn
     .   7.68D0,  8.16D0,  8.61D0,  9.01D0,  9.39D0,
C         Ga       Ge       As       Se       Br       Kr
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Rb       Sr
     .   1.56D0,  1.97D0,
C         Y        Zr       Nb       Mo       Tc
     .   2.41D0,  2.93D0,  3.45D0,  3.99D0,  4.54D0,
C         Ru       Rh       Pd       Ag       Cd
     .   5.12D0,  5.71D0,  6.33D0,  6.96D0,  7.60D0,
C         In       Sn       Sb       Te       I        Xe
     .   0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,  0.00D0,
C         Cs       Ba
     .   0.00D0,  0.00D0,
C         La
     .   5.516D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.609D0, 5.696D0, 5.777D0, 5.856D0, 5.932D0, 6.005D0, 6.074D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.140D0, 6.203D0, 6.264D0, 6.322D0, 6.377D0, 6.429D0, 6.478D0,
C                  Hf       Ta       W        Re
     x            5.936D0, 6.142D0, 6.348, 6.577D0,
C         Os       Ir       Pt       Au       Hg
     x   6.806D0, 7.033D0, 7.271D0, 7.5204D0, 0.00D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     X   25*0.00D0/
C     ----------------------------------------------------------
C     ----------------------------------------------------------
      DATA WSF/
C         H       He
C         Li      Be
C         B        C        N        O        F        Ne
C         Na       Mg
C         Al       Si       P        S        Cl       Ar
C         K        Ca
C         Sc       Ti       V        Cr       Mn
C         Fe       Co       Ni       Cu       Zn
C         Ga       Ge       As       Se       Br       Kr
C         Rb       Sr
C         Y        Zr       Nb       Mo       Tc
C         Ru       Rh       Pd       Ag       Cd
C         In       Sn       Sb       Te       I        Xe
C         Cs       Ba
     .   56*0.00D0,
C         La
     .   5.405D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   5.492D0, 5.580D0, 5.665D0, 5.748D0, 5.830D0, 5.911D0, 5.995D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   6.084D0, 6.172D0, 6.260D0, 6.350D0, 6.440D0, 6.532D0, 6.625D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   34*0.00D0/
C     ----------------------------------------------------------
      DATA WDF/
C         H       He
C         Li      Be
C         B        C        N        O        F        Ne
C         Na       Mg
C         Al       Si       P        S        Cl       Ar
C         K        Ca
C         Sc       Ti       V        Cr       Mn
C         Fe       Co       Ni       Cu       Zn
C         Ga       Ge       As       Se       Br       Kr
C         Rb       Sr
C         Y        Zr       Nb       Mo       Tc
C         Ru       Rh       Pd       Ag       Cd
C         In       Sn       Sb       Te       I        Xe
C         Cs       Ba
     .   56*0.00D0,
C         La
     .   8.674D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .   8.887D0, 9.100D0, 9.289D0, 9.459D0, 9.613D0, 9.754D0, 9.807D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .   9.860D0, 9.898D0, 9.916D0, 9.936D0, 9.935D0, 9.919D0, 9.888D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   34*0.00D0/
C     ----------------------------------------------------------
      DATA WFF/
C         H       He
C         Li      Be
C         B        C        N        O        F        Ne
C         Na       Mg
C         Al       Si       P        S        Cl       Ar
C         K        Ca
C         Sc       Ti       V        Cr       Mn
C         Fe       Co       Ni       Cu       Zn
C         Ga       Ge       As       Se       Br       Kr
C         Rb       Sr
C         Y        Zr       Nb       Mo       Tc
C         Ru       Rh       Pd       Ag       Cd
C         In       Sn       Sb       Te       I        Xe
C         Cs       Ba
     .   56*0.00D0,
C         La
     .  12.675D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     .  13.335D0,13.995D0,14.608D0,15.187D0,15.740D0,16.274D0,16.591D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .  17.093D0,17.582D0,18.060D0,18.528D0,18.988D0,19.440D0,19.885D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U
     .   34*0.00D0/
C     ----------------------------------------------------------
           DATA WSFP/105*0.0D0/
           DATA WPFP/105*0.0D0/
           DATA WDFP/105*0.0D0/
           DATA WFFP/105*0.0D0/
C     ----------------------------------------------------------
C
CJRY92
C PRINCIPLE QUANTUM NUMBER OF VALENCE SHELL ORBITALS
      DATA NPPP/2*1,8*2,8*3,18*4,18*5,32*6,19*7/
C FUNCTION STATEMENTS TO CALCULATE THE GAMMA FOR POLARIZATION FUNCTIONS
c     GAMPPP(K,NP,AMU)=WPP(K)*(2*NPPP(K)+1.D0)*AP(K,1)/((2*NP+1.D0)*AMU)
c     GAMDDP(K,NP,AMU)=WDD(K)*(2*NPPP(K)+1.D0)*AD(K,1)/((2*NP+1.D0)*AMU)
c     GAMFFP(K,NP,AMU)=WFF(K)*(2*NPPP(K)+1.D0)*AF(K,1)/((2*NP+1.D0)*AMU)
c     GAMSPP(K,NP,AMU)=WSP(K)*(2*NPPP(K)+1.D0)*AP(K,1)/((2*NP+1.D0)*AMU)
c     GAMSDP(K,NP,AMU)=WSD(K)*(2*NPPP(K)+1.D0)*AD(K,1)/((2*NP+1.D0)*AMU)
c     GAMSFP(K,NP,AMU)=WSF(K)*(2*NPPP(K)+1.D0)*AF(K,1)/((2*NP+1.D0)*AMU)
c     GAMDFP(K,NP,AMU)=WDF(K)*(2*NPPP(K)+1.D0)*AF(K,1)/((2*NP+1.D0)*AMU)
      COUL = .false.
      WOL = 0.0d0
      WOK = 0.0d0
CJRY92
C
c     WRITE (6,*)"#",M1,M2,K,L,R,LA,LB,AMU,BMU,NP1,NP2,"#"
C     WRITE (6,*) FCOUL(1,M1),FCOUL(1,M2)

      IZ=0
      FAC=1.2
C     FAC = SQRT(SQRT(ONE+0.40*R/BOHRS))
C     E=FAC*BOHRS*EVOLT   
      E=FAC*14.3986D0   
C
C     FLAGS KD,LD:  0 FOR S,P; 1 FOR D, 2 FOR F.
C
      IF (LA .LE. 0) THEN
         KD=0
      ELSEIF (LA .LE. 3) THEN
         KD = 1
      ELSEIF (LA .LE. 8) THEN
         KD = 2
      ELSE
         KD = 3
      ENDIF
      IF (LB .LE. 0) THEN
         LD = 0
      ELSEIF (LB .LE. 3) THEN
         LD = 1
      ELSEIF (LB .LE. 8) THEN
         LD = 2
      ELSE
         LD = 3
      ENDIF
      IK=0
      IL=0
C     THE ASSUMPTION HERE IS THAT EXTRA ORBITALS ARE RYDBERG LIKE -
C     NOT POLARIZATION LIKE.
      IF(K.NE.0) GO TO 10
    7 WOK=RRRR(IZ,IZ,IZ,NP1,IZ,IZ,AMU,NP1,IZ,IZ,AMU,NP1,IZ,IZ,AMU,
     1 NP1,IZ,IZ,AMU,FACT)*EVOLT
      IK=1
   10 CONTINUE
      IF(L.NE.0) GO TO 20
   17 WOL=RRRR(IZ,IZ,IZ,NP2,IZ,IZ,BMU,NP2,IZ,IZ,BMU,NP2,IZ,IZ,BMU,
     1 NP2,IZ,IZ,BMU,FACT)*EVOLT
      IL=1
   20 CONTINUE
      IF (ABS(R) .LT. 1.0D-2) GOTO 200
   90 IF (IK .EQ. 1) GOTO 110
      IF (KD .EQ. 0) THEN
c        IF(FCOUL(1,M1).LT.1.0D-8) THEN
           WOK = WSS(K)
c        ELSE
c          WOK=FCOUL(1,M1)
c        ENDIF
      ELSEIF (KD .EQ. 1) THEN
c        IF(FCOUL(2,M1).LT.1.0D-8) THEN
           WOK = WPP(K)
c        ELSE
c          WOK=FCOUL(2,M1)
c        ENDIF
C        IF(NEXTA.EQ.1) WOK = WPPP(K)
c        IF(NEXTA.EQ.1) WOK = GAMPPP(K,NP1,AMU)
      ELSEIF (KD .EQ. 2) THEN
c        IF(FCOUL(3,M1).LT.1.0D-8) THEN
           WOK = WDD(K)
c        ELSE
c           WOK = FCOUL(3,M1)
c        ENDIF
C        IF(NEXTA.EQ.1) WOK = WDDP(K)
c        IF(NEXTA.EQ.1) WOK = GAMDDP(K,NP1,AMU)
      ELSEIF (KD .EQ. 3) THEN
c        IF(FCOUL(4,M1).LT.1.0D-8) THEN
           WOK = WFF(K)
c        ELSE
c           WOK = FCOUL(4,M1)
c        ENDIF
C        IF(NEXTA.EQ.1) WOK = WFFP(K)
c        IF(NEXTA.EQ.1) WOK = GAMFFP(K,NP1,AMU)
      ENDIF
  110 IF (IL .EQ. 1) GOTO 300
      IF (LD .EQ. 0) THEN
c        IF(FCOUL(1,M2).LT.1.0D-8) THEN
           WOL = WSS(L)
c        ELSE
c          WOL=FCOUL(1,M2)
c        ENDIF
      ELSEIF (LD .EQ. 1) THEN
c        IF(FCOUL(2,M2).LT.1.0D-8) THEN
           WOL = WPP(L)
c        ELSE
c          WOL=FCOUL(2,M2)
c        ENDIF
C        IF(NEXTB.EQ.1) WOL = WPPP(L)
c        IF(NEXTB.EQ.1) WOL = GAMPPP(L,NP2,BMU)
      ELSEIF (LD .EQ. 2) THEN
c        IF(FCOUL(3,M2).LT.1.0D-8) THEN
           WOL = WDD(L)
c        ELSE
c           WOL = FCOUL(3,M2)
c        ENDIF
C        IF(NEXTB.EQ.1) WOL = WDDP(L)
c        IF(NEXTB.EQ.1) WOL = GAMDDP(L,NP2,BMU)
      ELSEIF (LD .EQ. 3) THEN
c        IF(FCOUL(4,M2).LT.1.0D-8) THEN
           WOL = WFF(L)
c        ELSE
c           WOL = FCOUL(4,M2)
c        ENDIF
C        IF(NEXTB.EQ.1) WOL = WFFP(L)
c        IF(NEXTB.EQ.1) WOL = GAMFFP(L,NP2,BMU)
      ENDIF
      GO TO 300
C
C   ONE CENTER INTEGRAL
C
  200 CONTINUE
      NTEST = MAX0(KD,LD)+1
      MTEST = MIN0(KD,LD)+1
      NTEST = (NTEST*(NTEST-1))/2 + MTEST
      IF (KD .EQ. LD) GOTO 90
      IF (NTEST .EQ. 2) THEN
         IF (IK .EQ. 0) THEN
           WOK = WSP(K)
c          IF(FCOUL(1,M1).GT.1.0D-8.AND.FCOUL(2,M1).GT.1.0D-8) THEN
c            WOK = SQRT (FCOUL(1,M1)*FCOUL(2,M1))
c          ELSEIF (FCOUL(1,M1).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(1,M1)*WPP(K))
c          ELSEIF (FCOUL(2,M1).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(2,M1)*WSS(K))
c          ENDIF
C          IF(NEXTA.EQ.1) WOK = WSPP(K)
c          IF(NEXTA.EQ.1) WOK = GAMSPP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WSP(L)
c          IF(FCOUL(1,M2).GT.1.0D-8.AND.FCOUL(2,M2).GT.1.0D-8) THEN
c            WOL = SQRT (FCOUL(1,M2)*FCOUL(2,M2))
c          ELSEIF (FCOUL(1,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(1,M2)*WPP(L))
c          ELSEIF (FCOUL(2,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(2,M2)*WSS(L))
c          ENDIF
C          IF(NEXTB.EQ.1) WOL = WSPP(L)
c          IF(NEXTB.EQ.1) WOL = GAMSPP(L,NP2,BMU)
         ENDIF
      ELSEIF (NTEST .EQ. 4) THEN
         IF (IK .EQ. 0) THEN
           WOK = WSD(K)
c          IF(FCOUL(1,M1).GT.1.0D-8.AND.FCOUL(3,M1).GT.1.0D-8) THEN
c            WOK = SQRT (FCOUL(1,M1)*FCOUL(3,M1))
c          ELSEIF (FCOUL(1,M2).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(1,M1)*WDD(K))
c          ELSEIF (FCOUL(2,M2).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(3,M1)*WSS(K))
c          ENDIF
C          IF(NEXTA.EQ.1) WOK = WSDP(K)
c          IF(NEXTA.EQ.1) WOK = GAMSDP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WSD(L)
c          IF(FCOUL(1,M2).GT.1.0D-8.AND.FCOUL(3,M2).GT.1.0D-8) THEN
c            WOL = SQRT (FCOUL(1,M2)*FCOUL(3,M2))
c          ELSEIF (FCOUL(1,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(1,M2)*WDD(L))
c          ELSEIF (FCOUL(2,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(3,M2)*WSS(L))
c          ENDIF
C          IF(NEXTB.EQ.1) WOL = WSDP(L)
c          IF(NEXTB.EQ.1) WOL = GAMSDP(L,NP2,BMU)
         ENDIF
      ELSEIF (NTEST .EQ. 5) THEN
         IF (IK .EQ. 0) THEN
           WOK = WSD(K)
c          IF(FCOUL(1,M1).GT.1.0D-8.AND.FCOUL(3,M1).GT.1.0D-8) THEN
c            WOK = SQRT (FCOUL(1,M1)*FCOUL(3,M1))
c          ELSEIF (FCOUL(1,M1).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(1,M1)*WDD(K))
c          ELSEIF (FCOUL(2,M1).GT.1.0D-8) THEN
c            WOK    = SQRT(FCOUL(3,M1)*WSS(K))
c          ENDIF
C          IF(NEXTA.EQ.1) WOK = WSDP(K)
c          IF(NEXTA.EQ.1) WOK = GAMSDP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WSD(L)
c          IF(FCOUL(1,M2).GT.1.0D-8.AND.FCOUL(3,M2).GT.1.0D-8) THEN
c            WOL = SQRT (FCOUL(1,M2)*FCOUL(3,M2))
c          ELSEIF (FCOUL(1,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(1,M2)*WDD(L))
c          ELSEIF (FCOUL(2,M2).GT.1.0D-8) THEN
c            WOL    = SQRT(FCOUL(3,M2)*WSS(L))
c          ENDIF
C          IF(NEXTB.EQ.1) WOL = WSDP(L)
c          IF(NEXTB.EQ.1) WOL = GAMSDP(L,NP2,BMU)
         ENDIF
      ELSEIF (NTEST .EQ. 7) THEN
         IF (IK .EQ. 0) THEN
           WOK = WSF(K)
C          IF(NEXTA.EQ.1) WOK = WSFP(K)
c          IF(NEXTA.EQ.1) WOK = GAMSFP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WSF(L)
C          IF(NEXTB.EQ.1) WOL = WSFP(L)
c          IF(NEXTB.EQ.1) WOL = GAMSFP(L,NP2,BMU)
         ENDIF
      ELSEIF (NTEST .EQ. 8) THEN
         IF (IK .EQ. 0) THEN
           WOK = WSF(K)
C          IF(NEXTA.EQ.1) WOK = WSFP(K)
c          IF(NEXTA.EQ.1) WOK = GAMSFP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WSF(L)
C          IF(NEXTB.EQ.1) WOL = WSFP(L)
c          IF(NEXTB.EQ.1) WOL = GAMSFP(L,NP2,BMU)
         ENDIF
      ELSEIF (NTEST .EQ. 9) THEN
         IF (IK .EQ. 0) THEN
           WOK = WDF(K)
C          IF(NEXTA.EQ.1) WOK = WDFP(K)
c          IF(NEXTA.EQ.1) WOK = GAMDFP(K,NP1,AMU)
         ENDIF
         IF (IL .EQ. 0) THEN
           WOL = WDF(L)
C          IF(NEXTB.EQ.1) WOL = WDFP(L)
c          IF(NEXTB.EQ.1) WOL = GAMDFP(L,NP2,BMU)
         ENDIF
      ENDIF
C
  300 CONTINUE
      IF ((WOK.LT.0.01).OR.(WOL.LT.1.0D-2)) THEN
         WRITE(6,1995) K,KD,L,LD
         CALL ZEXIT(3)
      ENDIF
 2000 IF(ISPP.EQ.5) THEN
C     FORM SUGGESTED BY WARSHEL FIXED BY MCZ
        RR=R!/BOHRS
c       write(6,*)"RR",RR
        F=1.0D0/((WOK+WOL)/54.42D0-0.05285D0*(AMU+BMU))
        S=0.05285D0*(AMU+BMU)*EXP(-0.03785D0*RR*(AMU+BMU))+1.D0/(F+RR)
        TAG=S
C     OHNO-KLOPPMAN DEWAR
      ELSEIF (ISPP.EQ.3) THEN
        RR=R/BOHRS
        SUM=(EVOLT/WOK+EVOLT/WOL)**2/FOUR
        TAG=ONE/SQRT(RR*RR+SUM)
      ELSEIF (ISPP.EQ.1) THEN
C     MATAGA-NISHIMOTO WEISS
c      write(6,*)"R",R*BOHRS,wok,wol
        RR = R*BOHRS
        SUM=(WOK+WOL)/2.0D0
c       write(6,*)"SUM",SUM
        Z=E/SUM
        S=E/(RR+Z)
        TAG=S/EVOLT
      ENDIF
c     write(6,*)"COUL",TAG
      RETURN
C
C   FORMAT
C
 1995 FORMAT('**** Integral in Coul missing for atomic no. ',I5,
     .       ', L=',I5,' or atomic no. ',I5,',L=',I5)
C
C
C
        END

