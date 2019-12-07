      SUBROUTINE PARA(NA,KAN,NTYP,EIS,EIP,EID,EIF,ISW1,ISW2,COFSQ2)
C
C     INSERTS THE DIAGONAL EIS, EIP, EID AND EIF.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      parameter (MPE =105,MXB=4000,ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/factor.cmn'
C
c     INCLUDE '../include/cf.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
C
      DIMENSION KAN(NA),NTYP(NA)
C    +     BETA(MPE,10), BETB(MPE, 4)
      DIMENSION
     +     EIS(MPE),   EIP(MPE),   EID(MPE),   EIF(MPE),
     +     EIS1(MPE),  EIP1(MPE),  EID1(MPE),  EIF1(MPE),
     +     EIS2(MPE),  EIP2(MPE),  EID2(MPE),  EIF2(MPE),
     +     EIS3(MPE),              EID3(MPE),  EIF3(MPE),
     +     EIS4(MPE),  EIF4(MPE),
     +     EISS1(MPE), EISS2(MPE), EIPP1(MPE) ,EIPP2(MPE),
     +     EIDD1(MPE), EIDD2(MPE), EIDD3(MPE),COFSQ2(MPE,5)
C
C     CNDO/1 AND INDO/1  : IONIZATION POTENTIALS
C
C
C     FROM D(N-1) S(1)
C
      DATA (EIS1(JJJ), JJJ = 1, 103)/
C         H         He        Li        Be
     S -13.06D0,   0.00D0,  -5.41D0,  -9.33D0,
C         B         C         N         O         F        Ne
     P -14.00D0, -19.42D0, -25.58D0, -32.49D0, -40.14D0,  0.00D0,
CJR  P -14.00D0, -19.42D0, -25.58D0, -30.39D0, -40.14D0,   0.0D0,
C         Na        Mg
     S  -4.86D0,  -8.11D0,
C         Al        Si        P         S         Cl       Ar
     P -11.42D0, -14.79D0, -18.23D0, -21.73D0, -25.29D0,  0.00D0,
C         K         Ca
     S  -4.34D0,  -5.13D0,
C         Sc        Ti        V         Cr        Mn
     D  -5.83D0,  -6.34D0,  -6.71D0,  -6.97D0,  -7.15D0,
C         Fe        Co        Ni        Cu        Zn
     D  -7.27D0,  -7.38D0,  -7.51D0,  -7.69D0,  -9.36D0,
C         Ga        Ge        As        Se        Br        Kr
     P -12.612D0,-1.00D2, -18.18D0, -20.95D0, -23.94D0,  -1.00D2,
C         Rb        Sr
     S  -4.43D0,  -5.19D0,
C         Y         Zr        Nb        Mo        Tc
     D  -5.85D0,  -6.40D0,  -6.84D0,  -7.18D0,  -7.42D0,
C         Ru        Rh        Pd        Ag        Cd
     D  -7.55D0,  -7.58D0,  -7.50D0,  -7.31D0,  -8.94D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2, -14.17D0,  -1.00D2,  -1.00D2, -20.84D0,  -1.00D2,
C         Cs        Ba
     S  -3.893D0,-4.646D0,
CCC**F    -5.9784D0,-5.6583D0,-5.6073D0,-5.5788D0,-5.5721D0,-5.5872D0,
CCC**F    -5.6240D0,-5.6833D0,-5.7628D0,-5.8648D0,-5.9886D0,-6.1341D0,
CCC**F    -6.3014D0,-6.4904D0,-7.0356D0,
C                                 LANTHAIDE IP FROM F(N-3) D S(2)
C         La
     F  -5.9784D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -5.93D0, -5.89D0, -5.93D0, -5.98D0, -6.04D0, -6.10D0, -6.17D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -6.25D0, -6.34D0, -6.42D0, -6.52D0, -6.62D0, -6.73D0, -6.77D0,
C         IP-S from d(n-1)S
C         Hf       Ta       W       Re       Os
     D -7.42D0,  -7.49D0, -7.66D0, -7.83D0, -8.09D0,
C         Ir       Pt       Au      Hg
     D -8.41D0,  -8.78D0, -9.22D0, -90.0D0,
     F   7* -1.00D2,
C              Ra        Ac       (ACTINIDE IP FROM 5F(N-3) 6D 7S(2))
     F       -1.00D2,  -5.86D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -6.00D0, -6.13D0, -6.26D0, -6.38D0, -6.48D0, -6.58D0, -6.67D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -6.75D0, -6.82D0, -6.88D0, -6.94D0, -6.98D0, -7.02D0, -7.05D0/
C
      DATA (EIP1(JJJ), JJJ = 1, 103)/
C         H         He        Li        Be
     +  -3.40D0,   0.00D0,  -3.61D0,  -5.88D0,
C         B         C         N         O         F         Ne
     P  -8.24D0, -10.70D0, -13.25D0, -15.88D0, -18.61D0,   0.00D0,
CJR  P  -8.24D0, -10.70D0, -13.25D0, -17.28D0, -18.61D0,   0.0D0,
C         Na        Mg
     S  -2.86D0,  -4.55D0,
C         Al        Si        P         S         Cl        Ar
     P  -6.29D0,  -8.10D0,  -9.98D0, -11.92D0, -13.93D0,   0.00D0,
C         K        Ca
     S  -2.73D0,  -2.99D0,
C     FROM D(N-1) P(1)
C         Sc        Ti        V         Cr        Mn
     D  -3.43D0,  -3.75D0,  -3.95D0,  -4.06D0,  -4.10D0,
C         Fe        Co        Ni        Cu        Zn
     D  -4.08D0,  -4.02D0,  -3.93D0,  -3.84D0,  -4.77D0,
C         Ga        Ge        As        Se        Br        Kr
     P  -5.93D0,  -1.00D2,  -9.19D0, -10.37D0, -12.44D0,  -1.00D2,
C         Rb        Sr
     S  -2.65D0,  -3.16D0,
C         Y         Zr        Nb        Mo        Tc
     D  -3.58D0,  -3.92D0,  -4.15D0,  -4.29D0,  -4.34D0,
C         Ru        Rh        Pd        Ag        Cd
     D  -4.30D0,  -4.17D0,  -3.94D0,  -3.62D0,  -4.75D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2,  -7.012D0, -1.00D2,  -1.00D2, -11.21D0,  -1.00D2,
C         Cs        Ba
     S  -2.461D0,  -3.434D0,
CCC**F-4.195D0,-3.30D0,-3.30D0,-3.20D0,-3.20D0,-3.20D0,-3.20D0,-3.30D0,
CCC**F-3.30D0,-3.30D0,-3.40D0,-3.40D0,-3.50D0,-3.60D0,-3.70D0,
C                                 LANTHAIDE IP FROM F(N-2) S P
C         La
     F  -4.195D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -3.69D0, -3.76D0, -3.82D0, -3.87D0, -3.91D0, -3.95D0, -3.99D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -4.03D0, -4.05D0, -4.08D0, -4.10D0, -4.12D0, -4.13D0, -3.70D0,
C         IP-P From d(N-1)P(1)
C         Hf       Ta       W       Re       Os
     D  -5.42D0, -5.67D0, -5.77D0, -5.72D0, -5.52D0,
C         Ir       Pt       Au      Hg
     D -5.18D0,  -4.69D0, -4.05D0, -90.0D0,
     F   7* -1.00D2,
C         Ra           Ac       (ACTINIDE IP FROM 5F(N-3) 7P 7S(2))
     F  -1.0D2,  -3.97D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -4.11D0, -4.25D0, -4.39D0, -4.51D0, -4.62D0, -4.73D0, -4.83D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -4.92D0, -5.00D0, -5.07D0, -5.13D0, -5.19D0, -5.23D0, -5.27D0/
C
      DATA (EID1(JJJ), JJJ = 1, 103)/
     S   -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,
     P   -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,
C         Na        Mg
     S  -1.51D0,  -1.76D0,
C         Al        Si        P         S         Cl        Ar
     P  -1.94D0,  -2.06D0,  -2.11D0,  -2.09D0,  -2.01D0,   0.00D0,
C         K         Ca
     S  -1.67D0,  -3.44D0,
C         Sc        Ti        V         Cr        Mn
     D  -4.85D0,  -5.93D0,  -6.77D0,  -7.43D0,  -7.99D0,
C         Fe       Co       Ni       Cu       Zn
     D  -8.53D0,  -9.10D0,  -9.79D0, -10.67D0, -17.57D0,
C         Ga        Ge        As        Se        Br        Kr
     P  -2.30D0,  -1.00D2,  -2.17D0,  -2.10D0,  -2.05D0,  -1.00D2,
C         Rb        Sr
     S  -1.00D2,  -3.66D0,
C         Y         Zr        Nb        Mo        Tc
     D  -4.74D0,  -5.79D0,  -6.82D0,  -7.82D0,  -8.80D0,
C         Ru        Rh        Pd        Ag        Cd
     D  -9.77D0, -10.70D0, -11.62D0, -12.51D0, -17.68D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -2.32D0,  -1.00D2,
C         Cs        Ba
     S  -2.0882D0,-3.982D0,
CCC**F     -6.4153D0,-6.9383D0,-7.0611D0,-7.1462D0,-7.1936D0,-7.2034D0,
CCC**F     -7.1755D0,-7.1104D0,-7.0066D0,-6.8656D0,-6.687D0,-6.4707D0,
CCC**F     -6.2167D0,-5.925D0,-5.2774D0,
C                                 LANTHAIDE IP FROM F(N-3) D S(2)
C         La
     F  -6.4153D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -6.74D0, -6.77D0, -6.77D0, -6.74D0, -6.73D0, -6.68D0, -6.61D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -6.54D0, -6.42D0, -6.30D0, -6.16D0, -6.10D0, -5.96D0, -5.31D0,
C        IP-d from D(N-1)S(1)
C         Hf       Ta       W       Re       Os
     D  -5.44D0, -6.47D0, -7.46D0, -8.41D0, -9.30D0,
C         Ir       Pt       Au      Hg
     D  -10.15D0,-10.97D0,-11.74D0,-90.0D0,
     F   7* -1.00D2,
C         Ra           Ac       (ACTINIDE IP FROM 5F(N-3) 6D 7S(2))
     F  -1.00D2,  -5.15D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -5.42D0, -5.46D0, -5.80D0, -5.90D0, -5.95D0, -5.93D0, -5.85D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -5.71D0, -5.51D0, -5.25D0, -4.93D0, -4.55D0, -4.12D0, -3.62D0/
C
      DATA (EIF1(JJJ), JJJ = 1, 103)/
     +  56 *   -1.00D2,
CCC**F-11.5064D0,-12.5819D0,-13.5711D0,-14.4739D0,-15.2905D0,-16.0207D0,
CCC**F-16.6646D0,-17.2222D0,-17.6935D0,-18.0785D0,-18.3771D0,-18.5895D0,
CCC**F-18.7155D0,-18.7552D0,-18.7086D0,
C                                 LANTHAIDE IP FROM F(N-3) D S(2)
C         La
     F -11.5064D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F -12.17D0,-12.79D0,-13.35D0,-13.85D0,-14.27D0,-14.62D0,-14.90D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F -15.13D0,-15.27D0,-15.35D0,-15.36D0,-15.30D0,-15.18D0,-18.60D0,
C         Hf       Ta       W       Re       Os
C        No F oORBITALS ON 3RD TRANSION SERIES
     D -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
C         Ir       Pt       Au      Hg
     D -100.0D0,-100.0D0,-100.0D0,-100.0D0,
     F   7* -1.00D2,
C          Ra          Ac       (ACTINIDE IP FROM 5F(N-3) 6D 7S(2)
     F   -1.00D2, -2.80D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -6.16D0, -7.36D0, -8.43D0, -9.40D0,-10.25D0,-11.00D0,-11.63D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F -12.15D0,-12.56D0,-12.86D0,-13.05D0,-13.12D0,-13.09D0,-12.94D0/
C
C     PARAMETERS FROM D(N-2) S(2)
C
      DATA (EIS2(JJJ), JJJ = 1, 103)/
C         H         He        Li        Be
     + -13.06D0,   0.00D0,  -5.41D0,  -9.33D0,
C         B         C         N         O         F         Ne
     P -14.00D0, -19.42D0, -25.58D0, -32.49D0, -40.14D0,   0.00D0,
CJR  P -14.00D0, -19.42D0, -25.58D0, -30.39D0, -40.14D0,   0.0D0,
C         Na        Mg
     S  -4.86D0,  -8.11D0,
C         Al        Si        P         S         Cl        Ar
     P -11.42D0, -14.79D0, -18.23D0, -21.73D0, -25.29D0,   0.0D0,
C         K         Ca
     S  -4.34D0,  -6.03D0,
C         Sc        Ti        V         Cr        Mn
     D  -6.72D0,  -7.28D0,  -7.73D0,  -8.07D0,  -8.35D0,
C         Fe        Co        Ni        Cu        Zn
     D  -8.57D0,  -8.76D0,  -8.94D0,  -9.13D0,  -9.36D0,
C         Ga        Ge        As        Se        Br        Kr
     P -12.612D0, -1.00D2, -18.18D0, -20.95D0, -23.94D0,  -1.00D2,
C         Rb        Sr
     S  -4.43D0,  -5.84D0,
C         Y         Zr        Nb        Mo        Tc
     D  -6.55D0,  -7.16D0,  -7.69D0,  -8.13D0,  -8.49D0,
C         Ru        Rh        Pd        Ag        Cd
     D  -8.75D0,  -8.93D0,  -9.02D0,  -9.03D0,  -8.94D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2, -14.17D0,  -1.00D2,  -1.00D2, -20.84D0,  -1.00D2,
C         Cs        Ba
     S  -3.893D0, -5.210D0,
CCC**F    -5.9784D0,-5.6583D0,-5.6073D0,-5.5788D0,-5.5721D0,-5.5872D0,
CCC**F    -5.6240D0,-5.6833D0,-5.7628D0,-5.8648D0,-5.9886D0,-6.1341D0,
CCC**F    -6.3014D0,-6.4904D0,-7.0356D0,
C                                 LANTHAIDE IP FROM F(N-3)D S(2)
C         La
     F  -5.9784D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -5.93D0, -5.89D0, -5.93D0, -5.98D0, -6.04D0, -6.10D0, -6.17D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -6.25D0, -6.34D0, -6.42D0, -6.52D0, -6.62D0, -6.73D0, -6.77D0,
C        IP-S FROM D(N-2)S(2)
C         Hf       Ta       W       Re       Os
     D  -7.93D0, -8.34D0, -8.72D0, -9.06D0, -9.39D0,
C         Ir       Pt       Au      Hg
     D  -9.67D0, -9.93D0, -10.15D0,-10.35D0,
     F   7* -1.00D2,
C         Ra           Ac       (ACTINIDE IP FROM 5F(N-3) 7P 7S(2))
     F   -1.00D2, -6.71D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -6.95D0, -7.16D0, -7.36D0, -7.54D0, -7.70D0, -7.84D0, -7.97D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -8.07D0, -8.16D0, -8.23D0, -8.28D0, -8.32D0, -8.33D0, -8.33D0/
C
      DATA (EIP2(JJJ), JJJ = 1, 103)/
C         H         He        Li        Be
     +  -3.40D0,   0.00D0,  -3.61D0,  -5.88D0,
C         B         C         N         O         F         Ne
     P  -8.24D0, -10.70D0, -13.25D0, -15.88D0, -18.61D0,   0.00D0,
CJR  P  -8.24D0, -10.70D0, -13.25D0, -17.28D0, -18.61D0,   0.0D0,
C         Na        Mg
     S  -2.86D0,  -4.55D0,
C         Al        Si        P         S         Cl        Ar
     P  -6.29D0,  -8.10D0,  -9.98D0, -11 .92D0, -13.93D0,   0.00D0,
C         K         Ca
     S  -2.73D0,  -3.96D0,
C     PARAMETERS FROM D(N-2) S(1) P(1)
C         Sc        Ti        V         Cr        Mn
     D  -4.20D0,  -4.48D0,  -4.77D0,  -5.04D0,  -5.27D0,
C         Fe        Co        Ni        Cu        Zn
     D  -5.42D0,  -5.48D0,  -5.41D0,  -5.18D0,  -4.77D0,
C         Ga        Ge        As        Se        Br        Kr
     P  -5.93D0,  -1.00D2,  -9.19D0, -10.37D0, -12.44D0,  -1.00D2,
C         Rb        Sr
     S  -2.65D0,  -3.76D0,
C         Y         Zr        Nb        Mo        Tc
     D  -4.13D0,  -4.43D0,  -4.67D0,  -4.87D0,  -5.00D0,
C         Ru        Rh        Pd        Ag        Cd
     D  -5.07D0,  -5.07D0,  -5.02D0,  -4.91D0,  -4.75D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2,  -7.012D0, -1.00D2,  -1.00D2, -11.21D0,  -1.00D2,
C         Cs        Ba
     S  -2.461D0, -2.891D0,
CCC**F-4.195D0,-3.30D0,-3.30D0,-3.20D0,-3.20D0,-3.20D0,-3.20D0,-3.30D0,
CCC**F-3.30D0,-3.30D0,-3.40D0,-3.40D0,-3.50D0,-3.60D0,-3.70D0,
C                                 LANTHAIDE IP FROM F(N-2) S P
C         La
     F  -4.195D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -3.69D0, -3.76D0, -3.82D0, -3.87D0, -3.91D0, -3.95D0, -3.99D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -4.03D0, -4.05D0, -4.08D0, -4.10D0, -4.12D0, -4.13D0, -3.70D0,
C        IP-P FROM D(N-2)S(1)P(1)
C         Hf       Ta       W       Re       Os
     D  -4.96D0, -5.47D0, -5.81D0, -6.01D0, -6.08D0,
C         Ir       Pt       Au      Hg
     D  -5.98D0, -5.70D0, -5.33D0, -4.05D0,
     F   7* -1.00D2,
C            Ra         Ac       (ACTINIDE IP FROM 5F(N-3) 7P 7S(2))
     F     -1.00D2,  -3.97D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -4.11D0, -4.25D0, -4.39D0, -4.51D0, -4.62D0, -4.73D0, -4.83D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -4.92D0, -5.00D0, -5.07D0, -5.13D0, -5.19D0, -5.23D0, -5.27D0/
C
CJiabo
C
c     CALL INDOS2(EIP2(8),BETA(8,1),EIS2(8))
CLi
      DATA (EID2(JJJ), JJJ = 1, 103)/
     S   -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,
     P   -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,
C         Na        Mg
     S  -1.51D0,  -1.76D0,
C         Al        Si        P         S         Cl        Ar
     P  -1.94D0,  -2.06D0,  -2.11D0,  -2.09D0,  -2.01D0,   0.00D0,
C         K         Ca
     S  -1.67D0,  -3.44D0,
C         Sc        Ti        V         Cr        Mn
     D  -8.16D0,  -9.07D0,  -9.89D0, -10.66D0, -11.45D0,
C         Fe        Co        Ni        Cu        Zn
     D -12.31D0, -13.30D0, -14.46D0, -15.87D0, -17.57D0,
C         Ga        Ge        As        Se        Br        Kr
     P  -2.30D0,  -1.00D2,  -2.17D0,  -2.10D0,  -2.05D0,  -1.00D2,
C         Rb        Sr
     S  -1.00D2,   -3.66D0,
C         Y         Zr        Nb        Mo        Tc
     D  -6.61D0,  -8.08D0,  -9.50D0, -10.85D0,  -12.13D0,
C         Ru        Rh        Pd        Ag       Cd
     D -13.38D0, -14.54D0, -15.65D0, -16.70D0, -17.68D0,
C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -2.32D0,  -1.00D2,
C         Cs        Ba
     S  -2.0882D0,-3.982D0,
CCC**F   -6.4153D0,-6.9383D0,-7.0611D0,-7.1462D0,-7.1936D0,-7.2034D0,
CCC**F   -7.1755D0,-7.1104D0,-7.0066D0,-6.8656D0,-6.6870D0,-6.4707D0,
CCC**F   -6.2167D0,-5.9250D0,-5.2774D0,
C                                 LANTHAIDE IP FROM F(N-3) D S(2)
C         La
     F  -6.4153D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -6.74D0, -6.77D0, -6.77D0, -6.74D0, -6.73D0, -6.68D0, -6.61D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -6.54D0, -6.42D0, -6.30D0, -6.16D0, -6.10D0, -5.96D0, -5.31D0,
C        IP-D FROM D(N-2)S(2)
C         Hf       Ta       W       Re       Os
     D  -6.67D0, -7.29D0, -8.06D0, -8.96D0, -10.02D0,
C         Ir       Pt       Au      Hg
     D  -11.21D0,-12.55D0,-14.02D0,-15.65D0,
     F   7* -1.00D2,
C         Ra           Ac       (ACTINIDE IP FROM 5F(N-4) 6D(2) 7S(2))
     F     -1.00D2,  -100D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -5.93D0, -6.45D0, -6.77D0, -6.88D0, -6.78D0, -6.47D0,  -100D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F   -100D0,  -100D0,  -100D0,  -100D0,  -100D0,  -100D0,  -100D0/
C
C   PARAMETERS FROM F(N-2)S(2)
C
      DATA (EIF2(JJJ), JJJ = 1, 103)/
     + 56 *   -1.00D2,
CCC**F-11.5064D0,-12.5819D0,-13.5711D0,-14.4739D0,-15.2905D0,-16.0207D0,
CCC**F-16.6646D0,-17.2222D0,-17.6935D0,-18.0785D0,-18.3771D0,-18.5895D0,
CCC**F-18.7155D0,-18.7552D0,-18.7086D0,
C                                 LANTHAIDE IP FROM F(N-2) S(2)
C         La
     F -11.5064D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     F  -7.24D0, -7.27D0, -7.29D0, -7.31D0, -7.35D0, -7.39D0,  -7.44D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     F  -7.49D0, -7.55D0, -7.61D0, -7.69D0, -7.76D0, -7.84D0, -18.60D0,
C         Hf       Ta       W       Re       Os
C      NO F ORBITALS FOR THIRD TRANSION SERIES
     D -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
C         Ir       Pt       Au      Hg
     D -100.0D0,-100.0D0,-100.0D0,-100.0D0,
     F   7* -1.00D2,
C                  (ACTINIDE IP FROM 5F(N-3) 7P 7S(2))
C         Ra           Ac (FROM F(N-2)S(2)
     F      -1.00D2, -2.80D2,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -7.09D0, -9.01D0,-10.65D0,-12.00D0,-13.07D0,-13.86D0,-14.37D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F -14.59D0,-14.54D0,-14.20D0,-14.73D0,-15.44D0,-16.15D0,-16.86D0/
C
C     PARAMETERS FROM D(N) S(0)
C
      DATA (EIS3(JJJ), JJJ = 1, 103)/
     . 88*-100.0D0,
C             Ac       (ACTINIDE IP FROM 5F(N-4) 6D(2) 7S(2))
     F   -100.0D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F   -6.26D0, -6.42D0, -6.50D0, -6.51D0, -6.45D0, -6.31D0, -100.D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F   -100.D0, -100.D0, -100.D0, -100.D0, -100.D0, -100.D0, -100.D0/
C
C
      DATA (EID3(JJJ), JJJ = 1, 103)/
     .   36 *   -1.00D2,
C         Rb        Sr
     D  -1.77D0,  -2.49D0,
C         Y         Zr        Nb        Mo        Tc
     D  -3.22D0,  -3.95D0,  -4.70D0,  -5.46D0,  -6.21D0,
C         Ru        Rh        Pd        Ag       Cd
     D  -6.98D0,  -7.75D0,  -8.53D0, -1.00D2,   -1.0D0,

C         In        Sn        Sb        Te        I         Xe
     P  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,  -1.00D2,
C         Cs        Ba
     S  -1.00D2   ,-2.964D0,
     F  15*-100D0,
C        IP-D FROM D(N)
C         Hf       Ta       W       Re       Os
     D  -5.59D0,  -5.47D0,  -5.58D0, -5.93D0,   -6.52D0,
C         Ir       Pt       Au      Hg
     D  -7.36D0,  -8.44D0,  -90.0D0, -90.0D0,
     F   8* -1.00D2,

C             Ac       (ACTINIDE IP FROM 5F(N-4) 6D(3) 7S
     F   -100.0D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -5.60D0, -5.13D0, -5.25D0, -5.34D0, -5.10D0, -5.40D0, -5.09D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -4.84D0, -5.48D0, -6.08D0, -6.35D0, -6.58D0, -6.89D0, -4.78D0/
C
C
      DATA (EIF3(JJJ), JJJ = 1, 103)/
     +   88 *-1.00D2,
C         Ac       (ACTINIDE IP FROM 5F(N-4) 6D(2) 7S(2))
     F  -1.00D2,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -4.80D0, -8.50D0,-11.33D0,-13.20D0,-14.19D0, -1.00D2, -15.51D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -1.00D2, -1.00D2, -1.00D2, -1.00D2, -1.00D2, -1.00D2, -1.00D2/
C
C      ONLY NEEDED FOR Lr
C      VALUE STORED AS EIF1
C
      DATA (EIS4(JJJ), JJJ = 1, 103)/
     . 88*-100.0D0,
C            Ac       (ACTINIDE IP FROM 5F(N-2) 7S(2)
     F    -5.12D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -5.33D0, -5.52D0, -5.69D0, -5.85D0, -5.99D0, -6.12D0, -6.23D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -6.33D0, -6.40D0, -6.47D0, -6.51D0, -6.54D0, -6.55D0, -1.00D2/
C
       DATA (EIF4(JJJ), JJJ = 1, 103)/
     F   88 *   -1.00D2,
C         Ac       (ACTINIDE IP FROM 5F(N-3) 6D 7S(2)
     F  -2.80D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     F  -3.65D0, -4.42D0, -5.12D0, -5.75D0, -6.30D0, -6.79D0, -7.20D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     F  -7.54D0, -7.81D0, -8.01D0, -8.13D0, -8.18D0, -8.17D0, -4.60D0/
C
C     CNDO/2 AND INDO/2 : 1/2(IONIZATION POTENTIAL + ELECTRON AFFINITY)
C
C     PARAMETERS FROM D(N-1) S(1)
C
C
      DATA (EISS1(JJJ), JJJ = 1, 105)/
C         H
     S -7.1761D0,                                            0.0D0,
C         Li        Be
     S -3.1055D0,-5.94557D0,
C          B             C         N           O          F
     P    -9.59407D0,-14.051D0,-19.31637D0,-25.39017D0,-32.2724D0,0.0D0,
C          Na      Mg
     S -2.805D0,-5.222D0,
C           Al        Si        P          S         Cl
     P    -8.288D0,-11.157D0,-13.551D0,-16.328D0,-19.841D0,0.0D0,
C          K       Ca
     S -2.385D0,-3.400D0,
     D          -4.345D0,-4.440D0,-5.240D0,-5.225D0,-5.335D0,
     D          -5.725D0,-5.120D0,-3.840D0,-3.985D0,-6.220D0,
     P    6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-3.930D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
      DATA (EIPP1(JJJ), JJJ = 1, 105)/
     . 0.0D0,                                                 0.0D0,
     S -1.258D0,-2.563D0,
     P          -4.001D0,-5.572D0,-7.275D0,-9.111D0,-11.080D0,0.0D0,
     S -1.565D0,-2.100D0,
     P          -2.950D0,-4.270D0,-6.080D0,-7.385D0,-9.380D0,0.0D0,
     S -1.265D0,-1.989D0,
     D          -2.399D0,-1.978D0,-2.542D0,-2.971D0,-2.612D0,
     D          -3.655D0,-2.550D0,-1.317D0,-1.165D0,-1.640D0,
     P          6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-0.710D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
      DATA (EIPP1(JJJ), JJJ = 1, 105)/
     . 0.0D0,                                                 0.0D0,
     S -1.258D0,-2.563D0,
     P          -4.001D0,-5.572D0,-7.275D0,-9.111D0,-11.080D0,0.0D0,
     S -1.565D0,-2.100D0,
     P          -2.950D0,-4.270D0,-6.080D0,-7.385D0,-9.380D0,0.0D0,
     S -1.265D0,-1.989D0,
     D          -2.399D0,-1.978D0,-2.542D0,-2.971D0,-2.612D0,
     D          -3.655D0,-2.550D0,-1.317D0,-1.165D0,-1.640D0,
     P          6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-0.710D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
      DATA (EIDD1(JJJ), JJJ = 1, 105)/
     .   10*  -1.00D2,
     S -0.445D0,-0.650D0,
     P    -0.890D0,+0.625D0,-0.025D0,-0.405D0,+1.900D0,+0.000D0,
     S -0.635D0,-0.770D0,
     D          -2.15D0,-2.70D0,-3.245D0,-3.545D0,-3.410D0,
     D          -3.735D0,-4.210D0,-3.975D0,-6.240D0,-13.850D0,
     P           6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-4.530D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
C
C     PARAMETERS FROM D(N-2) S(2)
C
      DATA (EISS2(JJJ), JJJ = 1, 105)/
     S -7.1761D0,                                             0.0D0,
     S -3.1055D0,-5.94557D0,
     P  -9.59407D0,-14.051D0,-19.31637D0,-25.39017D0,-32.2724D0,0.0D0,
     S -2.805D0,-5.222D0,
     P  -8.288D0,-11.157D0,-13.551D0,-16.328D0,-19.841D0,0.0D0,
     S -2.385D0,-3.400D0,
     D           -3.657D0,-3.770D0,-3.822D0,-3.909D0,-3.983D0,
     D           -4.120D0,-4.170D0,-4.306D0,-4.567D0,-6.220D0,
     P  6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-3.930D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
      DATA (EIPP2(JJJ), JJJ = 1, 105)/
     . 0.0D0,                                                 0.0D0,
     S -1.258D0,-2.563D0,
     P          -4.001D0,-5.572D0,-7.275D0,-9.111D0,-11.080D0,0.0D0,
     S -1.565D0,-2.100D0,
     P          -2.950D0,-4.270D0,-6.080D0,-7.385D0,-9.380D0,0.0D0,
     S -1.265D0,-1.989D0,
     D          -0.558D0,-0.690D0,-0.777D0,-0.876D0,-0.975D0,
     D                   -1.062D0,-1.160D0,-1.260D0,-1.347D0,-1.640D0,
     P           6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-0.710D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
      DATA (EIDD2(JJJ), JJJ = 1, 105)/
     .   10*  -1.00D2,
     S -0.445D0,-0.650D0,
     P          -0.890D0,+0.625D0,-0.025D0,-0.405D0,+1.900D0,+0.000D0,
     S -0.635D0,-0.770D0,
     D          -3.793D0,-4.140D0,-4.475D0,-4.822D0,-5.157D0,
     D                   -5.504D0,-5.839D0,-6.182D0,-6.520D0,-13.850D0,
     P           6*  -1.00D2,
     S 2*  -1.00D2,
     D          -100.0D0,-100.0D0,-100.0D0,-4.530D0,-100.0D0,
     D          -100.0D0,-100.0D0,-100.0D0,-100.0D0,-100.0D0,
     .   57*  -1.00D2/
C
C      PARAMETERS FROM D(N) S(0)
C
      DATA (EIDD3(JJJ), JJJ = 1, 105)/
     .   105*  -1.00D2/
C     ISW1=1 FOR CNDO/1     ISW1=2 FOR CNDO/2
C     ISW1=3 FOR INDO/1     ISW1=4 FOR INDO/2
C     ISW1=5 FOR NDDO/1
C
C     ISW2=4 FOR DEFAULT CONFIG MIXING
C
C     CNDO/1 AND INDO/1
C
c     write(6,*) "i am in param"
c     write(6,*)"i am in beta"
c     write(6,*)"i am in coflod"
c     DO  I = 1, MXB
c        NEXTR(I)= 0
c     END DO
      DO  I = 1, MPE
         COFSQ2(I,1) = ZERO
         COFSQ2(I,2) = ONE
         COFSQ2(I,3) = ZERO
         COFSQ2(I,4) = ZERO
         COFSQ2(I,5) = ZERO
       END DO
C
C     FIRST TRANSITION SERIES
      COFSQ2(20,2) = 0.9608D0
      COFSQ2(21,2) = 0.9399D0
      COFSQ2(22,2) = 0.9069D0
      COFSQ2(23,2) = 0.8395D0
      COFSQ2(24,2) = 0.7052D0
      COFSQ2(25,2) = 0.6652D0
      COFSQ2(26,2) = 0.3143D0
      COFSQ2(27,2) = 0.2065D0
      COFSQ2(28,2) = 0.1421D0
      COFSQ2(29,2) = 0.0956D0
      COFSQ2(30,2) = ONE
      COFSQ2(20,1) = 0.0392D0
      COFSQ2(21,1) = 0.0601D0
      COFSQ2(22,1) = 0.0931D0
      COFSQ2(23,1) = 0.1605D0
      COFSQ2(24,1) = 0.2948D0
      COFSQ2(25,1) = 0.3348D0
      COFSQ2(26,1) = 0.6857D0
      COFSQ2(27,1) = 0.7935D0
      COFSQ2(28,1) = 0.8579D0
      COFSQ2(29,1) = 0.9044D0
      COFSQ2(30,1) = ZERO
      CONTINUE
C     SECOND ROW TRANSITION SERIES
      COFSQ2(37,2) = ZERO
      COFSQ2(38,2) = 0.9390D0
      COFSQ2(39,2) = 0.9200D0
      COFSQ2(40,2) = 0.8840D0
      COFSQ2(41,2) = 0.8060D0
      COFSQ2(42,2) = 0.6160D0
      COFSQ2(43,2) = 0.2420D0
      COFSQ2(44,2) = 0.0160D0
      COFSQ2(45,2) = 0.0030D0
      COFSQ2(46,2) = 0.0070D0
      COFSQ2(47,2) = ZERO
      COFSQ2(48,2) = ONE
      COFSQ2(37,1) = ONE
      COFSQ2(38,1) = 0.0550D0
      COFSQ2(39,1) = 0.0730D0
      COFSQ2(40,1) = 0.1080D0
      COFSQ2(41,1) = 0.1870D0
      COFSQ2(42,1) = 0.3820D0
      COFSQ2(43,1) = 0.7360D0
      COFSQ2(44,1) = 0.7170D0
      COFSQ2(45,1) = 0.3190D0
      COFSQ2(46,1) = 0.1210D0
      COFSQ2(47,1) = ONE
      COFSQ2(48,1) = ZERO
      COFSQ2(37,3) = ZERO
      COFSQ2(38,3) = 0.0070D0
      COFSQ2(39,3) = 0.0070D0
      COFSQ2(40,3) = 0.0080D0
      COFSQ2(41,3) = 0.0070D0
      COFSQ2(42,3) = 0.0020D0
      COFSQ2(43,3) = 0.0220D0
      COFSQ2(44,3) = 0.2660D0
      COFSQ2(45,3) = 0.6780D0
      COFSQ2(46,3) = 0.8720D0
      COFSQ2(47,3) = ZERO
      COFSQ2(48,3) = ZERO
C
C   PUT IN THE MIXING COEFS FOR THE LANTHANIDES
C   PROCESS 1 = F(N-3)D S(2) : PROCESS 2 = F(N-2) S(2)
      COFSQ2(58,1) = 0.7558D0
      COFSQ2(58,2) = 0.2442D0
      COFSQ2(59,1) = 0.2764D0
      COFSQ2(59,2) = 0.7236D0
      COFSQ2(60,1) = 0.1772D0
      COFSQ2(60,2) = 0.8228D0
      COFSQ2(61,1) = 0.1465D0
      COFSQ2(61,2) = 0.8535D0
      COFSQ2(62,1) = 0.0557D0
      COFSQ2(62,2) = 0.9423D0
      COFSQ2(63,1) = 0.0236D0
      COFSQ2(63,2) = 0.9764D0
      COFSQ2(64,1) = 0.9037D0
      COFSQ2(64,2) = 0.0963D0
      COFSQ2(65,1) = 0.4821D0
      COFSQ2(65,2) = 0.5179D0
      COFSQ2(66,1) = 0.1564D0
      COFSQ2(66,2) = 0.8436D0
      COFSQ2(67,1) = 0.1533D0
      COFSQ2(67,2) = 0.8467D0
      COFSQ2(68,1) = 0.1661D0
      COFSQ2(68,2) = 0.8339D0
      COFSQ2(69,1) = 0.0731D0
      COFSQ2(69,2) = 0.9269D0
      COFSQ2(70,1) = 0.0273D0
      COFSQ2(70,2) = 0.9727D0
C     THESE VALUES ARE FIXED AND ARE AS GOOD AS ANY - JAN 1987 -
C
C          THIRD TRANSITION SERIES UNDER CONSTRUCTION
C          HEGE STROEMSNES
C
C     COFSQ2(ATOMIC NO., TYPE)  TYPE = 1, 2, 3, 4 FOR
C     D(N-1)S(1),  D(N-2)S(2), D(N), AND D(N-2)SP RESPECTIVELY.
C     D(N-1)P   The mixing coef is set to 5000cm-1
C
      COFSQ2(72,1) = 0.0590D+00
      COFSQ2(72,2) = 0.9411D+00
      COFSQ2(72,3) = 0.0000D+00
      COFSQ2(72,4) = 0.0000D+00
      COFSQ2(72,5) = 0.0000D+00
      COFSQ2(73,1) = 0.1240D+00
      COFSQ2(73,2) = 0.8760D+00
      COFSQ2(73,3) = 0.0000D+00
      COFSQ2(73,4) = 0.0000D+00
      COFSQ2(73,5) = 0.0000D+00
      COFSQ2(74,1) = 0.1890D+00
      COFSQ2(74,2) = 0.8110D+00
      COFSQ2(74,3) = 0.0000D+00
      COFSQ2(74,4) = 0.0000D+00
      COFSQ2(74,5) = 0.0000D+00
      COFSQ2(75,1) = 0.2550D+00
      COFSQ2(75,2) = 0.7450D+00
      COFSQ2(75,3) = 0.0000D+00
      COFSQ2(75,4) = 0.0000D+00
      COFSQ2(75,5) = 0.0000D+00
      COFSQ2(76,1) = 0.3200D+00
      COFSQ2(76,2) = 0.6800D+00
      COFSQ2(76,3) = 0.0000D+00
      COFSQ2(76,4) = 0.0000D+00
      COFSQ2(76,5) = 0.0000D+00
      COFSQ2(77,1) = 0.3850D+00
      COFSQ2(77,2) = 0.6150D+00
      COFSQ2(77,3) = 0.0000D+00
      COFSQ2(77,4) = 0.0000D+00
      COFSQ2(77,5) = 0.0000D+00
      COFSQ2(78,1) = 0.5860D+00
      COFSQ2(78,2) = 0.0000D+00
      COFSQ2(78,3) = 0.4140D+00
      COFSQ2(78,4) = 0.0000D+00
      COFSQ2(78,5) = 0.0000D+00
      COFSQ2(79,1) = 0.9100D+00
      COFSQ2(79,2) = 0.0900D+00
      COFSQ2(79,3) = 0.0000D+00
      COFSQ2(79,4) = 0.0000D+00
      COFSQ2(79,5) = 0.0000D+00
      COFSQ2(80,1) = 0.0000D+00
      COFSQ2(80,2) = 1.0000D+00
      COFSQ2(80,3) = 0.0000D+00
      COFSQ2(80,4) = 0.0000D+00
      COFSQ2(80,5) = 0.0000D+00

C PUT IN THE MIXING COEFS FOR THE ACTINIDES
C K=1, 5F**(N-3) 6d 7S**2: K=2, 5F**(N-3) 7P 7S**2:
C K=3, 5F**(N-4) 6d**2 7S**2: K=4, 5F**(N-4) 6d**3 7S:
C PUT IN THE MIXING COEFS FOR THE ACTINIDES
C K=1, 5F**(N-3) 6d 7S**2: K=2, 5F**(N-3) 7P 7S**2:
C K=3, 5F**(N-4) 6d**2 7S**2: K=4, 5F**(N-4) 6d**3 7S:
C K=5  5F**(N-2) 7S**2
C
C
C     MIXING COEF 0.5
      COFSQ2(89,1) = 0.7653D+00
      COFSQ2(89,2) = 0.2321D+00
      COFSQ2(89,5) = 0.2640D-02
      COFSQ2(90,1) = 0.2391D-01
      COFSQ2(90,2) = 0.9189D-02
      COFSQ2(90,3) = 0.8951D+00
      COFSQ2(90,4) = 0.7024D-01
      COFSQ2(90,5) = 0.1609D-02
      COFSQ2(91,1) = 0.1769D+00
      COFSQ2(91,2) = 0.8687D-02
      COFSQ2(91,3) = 0.7948D+00
      COFSQ2(91,4) = 0.1839D-01
      COFSQ2(91,5) = 0.1239D-02
      COFSQ2(92,1) = 0.8801D+00
      COFSQ2(92,2) = 0.3607D-01
      COFSQ2(92,3) = 0.6480D-01
      COFSQ2(92,4) = 0.1080D-01
      COFSQ2(92,5) = 0.8212D-02
      COFSQ2(93,1) = 0.8754D+00
      COFSQ2(93,2) = 0.5112D-01
      COFSQ2(93,3) = 0.1150D-01
      COFSQ2(93,5) = 0.6202D-01
      COFSQ2(94,5) = 0.4919D+00
      COFSQ2(94,2) = 0.0000D+00
      COFSQ2(94,3) = 0.5081D+00
      COFSQ2(94,1) = 0.0000D+00
      COFSQ2(95,1) = 0.8689D+00
      COFSQ2(95,2) = 0.1122D+00
      COFSQ2(95,5) = 0.1886D-01
      COFSQ2(96,1) = 0.8663D+00
      COFSQ2(96,2) = 0.8005D-01
      COFSQ2(96,3) = 0.1267D-01
      COFSQ2(96,5) = 0.4097D-01
      COFSQ2(97,1) = 0.8213D+00
      COFSQ2(97,2) = 0.1519D+00
      COFSQ2(97,5) = 0.1965D-01
      COFSQ2(98,1) = 0.1142D+00
      COFSQ2(98,2) = 0.4120D-01
      COFSQ2(98,5) = 0.8446D+00
      COFSQ2(99,1) = 0.3206D-01
      COFSQ2(99,2) = 0.2522D-01
      COFSQ2(99,5) = 0.9427D+00
      COFSQ2(100,1) = 0.2492D-01
      COFSQ2(100,2) = 0.2539D-01
      COFSQ2(100,5) = 0.9497D+00
      COFSQ2(101,1) = 0.1047D-01
      COFSQ2(101,2) = 0.1230D-01
      COFSQ2(101,5) = 0.9772D+00
      COFSQ2(102,1) = 0.1123D-01
      COFSQ2(102,2) = 0.1441D-01
      COFSQ2(102,5) = 0.9744D+00
      COFSQ2(103,1) = 0.4490D+00
      COFSQ2(103,2) = 0.5510D+00
      
c     write(6,*)"EIS1",EIS1
      IF (ISW1.EQ.2 .OR. ISW1.EQ.4) GO TO 60
      DO 10 I = 1, 88
         COF = COFSQ2(I,1) + COFSQ2(I,2) 
         EIS(I) = (EIS1(I)*COFSQ2(I,1)+EIS2(I)*COFSQ2(I,2))/COF
C        THIS ASSUMES WEIGHTING FOR d(N-2)s(2)/d(n-1)s(1) IS THE
C        SAME AS d(N-2)s(1)p(1)/d(N-1)p(1)

         EIP(I) = (EIP1(I)*COFSQ2(I,1)+EIP2(I)*COFSQ2(I,2))/COF
         EID(I) = (EID1(I)*COFSQ2(I,1)+EID2(I)*COFSQ2(I,2)+
     .            EID3(I)*COFSQ2(I,3))/ (COF+COFSQ2(I,3))
         EIF(I) = (EIF1(I)*COFSQ2(I,1)+EIF2(I)*COFSQ2(I,2))/COF
   10 CONTINUE
      DO 20 I = 89, 103
C
C        Actinides
C
         COF = COFSQ2(I,1) + COFSQ2(I,2) + COFSQ2(I,3)
         EIS(I) = (COFSQ2(I,1)*EIS1(I)+COFSQ2(I,2)*EIS2(I)+
     .            COFSQ2(I,3)*EIS3(I)+COFSQ2(I,5)*EIS4(I))/
     .            (COF+COFSQ2(I,5))
         EIP(I) = EIP1(I)
         EID(I) = (COFSQ2(I,1)*EID1(I)+COFSQ2(I,3)*EID2(I)+
     .            COFSQ2(I,4)*EID3(I))/ (COFSQ2(I,1)+COFSQ2(I,3)+
     .            COFSQ2(I,4))
         EIF(I) = (EIF1(I)*COFSQ2(I,1)+EIF2(I)*COFSQ2(I,2)+
     .            EIF3(I)*COFSQ2(I,3)+EIF4(I)*COFSQ2(I,5))/
     .            (COF+COFSQ2(I,5))
   20 CONTINUE
C     CORRECTION FOR Ac
      EIF(89) = EIF4(89)
      EIS(89) = (COFSQ2(I,1)*EIS1(I)+COFSQ2(I,2)*EIS2(I)+
     .          COFSQ2(I,5)*EIS4(I))/ (COFSQ2(I,1)+COFSQ2(I,2)+
     .          COFSQ2(I,5))
C     CORRECTION FOR Th
      EIF(90) = (EIF1(90)*COFSQ2(90,1)+EIF2(90)*COFSQ2(90,2)+
     .          EIF4(I)*COFSQ2(I,5))/ (COFSQ2(90,1)+COFSQ2(90,2)+
     .          COFSQ2(I,5))
C
      IF (ISW2.EQ.1) THEN
         DO 30 I = 1, MPE
            EIS(I) = EIS1(I)
            EIP(I) = EIP1(I)
            EID(I) = EID1(I)
            EIF(I) = EIF1(I)
   30    CONTINUE
      ELSE IF (ISW2.EQ.2) THEN
         DO 40 I = 1, MPE
            EIS(I) = EIS2(I)
            EIP(I) = EIP2(I)
            EID(I) = EID2(I)
            EIF(I) = EIF2(I)
   40    CONTINUE
      ELSE IF (ISW2.EQ.3) THEN
         DO 50 I = 1, MPE
            EID(I) = EID3(I)
            EIF(I) = EIF3(I)
   50    CONTINUE
      END IF
      GO TO 110
C
C      CNDO/2 AND INDO/2
C
   60 CONTINUE
      DO 70 I = 1, MPE
         COF = COFSQ2(I,1) + COFSQ2(I,2)
         EIS(I) = (EISS1(I)*COFSQ2(I,1)+EISS2(I)*COFSQ2(I,2))/COF
         EIP(I) = (EIPP1(I)*COFSQ2(I,1)+EIPP2(I)*COFSQ2(I,2))/COF
         EID(I) = (EIDD1(I)*COFSQ2(I,1)+EIDD2(I)*COFSQ2(I,2)+
     .            EIDD3(I)*COFSQ2(I,3))/ (COF+COFSQ2(I,3))
C        EIF(I) = EIFF1(I)
   70 CONTINUE
      IF (ISW2.EQ.1) THEN
         DO 80 I = 1, MPE
            EIS(I) = EISS1(I)
            EIP(I) = EIPP1(I)
            EID(I) = EIDD1(I)
C           EIF(I) = EIFF1(I)
   80    CONTINUE
      ELSE IF (ISW2.EQ.2) THEN
         DO 90 I = 1, MPE
            EIS(I) = EISS2(I)
            EIP(I) = EIPP2(I)
            EID(I) = EIDD2(I)
C           EIF(I) = EIFF1(I)
   90    CONTINUE
      ELSE IF (ISW2.EQ.3) THEN
         DO 100 I = 1, MPE
            EID(I) = EIDD3(I)
  100    CONTINUE
      END IF
C
C
  110 CONTINUE
      DO 120 I = 1, NA
         KK = KAN(I)
         IF ( KK      .LE.   0   ) GO TO 119
         IF ( EIS(KK) .LT. -99.0 )  THEN
          GO TO 130
         ENDIF
         IF ( NTYP(I) .EQ.   1   ) GO TO 119
         IF ( EIP(KK) .LT. -99.0 )  THEN
         GO TO 130
         ENDIF
         IF ( NTYP(I) .LT.   3   ) GO TO 119
         IF ( EID(KK) .LT. -99.0 )  THEN
         GO TO 130
         ENDIF
         IF ( NTYP(I) .LT.   5   ) GO TO 119
         IF ( EIF(KK) .LT. -99.0 )  THEN
         GO TO 130
         ENDIF
  119 CONTINUE
C     WRITE(6,*) ' PARAM', NFIRST(I), NEXTR(NFIRST(I))
  120 CONTINUE
C
      GO TO 140

  130 WRITE (6,9000) KK
      CALL ZEXIT(3)
CRAY  CALL FIN(' STOP IN PARAM 313      *PK* ')
  140 CONTINUE
C     WRITE(6,*) ' ##### RIDLEY IPs ##### '
 9000 FORMAT (' XXXXX NO PARAM FOR ATOM NO. = ',I5,' XXXXXXX')
C
c     write(6,*)"check", "EIS",EIS,"EIP",EIP,"EID",EID,"EIF",EIF
 5000 RETURN

      END
