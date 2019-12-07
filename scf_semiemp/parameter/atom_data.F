       subroutine atom_data(atnum,nat,sexpo,pexpo,dexpo,fexpo)
c
      implicit none
c
      integer atnum(*)
      integer nat
      double precision sexpo(*)
      double precision pexpo(*)
      double precision dexpo(*)
      double precision fexpo(*)
c
c     Get the atomic data
      call ATOMMN(atnum,nat,sexpo,pexpo,dexpo,fexpo)
c
      return
      end
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C*MODULE MNSOL  *DECK DATSOLMN
C
C     Source: http://comp.chem.umn.edu/solvation/mnsddb.pdf
C
      SUBROUTINE ATOMMN(ATOM_NUM,NAT,SEXPO,PEXPO,DEXPO,FEXPO)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
C     LOGICAL MASWRK
C
      CHARACTER*8 ATOM,ATOM_TYPE
      INTEGER AM,CM,JR,IBASS,FMASS
      double precision  AS,AP,AD,AF
      double precision CS,CP,CD,CF
      double precision COVRAD
      INTEGER ATOM_NUM,IBASIS
      double precision SEXPO,PEXPO,DEXPO,FEXPO
      INTEGER VDWR
      CHARACTER*8 AOS
      INTEGER NAT
C
      DIMENSION ATOM(105),AM(105),CM(105),JR(8,8),FMASS(105),
     $ AS(105,3),AP(105,3),AD(105,3),AF(105,3),CS(105,3),
     $CP(105,3),CD(105,3),CF(105,3),COVRAD(105),VDWR(54),
     $SEXPO(105,6),PEXPO(105,6),DEXPO(105,6),FEXPO(105,6),
     $AOS(16),IBASIS(8,6),ATOM_NUM(NAT)
C     ZINDO   UF   1990 - 1996


      DATA AOS/' S      ',' P(X)   ',' P(Y)   ',' P(Z)   ',' D(Z2)  ',
     1         'D(X2-Y2)',' D(XY)  ',' D(XZ)  ',' D(YZ)  ',' F(Z3)  ',
     2         ' F(XZ2) ',' F(YZ2) ','Z(X2-Y2)',' F(XYZ) ',' X3-3XY2',
     3         ' 3YX2-Y3'/
C
C
      DATA ATOM/' H','Ct','Li','Be',' B',' C',' N',' O',' F','Ne','Na',
     1 'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr',
     2 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     3 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     4 'Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     5 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re',
     6 'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     7 'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     8 'Md','No','Lr','XX','XX'/
C-----------------------------------------------------------------------
      DATA AM/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     .        0.0D0,0.831D0,0.831D0,1.7305D0,1.731D0,1.731D0,1.731D0,
     .        1.731D0,2.252D0,1.180D0,1.5097D0,1.6102D0,1.71016D0,
     .        1.81119D0,1.91367D0,2.01520D0,2.11093D0,2.21589D0,
     .        2.31682D0,2.41805D0,2.51974D0,75*0.0D0/
CMK92     .        2.31682D0,2.41805D0,2.51974D0,64*0.0D0/
C-----------------------------------------------------------------------
      DATA CM/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     .        0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,
     .        0.0D0,0.8105D0,0.8961D0,0.96428D0,1.04854D0,1.10712D0,
     .        1.1751D0,1.22250D0,1.23008D0,1.26242D0,1.29449D0,
     .        1.33244D0,75*0.0D0/
CMK92     .        1.33244D0,62*0.0D0/
C---------------------------------------------------------------------
      DATA JR /1,2,3,4,5,6,7,8,
     .         2,1,4,3,6,5,8,7,
     .         3,4,1,2,8,7,6,5,
     .         4,3,2,1,7,8,5,6,
     .         5,6,8,7,1,2,4,3,
     .         6,5,7,8,2,1,3,4,
     .         7,8,6,5,4,3,1,2,
     .         8,7,5,6,3,4,2,1/
C---------------------------------------------------------------------
      DATA IBASIS /2,1,2,3,3,4,7,1,
     .             2,1,2,4,3,3,6,1,
     .             1,2,2,1,2,2,4,1,
     .             2,2,1,4,4,4,8,1,
     .             2,2,1,3,4,3,5,1,
     .             1,1,1,2,1,2,3,1/
C---------------------------------------------------------------------
C     JR(I,J)=The direct product of irrep i with irrep J.
C     IBASIS(I,J)=In the I'th group the irrep of operator J.
C             J= 1  2  3  4    5    6
C                X  Y  Z R(X) R(Y) R(Z)
C     For sym. conventions see sub CORREL
C------------------------------------------------------------------
C     Table of atomic masses
C------------------------------------------------------------------
      DATA FMASS/1.0079D0,4.00260D0,6.941D0,9.01218D0,10.81D0,
     . 12.011D0,14.0067D0,15.9994D0,18.998403D0,20.179D0,22.98977D0,
     . 24.305D0,26.98154D0,28.0855D0,30.97376D0,32.06D0,35.453D0,
     . 39.948D0,39.0983D0,40.08D0,44.9559D0,47.90D0,50.9415D0,
     . 51.996D0,54.9380D0,55.847D0,58.9332D0,58.70D0,63.546D0,
     . 65.38D0,69.72D0,72.59D0,74.9216D0,78.96D0,79.904D0,83.80D0,
     . 85.4678D0,87.62D0,88.9059D0,91.22D0,92.9064D0,95.94D0,98.00D0,
     . 101.07D0,102.9055D0,106.4D0,107.868D0,112.41D0,114.82D0,
     . 118.69D0,121.75D0,127.60D0,126.9045D0,131.30D0,132.9054D0,
     . 137.33D0,138.9055D0,140.12D0,140.9077D0,144.24D0,145.0D0,
     . 150.4D0,151.96D0,157.25D0,158.9254D0,162.50D0,164.9304D0,
     . 167.26D0,168.9342D0,173.04D0,174.967D0,178.49D0,180.9479D0,
     . 183.85D0,186.207D0,190.2D0,192.22D0,195.09D0,196.9665D0,
     . 200.59D0,204.37D0,207.2D0,208.9804D0,209.0D0,210.0D0,222.0D0,
     . 223.0D0,226.0254D0,227.0278D0,232.0381D0,231.0359D0,238.029D0,
     . 237.0482D0,244.0D0,243.0D0,247.0D0,247.0D0,251.0D0,252.0D0,
     . 257.0D0,258.0D0,259.0D0,260.0D0,261.0D0,262.0D0/
CMK92     . 257.0D0,258.0D0,259.0D0,260.0D0,261.0D0,262.0D0,263.0D0/
C-----------------------------------------------------------------------
C      AS=S BASIS, AP=P BASIS, AD=D BASIS AND AF=F BASIS
C      DOUBLE ZETA FIRST ROW TRANSITION ELEMENTS ARE CLEMENTI ROETTI
C      ATOMIC AND NUC. DATA TABLES 14,177 (1974).
C-----------------------------------------------------------------------
      DATA (AS(I,1), I = 1, 105)/
C         H       He      Li      Be
     . 1.20D0, 1.70D0, 0.65D0, 0.975D0,
C         B        C         N        O        F        Ne
     . 1.30D0,  1.625D0,  1.95D0,  2.275D0,  2.60D0,  2.925D0,
C         Na       Mg
     . 0.8360D0, 1.103D0,
C         Al       Si       P        S        Cl       Ar
C     SI HAS BEEN CHANGED BY RAJIV BENDALE MAY 1994
     . 1.370D0, 1.600D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
CMZ  . 1.370D0, 1.520D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
C         K        Ca
     . 1.180D0, 1.210D0,
C         Sc       Ti       V        Cr       Mn
     . 1.2300D0, 1.300D0, 1.300D0, 1.320D0, 1.360D0,
C         Fe       Co       Ni       Cu       Zn
     . 1.3700D0, 1.423D0, 1.473D0, 1.482D0, 1.509D0,
C         Ga       Ge       As       Se       Br       Kr
     . 1.7667D0, 2.0109D0, 2.236D0, 2.439D0, 2.638D0, 0.000D0,
C         Rb       Sr
     . 0.000D0, 1.20866D0,
C         Y        Zr       Nb       Mo       Tc
     . 1.275D0, 1.337D0, 1.3930D0, 1.440D0, 1.461D0,
C         Ru       Rh       Pd       Ag       Cd
     . 1.470D0, 1.482D0, 1.4970D0, 1.507D0, 1.699D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 2.1198D0, 0.0000D0, 0.000D0, 2.671D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0,1.241682, 1.2550D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.548D0, 1.571D0, 1.593D0, 1.614D0, 1.635D0, 1.656D0, 1.678D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.700D0, 1.723D0, 1.745D0, 1.767D0, 1.789D0, 1.812D0, 1.834D0,
C                  Hf       Ta       W        Re
     .          1.7440D0, 1.8402D0, 1.9079D0, 1.9938D0, 
C         Os       Ir       Pt       Au       Hg
     . 2.0618D0, 2.1181D0, 2.1080D0, 2.1070D0, 2.3059D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 1.759D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.798D0, 1.831D0, 1.862D0, 1.891D0, 1.919D0, 1.946D0, 1.977D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 2.007D0, 2.037D0, 2.066D0, 2.096D0, 2.126D0, 2.157D0, 2.188D0,
     . 2 * 0.0D0/
C
      DATA (AS(I, 2), I = 1, 105)/
C         H       He       Li      Be
     . 1.200D0, 1.700D0, 0.650D0, 0.975D0,
C         B        C        N        O        F        Ne
     . 1.300D0, 1.625D0, 1.950D0, 2.275D0, 2.600D0, 2.925D0,
C         Na       Mg
     . 0.836D0, 1.103D0,
C         Al       Si       P        S        Cl       Ar
     . 1.370D0, 1.600D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
CMZ  . 1.370D0, 1.520D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
C         K        Ca
     . 1.180D0, 1.210D0,
C         Sc       Ti       V        Cr       Mn
     . 1.230D0, 1.300D0, 1.300D0, 1.320D0, 1.360D0,
C         Fe       Co       Ni       Cu       Zn
     . 1.370D0, 1.423D0, 1.473D0, 1.482D0, 1.509D0,
C         Ga       Ge       As       Se       Br       Kr
     . 1.7667D0,2.0109D0,2.921D0, 2.439D0, 2.638D0, 0.000D0,
C         Rb       Sr
     + 0.000D0, 1.63036,
C         Y        Zr       Nb       Mo       Tc
     . 1.275D0, 1.337D0, 1.393D0, 1.440D0, 1.461D0,
C         Ru       Rh       Pd       Ag       Cd
     . 1.470D0, 1.482D0, 1.497D0, 1.507D0, 1.699D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 2.7626D0, 0.000D0, 0.000D0, 3.341D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0,1.1024D0,1.255D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.548D0, 1.571D0, 1.593D0, 1.614D0, 1.635D0, 1.656D0, 1.678D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.700D0, 1.723D0, 1.745D0, 1.767D0, 1.789D0, 1.812D0, 1.834D0,
C                  Hf       Ta       W        Re
     .          1.5000D0,1.5500D0,1.6000D0,1.6500D0,
C         Os       Ir       Pt       Au       Hg
     . 1.7000D0,1.7500D0,1.6500D0,1.7500D0,1.8500D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 1.759D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.798D0, 1.831D0, 1.862D0, 1.891D0, 1.919D0, 1.946D0, 1.977D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 2.007D0, 2.037D0, 2.066D0, 2.096D0, 2.126D0, 2.157D0, 2.188D0,
     . 2 * 0.0D0/
C
      DATA (AS(I, 3), I = 1, 105)/
C         H       He
C         Li      Be
C         B        C        N        O        F        Ne
C         Na       Mg
C         Al       Si       P        S        Cl       Ar
C         K        Ca
C         Sc       Ti       V        Cr       Mn
C         Fe       Co       Ni       Cu       Zn
C         Ga       Ge       As       Se       Br       Kr
     . 32  *0.0D0,       1.7361D0, 0.000D0, 0.000D0, 0.000D0,
C
C         Rb       Sr
     + 0.000000, 0.96129D0,
C         Y        Zr       Nb       Mo       Tc
C         Ru       Rh       Pd       Ag       Cd
C         In
     + 11 * 0.000D0,
C         Sn
     . 1.625D0,
C         Sb        Te       I        Xe
     . 0.000D0, 0.000D0, 2.04555D0, 0.000D0,
C         Cs       Ba 
     . 0.0000D0,2.8012D0,
C         La
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     + 15*0.0D0,
C                  Hf       Ta       W        Re
     .         4.3000D0,4.3000D0,4.0500D0,4.4000D0,
C         Os       Ir       Pt       Au       Hg
     . 4.6000D0,4.7000D0,4.1000D0,4.6000D0,4.6000D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U        Np       Pu       Am       Cm
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 25*0.0D0/
C
      DATA (AP(I, 1), I = 1, 105)/
C         H       He        Li      Be
     . 1.200D0, 1.700D0, 0.650D0, 0.975D0,
C         B        C        N        O        F        Ne
     . 1.300D0, 1.625D0, 1.950D0, 2.275D0, 2.600D0, 2.925D0,
C         Na       Mg
     . 0.836D0, 1.103D0,
C         Al       Si       P        S        Cl       Ar
     . 1.370D0, 1.600D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
CMZ  . 1.370D0, 1.520D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
C         K        Ca
     . 1.180D0, 1.210D0,
C         Sc       Ti       V        Cr       Mn
     . 1.230D0, 1.300D0, 1.300D0, 1.320D0, 1.360D0,
C         Fe       Co       Ni       Cu       Zn
     . 1.370D0, 1.423D0, 1.473D0, 1.482D0, 1.509D0,
C         Ga       Ge       As       Se       Br       Kr
     . 1.555D0, 1.695D0, 1.862D0, 2.072D0, 2.257D0, 0.000D0,
C         Rb       Sr
     + 0.000D0, 1.20866D0,
C         Y        Zr       Nb       Mo       Tc
     . 1.275D0, 1.337D0, 1.393D0, 1.440D0, 1.461D0,
C         Ru       Rh       Pd       Ag       Cd
     . 1.470D0, 1.482D0, 1.497D0, 1.507D0, 1.699D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 1.69197D0, 0.000D0,0.000D0, 2.195D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0,1.2416820,0.949D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.361D0, 1.380D0, 1.398D0, 1.415D0, 1.432D0, 1.448D0, 1.463D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.478D0, 1.492D0, 1.507D0, 1.521D0, 1.535D0, 1.549D0, 1.5620D0,
C                  Hf       Ta       W        Re
     .          1.7440D0, 1.8402D0, 1.9079D0, 1.9938D0, 
C         Os       Ir       Pt       Au       Hg
     . 2.0618D0, 2.1181D0, 2.1080D0, 2.1070D0, 2.3059D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 1.504D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.538D0, 1.567D0, 1.593D0, 1.618D0, 1.642D0, 1.665D0, 1.686D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 1.707D0, 1.727D0, 1.747D0, 1.766D0, 1.786D0, 1.805D0, 1.825D0,
     . 2 * 0.0D0/
C
      DATA (AP(I, 2), I = 1, 105)/
C         H       He        Li      Be
     . 1.200D0, 1.700D0, 0.650D0, 0.975D0,
C         B        C        N        O        F        Ne
     . 1.300D0, 1.625D0, 1.950D0, 2.275D0, 2.600D0, 2.925D0,
C         Na       Mg
     . 0.836D0, 1.103D0,
C         Al       Si       P        S        Cl       Ar
     . 1.370D0, 1.600D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
CMZ  . 1.370D0, 1.520D0, 1.730D0, 1.925D0, 2.130D0, 2.365D0,
C        K        Ca
     . 1.180D0, 1.210D0,
C         Sc       Ti       V        Cr       Mn
     . 1.230D0, 1.300D0, 1.300D0, 1.320D0, 1.360D0,
C         Fe       Co       Ni       Cu       Zn
     . 1.370D0, 1.423D0, 1.473D0, 1.482D0, 1.509D0,
C         Ga       Ge       As       Se       Br       Kr
     . 1.5554D0,1.695D0, 2.497D0, 2.072D0, 2.257D0, 0.000D0,
C         Rb       Sr
     + 0.000D0, 1.63036,
C         Y        Zr       Nb       Mo       Tc
     . 1.275D0, 1.337D0, 1.393D0, 1.440D0, 1.461D0,
C         Ru       Rh       Pd       Ag       Cd
     . 1.470D0, 1.482D0, 1.497D0, 1.507D0, 1.699D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 2.35868D0, 0.000D0,0.000D0, 2.920D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0, 1.1024D0,0.949D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.361D0, 1.380D0, 1.398D0, 1.415D0, 1.432D0, 1.448D0, 1.463D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.478D0, 1.492D0, 1.507D0, 1.521D0, 1.535D0, 1.549D0, 1.562D0,
C                  Hf       Ta       W        Re
     .          1.5000D0,1.5500D0,1.6000D0,1.6500D0,
C         Os       Ir       Pt       Au       Hg
     . 1.7000D0,1.7500D0,1.6500D0,1.7500D0,1.8500D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 1.504D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.538D0, 1.567D0, 1.593D0, 1.618D0, 1.642D0, 1.665D0, 1.686D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 1.707D0, 1.727D0, 1.747D0, 1.766D0, 1.786D0, 1.805D0, 1.825D0,
     . 2 * 0.0D0/
C
      DATA (AP(I, 3), I = 1, 105)/
C                   As                 Sr                 Sn
     . 32*0.0D0, 1.4078D0, 4*0.0000, 0.96129D0, 11*0.0D0, 1.32009D0,
C         Sb        Te        I        Xe
     . 0.000D0,  0.000D0, 1.67111D0, 0.0D0,
C         Cs       Ba 
     . 0.0000D0,2.8012D0,
C         La
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     + 15 * 0.000D0,
C                  Hf       Ta       W        Re
     .         4.3000D0,4.3000D0,4.0500D0,4.4000D0,
C         Os       Ir       Pt       Au       Hg
     . 4.6000D0,4.7000D0,4.1000D0,4.6000D0,4.6000D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C         Ac
C         Th       Pa       U        Np       Pu       Am       Cm
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 25*0.0D0/
C
C----------------------------------------------------------------
C AP Data set for P orbitals included. presently data is equal to
C S basis (earlier AS(*) was used for P orbitals too).
C
C The D orbitals for the first rowd0, LI-FD0, are assumed to be
C polarization with sd exponent = 2S
C----------------------------------------------------------------
C
      DATA (AD(I, 1), I = 1, 105)/
C         H       He        Li      Be
     . 1.200D0, 0.000D0, 0.650D0, 0.975D0,
C         B        C        N        O        F        Ne
     . 1.300D0, 1.625D0, 1.950D0, 2.275D0, 2.600D0, 0.000D0,
C         Na       Mg
     . 0.831D0, 0.831D0,
C         Al       Si       P        S        Cl       Ar
     . 1.731D0, 1.731D0, 1.731D0, 1.731D0, 1.731D0, 2.252D0,
C         K        Ca
     . 1.180D0, 1.850D0,
C         Sc       Ti       V        Cr       Mn
     . 2.020D0, 2.497D0, 2.738D0, 2.966D0, 3.195D0,
C         Fe       Co       Ni       Cu       Zn
     . 3.383D0, 3.580D0, 3.776D0, 3.972D0, 4.173D0,
C         Ga       Ge       As       Se       Br       Kr
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Rb       Sr
     . 0.000D0, 0.000D0,
C         Y        Zr       Nb       Mo       Tc
CMCZ   WALTER THIEL Zr = 1.658
     . 1.878D0, 2.140D0, 2.355D0, 2.532D0, 2.660D0,
C         Ru       Rh       Pd       Ag       Cd
     . 2.748D0, 2.872D0, 2.998D0, 3.321D0, 3.617D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 2.671D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0,1.241682, 1.869D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 2.061D0, 2.102D0, 2.138D0, 2.171D0, 2.201D0, 2.229D0, 2.243D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 2.254D0, 2.262D0, 2.269D0, 2.274D0, 2.277D0, 2.279D0, 2.278D0,
C                  Hf       Ta       W        Re
     .          2.5197D0, 2.7140D0, 2.9326D0, 3.1148D0, 
C         Os       Ir       Pt       Au       Hg
     . 3.2544D0, 3.4062D0, 3.5519D0, 3.6916D0, 3.8434D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     . 8 *0.0D0, 2.049D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 2.136D0, 2.204D0, 2.261D0, 2.311D0, 2.356D0, 2.395D0, 2.413D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 2.427D0, 2.438D0, 2.446D0, 2.451D0, 2.453D0, 2.453D0, 2.450D0,
     . 2 * 0.0D0/
C
      DATA (AD(I, 2), I = 1, 105)/
C         H       He        Li      Be
     . 1.200D0, 0.000D0, 0.650D0, 0.975D0,
C         B        C        N        O        F        Ne
     . 1.300D0, 1.625D0, 1.950D0, 2.275D0, 2.600D0, 0.000D0,
C         Na       Mg
     . 0.831D0, 0.831D0,
C         Al       Si       P        S        Cl       Ar
     . 1.731D0, 1.731D0, 1.731D0, 1.731D0, 1.731D0, 2.252D0,
C         K        Ca
     . 1.180D0, 1.850D0,
C         Sc       Ti       V        Cr       Mn
     . 4.2224D0, 4.670D0, 5.052D0, 5.138D0, 5.767D0,
C         Fe       Co       Ni       Cu       Zn
     . 6.06828D0,  6.386D0, 6.706D0, 6.795D0, 7.349D0,
C         Ga       Ge       As       Se       Br       Kr
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Rb       Sr
     . 0.000D0, 0.000D0,
C         Y        Zr       Nb       Mo       Tc
     . 3.837D0, 3.639D0, 3.7741D0, 3.954D0, 4.124D0,
C         Ru       Rh       Pd       Ag       Cd
     . 4.259D0, 4.485D0, 4.723D0,  4.989D0, 5.213D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 2.671D0, 0.000D0,
C         Cs       Ba       La
     . 0.000D0,1.241682, 2.987D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 3.077D0, 3.126D0, 3.169D0, 3.208D0, 3.244D0, 3.278D0, 3.288D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 3.298D0, 3.307D0 ,3.316D0, 3.324D0, 3.332D0, 3.338D0, 3.343D0,
C                  Hf       Ta       W        Re
     .          1.60D0,   1.80D0,  2.00D0,  2.200D0,
C         Os       Ir       Pt       Au       Hg
     .  2.35D0,  2.55D0  , 2.40D0,  2.60D0,  2.90D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 3.059D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 3.191D0, 3.290D0, 3.370D0, 3.437D0, 3.494D0, 3.545D0, 3.552D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 3.558D0, 3.562D0, 3.566D0, 3.568D0, 3.569D0, 3.569D0, 3.567D0,
     . 2 * 0.0D0/
C
      DATA (AD(I, 3), I = 1, 105)/
     . 20*0.0D0,
C         Sc       Ti       V        Cr       Mn
     . 1.7465D0, 1.986D0, 2.173D0, 2.077D0, 2.510D0,
C         Fe       Co       Ni       Cu       Zn
     . 2.61836D0,  2.745D0, 2.874D0, 2.765D0, 3.139D0,
C         Ga       Ge       As       Se       Br       Kr
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Rb       Sr
     . 0.000D0, 0.000D0,
C         Y        Zr       Nb       Mo       Tc
     . 1.739D0, 1.804D0, 1.9251D0, 2.047D0, 2.155D0,
C         Ru       Rh       Pd       Ag       Cd
     . 2.094D0, 2.217D0, 2.3420D0, 2.584D0, 2.706D0,
C         In       Sn       Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Cs       Ba
     . 0.000D0, 0.000D0,
C         La
     . 1.589D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.581D0, 1.598D0, 1.612D0, 1.623D0, 1.633D0, 1.640D0, 1.634D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.628D0, 1.621D0, 1.613D0, 1.605D0, 1.596D0, 1.587D0, 1.577D0,
C                  Hf       Ta       W        Re
     .          4.05D0,  4.60D0,  5.05D0,  5.75D0, 
C         Os       Ir       Pt       Au       Hg
     .  6.05D0,  6.80D0,  6.20D0, 6.85D0,  7.45D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 1.569D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.630D0, 1.674D0, 1.708D0, 1.734D0, 1.755D0, 1.773D0, 1.767D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 1.759D0, 1.751D0, 1.741D0, 1.730D0, 1.718D0, 1.706D0, 1.693D0,
     . 2 * 0.0D0/
C
C       DATA AF/10*0.0D0,
C     X 26*0.0D0,
C     X 18*0.0D0,
C     X 0.000D0,0.000D0,3.940D0,4.253D0,4.509D0,4.712D0,4.902D0,
C     . 5.085D0,5.269D0,5.396D0,
C     X 5.535D0,5.682D0,5.829D0,5.968D0,6.106D0,6.259D0,6.667D0,
C     X 21*0.00D0,
C     X 10*0.0D0,
C     X 26*0.0D0,
C     X 18*0.0D0,
C     X 0.000D0,0.000D0,3.940D0,4.253D0,4.509D0,4.712D0,4.902D0,
C     . 5.085D0,5.269D0,5.396D0,
C     X 5.535D0,5.682D0,5.829D0,5.968D0,6.106D0,6.259D0,6.667D0,
C     X 21*0.00D0,
C     X 10*0.0D0,
C     X 26*0.0D0,
C     X 18*0.0D0,
C     X 17*0.0D0,
C     X 21*0.00D0/
C
      DATA (AF(I, 1), I = 1, 105)/
     . 10*0.00D0,
     . 40*0.00D0,
C         Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Cs       Ba
     . 0.000D0, 0.000D0,
C         La
     . 4.022D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 4.439D0, 4.657D0, 4.681D0, 5.053D0, 5.236D0, 5.414D0, 5.565D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 5.717D0, 5.869D0, 6.019D0, 6.168D0, 6.316D0 ,6.426D0, 6.607D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     . 17 *0.0D0, 2.634D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 3.395D0, 3.667D0, 3.900D0, 4.110D0, 4.304D0, 4.485D0, 4.580D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 4.749D0, 4.911D0, 5.067D0, 5.219D0, 5.367D0, 5.512D0, 5.653D0,
     . 2 * 0.0D0/
C
      DATA (AF(I, 2), I = 1, 105)/
     . 10*0.0D0,
     . 40*0.00D0,
C         Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Cs       Ba
     . 0.000D0, 0.000D0,
C         La
     . 7.802D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 6.118D0 ,6.393D0, 6.648D0, 6.889D0, 7.120D0, 7.342D0, 7.523D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 7.707D0, 7.892D0, 8.076D0, 8.258D0, 8.438D0, 8.616D0, 8.791D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     . 17 *0.0D0, 3.718D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 4.692D0, 5.058D0, 5.381D0, 5.678D0, 5.956D0, 6.212D0, 6.309D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 6.554D0, 6.790D0, 7.017D0, 7.238D0, 7.452D0, 7.661D0, 7.859D0,
     . 2 * 0.0D0/
C
      DATA (AF(I, 3), I = 1, 105)/
     . 10*0.0D0,
     . 40*0.0D0,
C         Sb       Te       I        Xe
     . 0.000D0, 0.000D0, 0.000D0, 0.000D0,
C         Cs       Ba
     . 0.000D0, 0.000D0,
C         La
     . 3.371D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 2.522D0, 2.646D0, 2.757D0, 2.858D0, 2.951D0, 3.038D0, 3.091D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 3.148D0, 3.206D0, 2.264D0, 3.320D0, 3.375D0, 3.428D0, 3.480D0,
C                  Hf       Ta       W        Re
C         Os       Ir       Pt       Au       Hg
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     . 17 *0.0D0, 1.329D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.988D0, 2.201D0, 2.384D0, 2.549D0, 2.702D0, 2.842D0, 2.862D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 2.990D0, 3.113D0, 3.229D0, 3.340D0, 3.447D0, 3.551D0, 3.649D0,
     . 2 * 0.0D0/

      DATA CS/
C     1ST SET
     + 105*1.0D0,
C
C     2ND SET
C                    As                   Sr                 Sn
     + 32*1.0D0, 0.566645D0, 4*1.00D0 ,0.495523D0,11*1.0D0, 0.592426D0,
C                    I              Xe
     .  2*1.0D0, 0.618030D0,   1.0000D0,
C         Cs       Ba
     . 1.000D0, 0.8442D0,
C         La
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .  15*1.00D0,
C                  Hf       Ta       W        Re
     .          0.8600D0,0.8400D0,0.8200D0,0.8200D0,
C         Os       Ir       Pt       Au       Hg
     . 0.8200D0,0.8200D0,0.7600D0,0.8200D0,0.7800D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *1.0D0, 1.000D0,
C         Th       Pa       U        Np       Pu       Am       Cm
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 16*1.0D0,
C     3RD SET
C                    As                    Sr                  Sn
     + 32*0.0D0, 0.505256D0, 4*0.0D0, 0.593348D0, 11*0.0D0, 0.497235D0,
C                    I           Xe
     .  2*0.0D0, 0.458166D0,  0.00D0,
C         Cs       Ba
     . 0.000D0, 0.362138D0,
C         La
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     .  15*0.00D0,
C                  Hf       Ta       W        Re
     .          0.3795D0,0.4016D0,0.3985D0,0.4185D0,
C         Os       Ir       Pt       Au       Hg
     . 0.4239D0,0.4211D0,0.4753D0,0.4131D0,0.4489D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
C                  Ac
     .  8 *0.0D0, 0.000D0,
C         Th       Pa       U        Np       Pu       Am       Cm
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 16*0.0D0/
C
      DATA CP/
C     1ST SET
     + 105*1.0D0,
C
C     2ND SET
C                    As                    Sr                  Sn
     + 32*1.0D0, 0.503682D0, 4*1.0D0, 0.495523D0, 11*1.0D0, 0.522209D0,
C                    I          Xe
     .  2*1.0D0, 0.591824D0, 1.000D0,
C         Cs       Ba
     . 1.000D0, 0.8442D0,
     + 15*1.00D0,
C                  Hf       Ta       W        Re
     .          0.8600D0,0.8400D0,0.8200D0,0.8200D0,
C         Os       Ir       Pt       Au       Hg
     . 0.8200D0,0.8200D0,0.7600D0,0.8200D0,0.7800D0,
     . 25*1.0D0,
C     3RD SET
C                    As                    Sr                  Sn
     + 32*0.0D0, 0.582205D0, 4*0.0D0, 0.593348D0, 11*0.0D0, 0.583986D0,
C                    I           Xe
     .  2*0.0D0, 0.506777D0,  0.0000D0,
C         Cs       Ba
     . 0.000D0, 0.362138D0,
     + 15*0.00D0,
C                  Hf       Ta       W        Re
     .          0.3795D0,0.4016D0,0.3985D0,0.4185D0,
C         Os       Ir       Pt       Au       Hg
     . 0.4239D0,0.4211D0,0.4753D0,0.4131D0,0.4489D0,
     . 25*0.0D0/
C
      DATA (CD(I, 1), I = 1, 105)/
     . 105*1.0D0/
C
      DATA (CD(I, 2), I = 1, 105)/
     . 20*1.0D0,
C          Sc         Ti         V          Cr         Mn
     . 0.35922D0, 0.36461D0, 0.37378D0, 0.40714D0, 0.38984D0,
C         Fe       Co       Ni       Cu       Zn
     . 0.40379D0, 0.41333D0, 0.42120D0, 0.44729D0, 0.43305D0,
     . 8*1.000D0,
C          Y          Zr         Nb         Mo         Tc
     . 0.28655D0, 0.39924D0, 0.44775D0, 0.48013D0, 0.51241D0,
C          Ru         Rh         Pd         Ag         Cd
     . 0.53423D0, 0.54409D0, 0.54902D0, 0.55763D0, 0.56165D0,
     . 8*1.0D0,
C          La         Ce         Pr         Nd         Pm
     . 0.50392D0, 0.50030D0, 0.51040D0, 0.52020D0, 0.52950D0,
C          Sm         Eu         Gd         Tb         Dy
     . 0.53830D0, 0.54660D0, 0.55620D0, 0.56480D0, 0.57240D0,
C          Ho         Er         Tm         Yb         Lu
     . 0.57910D0, 0.58520D0, 0.59050D0, 0.59520D0, 0.59940D0,
C                     Hf         Ta         W          Re
     .            0.6400D0,  0.6600D0,  0.6800D0,  0.7200D0,
C         Os          Ir         Pt         Au         Hg
     .  0.7200D0,  0.7600D0,  0.6800D0,  0.7200D0,  0.7400D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
     .  8*1.0D0,
C          Ac         Th         Pa         U          Np
     . 0.52280D0, 0.52590D0, 0.53150D0, 0.53820D0, 0.54550D0,
C          Pu         Am         Cm         Bk         Cf
     . 0.55302D0, 0.56050D0, 0.57200D0, 0.58230D0, 0.59160D0,
C          Es         Fm         Md         No         Lr
     . 0.59980D0, 0.60720D0, 0.61370D0, 0.61940D0, 0.62440D0,
     . 2 * 0.0D0/
C
      DATA (CD(I, 3), I = 1, 105)/
     . 20*0.0D0,
C          Sc         Ti         V          Cr         Mn
     . 0.76601D0, 0.75561D0, 0.74564D0, 0.73242D0, 0.72965D0,
C         Fe       Co       Ni       Cu       Zn
     . 0.71984D0, 0.71262D0, 0.70658D0, 0.69683D0, 0.69712D0,
     . 8*0.000D0,
C          Y          Zr         Nb         Mo         Tc
     . 0.82468D0, 0.71381D0, 0.66301D0, 0.62892D0, 0.59542D0,
C         Ru       Rh       Pd       Ag       Cd
     . 0.59271D0, 0.58137D0, 0.57561D0, 0.55359D0, 0.54887D0,
     . 8*0.0D0,
C          La         Ce         Pr         Nd         Pm
     . 0.61851D0, 0.63340D0, 0.62570D0, 0.61820D0, 0.61120D0,
C          Sm         Eu         Gd         Tb         Dy
     . 0.60460D0, 0.59860D0, 0.59160D0, 0.58570D0, 0.58080D0,
C          Ho         Er         Tm         Yb         Lu
     . 0.57680D0, 0.57350D0, 0.57090D0, 0.56900D0, 0.56760D0,
C                     Hf         Ta         W          Re
     .            0.5913D0,  0.5736D0,  0.5473D0,  0.5137D0,
C         Os          Ir         Pt         Au         Hg
     .  0.5081D0,  0.4711D0,  0.5557D0,  0.5166D0,  0.4828D0,
C         Tl       Pb       Bi       Po       At       Rn
C         Fr       Ra
     .  8*0.0D0,
C          Ac         Th         Pa         U          Np
     . 0.63360D0, 0.63200D0, 0.62840D0, 0.62370D0, 0.61840D0,
C          Pu         Am         Cm         Bk         Cf
     . 0.61300D0, 0.60750D0, 0.59830D0, 0.59050D0, 0.58370D0,
C          Es         Fm         Md         No         Lr
     . 0.57810D0, 0.57340D0, 0.56960D0, 0.56660D0, 0.56440D0,
     . 2 * 0.0D0/
C
      DATA (CF(I, 1), I = 1, 105)/
     . 105*1.0D0/
C
      DATA (CF(I, 2), I = 1, 105)/
     . 56*1.0D0,
C          La         Ce         Pr         Nd         Pm
     . 0.41259D0, 0.71590D0, 0.71750D0, 0.71960D0, 0.72200D0,
C          Sm         Eu         Gd         Tb         Dy
     . 0.72470D0, 0.72750D0, 0.73320D0, 0.73790D0, 0.74210D0,
C          Ho         Er         Tm         Yb         Lu
     . 0.74610D0, 0.74980D0, 0.75340D0, 0.75700D0, 0.76050D0,
     . 17*1.0D0,
C          Ac         Th         Pa         U          Np
     . 0.77830D0, 0.72410D0, 0.71260D0, 0.70290D0, 0.69440D0,
C          Pu         Am         Cm         Bk         Cf
     . 0.68700D0, 0.68140D0, 0.69190D0, 0.68660D0, 0.88210D0,
C          Es         Fm         Md         No         Lr
     . 0.67840D0, 0.67520D0, 0.67260D0, 0.67050D0, 0.66920D0,
     . 2 * 0.0D0/
C
      DATA (CF(I, 3), I = 1, 105)/
     . 56*0.0D0,
C          La         Ce         Pr         Nd         Pm
     . 0.73974D0, 0.45750D0, 0.45420D0, 0.45100D0, 0.44790D0,
C          Sm         Eu         Gd         Tb         Dy
     . 0.44500D0, 0.44210D0, 0.43750D0, 0.43370D0, 0.43040D0,
C          Ho         Er         Tm         Yb         Lu
     . 0.42730D0, 0.42430D0, 0.42140D0, 0.41860D0, 0.41590D0,
     . 17*1.0D0,
C          Ac         Th         Pa         U          Np
     . 0.46400D0, 0.47030D0, 0.47390D0, 0.47840D0, 0.48320D0,
C          Pu         Am         Cm         Bk         Cf
     . 0.48780D0, 0.49110D0, 0.48230D0, 0.48610D0, 0.48940D0,
C          Es         Fm         Md         No         Lr
     . 0.49220D0, 0.49460D0, 0.49670D0, 0.49840D0, 0.49910D0,
     . 2 * 0.0D0/
C
C-----------------------------------------------------------------------
C D orbital exponents for Na, Mg, Al, Cl from...
C
C              Collins et al.
C              J. Chem. Phys. 64,5142(1976)
C-----------------------------------------------------------------------
      DATA (COVRAD(I), I = 1, 88) /
C         H       He        Li      Be
     . 0.320D0, 0.930D0, 1.230D0, 0.900D0,
C        B        C        N        O        F        Ne
     . 0.820D0, 0.770D0, 0.750D0, 0.730D0, 0.720D0, 0.710D0,
C         Na       Mg
     . 1.540D0, 1.360D0,
C         Al       Si       P        S        Cl       Ar
     . 1.180D0, 1.110D0, 1.060D0, 1.020D0, 0.990D0, 0.980D0,
C         K        Ca
     . 2.030D0, 1.740D0,
C         Sc       Ti       V        Cr       Mn
     . 1.440D0, 1.320D0, 1.220D0, 1.180D0, 1.170D0,
C         Fe       Co       Ni       Cu       Zn
     . 1.170D0, 1.160D0, 1.150D0, 1.170D0, 1.250D0,
C         Ga       Ge       As       Se       Br       Kr
     . 1.260D0, 1.220D0, 1.200D0, 1.160D0, 1.140D0, 1.120D0,
C         Rb       Sr
     . 2.160D0, 1.910D0,
C         Y        Zr       Nb       Mo       Tc
     . 1.620D0, 1.450D0, 1.340D0, 1.300D0, 1.270D0,
C         Ru       Rh       Pd       Ag       Cd
     . 1.250D0, 1.250D0, 1.280D0, 1.340D0, 1.480D0,
C         In       Sn       Sb       Te       I        Xe
     . 1.440D0, 1.410D0, 1.400D0, 1.360D0, 1.330D0, 1.310D0,
C         Cs       Ba       La
     . 2.350D0, 1.980D0, 1.690D0,
C         Ce       Pr       Nd       Pm       Sm       Eu       Gd
     . 1.650D0, 1.650D0, 1.640D0, 1.630D0, 1.620D0, 1.850D0, 1.610D0,
C         Tb       Dy       Ho       Er       Tm       Yb       Lu
     . 1.590D0, 1.590D0, 1.580D0, 1.570D0, 1.560D0, 1.560D0, 1.560D0,
C                  Hf       Ta       W        Re
     .          1.440D0, 1.340D0, 1.300D0, 1.280D0,
C         Os       Ir       Pt       Au       Hg
     . 1.260D0, 1.270D0, 1.300D0, 1.340D0, 1.490D0,
C         Tl       Pb       Bi       Po       At       Rn
     . 1.480D0, 1.470D0, 1.460D0, 1.460D0, 1.450D0, 0.000D0,
C         Fr       Ra
     . 0.000D0, 0.000D0/
C
      DATA (COVRAD(I), I = 89, 105) /
C
C Using the Lanthanides' values is probably the best approximation
C         Ac
     . 1.690D0,
C         Th       Pa       U        Np       Pu       Am       Cm
     . 1.650D0, 1.650D0, 1.640D0, 1.630D0, 1.620D0, 1.850D0, 1.610D0,
C         Bk       Cf       Es       Fm       Md       No       Lr
     . 1.590D0, 1.590D0, 1.580D0, 1.570D0, 1.560D0, 1.560D0, 1.560D0,
     . 2 * 0.000D0/
C
C
c-----------------------------------------------------------------
c
c    Van der Waals radii from A. Bondi, J.Phys.Chem. 68 (1964) 441
c
c
c    (missing radii are represented by a value of 1.99 A)
c   
c-----------------------------------------------------------------
c                 
      DATA VDWR/
c
c       H   He   
     .1.20,1.40,
c
c      Li   Be    B    C    N    O    F   Ne 
     .1.99,1.99,1.99,1.70,1.55,1.52,1.47,1.54,
c
c      Na   Mg   Al   Si    P    S   Cl   Ar 
     .1.99,1.99,1.99,2.10,1.80,1.80,1.75,1.88, 
c
c       K   Ca   Sc   Ti    V   Cr   Mn   Fe   Co   Ni   Cu   Zn            
     .1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,
c
c      Ga   Ge   As   Se   Br   Kr 
     .1.99,1.99,1.85,1.90,1.85,2.02, 
c
c      Rb   Sr    Y   Zr   Nb   Mo   Tc   Ru   Rh   Pd   Ag   Cd 
     .1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,1.99,
c
c      In   Sn   Sb   Te    I   Xe 
     .1.99,1.99,1.99,2.06,1.98,2.16/ 
c
c
c-----------------------------------------------
c
c   Van der Waals radii from MOLDEN.2.4 (fishy?)
c
c-----------------------------------------------
c
c      data vdwr/0.430d0,0.200d0,0.880d0,0.550d0,1.030d0,0.900d0,0.880d0,
c     1          0.880d0,0.840d0,0.200d0,1.170d0,1.300d0,1.550d0,1.400d0,
c     2          1.250d0,1.220d0,1.190d0,0.200d0,1.530d0,1.190d0,1.640d0,
c     3          1.670d0,1.530d0,1.550d0,1.550d0,1.540d0,1.430d0,1.700d0,
c     4          1.720d0,1.650d0,1.420d0,1.370d0,1.410d0,1.420d0,1.410d0,
c     5          0.200d0,1.670d0,1.320d0,1.980d0,1.760d0,1.680d0,1.670d0,
c     6          1.550d0,1.600d0,1.650d0,1.700d0,1.790d0,1.890d0,1.830d0,
c     7          1.660d0,1.660d0,1.670d0,1.600d0,0.200d0,1.870d0,1.540d0,
c     8          2.070d0,2.030d0,2.020d0,2.010d0,2.000d0,2.000d0,2.190d0,
c     9          1.990d0,1.960d0,1.950d0,1.940d0,1.930d0,1.920d0,2.140d0,
c     &          1.920d0,1.770d0,1.630d0,1.570d0,1.550d0,1.570d0,1.520d0,
c     &          1.700d0,1.700d0,1.900d0,1.750d0,1.740d0,1.740d0,1.880d0,
c     &          0.200d0,0.200d0,0.200d0,2.100d0,2.080d0,1.990d0,1.810d0,
c     &          1.780d0,1.750d0,0.200d0,1.710d0,0.200d0,0.200d0,1.730d0/
c
c        do j = 1,nat
c          write(6,*) j,atom_num(j)
c        end do
c        stop
c
        DO J=1, NAT
         DO I=1,105
          IF (ATOM_NUM(J).EQ.I) THEN
           SEXPO(I,1)=AS(I,1)
           SEXPO(I,2)=AS(I,2)
           SEXPO(I,3)=AS(I,3)
           SEXPO(I,4)=CS(I,1)
           SEXPO(I,5)=CS(I,2)
           SEXPO(I,6)=CS(I,3)
           PEXPO(I,1)=AP(I,1)
           PEXPO(I,2)=AP(I,2)
           PEXPO(I,3)=AP(I,3)
           PEXPO(I,4)=CP(I,1)
           PEXPO(I,5)=CP(I,2)
           PEXPO(I,6)=CP(I,3)
           DEXPO(I,1)=AD(I,1)
           DEXPO(I,2)=AD(I,2)
           DEXPO(I,3)=AD(I,3)
           DEXPO(I,4)=CD(I,1)
           DEXPO(I,5)=CD(I,2)
           DEXPO(I,6)=CD(I,3)
           FEXPO(I,1)=AF(I,1)
           FEXPO(I,2)=AF(I,2)
           FEXPO(I,3)=AF(I,3)
           FEXPO(I,4)=CF(I,1)
           FEXPO(I,5)=CF(I,2)
           FEXPO(I,6)=CF(I,3)
          ENDIF
         ENDDO  ! I LOOP
        ENDDO ! J LOOP
c
        RETURN
        END                                          
