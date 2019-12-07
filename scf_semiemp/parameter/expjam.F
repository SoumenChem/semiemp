C---------------------------------------------------------------------

      SUBROUTINE expjam(KAN,nat,TNB,NB,sexpo,pexpo,dexpo,fexpo
     1   ,ALPHA,NW,NU,NP,PQN,BAN,LQNO,NFIRST,ISW1,NTYP,orbtype,
     1 atom_type,atnum)

C-----------------------------------------------------------------------
C     THIS SUBROUTINE FILLS IN THE ORBITAL EXPONENTS, ALPHA(I,J).
C     EXCEPT WHERE NOTED THESE ARE SLATER EXPONENTS.
C-----------------------------------------------------------------------

C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C     PARAMETER  (ZERO = 0.0D0, ONE = 1.0D0)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/cofexp.cmn'
c     INCLUDE '../include/z.cmn'
c     INCLUDE '../include/switch.cmn'
C----------------------------------------------------------------------
C     VALENCE BASIS SET
c      implicit none
c
      integer nat !number of atom
      integer TNB !total number of basis function
      character*8 orbtype(TNB)
      character*8 atom_type(TNB)
      character*8 AOS(16)
      character*8 ATOM(105)
      integer KAN(nat) !atomic number array
      integer BAN(TNB) !atomic number for basis set
      integer NB(nat) !number of basis function in each atom
      integer NP(nat) !priciple quantum number for each atom
      integer NFIRST(NAT) !THE INDEX OF THE FIRST AO OF ATOM
      integer NTYP(*)
      integer atnum(*)
      double precision sexpo(105,6)
      double precision pexpo(105,6)
      double precision dexpo(105,6)
      double precision fexpo(105,6)
      double precision ALPHA(TNB,6)
      integer Basis
      integer I,J,K,IJ,NT,NLQN
      PARAMETER  (ZERO = 0.0D0, ONE = 1.0D0)
      integer NW(TNB)!type of basis array
      integer NU(TNB) !type of atom of the basis array
      integer PQN(TNB) !Principle quantum number array
      integer NL(16) !type of basis function
      integer LQNO(*) !L quantum number array
      dimension nlqn(16)
      data NL/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/
      data nlqn / 0,3*1,5*2,7*3/
      DATA AOS/'S', 'P(X)', 'P(Y)', 'P(Z)',
     X   'D(Z2)', 'D(X2-Y2)', 'D(XY)', 'D(XZ)',
     X   'D(YZ)', 'F(Z3)', 'F(XZ2)', 'F(YZ2)',
     X   'Z(X2-Y2)', 'F(XYZ)', 'X3-3XY2', '3YX2-Y3'/
      DATA ATOM/' H','Ct','Li','Be',' B',' C',' N',' O',' F','Ne','Na',
     1 'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr',
     2 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     3 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     4 'Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     5 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re',
     6 'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     7 'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     8 'Md','No','Lr','XX','XX'/

c     write(6,*) "i am in expjam"
        Basis= 0
      LL = 1
      NFIRST(1) = 1
      DO I = 1, nat
       NFIRST(I) = LL
       LL = LL + NB(I) 
      END DO
      DO 40 I = 1,nat
         J=NB(I)
c       write(6,*)"J",J
         K = KAN(I)
        DO IJ= 1,J
            Basis=Basis+1
            NU(Basis)=I
            PQN(Basis)=NP(I)
            BAN(Basis)=K

c        IF (ALPHA(Basis,1).GT.ZERO) THEN
c           write(6,*)"NT1", NT
c           ALPHA(Basis,4) = ONE
c           ALPHA(Basis,3) = ZERO
c           ALPHA(Basis,2) = ALPHA(Basis,1)
c           ALPHA(Basis,5) = ONE
c           ALPHA(Basis,6) = ZERO
c           GO TO 40

c        ELSE
            IF (K.EQ.0) THEN
               WRITE (6,9000) K
               Write(6,*)"dummy atom"
            END IF
            NT = NL(IJ)
c           NLQN = nlqn(IJ)
            NW(Basis)=NT
            LQNO(Basis) = NLQN(IJ)
            atom_type(Basis) = ATOM(K)
            orbtype(Basis) = AOS(IJ)
            atnum(basis) = I
c           write(6,*) "ORBTYPE =",orbtype(Basis)
c           write(6,*) "ATOMTYPE =", atom_type(basis)
C        CORRECT PRINC Q. NO.
C        FOR 4S 4P 3D (NTYP=4): nS, nP, (n-1)d
C        FOR 6S 6P 5D 4F (NTYP=5): ns, np, (n-1)d, (n-2)f
C

         IF (NTYP(I).GE.4 .AND. NW(BASIS).GT.3) PQN(BASIS) = NP(I) - 1
         IF (NTYP(I).EQ.5 .AND. NW(BASIS).GT.8) PQN(BASIS) = NP(I) - 2
c           write(6,*)"NLQN",NLQN(IJ)
c           write(6,*)"PQN",PQN(Basis)
            IF (ISW1.LT.50) THEN
            IF  (NT.GT.0.AND.NT.LE.3) THEN
c           write(6,*)"exponent watson"
C----------------------------------------------------------------------
C             3D EXPONENT FOR R=0 FROM WATSON
C----------------------------------------------------------------------
   10          ALPHA(Basis,1) = pexpo(K,1)
               ALPHA(Basis,2) = pexpo(K,2)
               ALPHA(Basis,3) = pexpo(K,3)
               ALPHA(Basis,4) = pexpo(K,4)
               ALPHA(Basis,5) = pexpo(K,5)
               ALPHA(Basis,6) = pexpo(K,6)
            ELSE IF (NT.EQ.0) THEN
c           write(6,*)"exponent zener"
C----------------------------------------------------------------------
C             HYDROGEN EXPONENT SET TO 1.2
C             ZERNER 'BEST SINGLE' EXPONENT FOR 4S, 4P
C----------------------------------------------------------------------
   20          ALPHA(Basis,1) = sexpo(K,1)
               ALPHA(Basis,2) = sexpo(K,2)
               ALPHA(Basis,3) = sexpo(K,3)
               ALPHA(Basis,4) = sexpo(K,4)
               ALPHA(Basis,5) = sexpo(K,5)
               ALPHA(Basis,6) = sexpo(K,6)
            ELSE IF (NT.LE.8) THEN
c              write(6,*)"exponent watson 2"
C----------------------------------------------------------------------
C             3D EXPONENT FOR R=0 FROM WATSON
C----------------------------------------------------------------------
   30          ALPHA(Basis,1) = dexpo(K,1)
               ALPHA(Basis,2) = dexpo(K,2)
               ALPHA(Basis,3) = dexpo(K,3)
               ALPHA(Basis,4) = dexpo(K,4)
               ALPHA(Basis,5) = dexpo(K,5)
               ALPHA(Basis,6) = dexpo(K,6)
            ELSE
C----------------------------------------------------------------------
C             F ORBITAL EXPONENTS
C----------------------------------------------------------------------
               ALPHA(I,1) = fexpo(K,1)
               ALPHA(I,2) = fexpo(K,2)
               ALPHA(I,3) = fexpo(K,3)
               ALPHA(I,4) = fexpo(K,4)
               ALPHA(I,5) = fexpo(K,5)
               ALPHA(I,6) = fexpo(K,6)
            END IF
         END IF
       END DO
   40 CONTINUE
      RETURN
C----------------------------------------------------------------------
 9000 FORMAT (' **** NO EXPONENTS FOR CENTER',I5,' ****')

      END
