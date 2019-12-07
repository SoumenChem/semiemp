      SUBROUTINE FETCH(NB,NW,NU,NP,KAN,ALPHA,NTYP,SCG1SP,SCF2PP,
     1 SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,SCG3SF,SCG2PF,SCG4PF,
     2 SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,SCF2FF,SCF4FF,SCF6FF
     3 ,R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,R2PFFF,R2PPPF,R2SDFF,
     4 R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF)

C
C     MODIFIED BY N. ROESCH GAINESVILLE        VERSION : 29.10.84
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/ftch.cmn'
      PARAMETER (MXA=1200,ISPEC = 0,MXB = 4000)
      DIMENSION KAN(NB),NW(NB),NU(NB),ALPHA(MXB,6),NP(NB),NTYP(*)
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
      double precision FACT(35)
c     logical RRRR
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
C****************************************************************
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     .   ATOM(100),NIN(MXB),
C***     .   NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     .   NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     .   Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
C
c     write(6,*)"i am in fetch"!,(NTYP(I),I=1,10)
      IZ = 0
      I1 = 1
      I2 = 2
      I3 = 3
      I4 = 4
      IATOM = 0
      DO 50 J = 1, NB
         L = NW(J)
         I = NU(J)
         IF (IATOM.EQ.I) GO TO 50
         IATOM = I
c     IF (L .NE. 0) GOTO 15
         NP1 = NP(J)
         NP2 = NP(J+1)
         AMU = ALPHA(J,1)
         BMU = ALPHA(J+1,1)
         K = NTYP(I) + 1
         IF(ISPEC.GE.0) THEN
           KAA = KAN(I)
C          THIS USES AB-INITIO SLATER-CONDON FACTORS FOR 2ND TRANSITION
C          SERIES
           IF(KAA.GT.36.AND.KAA.LE.48) KAA=0
         ELSE
           KAA = 0
         ENDIF
C     ******
C     CONVENTIONAL PPP DOES NOT CONTAIN THE HIGHER SLATER CONDON
C     FACTORS. IF YOU WANT THEM ACTIVATE NEXT CARD.
C      IF(ISW1.EQ.7) GO TO 2
C     ******
c        write(6,*)"i am before k",KAA
         GO TO (40,40,10,20,20,30), K

   10    CONTINUE
         SCG1SP(I) = RRRR(1,1,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP1,IZ,IZ,AMU,FACT)/3.0D0
         SCF2PP(I) = RRRR(2,2,KAA,NP2,I1,IZ,BMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP2,I1,IZ,BMU,FACT)/25.0D0
C
C     WRITE(6,*) I, SCG1SP(I), SCF2PP(I)
C
         GO TO 40

   20    CONTINUE
         NP3 = NP(J+4)
         CMU = ALPHA(J+4,1)
         SCG1SP(I) = RRRR(1,1,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP1,IZ,IZ,AMU,FACT)/3.0D0
         SCF2PP(I) = RRRR(2,2,KAA,NP2,I1,IZ,BMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP2,I1,IZ,BMU,FACT)/25.0D0
         SCG2SD(I) = RRRR(2,3,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP1,IZ,IZ,AMU,FACT)/5.0D0
         SCG1PD(I) = RRRR(1,4,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP2,I1,IZ,BMU,FACT)/15.0D0
         SCF2PD(I) = RRRR(2,5,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP2,I1,IZ,
     .               BMU,NP3,I2,IZ,CMU,FACT)/35.0D0
         SCG3PD(I) = RRRR(3,6,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP2,I1,IZ,BMU,FACT)/245.0D0
         SCF2DD(I) = RRRR(2,7,KAA,NP3,I2,IZ,CMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP3,I2,IZ,CMU,FACT)/49.0D0
         SCF4DD(I) = RRRR(4,8,KAA,NP3,I2,IZ,CMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP3,I2,IZ,CMU,FACT)/441.0D0
C
C   S, P, D Basis R-integrals
C
         R1SPPD(I) = RRRR(1,21,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP2,I1,
     .               IZ,BMU,NP3,I2,IZ,CMU,FACT)
         R2SDDD(I) = RRRR(2,22,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP3,I2,
     .               IZ,CMU,NP3,I2,IZ,CMU,FACT)
         R2SDPP(I) = RRRR(2,23,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP3,I2,
     .               IZ,CMU,NP2,I1,IZ,BMU,FACT)
C
c     WRITE(6,610) I, SCG1SP(I), SCF2PP(I), SCG2SD(I), SCG1PD(I),
c    1                SCF2PD(I), SCG3PD(I), SCF2DD(I), SCF4DD(I),
c    2             I, R1SPPD(I), R2SDDD(I), R2SDPP(I)
C
         GO TO 40

   30    CONTINUE
C     S, P, D, F, BASIS
         NP3 = NP(J+4)
         DMU = ALPHA(J+9,1)
         NP4 = NP(J+9)
         CMU = ALPHA(J+4,1)
         SCG1SP(I) = RRRR(1,1,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP1,IZ,IZ,AMU,FACT)/3.0D0
         SCF2PP(I) = RRRR(2,2,KAA,NP2,I1,IZ,BMU,NP2,I1,IZ,BMU,NP2,I1,IZ,
     .               BMU,NP2,I1,IZ,BMU,FACT)/25.0D0
         SCG2SD(I) = RRRR(2,3,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP1,IZ,IZ,AMU,FACT)/5.0D0
         SCG1PD(I) = RRRR(1,4,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP2,I1,IZ,BMU,FACT)/15.0D0
         SCF2PD(I) = RRRR(2,5,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP2,I1,IZ,
     .               BMU,NP3,I2,IZ,CMU,FACT)/35.0D0
         SCG3PD(I) = RRRR(3,6,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP2,I1,IZ,BMU,FACT)/245.0D0
         SCF2DD(I) = RRRR(2,7,KAA,NP3,I2,IZ,CMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP3,I2,IZ,CMU,FACT)/49.0D0
         SCF4DD(I) = RRRR(4,8,KAA,NP3,I2,IZ,CMU,NP3,I2,IZ,CMU,NP3,I2,IZ,
     .               CMU,NP3,I2,IZ,CMU,FACT)/441.0D0
C     F ORBITAL ELEMENTS
C
CMK90
C     IMPROVED G3SF INTEGRALS
CMK90
         SCG3SF(I) = RRRR(3,9,KAA,NP1,IZ,IZ,AMU,NP4,I3,IZ,DMU,NP4,I3,IZ,
     .               DMU,NP1,IZ,IZ,AMU,FACT)/7.0D0
CMK89
         FACTOR = ((KAA-58)*0.165D0) + 1.0D0
         FACTOR = FACTOR*5.5D0
         SCG3SF(I) = SCG3SF(I)*FACTOR
CMK89
         SCG2PF(I) = RRRR(2,10,KAA,NP2,I1,IZ,BMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,
CMK89     1 NP2,I1,IZ,BMU,FACT)/175.0D0
     .               NP2,I1,IZ,BMU,FACT)/175.0D0*6.0D0
CMK89
         SCG4PF(I) = RRRR(4,11,KAA,NP2,I1,IZ,BMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,
CMK89     1 NP2,I1,IZ,BMU,FACT)/189.0D0
     .               NP2,I1,IZ,BMU,FACT)/189.0D0*6.0D0
CMK89
         SCG1DF(I) = RRRR(1,12,KAA,NP3,I2,IZ,CMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,
CMK89     1 NP3,I2,IZ,CMU,FACT)/35.0D0
CMK89     1 NP3,I2,IZ,CMU,FACT)/35.0D0 * 1.57
     .               NP3,I2,IZ,CMU,FACT)/35.0D0*1.4
         SCG3DF(I) = RRRR(3,13,KAA,NP3,I2,IZ,CMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,
CMK89     1 NP3,I2,IZ,CMU,FACT)/315.0D0
CMK89     1 NP3,I2,IZ,CMU,FACT)/315.0D0 * 1.93
     .               NP3,I2,IZ,CMU,FACT)/315.0D0*1.4
         SCG5DF(I) = RRRR(5,14,KAA,NP3,I2,IZ,CMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,
CMK89     1 NP3,I2,IZ,CMU,FACT)/15246.0D0
CMK89     1 NP3,I2,IZ,CMU,FACT)/15246.0D0 * 1.93
     .               NP3,I2,IZ,CMU,FACT)/15246.0D0*1.4
CMK89
         SCF2PF(I) = RRRR(2,15,KAA,NP2,I1,IZ,BMU,NP4,I3,IZ,DMU,NP2,I1,
     .               IZ,BMU,NP4,I3,IZ,DMU,FACT)/75.0D0
         SCF2DF(I) = RRRR(2,16,KAA,NP3,I2,IZ,CMU,NP4,I3,IZ,DMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)/105.0D0
         SCF4DF(I) = RRRR(4,17,KAA,NP3,I2,IZ,CMU,NP4,I3,IZ,DMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)/693.0D0
         SCF2FF(I) = RRRR(2,18,KAA,NP4,I3,IZ,DMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)/225.0D0
         SCF4FF(I) = RRRR(4,19,KAA,NP4,I3,IZ,DMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)/1089.0D0
         SCF6FF(I) = RRRR(6,20,KAA,NP4,I3,IZ,DMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)/736164.0D0
C
C   S, P, D Basis R-integrals
C
         R1SPPD(I) = RRRR(1,21,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP2,I1,
     .               IZ,BMU,NP3,I2,IZ,CMU,FACT)
         R2SDDD(I) = RRRR(2,22,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP3,I2,
     .               IZ,CMU,NP3,I2,IZ,CMU,FACT)
         R2SDPP(I) = RRRR(2,23,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP3,I2,
     .               IZ,CMU,NP2,I1,IZ,BMU,FACT)
C
C   R-integrals from the F orbitals
C
         R1PDDF(I) = RRRR(1,23,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)
         R1SPDF(I) = RRRR(1,24,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP2,I1,
     .               IZ,BMU,NP4,I3,IZ,DMU,FACT)
         R2PFDD(I) = RRRR(2,25,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP4,I3,
     .               IZ,DMU,NP3,I2,IZ,CMU,FACT)
         R2PFFF(I) = RRRR(2,26,KAA,NP2,I1,IZ,BMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)
         R2PPPF(I) = RRRR(2,27,KAA,NP2,I1,IZ,BMU,NP2,I1,IZ,BMU,NP2,I1,
     .               IZ,BMU,NP4,I3,IZ,DMU,FACT)
         R2SDFF(I) = RRRR(2,28,KAA,NP1,IZ,IZ,AMU,NP4,I3,IZ,DMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)
         R2SDPF(I) = RRRR(2,29,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)
         R3PDDF(I) = RRRR(3,30,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP3,I2,
     .               IZ,CMU,NP4,I3,IZ,DMU,FACT)
         R3SFDF(I) = RRRR(3,31,KAA,NP1,IZ,IZ,AMU,NP3,I2,IZ,CMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)
         R3SFPD(I) = RRRR(3,32,KAA,NP1,IZ,IZ,AMU,NP2,I1,IZ,BMU,NP4,I3,
     .               IZ,DMU,NP3,I2,IZ,CMU,FACT)
         R4PFDD(I) = RRRR(4,33,KAA,NP2,I1,IZ,BMU,NP3,I2,IZ,CMU,NP4,I3,
     .               IZ,DMU,NP3,I2,IZ,CMU,FACT)
         R4PFFF(I) = RRRR(4,34,KAA,NP2,I1,IZ,BMU,NP4,I3,IZ,DMU,NP4,I3,
     .               IZ,DMU,NP4,I3,IZ,DMU,FACT)
C
C     WRITE(6,620) I,SCG1SP(I),SCF2PP(I),SCG2SD(I),SCG1PD(I),SCF2PD(I),
C    1               SCG3PD(I),SCF2DD(I),SCF4DD(I),SCG3SF(I),SCG2PF(I),
C    2               SCG4PF(I),SCG1DF(I),SCG3DF(I),SCG5DF(I),SCF2PF(I),
C    3               SCF2DF(I),SCF4DF(I),SCF2FF(I),SCF4FF(I),SCF6FF(I),
C    4             I,R1SPPD(I),R2SDDD(I),R2SDPP(I),
C    5               R1PDDF(I),R1SPDF(I),R2PFDD(I),R2PFFF(I),R2PPPF(I),
C    6               R2SDFF(I),R2SDPF(I),R3PDDF(I),R3SFDF(I),R3SFPD(I),
C    7               R4PFDD(I),R4PFFF(I)
C
   40    CONTINUE
   50 CONTINUE
C
      RETURN
C
C 600 FORMAT('0 >>>>>  F AND G INTEGRALS FOR ATOM NO.',I5,/
C    O       ,12X, 'G1SP', 8X, 'F2PP', /, 4X,2D12.4)
c 610 FORMAT('0 >>>>>  F AND G INTEGRALS FOR ATOM NO.',I5,/,
c    O    12X, 'G1SP', 8X, 'F2PP', 8X, 'G2SD', 8X, 'G1PD', 8X, 'F2PD',
c    1     8X, 'G3PD', 8X, 'F2DD', 8X, 'F4DD', 8X, /, 4X, 8D12.4,/,
c    2       '0 >>>>>     R    INTEGRALS FOR ATOM NO.',I5,/,
c    3    11X,'R1SPPD',6X,'R2SDDD',6X,'R2SDPP',/,4X,3D12.4)
C 620 FORMAT('0 >>>>>  F AND G INTEGRALS FOR ATOM NO.',I5,/,
C    O    12X, 'G1SP', 8X, 'F2PP', 8X, 'G2SD', 8X, 'G1PD', 8X, 'F2PD',
C    1     8X, 'G3PD' ,8X, 'F2DD', 8X, 'F4DD', 8X, 'G3SF', 8X, 'G2PF'/
C    $     4X, 10D12.4,/,
C    2    12X, 'G4PF', 8X, 'G1DF', 8X, 'G3DF', 8X, 'G5DF', 8X, 'F2PF',
C    3     8X, 'F2DF', 8X, 'F4DF', 8X, 'F2FF', 8X, 'F4FF', 8X, 'F6FF',/,
C    $     4X,10D12.4,/,
C    4       '0 >>>>>     R    INTEGRALS FOR ATOM NO.',I5,/,
C    5    11X,'R1SPPD',6X,'R2SDDD',6X,'R2SDPP',/,4X,3D12.4,/,
C    6    11X,'R1PDDF',6X,'R1SPDF',6X,'R2PFDD',6X,'R2PFFF',6X,'R2PPPF',
C    7     6X,'R2SDFF',6X,'R2SDPF',6X,'R3PDDF',6X,'R3SFDF',6X,'R3SFPD',/
C    $     ,4X,10D12.4,/,
C    8    11X,'R4PFDD',6X,'R4PFFF',/,4X,2D12.4)
C
      END
c-----------------------------------------------------------------------------
c----------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION RRRR(K,KI,NZ,N1,L1,M1,
     1       A1,N2,L2,M2,A2,N3,L3,M3,A3,N4,L4,M4,A4,FACT)
C     COPYRIGHTED 1982, 1984, 1990
C     M. C. ZERNER UPPSALA
C     M. C. ZERNER FLORIDA
C     +++++++  ZINDO +++++
C     MODIFIED BY N. ROESCH GAINESVILLE        VERSION : 29.10.84
C     EVALUATES SLATER CONDON RADIAL INTEGRALS R(K)(1,2/1,2) BETWEEN
C     S.T.O.'S CHARACTERIZED BY Q. NO.'S N, L, M AND EXPONENT A.
C     IF NZ THE ATOMIC NUMBER IS POSITIVE, EMPIRICAL VALUES ARE USED.
C
CNR   IF NZ IS NEGATIVE, THE TABLE LOOK UP IS GOVERNED BY THE POINTER
CNR   ARRAY INDX. THE INDEX (KI ,K,L1,  L4)    POINTS TO THE
CNR   SLATER-CONDON INTEGRAL (F/G)K(L1+1,L4+1)
C     KI = THE INDEX OF STORAGE OF SPK, SEE BELOW
C  
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      parameter (MPE = 105,ONE=1.0d0,ZERO =0.0d0,TWO=2.0d0,THREE=3.0d0)
      parameter ( ENUBAR   = 2.1947463067D5)
      DIMENSION SPK(20*MPE),SCR(20,MPE),FACT(*),INDX(2,6,4,4)
      EQUIVALENCE (SPK(1),SCR(1,1))
C     SCR(I,M) HOLDS INTEGRALS FOR ATOM NO. M
C     STORAGE OF SPK,G1(S,P),F2(P,P),G2(S,D),G1(P,D),F2(P,D),G3(P,D),
C     F2(D,D),F4(D,D),
C     G3(S,F),G2(P,F),G4(P,F),G1(D,F),G3(D,F),G5(D,F),
C     F2(P,F),F2(D,F),F4(D,F),F2(F,F),F4(F,F),F6(F,F).
C     TWENTY VALUES EACH ATOM, THEN NEXT ATOM.
C     SLATER-CONDON FACTORS FOR THE TRANSITION METALS ARE THOSE FOR THE
C     NEUTRAL ATOMS.
C     VALUES FOR NA-K ARE ZERNER'S
C     VALUES FOR CA-CU ARE FROM ANNO + TERUYA J.CHEM.PHYS.,52,2840(1970)
C     VALUES FOR ZN ARE FROM THE FORMULA OF ANNO + TERUYA
C     VALUES OF F2PP FOR CA-ZN ARE FROM HINZE + JAFFE
C========== 810520 - HJAAJ
C  FIRST ATOM IS Z = 3 (LI).  THE UNIT IS INVERSE CM.
C
C **** THE UNITS ARE CM-1 ****  THE UNITS ARE CM-1 *****
C
C  REF. ABOVE: HINZE + JAFFE, J.CHEM.PHYS.,38,1834(1963)
C==========
      SAVE INDX
      DATA (SCR(I, 1),I=1,20)/20*0.D0/, (SCR(I, 2),I=1,20)/20*0.D0/,
     $     (SCR(I, 3),I=1,20)/20194.D0,10944.D0,18*0.D0/,
     $     (SCR(I, 4),I=1,20)/30876.D0,21425.D0,18*0.D0/,
     $     (SCR(I, 5),I=1,20)/43566.D0,28075.D0,18*0.D0/,
     $     (SCR(I, 6),I=1,20)/55635.D0,36375.D0,18*0.D0/,
     $     (SCR(I, 7),I=1,20)/72255.D0,52100.D0,18*0.D0/,
     $     (SCR(I, 8),I=1,20)/95298.D0,55675.D0,18*0.D0/,
     $     (SCR(I, 9),I=1,20)/116825.D0,69309.D0,18*0.D0/,
     $     (SCR(I,10),I=1,20)/20*0.D0/,
     $     (SCR(I,11),I=1,20)/13450.D0,6000.D0,3042.D0,9224.D0,
     $               11000.D0,5506.D0,7477.D0,4876.D0,12*0.D0/,
     $     (SCR(I,12),I=1,20)/19974.D0,26400.D0,3885.D0,3070.D0,
     $               4416.D0,1833.D0,7477.D0,4876.D0,12*0.D0/,
     $     (SCR(I,13),I=1,20)/27093.D0,12925.D0,1426.D0,1896.D0,
     $               3461.D0,1130.D0,7477.D0,4876.D0,12*0.D0/,
     $     (SCR(I,14),I=1,20)/38814.D0,18250.D0,14496.D0,2506.D0,
     $               19675.D0,11828.D0,19438.D0,12677.D0,12*0.D0/,
     $     (SCR(I,15),I=1,20)/8451.D0,23775.D0,19081.D0,25334.D0,
     $               23176.D0,15124.D0,23925.D0,15603.D0,12*0.D0/
      DATA (SCR(I,16),I=1,20)/24807.D0,36600.D0,25972.D0,34486.D0,
     $               29173.D0,20587.D0,28411.D0,18529.D0,12*0.D0/,
     $     (SCR(I,17),I=1,20)/71000.D0,52000.D0,17131.D0,2274.D0,
     $               5102.D0,1358.D0,12860.D0,8387.D0,12*0.D0/,
     $     (SCR(I,18),I=1,20)/20*0.D0/,
     $     (SCR(I,19),I=1,20)/8960.D0,4000.D0,18*0.D0/,
     $     (SCR(I,20),I=1,20)/12600.D0,2325.D0,3730.D0,5890.D0,
     $               4480.D0,-238.D0,18900.D0,9500.D0,12*0.D0/,
     $     (SCR(I,21),I=1,20)/12100.D0,5000.D0,5870.D0,5650.D0,
     $               11000.D0,2210.D0,29500.D0,14600.D0,12*0.D0/,
     $     (SCR(I,22),I=1,20)/13100.D0,5500.D0,6200.D0,7320.D0,
     $               13700.D0,10300.D0,44900.D0,29700.D0,12*0.D0/,
     $     (SCR(I,23),I=1,20)/15100.D0,6000.D0,6240.D0,5180.D0,
     $               11200.D0,1710.D0,50800.D0,35400.D0,12*0.D0/,
     $     (SCR(I,24),I=1,20)/14400.D0,6500.D0,5220.D0,5580.D0,
     $               11400.D0,297.D0,63500.D0,36800.D0,12*0.D0/,
     $     (SCR(I,25),I=1,20)/18900.D0,7000.D0,6110.D0,1240.D0,
     $               8010.D0,4970.D0,66000.D0,37900.D0,12*0.D0/,
     $     (SCR(I,26),I=1,20)/16300.D0,7500.D0,6640.D0,2450.D0,
     $               5020.D0,3520.D0,61000.D0,38400.D0,12*0.D0/
      DATA (SCR(I,27),I=1,20)/22700.D0,8000.D0,6340.D0,3170.D0,
     $               6290.D0,2260.D0,64500.D0,48100.D0,12*0.D0/,
     $     (SCR(I,28),I=1,20)/19400.D0,8500.D0,6700.D0,3010.D0,
     $               6050.D0,3250.D0,79800.D0,53300.D0,12*0.D0/,
     $     (SCR(I,29),I=1,20)/20700.D0,9000.D0,4460.D0,5620.D0,
     $               10700.D0,6930.D0,85980.D0,57965.D0,12*0.D0/,
     $     (SCR(I,30),I=1,20)/20400.D0,9500.D0,8890.D0,4070.D0,
     $               14170.D0,7220.D0,91100.D0,61950.D0,12*0.D0/,
     $     ((SCR(I,M),I=1,20),M=31,43)/260*0.D0/
     $     (SCR(I,44),I=1,20)/20000.D0,9034.D0,7985.D0,2662.D0,
     $               11615.D0,1452.D0,49767.D0,37345.D0,12*0.0D0/,
C    $     (SCR(I,44),I=1,20)/20*0.0D0/,
     $     ((SCR(I,M),I=1,20),M=45,105)/1220*0.D0/
C=== 810521 === HJAAJ === THIS IS THE LAST CONTINUATION CARD ALLOWED]]]]
CNR
      DATA
     $ INDX(1,1,1,1)/0/, INDX(1,1,1,2)/0/, INDX(1,1,1,3)/0/,
     $ INDX(1,1,1,4)/0/,
     $ INDX(1,1,2,1)/0/, INDX(1,1,2,2)/0/, INDX(1,1,2,3)/0/,
     $ INDX(1,1,2,4)/0/,
     $ INDX(1,1,3,1)/0/, INDX(1,1,3,2)/0/, INDX(1,1,3,3)/0/,
     $ INDX(1,1,3,4)/0/,
     $ INDX(1,1,4,1)/0/, INDX(1,1,4,2)/0/, INDX(1,1,4,3)/0/,
     $ INDX(1,1,4,4)/0/,
     $ INDX(1,2,1,1)/0/, INDX(1,2,1,2)/0/, INDX(1,2,1,3)/0/,
     $ INDX(1,2,1,4)/0/
      DATA
     $ INDX(1,2,2,1)/0/, INDX(1,2,2,2)/2/, INDX(1,2,2,3)/5/,
     $ INDX(1,2,2,4)/15/,
     $ INDX(1,2,3,1)/0/, INDX(1,2,3,2)/0/, INDX(1,2,3,3)/7/,
     $ INDX(1,2,3,4)/16/,
     $ INDX(1,2,4,1)/0/, INDX(1,2,4,2)/0/, INDX(1,2,4,3)/0/,
     $ INDX(1,2,4,4)/18/,
     $ INDX(1,3,1,1)/0/, INDX(1,3,1,2)/0/, INDX(1,3,1,3)/0/,
     $ INDX(1,3,1,4)/0/,
     $ INDX(1,3,2,1)/0/, INDX(1,3,2,2)/0/, INDX(1,3,2,3)/0/,
     $ INDX(1,3,2,4)/0/
      DATA
     $ INDX(1,3,3,1)/0/, INDX(1,3,3,2)/0/, INDX(1,3,3,3)/0/,
     $ INDX(1,3,3,4)/0/,
     $ INDX(1,3,4,1)/0/, INDX(1,3,4,2)/0/, INDX(1,3,4,3)/0/,
     $ INDX(1,3,4,4)/0/,
     $ INDX(1,4,1,1)/0/, INDX(1,4,1,2)/0/, INDX(1,4,1,3)/0/,
     $ INDX(1,4,1,4)/0/,
     $ INDX(1,4,2,1)/0/, INDX(1,4,2,2)/0/, INDX(1,4,2,3)/0/,
     $ INDX(1,4,2,4)/0/,
     $ INDX(1,4,3,1)/0/, INDX(1,4,3,2)/0/, INDX(1,4,3,3)/8/,
     $ INDX(1,4,3,4)/17/
      DATA
     $ INDX(1,4,4,1)/0/, INDX(1,4,4,2)/0/, INDX(1,4,4,3)/0/,
     $ INDX(1,4,4,4)/19/,
     $ INDX(1,5,1,1)/0/, INDX(1,5,1,2)/0/, INDX(1,5,1,3)/0/,
     $ INDX(1,5,1,4)/0/,
     $ INDX(1,5,2,1)/0/, INDX(1,5,2,2)/0/, INDX(1,5,2,3)/0/,
     $ INDX(1,5,2,4)/0/,
     $ INDX(1,5,3,1)/0/, INDX(1,5,3,2)/0/, INDX(1,5,3,3)/0/,
     $ INDX(1,5,3,4)/0/,
     $ INDX(1,5,4,1)/0/, INDX(1,5,4,2)/0/, INDX(1,5,4,3)/0/,
     $ INDX(1,5,4,4)/0/
      DATA
     $ INDX(1,6,1,1)/0/, INDX(1,6,1,2)/0/, INDX(1,6,1,3)/0/,
     $ INDX(1,6,1,4)/0/,
     $ INDX(1,6,2,1)/0/, INDX(1,6,2,2)/0/, INDX(1,6,2,3)/0/,
     $ INDX(1,6,2,4)/0/,
     $ INDX(1,6,3,1)/0/, INDX(1,6,3,2)/0/, INDX(1,6,3,3)/0/,
     $ INDX(1,6,3,4)/0/,
     $ INDX(1,6,4,1)/0/, INDX(1,6,4,2)/0/, INDX(1,6,4,3)/0/,
     $ INDX(1,6,4,4)/20/,
     $ INDX(2,1,1,1)/0/, INDX(2,1,1,2)/1/, INDX(2,1,1,3)/0/,
     $ INDX(2,1,1,4)/0/
      DATA
     $ INDX(2,1,2,1)/0/, INDX(2,1,2,2)/0/, INDX(2,1,2,3)/4/,
     $ INDX(2,1,2,4)/0/,
     $ INDX(2,1,3,1)/0/, INDX(2,1,3,2)/0/, INDX(2,1,3,3)/0/,
     $ INDX(2,1,3,4)/12/,
     $ INDX(2,1,4,1)/0/, INDX(2,1,4,2)/0/, INDX(2,1,4,3)/0/,
     $ INDX(2,1,4,4)/0/,
     $ INDX(2,2,1,1)/0/, INDX(2,2,1,2)/0/, INDX(2,2,1,3)/3/,
     $ INDX(2,2,1,4)/0/,
     $ INDX(2,2,2,1)/0/, INDX(2,2,2,2)/0/, INDX(2,2,2,3)/0/,
     $ INDX(2,2,2,4)/10/,
     $ INDX(2,2,3,1)/0/, INDX(2,2,3,2)/0/, INDX(2,2,3,3)/0/,
     $ INDX(2,2,3,4)/0/
      DATA
     $ INDX(2,2,4,1)/0/, INDX(2,2,4,2)/0/, INDX(2,2,4,3)/0/,
     $ INDX(2,2,4,4)/0/,
     $ INDX(2,3,1,1)/0/, INDX(2,3,1,2)/0/, INDX(2,3,1,3)/0/,
     $ INDX(2,3,1,4)/9/,
     $ INDX(2,3,2,1)/0/, INDX(2,3,2,2)/0/, INDX(2,3,2,3)/6/,
     $ INDX(2,3,2,4)/0/,
     $ INDX(2,3,3,1)/0/, INDX(2,3,3,2)/0/, INDX(2,3,3,3)/0/,
     $ INDX(2,3,3,4)/13/,
     $ INDX(2,3,4,1)/0/, INDX(2,3,4,2)/0/, INDX(2,3,4,3)/0/,
     $ INDX(2,3,4,4)/0/
      DATA
     $ INDX(2,4,1,1)/0/, INDX(2,4,1,2)/0/, INDX(2,4,1,3)/0/,
     $ INDX(2,4,1,4)/0/,
     $ INDX(2,4,2,1)/0/, INDX(2,4,2,2)/0/, INDX(2,4,2,3)/0/,
     $ INDX(2,4,2,4)/11/,
     $ INDX(2,4,3,1)/0/, INDX(2,4,3,2)/0/, INDX(2,4,3,3)/0/,
     $ INDX(2,4,3,4)/0/,
     $ INDX(2,4,4,1)/0/, INDX(2,4,4,2)/0/, INDX(2,4,4,3)/0/,
     $ INDX(2,4,4,4)/0/,
     $ INDX(2,5,1,1)/0/, INDX(2,5,1,2)/0/, INDX(2,5,1,3)/0/,
     $ INDX(2,5,1,4)/0/
      DATA
     $ INDX(2,5,2,1)/0/, INDX(2,5,2,2)/0/, INDX(2,5,2,3)/0/,
     $ INDX(2,5,2,4)/0/,
     $ INDX(2,5,3,1)/0/, INDX(2,5,3,2)/0/, INDX(2,5,3,3)/0/,
     $ INDX(2,5,3,4)/14/,
     $ INDX(2,5,4,1)/0/, INDX(2,5,4,2)/0/, INDX(2,5,4,3)/0/,
     $ INDX(2,5,4,4)/0/,
     $ INDX(2,6,1,1)/0/, INDX(2,6,1,2)/0/, INDX(2,6,1,3)/0/,
     $ INDX(2,6,1,4)/0/,
     $ INDX(2,6,2,1)/0/, INDX(2,6,2,2)/0/, INDX(2,6,2,3)/0/,
     $ INDX(2,6,2,4)/0/
      DATA
     $ INDX(2,6,3,1)/0/, INDX(2,6,3,2)/0/, INDX(2,6,3,3)/0/,
     $ INDX(2,6,3,4)/0/,
     $ INDX(2,6,4,1)/0/, INDX(2,6,4,2)/0/, INDX(2,6,4,3)/0/,
     $ INDX(2,6,4,4)/0/
C
C    $FOR THE MOMENT*** CALC. ALL INTEGRALS INVOLVING F ORBITALS
C    $SCALED BY 0.666*****
C
C     LMAX=MAX0(L1,L2,L3,L4)
C     IF (LMAX.EQ.3) GO TO 45
C
c     write(6,*)"i am in rrrr"
c     write(6,*)"A1,A2,A3,A4,N1,N2,N3,N4",A1,A2,A3,A4,N1,N2,N3,N4
      KIND = KI
      IZ = IABS(NZ)
c     write(6,*)"NZ",NZ
      IF (NZ.GT.0) GO TO 40
      IF (NZ.EQ.0) GO TO 10
CNR   NZ.LT.0
c     write(6,*)"i am in 1"
      KIND = INDX(KI,K,L1+1,L4+1)
      GO TO 40
CNR
   10 CONTINUE
      FACTOR = ONE
   20 CONTINUE
      F = ZERO
C     GENERAL SLATER CONDON RADIAL INTEGRAL,
      FN = SQRT((2.0D0*A1)** (2*N1+1))
      FN1 = SQRT((2.0D0*A2)** (2*N2+1))
      FN2 = SQRT((2.0D0*A3)** (2*N3+1))
      FN3 = SQRT((2.0D0*A4)** (2*N4+1))
      FN = FN*FN1*FN2*FN3/SQRT(FACT(2*N1+1)*FACT(2*N2+1)*FACT(2*N3+1)*
     .     FACT(2*N4+1))
      N = N1 + N3
      M = N2 + N4
      A = A1 + A3
      B = A2 + A4
      NL = N + K + 1
      ML = N - K
      MN = M + N + 1
      FA = FACT(NL)*FACT(M-K)/ (A**NL*B** (M-K))
      AA = ONE
      BB = (A+B)**MN
      DO 30 I = 1, NL
         AA = AA*A
         BB = BB/ (A+B)
         F = FACT(MN-I)/ (AA*BB)
         FA = FA - FACT(NL)*F/FACT(NL-I+1)
         IF (I.LE.ML) FA = FA + FACT(ML)*F/FACT(ML-I+1)
   30 CONTINUE
C
C     FACTOR =1.0 FOR AB INITIO AND RYDBERG ORBITALS
C     FACTOR =2/3 FOR SEMI EMPIRICAL
C
c     write(6,*)"FN,FA,FACTOR",FN,FA,FACTOR
      F = FN*FA*FACTOR
C
      GO TO 60

   40 CONTINUE
C
C     DO LOOK UP
C     IF ZERO FOUND SET FACTOR TO 2/3 AND CALCULATE
C
C     KKK=KI+8*(NZ-3)
C
C   IF (KIND .GT. 20) TABLE LOOKUP WILL NOT WORK.
C                    BRANCH TO SEMI-EMPIRICAL
C
      IF (KIND.GT.20) GO TO 50
C
      KKK = KIND + 20* (IZ-1)
      F = SPK(KKK)/ENUBAR
      IF (ABS(F).GT.1.0D-8) GO TO 60
   50 CONTINUE
      FACTOR = TWO/THREE
C     IF(KI.GT.20) FACTOR = ONE/THREE
C     WRITE(6,600)
      GO TO 20

   60 RRRR = F
c     write(6,*)"RRRR",RRRR
      RETURN

 9000 FORMAT ('   >>>F AND G NOT IN TABLE-SUBROUTINE RRRR<<<')

      END

