      SUBROUTINE SING(HA,HB,F,H,FA,FZ,SVEC,GAMMA,SC,RAD,E1,BH,
     X IT,ITGL,ITMAX,ISW1,ISWANN,NDN,NB,NA,ISW3,NTYP,NU,W,JSLIT,FOP
     X ,SUMELEL,NG,NRR,Q,REP,D12,NIN,PNOW,RA,
     X RB,PRA,PRB,XBACK1,XBACK2,PBACK1,slowscf,LL)
C
C     ..................................................................
C
C     THIS SUBROUTINE CALCULATES THE ENERGY E1
C     BB IS THE FOCK MATRIX WHICH MUST BE CALCULATED BY SUBROUTINE FOCK
C     F IS THE BOND ORDER MATRIX
C     H IS THE ORIGINAL HAMILTONIAN MATRIX
C     GAMMA IS THE COULOMB MATRIX
C     FZ IS USED FOR TEMPORARY STORAGE IN SUBROUTINE FOCK
C
C     ..................................................................
C

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I,N)
#include "global.fh"
c     Integer rtdb
C     WARNING   RPA PGM MUST HAVE THESE DIMENSIONS
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10,IDD1=0, IDD2 = 0)
      LOGICAL      DIIS, DIISCHK, DDFLAG, DRYRUN
      COMMON       /DIISX/DIIS, DIISCHK, DDFLAG
C
      DIMENSION    HA(*),HB(*),BH(*),F(*),NTYP(*),NU(*),NG(*)
      DIMENSION    H(*),FA(*),FZ(*),GAMMA(*),RAD(*),DIFF(10)
     . ,W(*),FOP(*),SUMELEL(*),RA(*),RB(*),
     . D12(NA), Q(*),NIN(*),PNOW(*),PRA(*),PRB(*),
     . XBACK1(*), XBACK2(*), PBACK1(*)
      INTEGER JSLIT(4)
      logical slowscf
      SAVE         ARF1, BRF1, E11
      DATA         ARF1/0.0D0/, BRF1/0.0D0/

c     write(6,*)"i am in sing"
c     write(6,*)"HB",(HB(I),I=1,nrr)
c     write(6,*)"JSLIT in sing",JSLIT
      IRFAC=0
c     if (ga_nodeid().eq.0) then
c     IF(DIIS .AND. IDD2.EQ.1) THEN
C
C     IF DIIS IS ON THEN THE DENSITY EXTRAPOLATION IS NOT
C     DONE AND ALL STUFF ON UNIT 14 IS PUT ON ITAPE.
C
C     REWIND ITAPE
c      CALL ZREWIND (ITAPE, 'SING')
C     REWIND 14
c      CALL ZREWIND (IO_VEC_DEN, 'SING')
c      DO 10 ID = 1, NDN
C     READ(14)
c       CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'SING')
C     READ(14) (F(I), I=1,NRR)
c       CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'SING')
C     WRITE(ITAPE) (F(I),I=1,NRR)
c       CALL ZRWRITE (ITAPE, F, 1, NRR, 'SING')
C     SCFCHK IS CALLED HERE ONLY TO UPDATE THE DENSITY DIFFERENCE.
c       CALL SCFCHK(F,FA,(IT-1),ITGL,ID,.TRUE.)
c  10  CONTINUE
C
c     ELSE
C
C       DO THE DENSITY EXTRAPOLATION (DYNAMIC DAMPING) IF NOT DIIS.
C       PICK UP THE CURRENT DENSITIES, CALL SCFCHK AND THEN FORM
C       EXTRAPOLATED TOTAL DENSITY.
c      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
c      if (IT.eq.1) then
c      DO  I = 1, NRR
c       write(6,*)"RA(I),RB(I)",RA(I),RB(I)
c       F(I) = RA(I) 
c       FA(I) = RB(I)
c      END DO
c      else
c      DO  I = 1, NRR
c       write(6,*)"RA(I),RB(I)",RA(I),RB(I)
c       F(I) = RB(I)
c       FA(I) = RA(I)
c      END DO
c      end if
c      END IF
C
        DIIS = .FALSE.
        DO 20 ID = 1, NDN
        IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
C       REWIND 14
c       CALL ZREWIND (IO_VEC_DEN, 'SING')
C        IF(ISAVE.EQ.2.AND.IT.EQ.1) THEN
C          CALL DECALL(BB,F,ID)
C        ELSE
C          READ(14)
c          CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'SING')
C          READ(14) (F(I),I=1,NRR)
c          CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'SING')
c        ENDIF
          if (ID.eq.1) then
c        write(6,*)"RB before SCFCHK",(RB(I),I=1,15)
c        write(6,*)"FA before SCFCHK",(FA(I),I=1,15)
         CALL SCFCHK(RA,FA,(IT-1),ITGL,ID,.FALSE.,DE1,NDN,NB,NA,ISW1,
     1   ISW3,NTYP,NU,JSLIT,NRR,PNOW,ISWANN,PRA,PRB,RA,RB,XBACK1,XBACK2,
     2   PBACK1,slowscf)
c        write(6,*)"RB after SCFCHK",(RB(I),I=1,15)
c        write(6,*)"FA after SCFCHK",(FA(I),I=1,15)
          else
c        write(6,*)"RA before SCFCHK",(RA(I),I=1,21)
c        write(6,*)"FA before SCFCHK",(RB(I),I=1,21)
         CALL SCFCHK(RB,RA,(IT-1),ITGL,ID,.FALSE.,DE1,NDN,NB,NA,ISW1,
     1   ISW3,NTYP,NU,JSLIT,NRR,PNOW,ISWANN,PRA,PRB,RA,RB,XBACK1,XBACK2,
     2   PBACK1,slowscf)
c        write(6,*)"RA after SCFCHK",(RA(I),I=1,21)
c        write(6,*)"FA after SCFCHK",(FA(I),I=1,21)
          ENDIF
         else
c        write(6,*)"FA after SCFCHK",(FA(I),I=1,21)
c        write(6,*)"F  after SCFCHK",(F(I),I=1,21)
        CALL SCFCHK(F,FA,(IT-1),ITGL,ID,.FALSE.,DE1,NDN,NB,NA,ISW1,
     1  ISW3,NTYP,NU,JSLIT,NRR,PNOW,ISWANN,PRA,PRB,RA,RB,XBACK1,XBACK2,
     2   PBACK1,slowscf)
        ENDIF
   20   CONTINUE
C
C     SET UP FACTOR TEST - FACTOR = 0 TWICE
C     IN A ROW BEFORE DIIS IS USED.
C
       IF(IDD2.EQ.1) THEN
        RFAC = ZERO
        DO 30 IR = 1, NDN
         RFAC = RFAC + FAC(IR)
   30   CONTINUE
        ARF1 = BRF1
        BRF1 = RFAC/NDN
        RFAC = ARF1 + BRF1
c     IF(RFAC.LE.0.1 .AND. DIISCHK) THEN
        IF(RFAC.LE.0.1) THEN
         IRFAC = IRFAC + 1
        ELSE
         IRFAC = 0
        END IF
       END IF
c      END IF
C     WRITE(6,67) IRFAC, RFAC
C     67    FORMAT(5X,'IRFAC = ',I4,'  RFAC = ',F11.6)
C
c     END IF
C
C     FORM TOTAL DENSITY.
C
C     370  REWIND ITAPE


   40 CONTINUE
c     CALL ZREWIND (ITAPE, 'SING')
c     IF(ISAVE.EQ.2.AND.IT.EQ.0) GO TO 70
C     READ (ITAPE) (F(I),I=1,NRR)
c     CALL ZRREAD (ITAPE, F, 1, NRR, 'SING')
      IF(NDN.EQ.1) GO TO 70
c     DO 60 J = 2, NDN
C     READ(ITAPE) (FA(I),I=1,NRR)
c      CALL ZRREAD (ITAPE, FA, 1, NRR, 'SING')
c      if (IT.eq.1) then 
       DO 50 I = 1, NRR
c       write(6,*)"RA(I),RB(I)",RA(I),RB(I)
        F(I) = RA(I) + RB(I)
        FA(I) = F(I) - RA(I)
   50  CONTINUE
c      else
c      DO  I = 1, NRR
c       write(6,*)"RA(I),RB(I)",RA(I),RB(I)
c       F(I) = RA(I) + RB(I)
c       FA(I) = F(I) - RA(I)
c      END DO
c      end if
   60 CONTINUE
   70 CONTINUE
C     WRITE(ITAPE) (F(I),I=1,NRR)
c     CALL ZRWRITE (ITAPE, F, 1, NRR, 'SING')
C
C     CALCULATE THE ATOMIC DENSITIES, FOR SCONS.
C
      L = 0
      DO 80 I = 1, NA
       D12(I) = ZERO
   80 CONTINUE
c     write(6,*)"F(I)",(F(I),I=1,nrr)
      DO 90 I = 1, NB
       K = NU(I)
       L = L + I
c      write(6,*)"D12,F",D12(K),F(L)
       D12(K) = D12(K) + F(L)
c     write(6,*)"D12",D12(K)
   90 CONTINUE
c     write(6,*)"D12(K)",(D12(K),K=1,5)
c     IF(JSW(6).EQ.50) THEN
CJiabo
C
C     calculate bond order for CM2 model
C
c       LPS=NB*(NB+1)/2
c       NAR = (NA*NA+NA)/2
c       DO 15 I = 1, NAR
c       BOND(I) = ZERO
c  15   CONTINUE
c       DO 35 I = 1, NB
c        NI = NU(I)
c        DO 25 J = 1, I
c           NJ = NU(J)
c           IJ = NIN(I) + J
c           NINJ = NIN(NI) + NJ
c           BOND(NINJ) = BOND(NINJ) + F(IJ)**2
c           W(LPS+I+(J-1)*NB)=F(IJ)
c           W(LPS+J+(I-1)*NB)=F(IJ)
c           DENG(IJ)=F(IJ)
c  25    CONTINUE
c  35   CONTINUE
C       IF(JSW(6).LT.0) CALL AMSOLE(F,NB,ES)
C       WRITE(*,*)"ES=",ES
c     ENDIF
CLi
C
C
C     CALL FOCK TO FORM THE FOCK MATRICES BB=HA AND HB
C
c     write(6,*)"before fock HA",(HA(I),I=1,nrr)
c     write(6,*)"before fock HB",(HB(I),I=1,nrr)
c     write(6,*)"before fock F",(F(I),I=1,nrr)
c     write(6,*)"before fock FA",(FA(I),I=1,nrr)
c     write(6,*)"before fock H",(H(I),I=1,nrr)
      CALL FOCK(HA,HB,F,H,FA,Q,FZ,SVEC,GAMMA,SC,BH,
     x IT,ITMAX,RAD,ISW1,ISWANN,FOP,SUMELEL,NB,NG,NRR,NIN,LL)
c     if (ga_nodeid().eq.0) write(6,*)"after fock HA",(HA(I),I=1,nrr)
c     if (ga_nodeid().eq.0) write(6,*)"after fock HB",(HB(I),I=1,nrr)
c     if (ga_nodeid().eq.0) write(6,*)"after fock F",(F(I),I=1,nrr)
c     if (ga_nodeid().eq.0) write(6,*)"after fock FA",(FA(I),I=1,nrr)
      DO I = 1,NRR
        H(I) = Q(I)
      END DO
c     end if
C
C
C     March 17, 1998
C     Jiabo Li
C
c     LPS=NB*(NB+1)/2
c     DO I=1,NB*(NB+1)/2
c     W(LPS+NB*NB+I)=0.0D0
c     END DO
C
C     CALCULATE ENERGY OF GASPHASE HAMILTONIAN USING SOLVATED DENSITY
C     E1-SOLATIN CONTRIBUTION

c     CALL GMSCRF(NB,NU,KAN,W(LPS+NB*NB+1))
c     NADD=LPS+NB*NB
c     DO I=1,NB*(NB+1)/2
c     H(I)=H(I)-W(NADD+I)
c     HA(I)=HA(I)-W(NADD+I)
c     HB(I)=HB(I)-W(NADD+I)
c     END DO

C
CLi
C     WRITE(6,*) ' IN SING FOCK'
C     CALL PRNT (HA,NB,NRR)
C
      ENERGY = ZERO
      IJ = 0
C     ELECTRONIC ENERGY, OFF DIAGONAL CONTRIBUTION.
c     IF(ISW1.EQ.8.OR.ISW1.EQ.9) REWIND (IO_XINDO)
c     if (ga_nodeid().eq.0) then
      IF(ISWANN.NE.0) GO TO 120
C
C     CLOSED SHELL RHF
C
c     IF(ISW1.EQ.8) THEN
C       XINDO REQUIRES THE NON-ORTHOGONAL DENSITY FOR ENERGY
C
C
C       REWIND 16
c       CALL ZREWIND (IO_SCF_MATRICES, 'SING')
C       READ (16)
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'SING')
C       READ (16)
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'SING')
c         CALL PRNT( HA, NB, NRR)
c       ENDIF
c       IF(JSLIT(1).EQ.2) WRITE(IO_XINDO) (F(I) , I = 1,NRR)
c     ENDIF
C     WRITE(6,*) ' IN SING, DENSITY MATRICES'
C     CALL PRNT(F,NB,NRR)
c     write(6,*)"F(1),H(1),HA(1)",F(1),H(1),HA(1)
      ENERGY = HALF*F(1)* (H(1)+HA(1))
      IJ = 2
      DO 110 I = 2, NB
       II = I - 1
       DO 100 J = 1, II
         ENERGY = ENERGY + F(IJ)* (H(IJ)+HA(IJ))
        IJ = IJ + 1
  100  CONTINUE
       ENERGY = ENERGY + HALF*F(IJ)* (H(IJ)+HA(IJ))
       IJ = IJ + 1
  110 CONTINUE
c     write(6,*)"Diagonal energy",ENERGY
C
C     ELECTRONIC ENERGY, DIAGONAL CONTRIBUTION.
C
C     II=0
C     DO 12 I=1,NB
C     II=II+I
C     ENERGY=ENERGY+HALF *F(II)*(HA(II)+H(II))
C     12 CONTINUE
      call ga_sync
      GO TO 210

  120 CONTINUE
      IF(ISWANN.EQ.3) GO TO 160
c     write(6,*)"F total",(F(I),I=1,nrr)
c     write(6,*)"FA",(FA(I),I=1,nrr)
c     write(6,*)"H(I)",(H(I),I=1,nrr)
c     write(6,*)"HA(I)",(HA(I),I=1,nrr)
c     write(6,*)"HB(I)",(HB(I),I=1,nrr)
C
C     UHF,F=TOTAL DENSITY, FA=ALPHA DENSITY
C
c     IF(ISW1.EQ.8) THEN
C     XINDO REQUIRES THE NON-ORTHOGONAL DENSITY FOR ENERGY
C
C     REWIND 16
c      CALL ZREWIND (IO_SCF_MATRICES, 'SING')
C     READ (16)
c      CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'SING')
C     READ (16)
c      CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'SING')
C
C     GET FOCK TRANSFORMATION MATRIX
C
C     READ (16) (FZ(I),I=1,NRR)
c      CALL ZRREAD (IO_SCF_MATRICES, FZ, 1, NRR, 'SING')
c      CALL FOCKTR(F,FZ,NB,DRYRUN)
c      CALL FOCKTR(FA,FZ,NB,DRYRUN)
c      IF(JSLIT(1).EQ.2)  THEN
c       WRITE(IO_XINDO) (F(I) , I = 1,NRR)
c       WRITE(IO_XINDO) (FA(I), I = 1,NRR)
c      ENDIF
c     ENDIF
C     ELECTRONIC ENERGY, OFF DIAGONAL CONTRIBUTION
      DO 140 I = 2, NB
       IJ = IJ + 1
       II = I - 1
       DO 130 J = 1, II
        IJ = IJ + 1
        ENERGY = ENERGY + F(IJ)*H(IJ) + FA(IJ)*HA(IJ) +
     .     (F(IJ)-FA(IJ))*HB(IJ)
  130  CONTINUE
  140 CONTINUE
C     ELECTRONIC ENERGY , DIAGONAL CONTRIBUTION
      II = 0
C     WRITE(6,*) ' IN SING, DENSITY MATRICES'
C     CALL PRNT(F,NB,NRR)
C     CALL PRNT(FA,NB,NRR)

      DO 150 I = 1, NB
       II = II + I
       ENERGY = ENERGY + HALF* (F(II)*H(II)+FA(II)*HA(II)+
     .    (F(II)-FA(II))*HB(II))
  150 CONTINUE
C
      call ga_sync()
      GO TO 210
C

C     OPEN SHELL RHF
C
C     F=TOTAL DENSITY   FA=OPEN SHELL DENSITY   HA=P MATRIX  HB=Q MATRIX
  160 CONTINUE
C     REWIND ITAPE
c     CALL ZREWIND(ITAPE, 'SING')
C     REWIND 16
c     CALL ZREWIND(IO_SCF_MATRICES, 'SING')
c     READ (16)  OVERLAP
c     CALL ZIREAD(IO_SCF_MATRICES, IDU, 1, 1, 'SING')
C     READ (16) WEIGHTED OVERLAP
c     CALL ZIREAD(IO_SCF_MATRICES, IDU, 1, 1, 'SING')
C     READ (16)  S**(-1/2)
c     CALL ZRREAD(IO_SCF_MATRICES, FZ, 1, NRR, 'SING')
C     READ (16)  SKIP ONE ELECTRON MATRIX
c     CALL ZIREAD(IO_SCF_MATRICES, IDU, 1, 1, 'SING')
c     DO 200 K = 1, NDN
c      IJ = 0
C     READ (16) (HA(I),I=1,NRR)
c      CALL ZRREAD(IO_SCF_MATRICES,HA,1,NRR,'SING')
c     READ (ITAPE) (FA(I),I=1,NRR)
c      CALL ZRREAD(ITAPE,FA,1,NRR,'SING')
c      IF(ISW1 .EQ. 8) THEN
C     XINDO REQUIRES THE NON-ORTHOGONAL DENSITY FOR ENERGY
C
C     WRITE(6,*) ' SING F OP = ', K
C     CALL PRNT(HA,NB,NRR)
c       CALL FOCKTR(FA,FZ,NB,DRYRUN)
c       IF(JSLIT(1) .EQ. 2) WRITE(IO_XINDO) (FA(I), I = 1,NRR)
c      ENDIF
c      DO 180 I = 2, NB
c       IJ = IJ + 1
c       II = I - 1
c       DO 170 J = 1, II
c        IJ = IJ + 1
c        ENERGY = ENERGY + FA(IJ) * (H(IJ) + HA(IJ))
c 170   CONTINUE
c 180  CONTINUE
c      II = 0
c      DO 190 I = 1, NB
c       II = II + I
c       ENERGY = ENERGY + FA(II) * (H(II) + HA(II)) / TWO
c 190  CONTINUE
c 200 CONTINUE
C     READ(ITAPE) (F(I),I=1,NRR)
c     CALL ZRREAD(ITAPE, F, 1, NRR, 'SING')
  210 CONTINUE
C
C     NOW THE NUCLEAR REPULSION. PT. CHARGE PT. CHARGE TERM IS NOT INCLU
C     DED
 9000 FORMAT (//' CORE REPULSION ENERGY =',F13.6,' au',//)
C
C     Here the TOTAL energy E1 is calculated. When adding new terms to
C     it, please use new variables in ennuc.cmn. Do NOT add your energy
C     to any of the existing variables (REP in particular).
c
c     write(6,*)"ENERGY1,REP",ENERGY,REP
       E1 = ENERGY + REP ! + EPTCHG
c     write(6,*)"E1",E1
C     Additions due to reaction field:
C
      E11 = ZERO
c     end if
CTT   I ASSUME THAT THE 20 WAS FOR DRF - NOW IT IS 40.
c     IF((JSW(6) .NE. 0) .AND. (JSW(6) .NE. 4) .AND.
c    x   (JSW(6) .NE. 40) . AND.(JSW(6).NE.50) .AND.
c    x   (JSW(6).NE.30)) THEN
c        CHRGE = ZERO
c        DO 220 I = 1, NA
c         CHRGE = CHRGE + Z12(I)
c 220    CONTINUE
c        CHRGE = CHRGE - FLOAT(NEL)
c        E22   = ENERGY + REP
c        CALL ESCRF(E22,E11,Z12,CO,CHRGE,IT,NA,IXS)

C       WRITE(6,*) 'sing: energy,rep', ENERGY,REP
C       WRITE(6,*) 'sing: e11,e22', E11,E22
C       E11 RETURNS AS THE TOTAL SOLVATION ENERGY IN THE SCF.
C       FOR GEOMETRY OPTIMISATION THIS IS THE ENERGY E1 THAT
C       IS VARIED.
c       E_SOLV = E11
c       E1     = E1 + E11
c     ELSE IF(JSW(6) .EQ. 4) THEN

C       PCM model: add 1/2 of M-S energy, similar to SCRF model A.
C       See also the related printouts in PTCHG and ESCRF

c       E1 =     E1-E_MS/TWO
c     ELSE IF(JSW(6) .EQ. 40) THEN
c       DO PIETS THING
c       IF(USEUNIT0) write(6,*) 'zindrf5 in'
c       EDRFN=ZERO
c       CALL ZINDRF5(EDRFN)
c       E1 = E1 + EDRFN
c       IF(USEUNIT0) write(6,*) 'zindrf5 out'
c     END IF
C     Extended Huckel:
c     IF(ISW1 .EQ. 0) E1 = E1 - REP + E11
C
C     DDFLAG=TRUE => USE DIIS ON DENSITIES.
C
C     IF(DDFLAG) GO TO 360 THIS MUST BE CHECKED.
C
C     NOW THAT WE HAVE THE ENERGY, CHECK FOR CREEP
C     AND TURN ON DIIS IF ALRIGHT.
C
c     IF(IDD2 .EQ. 1) THEN
c      IF(IT .GT. 10) THEN
c       DO 230 I = 2, 10
c        DIFF(I-1) = DIFF(I)
c 230   CONTINUE
c       DIFF(10) = DE1
c      ELSE
c       DIFF(IT) = DE1
c      END IF
C     WRITE(6,66) ABS(DIFF(IA)/DIFF(IA-1)),
C     .              ABS(DIFF(IA-1)/DIFF(IA-2))
C     66    FORMAT(5X,'D/D-1 = ',F11.6,'   D-1/D-2 =',F11.6)
C
c      IA    = 10
c      IDIIS = - 7
c      CALL GET_INT('CONTRL','DIIS',IDIIS,1)
c      RDIFF = 10.0D0**(IDIIS)
c      IF(DIIS) THEN
c       CONTINUE
c      ELSE IF(ABS(DIFF(IA-1)) .LT. 1.0D-15) THEN
c       CONTINUE
c      ELSE IF(ABS(DIFF(IA-2)) .LT. 1.0D-15) THEN
c       CONTINUE
c      ELSE IF(ABS(DIFF(IA)/DIFF(IA-1)) .LT. HALF .OR.
c    .     ABS(DIFF(IA)/DIFF(IA-1)) .GT. 1.08D00 .OR.
c    .     ABS(DIFF(IA-1)/DIFF(IA-2)) .LT. HALF .OR.
c    .     ABS(DIFF(IA-1)/DIFF(IA-2)) .GT. 1.08D00 .OR.
c    .     DE1 .GE. RDIFF) THEN
C
C     DE1 ABOVE IS DIFF, CURRENTLY IF DIFF < 1.0D-07 THEN
C     DIIS WILL TURN ON.
C
c       CONTINUE
c      ELSE
c       IF(IRFAC .GE. 2) DIIS = .TRUE.
c      END IF
C
C
c      IF(((IDD2 .EQ. 1) .AND. (ABS(DE1) .LE. RDIFF)) .AND.
c    .    ((.NOT. DIIS) .AND. (IT .GT. 4))) THEN
CCC   write(6,*)' ## Turning on DIIS now!'
CCC   write(6,*)' ## DE1=',de1,' RDIFF=',rdiff
c       DIIS = .TRUE.
c      END IF
c      IF(DIIS) THEN
C     RDB     IF(DIIS .AND. IRFAC .GE. 2) THEN
C     ***** USE DIIS *****
C
c       ISTCHK = ISTCHK + 1
C     WRITE(6,*) ' SING IDD2 DIIS IRFAC',IDD2,DIIS,IRFAC
C
c       CALL UDIIS(IT,IERR,F,FA,H,HB)
C
C     PICK UP TOTAL DENSITY TO FORM D13 IN SCFIT.
C
C     REWIND 14
C     READ(14)
C     READ(14) (F(I),I=1,NRR)
C     READ(14) (F(I),I=1,NRR)

c       CALL ZREWIND(IO_VEC_DEN, 'SING')
c       DO 240 IL=1,NDN
c        CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'SING')
c        CALL ZRREAD(IO_VEC_DEN, F, 1, NRR, 'SING')
c 240   CONTINUE
c       IF(ISWANN .NE. IZERO) THEN
c        CALL ZRREAD(IO_VEC_DEN, F, 1, NRR, 'SING')
c       END IF
c      END IF
c     END IF
 5000 RETURN
      END

