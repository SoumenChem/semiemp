       SUBROUTINE FOCK(BB,HB,F,H,FA,Q,FZ,SVEC,GAMMA,SC,BH,
     x IT,ITMAX,RAD,ISW1,ISWANN,FOP,SUMELEL,NB,NG,QQ,NRR,LABI,LABJ,
     x LABK,LABL,LL,NIN)
C
C     ..................................................................
C
C     ALPHA FOCK MATRIX ELEMENTS CALCULATED AND STORED IN ARRAY BB
C     BETA  FOCK MATRIX ELEMENTS CALCULATED AND STORED IN ARRAY HB
C     BOND ORDER MATRICES ARE FA(ALPHA) AND F(TOTAL)
C     COULOMB MATRIX IS GAMMA IS F(0)
C     Q IS CORE HAMILTONIAN, POSSIBLY MODIFIED FOR SCRF IN FILLF(U)
C     FZ IS USED AS SCRATCH SPACE IN FILLFR
C
C     NP(I)  PRINCIPAL QUANTUM NUMBER OF ORBITAL I
C     D12(I) TOTAL ATOMIC CHARGE DENSITY OF S,P+D ORBITALS
C     Z12(I) NUMBER OF VALENCE ELECTRONS ON ATOM I,OR THE CORE CHARGE
C     M.O. TYPE : 0=S , 1=PX , 2=PY , 3=PZ , 4=D(ZZ) , 5=D(XX-YY) ,
C     6=D(XY) , 7=D(XZ) , 8=D(YZ)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      double precision mult,FOP(*)
c     Integer rtdb
C     WARNING   RPA PGM MUST HAVE THESE DIMENSIONS
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
      dimension  F(*),QQ(*),FA(*),GAMMA(*),RAD(*),FZ(*), SVEC(*),
     .    SC(NB,33),BB(*),HB(*),SUMELEL(*),NG(*),Q(*),LABI(*),
     . LABJ(*),LABK(*),LABL(*),NIN(*)
C
      INS = ISWANN + 1

c     write(6,*)"i am in fock"
c     write(6,*)"Q(I)",(Q(I),I=1,20)

C
C     Extended Huckel section
C
c     IF (ISW1.NE.0) GO TO 10
c     CALL FILLEX(GAMMA,BB,Q,F,FA,SC,IT,ITMAX,DRYRUN)
c     GO TO 50

   10 CONTINUE
      DO 15 I = 1, NB
       SUMELEL(I) = ZERO
   15 CONTINUE
C
        GO TO (20,30,30,40), INS
C
C     Restricted Hartree-Fock for closed shell
C
   20 CONTINUE
c     write(6,*)"before fillf BB",(BB(I),I=1,10)
c     write(6,*)"before fillf HB",(HB(I),I=1,10)
c     write(6,*)"before fillf F",(F(I),I=1,10)
c     write(6,*)"before fillf FA",(FA(I),I=1,10)
       CALL FILLF(BB,FA,F,Q,GAMMA,SUMELEL,NG,QQ,NRR,LABI,LABJ
     1 ,LABK,LABL,LL,NIN,1)
C     WRITE(6,*) ' ONE ELECTRON MATRIX'
C     CALL PRNT(Q,NB,NRR)
C
c     write(6,*)"after fillf BB",(BB(I),I=1,nrr)
c     write(6,*)"after fillf HB",(HB(I),I=1,nrr)
c     write(6,*)"after fillf F",(F(I),I=1,nrr)
c     write(6,*)"after fillf FA",(FA(I),I=1,nrr)

c     IF(ISW1.EQ.8) THEN
C     WRITE(6,*) ' IN FOCK FOCK MATRIX'
C     CALL PRNT(BB,NB,NRR)
c      WRITE(6,*) ' IN FOCK SUMELEL'
c      WRITE(6,*) (SUMELEL(I), I = 1,NB)
C     WRITE(6,*) (SUMNUCEL(I), I = 1,NB)
c     ENDIF
      GO TO 50
C
C     Unrestricted Hartree-Fock (UHF) calculation
C
   30 CONTINUE
      CALL FILLFU(BB,HB,FA,F,Q,GAMMA,SUMELEL,NG,QQ,NRR,LABI,
     1 LABJ,LABK,LABL,LL,NIN,1)
C     IF XINDO , CONSIDER THIS!!
c     IF(ISW1.EQ.8) THEN
c      WRITE(6,*) ' IN FOCK SUMELEL'
c      WRITE(6,*) (SUMELEL(I), I = 1,NB)
C      WRITE(6,*) (SUMNUCEL(I), I = 1,NB)
C      WRITE(6,*) ' FOCK '
C      WRITE(6,*) ' FOCK MATRIX'
C      CALL PRNT(BB,NB,NRR)
C      WRITE(6,*) ' FOCK MATRIX'
C      CALL PRNT(HB,NB,NRR)
c     ENDIF
      GO TO 50
C
C     Restricted open-shell Hartree-Fock (ROHF)
C
   40 CONTINUE
c     CALL FILLF(BB,HB,FA,F,Q,RAD,GAMMA,DRYRUN)
C
C     BB IS NOW H+P
C
c     CALL FILLFR(BB,HB,FA,F,FZ,RAD,GAMMA,DRYRUN)
   50 CONTINUE

C     THIS IS A COSMO SECTION
c     IF(JSW(6) .EQ. 30) THEN

c      NFM1 = 10 * NA
c      NFM2 = 34
c      I_FM = MEM_GET(NFM1*NFM2,      'real',    'GSN1-FM')

c      IF(.NOT.DRYRUN) THEN
c        CALL ZREWIND(IO_VEC_DEN, 'FOCK')
c        IF(ISWANN .EQ. 0) THEN
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZRREAD(IO_VEC_DEN, F, 1, NRR, 'FOCK')
c        ELSEIF(ISWANN .LT. 3) THEN
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZRREAD(IO_VEC_DEN, F, 1, NRR, 'FOCK')
c        ELSEIF(ISWANN.EQ. 3) THEN
c                   ROHF
c          DO  I = 1, NDN
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c          ENDDO
c          CALL ZRREAD (IO_VEC_DEN, F , 1, NRR, 'FOCK')
C                   READ (IO_VEC_DEN)
C                   READ (IO_VEC_DEN) (BB(I),I=1,NRR)
c          CALL ZIREAD(IO_VEC_DEN, IDU, 1, 1, 'FOCK')
c        ELSE
c          WRITE(6,2000)
c          CALL ZEXIT(3)
c        ENDIF
c        CALL REFMOM(BIGMEM(I_FM),FZ,NFM1,NFM2)
C        THIS FILE WITH REFMOM
c        WRITE(6,*) ' FOCK NFM1, NFM2', NFM1,NFM2
c        CALL ATMOM(BIGMEM(I_FM),F,FZ,NFM1,NFM2)
c        CALL ADDFCK(BB,F,AMAT,BH,BIGMEM(I_FM),NFM1,NFM2)
c      ENDIF
c      CALL MEM_FREE(NFM1*NFM2,     'real',    'GSN1-FM')
c     ENDIF
C

 2000 FORMAT(' **** ERROR IN FOCK, WRONG SCFTYP IN $CONTRL ****')
 5000 RETURN
      END
             

