      SUBROUTINE DIAGU(rtdb,BB,F,H,FA,FZ,RAD,ITMAX,IT,PREVS,SCFTOL,MIM,
     1  ISW3,NB,EIG,FOP,ISWANN,NEL,MULT,W,JSLIT,NRR,Q,NIN,NU,NTYP,NDN,RA
     2  ,RB,NCORE,NASS)
C
C     ..................................................................
C
C     THIS ROUTINE DIAG. THE FOCK MATRIX AND CALCULATES THE BOND ORDER
C     MATRIX, WHICH IS THEN CHECKED BY SUBROUTINE SCFCHK FOR OSCILLATION
C     .... UNRESTRICTED HARTREE FOCK .....  UHF ....
C
C     ..................................................................
C
C     NR= DIMENSION OF THE SQUARE MATRIX THAT EIGENVALUES ARE RETURNED
C     NRR= DIMENSION OF THE TRIANGULAR MATRIX THAT IS TO BE DIAGONALIZED
C     NB = DIMENSION OF THE CORRESPONDING SQUARE MATRICES
C     NEL = THE NUMBER OF VALENCE ELECTRONS
C     MULT=SPIN MULTIPLICITY=2S+1
C     W AND WA - OCCUP OF THE BETA AND ALPHA ELECTRON MO'S
C     MATRIX IS DIAGONALIZED BY THE SUBROUTINE GIVENS
C
C     ..................................................................
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
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
#include "inp.fh"
      integer rtdb
      double precision MULT,FOP(*)
      LOGICAL      IBGIV
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
      dimension  F(*),H(*),FA(*),RAD(*),W(*),Q(*),D(NB),E(NB),E2(NB)
     .           ,NASS(7),MIM(*),BB(*),EIG(*),NIN(*),A(NRR),EIGA(NB),
     .            WA(NB),NU(*),NTYP(*),RA(*),RB(*)
      DIMENSION    IIMO(NB),JSLIT(4)
      SAVE         ISWAN1, IL, NXT

      integer g_A,g_F,g_alpha
C
C****************************************************************
C
CGP      DATA         IZERO, IONE, ITWO, ITHREE/0, 1, 2, 3/
CGP      DATA         IFOUR, IFIVE, IT12/4, 5, 12/
      DATA IT12/12/
C
C
c     write(6,*)"i am in diagu",MIM
         
      IF (IT.NE.0) GO TO 20
C
C     FIRST ITERATION ONLY!
C     
      IPUN = 0
      ISWAN1 = (ISWANN+1)/2 + 1
      IL = 0
      NXT = 0
      IF (ISW3.EQ.0) GO TO 10
      GO TO 220

   10 CONTINUE
      IPUN1 = IPUN + 1
c     write(6,*)"IPUN",IPUN
      GO TO (20,20,220,220,20,70,70,70,70,
     .       20,70), IPUN1
C
C     FOR ALL CYCLES (JUST PLAIN START)
C
   20 CONTINUE
      IL = IL + 1
C
C     CLEAN UP MATRIX TO AVOID PROBLEMS IN DIAGONALIZATION
C
C   REMOVED BY JDHD JAN 17 84
C   IF ZEROING F(I) USE THRESH OF ORDER 1.D-11
C      DO 60 I=1,NRR
C        IF(ABS (F(I)).LE.1.0D-8)  F(I)=ZERO
C   60 CONTINUE
      NROOT = MIM(2) + 1
      IF (NASS(1).LT.0) NROOT = NROOT + 22
      IF (NASS(1).GT.0) NROOT = NROOT + NASS(1)
      IF (NROOT.GT.NB) NROOT = NB
      CALL SLITET(ITWO,J,JSLIT)
C     IF CALL FOR NO ITERATIONS QUIT
      IF (ITMAX.EQ.0) GO TO 30
      IF (J.EQ.2) GO TO 40
C
C     LAST TIME THROUGH--CONVERGED!! (MAYBE)
C
   30 CONTINUE
      NXT = 1
      NROOT = NB
      CALL SLITE(ITWO,JSLIT)
   40 CONTINUE
C
C     BETA UHF
C
CVAX      STATPR=SYS$SETPRN('GIVENS')
      IBGIV = .TRUE.
C     IS LEVEL SHIFTING ADVANTAGEOUS?
C     IF(IDD2.NE.0.AND.(EIG(NEL/2+1)-EIG(NEL/2)).LT.0.05) THEN
C       IF(IT.NE.0) IDD2 = 2
C     ENDIF

c     IF (IDD2.EQ.20 .AND. IT.GE.10) IDD2 = 0
c     IF (IDD2.EQ.21 .AND. IT.GE.10) IDD2 = 1
c     IF ((IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1 .AND. IT.NE.0) THEN
c        FLAM = TWO/ (FLOAT(IT))
c        CALL LEVELSHF(BB(1),F,BB(NRR+1),FLAM,NB,IONE)
C     REWIND 14
c        CALL ZREWIND (IO_VEC_DEN, 'DIAGU')
C     END IF
      IF(ISW1.EQ.8.OR.ISW1.EQ.9) THEN
C
C        REWIND 16
c        CALL ZREWIND (IO_SCF_MATRICES, 'DIAGU')
C        READ (16) OVERLAP
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGU')
C        READ (16) WEIGHTED OVERLAP
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGU')
C
C        GET S**(-1/2)                   
C
C        READ (16) (H(I),I=1,NRR)
c        CALL ZRREAD (IO_SCF_MATRICES, H, 1, NRR, 'DIAGU')
c        CALL FOCKTR(F,H,NB,DRYRUN)
C        READ (16)  SKIP ONE ELECTRON MATRIX
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGU')
C        READ (16)  SKIP BETA FOCK MATRIX F
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGU')
      ENDIF

c     WRITE(6,*) ' IN DIAGU F'
C     CALL PRNT(F,NB,NRR)
c     write(6,*)"FA before givens",(FA(I),I=1,nrr)
c     CALL GIVENS(F,BB,EIG,NB,NROOT,IBGIV)
c     write(6,*)"F after givens", (F(I),I=1,nrr)
c     write(6,*)"BB after givens", (BB(I),I=1,nrr)
c     write(6,*)"EIG after givens", (EIG(I),I=1,nb)
c     WRITE(6,*) ' DIAGR BEFORE GIVENS H ', (H(I),I=1,40)
c     call mat_sym(F,nb)
c     CALL GIVENS(F,BB,EIG,NB,NROOT,IBGIV,QQ,NIN)

       call diag_para(F,NRR,NB,BB,EIG) 

c      call ga_sync()
c     IF ((IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1 .AND. IT.NE.0) THEN
c        DO 50 I = MIM(2) + 1, NROOT
c           EIG(I) = EIG(I) - FLAM
c  50    CONTINUE
c     END IF
CVAX      STATPR=SYS$SETPRN('EXDIAG')
C
C     CLEAN UP VECTORS
C  DO NOT ZERO VECTORS IF RE-INCLUDING ZEROING
C  TENDS TO CREATE MORE PROBLEMS THAN ITS WORTH
C      DO 90 I=1,NR
C   90 IF(ABS(BB(I)).LE.1.0D-8) BB(I)=ZERO
C

C
      CALL ASSIGN(BB,IZERO,NROOT,IONE,IT,NASS,MIM,FOP,ISWANN,ISW3,IPUN
     1 ,NB,EIG,NEL,MULT,W,NCORE)

c     IF (ISW1.EQ.0 .OR. ISW1.EQ.9) CALL EXOUT(RAD,W,BB,H,FA,F,EIG,NB,
c    .    NXT,DRYRUN)
c     IF (SCFTOL.LT.ZERO) CALL BACKTR(BB,H,FA,NB)
c     IF (ISW1.EQ.0.OR.ISW1.EQ.9) CALL BACKTR(BB,H,FA,NB)
c     IF (NXT.EQ.1) GO TO 60
c     IF (ISW3.EQ.0 .AND. IPRINT.EQ.1.OR.(ISW3.EQ.0.AND.
c    .  INDEX(CONTRLO,'SCF_IT').NE.0.)) THEN 
c        WRITE (6,9050)
c        WRITE (6,9000) IL
c        IF (NB.GT.999) THEN
c           WRITE (6,9020) (I,EIG(I)+EIGSHF,I=1,NROOT)
c        ELSE
c           WRITE (6,9010) (I,EIG(I)+EIGSHF,I=1,NROOT)
c        END IF
c     END IF
   60 CALL REAR(BB,EIG,W,NB,NROOT,IFLAG)
c     IF (NXT.EQ.1) CALL SYMORB(BB,EIG,NU,NW,IIMO,NB,IONE,IFLAG)
c     IAB = 1
c     IF (IPUN.GT.8) CALL BACKUP(BB,IT,NR,IAB)

        if (.not.rtdb_put(rtdb,'semiemp:eiga',mt_dbl,nb,EIG))
     1  call errquit('diagu: failed to write eig',0,
     1  RTDB_ERR)

         if (.not.rtdb_put(rtdb,'semiemp:acoefficient',mt_dbl,nb*nb,BB))
     1  call errquit('diagu: failed to write alpha coefficient',0,
     1  RTDB_ERR)

       if (.not.rtdb_put(rtdb,'semiemp:aoccupation',mt_dbl,nb,W))
     1  call errquit('diagu: failed to write alpha occupation',0,
     1  RTDB_ERR)

       
   70 CONTINUE
c ***********
C **********
C     WRITE(6,*) ' DIAGU BB'
C     CALL PRNT(BB,NB,NR)
C     WRITE (14) (BB(I),I=1,NR)
c     CALL ZRWRITE (IO_VEC_DEN, BB, 1, NR, 'DIAGU')
C
C      CLOSED SHELL RHF DENSITY OR BETA UHF DENSITY
C
      NROOT1 = NB
c     write(6,*)"BB after densit",(BB(I),I=1,nrr)
c     write(6,*)"FA after densit",(F(I),I=1,nrr)
c     write(6,*)"W after densit",(W(I),I=1,nb)

      CALL DENSIT(BB,F,W,NB,NROOT1,NRR)
c     IF (MULT.LT.3) THEN
         call ga_dgop(Msg_diagu,F,NRR,'+')

c     END IF
c     write(6,*)"BB after densit",(BB(I),I=1,nrr)
c     write(6,*)"FA after densit",(F(I),I=1,nrr)
c     write(6,*)"FA after densit",(W(I),I=1,nb)
c     write(6,*)"MULT",MULT
c     write(6,*)"F",(F(I),I=1,nrr)
C     CALL DROPDEN(F,RAD,IT)

      IF (IT.NE.0) GO TO 90
      IF (MULT.NE.1) GO TO 90
c     IF (IPUN.EQ.8 .OR. IPUN.EQ.10) GO TO 90
C
C     CASE OF SINGLET STATE UHF. IS THERE A CRACKED SYM. SLN?
C

      L = 0
      DO 80 I = 1, NB
         L = L + I
         KA = NU(I)
         DENOM = DFLOAT(NTYP(KA)**2)
         IF (DENOM.GT.9.0D0) DENOM = 9.0D0
         FAC = 1.0D-3/DENOM
         F(L) = F(L) - ((-ONE)**KA)*FAC
c      write(6,*)"F(L)",F(L)
   80 CONTINUE

C
C      F IS NOW THE  CHARGE AND BOND ORDER MATRIX
C
   90 CONTINUE

C     WRITE (14) (F(I),I=1,NRR)
c     CALL ZRWRITE (IO_VEC_DEN, F, 1, NRR, 'DIAGU')
c     WRITE(6,*) ' DIAGU F'
C     CALL PRNT(F,NB,NRR)
C
C
C     ALPHA UHF
C
C
C

      DO 100 I = 1, NB
         EIGA(I) = EIG(I)
         WA(I) = W(I)
c     write(6,*)"WA,EIG",WA(I),EIG(I)
  100 CONTINUE
C
      NROOT = MIM(2) + 1
      IF (NASS(1).NE.0) NROOT = NROOT + 22
      IF (NROOT.GT.NB) NROOT = NB
      NC = 0
      NORBS = NB - NC
      NSTART = NC*NB + 1
c     write(6,*)"NDN",NDN
      DO 210 JJ = 2, NDN
c     write(6,*)"IT",IT
         IF (IT.GT.0) GO TO 120
c     write(6,*)"HERE"
C
C        START WITH VECTORS
C
c        IF (IPUN.NE.8 .AND. IPUN.NE.10) GO TO 110
c        CALL RECOVR('        ',BB,IT,NR,JJ)
  110    CONTINUE
         NROOT = NORBS
c        write(6,*)"nroot,MIM(1),MIM(2),NEL,MULT",nroot,MIM(1),MIM(2),
c    1   NEL,MULT
         CALL ASSIGN(BB,NC,NROOT,JJ,IT,NASS,MIM,FOP,ISWANN,
     1        ISW3,IPUN,NB,EIG,NEL,MULT,W)

         CALL REAR(BB(NSTART),EIG(NC+1),W(NC+1),NB,NROOT,IFLAG)
         GO TO 180
C
C
C
  120    CONTINUE
C        READ (16) (FA(I),I=1,NRR)
C        WRITE(6,*) ' DIAGU  FOCK BETA'
C        CALL PRNT(FA,NB,NRR)
c        CALL ZRREAD (IO_SCF_MATRICES, FA, 1, NRR, 'DIAGU')
c        IF (SCFTOL.LT.ZERO) CALL FOCKTR(FA,H,NB,DRYRUN)
c        IF (ISW1.EQ.0.OR.ISW1.EQ.9) CALL FOCKTR(FA,H,NB,DRYRUN)
C
C        CLEAN UP FOCK MATRIX
C
C  DO YOU NEED THIS???
c        DO 130 I = 1, NRR
c           IF (ABS(FA(I)).LE.1.0D-8) FA(I) = ZERO
  130    CONTINUE
C
C
C
         NROOT = MIM(JJ) + 1
         IF (NASS(1).NE.0) NROOT = NROOT + 22
         IF (NXT.EQ.1) NROOT = NB
         IF (NROOT.GT.NORBS) NROOT = NORBS
C
c 140    IF (IT.NE.0) THEN
c           IBGIV = .FALSE.
c           IF ((IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1) THEN
c              FLAM = TWO/SQRT(FLOAT(IL))
c              CALL LEVELSHF(BB(1),FA,BB(NRR+1),FLAM,NB,ITWO)
c            REWIND 14
c              CALL ZREWIND (IO_VEC_DEN, 'DIAGU')
c            READ (14)
c              CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'DIAGU')
C            READ (14)
c              CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'DIAGU')
c           END IF
c        IF(ISW1.EQ.8.OR.ISW1.EQ.9.OR.JSW(16).EQ.1) THEN
C
C           REWIND 16
C           READ (16)
C           READ (16)
C
C           GET FOCK TRANSFORMATION MATRIX
C
C           READ (16) (H(I),I=1,NRR)
C           READ (16)  SKIP BETA SPIN FOCK MATRIX
C           READ (16)  READ ALPHA SPIN FOCK MATRIX
C           WRITE(6,*) ' DIAGU  FOCK ALPHA'
C           CALL PRNT(FA,NB,NRR)
c           CALL FOCKTR(FA,H,NB,DRYRUN)
C           WRITE(6,*) ' DIAGU  FOCK ALPHA'
C           CALL PRNT(FA,NB,NRR)
C        ENDIF

C        WRITE(6,*) ' IN DIAGU FA'
C        CALL PRNT(FA,NB,NRR)
       call diag_para(FA,NRR,NB,BB,EIG(NC+1))

c      CALL GIVENS(FA,BB,EIG(NC+1),NB,NROOT,IBGIV,QQ,NIN)
c        write(6,*)"FA(I)",(FA(I),I=1,nrr)
c        write(6,*)"AA",(BB(nrr+I),I=1,nrr)
c           IF ((IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1) THEN
c              DO 150 I = MIM(1) + 1, NROOT
c                 EIG(I+NC) = EIG(I+NC) - FLAM
c 150          CONTINUE
c           END IF
c        END IF
C
C        CLEAN UP VECTORS
C

C   DO YOU NEED THIS???
         DO 160 I = 1, NRR
            IF (ABS(BB(I)).LE.1.0D-8) BB(I) = ZERO
c           FA(I) = BB(I)
  160    CONTINUE

C
C
         CALL ASSIGN(BB,NC,NROOT,JJ,IT,NASS,MIM,FOP,ISWANN,
     1        ISW3,IPUN,NB,EIG(NC+1),NEL,MULT,W(NC+1))

c        IF (ISW1.EQ.0 .OR. ISW1.EQ.9) CALL EXOUT(RAD,W,BB,H,FA,F,
c    .       EIG(NC+1),NB,NXT,DRYRUN)
c        IF (SCFTOL.LT.ZERO) CALL BACKTR(BB,H,F,NB)
c        IF (ISW1.EQ.0.OR.ISW1.EQ.9) CALL BACKTR(BB,H,F,NB)
c        IF (NXT.NE.0) GO TO 170
c        IF (ISW3.EQ.0 .AND. IPRINT.EQ.1.OR.(ISW3.EQ.0.AND.
c    .     INDEX(CONTRLO,'SCF_IT').NE.0.)) THEN 
c           WRITE (6,9040)
c           WRITE (6,9030)
c           IF (NB.GT.999) THEN
c              WRITE (6,9020) (I,EIG(I)+EIGSHF,I=NC+1,NC+NORBS)
c           ELSE
c              WRITE (6,9010) (I,EIG(I)+EIGSHF,I=NC+1,NC+NORBS)
c           END IF
c        END IF
c        write(6,*)"here 3",NC
c        write(6,*)"EIG",(EIG(I),I=1,6)
  170    CALL REAR(BB,EIG(NC+1),W(NC+1),NB,NROOT,IFLAG)
c        write(6,*)"EIG",(EIG(I),I=1,6)
c        IF (NXT.EQ.1) CALL SYMORB(BB,EIGA,NU,NW,IIMO,NB,ITWO,IFLAG)
c        IAB = 2
c        IF (IPUN.GT.8) CALL BACKUP(BB,IT,NR,IAB)
c
C        OPEN SHELL RHF OR ALPHA UHF DENSITY
C
        if (.not.rtdb_put(rtdb,'semiemp:bcoefficient',mt_dbl,nb*nb,BB))
     1  call errquit('diagu: failed to write beta coefficient',0,
     1  RTDB_ERR)

       if (.not.rtdb_put(rtdb,'semiemp:boccupation',mt_dbl,nb,W(NC+1)))
     1 call errquit('diagu: failed to write beta ocuupation',0,RTDB_ERR)

        if (.not.rtdb_put(rtdb,'semiemp:eigb',mt_dbl,nb,EIG))
     1  call errquit('diagu: failed to write eig',0,
     1  RTDB_ERR)
       
  180    CONTINUE
C        WRITE (14) (BB(I),I=1,NR)
c        CALL ZRWRITE (IO_VEC_DEN, BB, 1, NR, 'DIAGU')
C
C
C
         CALL DENSIT(BB,FA,W(NC+1),NB,NROOT,NRR)
         call ga_dgop(Msg_diagu,FA,NRR,'+')

         
         NC = NC + MIM(JJ)
         NORBS = NB - NC
         NSTART = NC*NB + 1
         IF (IT.NE.0) GO TO 200
         IF (MULT.NE.1) GO TO 200
         IF (IPUN.EQ.8 .OR. IPUN.EQ.10) GO TO 200
C
C        CASE SINGLET STATE BUT UHF. IS THERE A CRACKED SYM. SLN.?
C
         L = 0
         DO 190 I = 1, NB
            L = L + I
            KA = NU(I)
            DENOM = DFLOAT(NTYP(KA)**2)
            IF (DENOM.GT.9.0D0) DENOM = 9.0D0
            FAC = 1.0D-3/DENOM
            FA(L) = FA(L) + ((-ONE)**KA)*FAC
  190    CONTINUE

C
C
  200    CONTINUE


C        CALL DROPDEN(FA,RAD,IT)

C        WRITE (14) (FA(I),I=1,NRR)
c        CALL ZRWRITE (IO_VEC_DEN, FA, 1, NRR, 'DIAGU')
  210 CONTINUE
C
C     FORM NEW TOTAL DENSITY AND UHF SPIN DENSITY
C
  220 CONTINUE
C     REWIND 14
c     CALL ZREWIND (IO_VEC_DEN, 'DIAGU')
C     READ (14)
c     CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'DIAGU')
C     READ (14) (F(J),J=1,NRR)
c     CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'DIAGU')
C        READ (14)
c     CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'DIAGU')
C        READ (14) (FA(J),J=1,NRR)
c     CALL ZRREAD (IO_VEC_DEN, FA, 1, NRR, 'DIAGU')
c     write(6,*)"F(J)",(F(J),J=1,nrr)
c     write(6,*)"FA(J)",(FA(J),J=1,nrr)
c     if (IT.ne.0) then
c     DO  J = 1, NRR
c     write(6,*)"F(J)",F(J)
c     write(6,*)"FA(J)",FA(J)
c        RB(J) = FA(J)
c        RA(J) = F(J)
c     write(6,*)"RA,RB",RA(J),RB(J)
c     write(6,*)"RB(J)",RB(J)
c        ADD = F(J)
c        F(J) = FA(J) + ADD
c        FA(J) = FA(J) - ADD
c     write(6,*)"FA,F",FA(J),F(J)
c     END DO
c     else
      DO 230 J = 1, NRR
c     write(6,*)"F(J)",F(J)
c     write(6,*)"FA(J)",FA(J)
         RA(J) = FA(J)
         RB(J) = F(J)
c     write(6,*)"RA,RB",RA(J),RB(J)
c     write(6,*)"RB(J)",RB(J)
         ADD = F(J)
         F(J) = FA(J) + ADD
         FA(J) = FA(J) - ADD
c     write(6,*)"FA,F",FA(J),F(J)
  230 CONTINUE
c     end if
  240 CONTINUE
C
CAC       F = TOTAL DENSITY
C
C     WRITE (14) (F(I),I=1,NRR)
c     CALL ZRWRITE (IO_VEC_DEN, F, 1, NRR, 'DIAGU')
C
CAC       FA = SPIN DENSITY
C
C     WRITE (14) (FA(I),I=1,NRR)
c     CALL ZRWRITE (IO_VEC_DEN, FA, 1, NRR, 'DIAGU')
c 250 CONTINUE
C
C
c     write(6,*)"FA last",(FA(I),I=1,nrr)
  260 IF (NXT.LE.0) GO TO 310
      DO 270 I = 1, NB
         TEMP = EIGA(I)
         EIGA(I) = EIG(I)
         EIG(I) = TEMP
  270 CONTINUE
  280 CONTINUE
C
C  VECTOR AND DENSITY OUTPUT SECTION
C    NOTE:  F MUST CONTAIN TOTAL DENSITY WHEN RETURNING.
C     See also subr diagr, same comment - TT
C
  290 CONTINUE
c     IF(JSW(12).EQ.1) THEN
c       REWIND(IO_VCORE)
c       READ (IO_VCORE) (BB(I), I = 1,NRR)
c       TRACEV = ZERO
c       DO 295 I = 1, NRR
c       TRACEV = TRACEV + BB(I)*F(I)
c 295   CONTINUE
c       ID = 0
c       TRACEV = TWO*TRACEV
c       DO 298 I = 1,NB
c       ID = ID + I
c       TRACEV = TRACEV - BB(ID)*F(ID)
c 298   CONTINUE
c       WRITE(6,2022) TRACEV
c2022   FORMAT(/,16X,'     VCORE ENERGY =', F12.8,' au',//)
C     ENDIF                                                  
c     CALL SCFOUT(H,IIMO,IT,IXS,DRYRUN)
      IL = 0
      NXT = 0
      CALL SLITE(IONE,JSLIT)
C      if(nxt.eq.0) go to 5000
      DO 300 I = 1, NRR
         BB(I) = ZERO
  300 CONTINUE
C      GO TO 5000
C 4900 WRITE(6,2045)
C      CALL ZEXIT(3)
  310  CONTINUE
       RETURN
C
C     Format Statements
C
 9000 FORMAT (/,' EIGENVALUES',10X,'DIAG. NO. ',I3)
C 2045 FORMAT(' ***** ERROR - AT LEAST ONE SCF CALC. FOR C.I. CALC IS
C     1 NEEDED *****')
 9010 FORMAT (6 (I4,F9.5))
 9020 FORMAT (6 (I5,F8.4))
 9030 FORMAT (' EIGENVALUES')
 9040 FORMAT (' ALPHA SPIN')
 9050 FORMAT (' BETA SPIN')

      END

