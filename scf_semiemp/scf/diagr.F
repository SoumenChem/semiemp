       SUBROUTINE DIAGR(rtdb,BB,F,H,FA,RAD,ITMAX,IT,PREVS,SCFTOL,MIM,
     1  ISW3,NB,EIG,FOP,ISWANN,NEL,MULT,W,JSLIT,NRR,Q,NIN,NCORE,NASS)
C
C     ..................................................................
C
C     THIS ROUTINE DIAG. THE FOCK MATRIX AND CALCULATES THE BOND ORDER
C     MATRIX, WHICH IS THEN CHECKED BY SUBROUTINE SCFCHK FOR OSCILLATION
C     .... FOR RHF CLOSED SHELL ....
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

      IMPLICIT DOUBLE PRECISION(A-F,O-Z)
      IMPLICIT INTEGER (G-N)
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
      LOGICAL      IBGIV, DRYRUN
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
      dimension  F(*),H(*),FA(*),RAD(*),W(NB),Q(*)
     .           ,NASS(7),MIM(31),BB(*),EIG(*),NIN(*)
      DIMENSION    IIMO(NB),JSLIT(4)
      SAVE         ISWAN1, IL, NXT

      integer g_A,g_F,g_alpha 

      IF (IT.NE.0) GO TO 20
C
C     FIRST ITERATION ONLY!
c      if (ga_nodeid().eq.0)
c    . write(6,*)"i am in diagr"!,(BB(I),I=1,15)
      ISWAN1 = (ISWANN+1)/2 + 1
      IL = 0
      NXT = 0
C     STARTING VECTORS FROM SUBROUTINE DIAGST?

c     IF (LINFLAG.AND.ITGL.GT.ITWO) THEN
c      CALL GLOC(F,NOWF,IT,DRYRUN)
C       CALL PRNT(F,NB,NRR)
c      CALL ZREWIND(IO_VEC_DEN,'DIAGR')
c      CALL ZRWRITE(IO_VEC_DEN,BB,1,NR,'DIAGR')
c      CALL ZRWRITE(IO_VEC_DEN,F,1,NRR,'DIAGR')
c      GO TO 80
c     END IF

      IF (ISW3.EQ.0) GO TO 10
c     GO TO 80

C     STARTING WITH DENSITY
c     IF(ISAVE.EQ.2) GO TO 120

   10 CONTINUE
      IPUN1 = IPUN + 1
c     GO TO (20,20,80,80,20,70,70,70,70,
c    .       20,70), IPUN1

C      STARTING ITERATION WITH PREVIOUSLY CALCULATED DENSITY
C
C
C     FOR ALL CYCLES (JUST PLAIN START)
C
   20 CONTINUE


      IL = IL + 1
C     WRITE(6,*) ' DIAGR START ', IPUN, NXT
C
C     CLEAN UP MATRIX TO AVOID PROBLEMS IN DIAGONALIZATION
C
C   REMOVED BY JDHD JAN 17 84
C   IF ZEROING F(I) USE THRESH OF ORDER 1.D-11
C      DO 60 I=1,NRR
C        IF(ABS (F(I)).LE.1.0D-8)  F(I)=ZERO
C   60 CONTINUE
c     write(6,*)"MIM",MIM
      NROOT = MIM(1) + 5
      IF (NASS(1).LT.0) NROOT = NROOT + 22
      IF (NASS(1).GT.0) NROOT = NROOT + NASS(1)
      IF (NROOT.GT.NB) NROOT = NB
      IF (ISWANN.EQ.3) NROOT = NB
      CALL SLITET(2,J,JSLIT)
C     IF CALL FOR NO ITERATIONS QUIT
      IF (ITMAX.EQ.0.OR.J.EQ.1) THEN
C
C     LAST TIME THROUGH--CONVERGED!! (MAYBE)
C
        NXT = 1
        NROOT = NB
        CALL SLITE(2,JSLIT)
        CALL SLITE(1,JSLIT)
      ELSE
        NXT = 0
      ENDIF
C
C     CLOSED SHELL RHF
C
      IBGIV = .TRUE.
C     IS LEVEL SHIFTING ADVANTAGEOUS?
C     IF(IDD2.NE.0.AND.(EIG(NEL/2+1)-EIG(NEL/2)).LT.0.05) THEN
C       IF(IT.NE.0) IDD2 = 2
C     ENDIF
c     write(6,*)"IDD2 in diagr",IDD2
c     IF (IT.GT.10 .AND. IDD2.EQ.21) IDD2 = 1
c     IF (IT.GT.10 .AND. IDD2.EQ.20) IDD2 = 0
c     IF ( (IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1 .AND. IT.NE.0) THEN
c        FLAM = ONE/ (FLOAT(IT))
c        CALL LEVELSHF(BB(1),F,BB(NRR+1),FLAM,NB,IONE)
c      REWIND THIS UNIT FOR WRITE
c        CALL ZREWIND (IO_VEC_DEN, 'DIAGR')
c     END IF
C     WRITE(6,*) ' DIAGR PAST ALL SHIFTING', IDD2, NXT
c     IF(ISW1.EQ.8.OR.ISW1.EQ.9.OR.JSW(16).EQ.1) THEN
C
C        REWIND 16
c        CALL ZREWIND (IO_SCF_MATRICES, 'DIAGR')
C        READ (16)
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGR')
C        READ (16)
c        CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, 1, 'DIAGR')
C
C        GET FOCK TRANSFORMATION MATRIX
C
C        READ (16) (H(I),I=1,NRR)
c        CALL ZRREAD (IO_SCF_MATRICES, H, 1, NRR, 'DIAGR')
c        IF (IT.NE.0) CALL FOCKTR(F,H,NB,DRYRUN)
c     ENDIF
c     call mat_sym(F,NB)
c     WRITE(6,*) ' DIAGR BEFORE GIVENS F ',(F(I),I=1,nrr)
c     WRITE(6,*) ' DIAGR BEFORE GIVENS BB',(BB(I),I=1,20)
c     WRITE(6,*) ' DIAGR BEFORE GIVENS EIG',(EIG(I),I=1,20)
c      call ga_sync()
C----------------------------------------------------------------------------------------
C ORTHOGONAL SIMILARTIY TRANSFORMATION OF OVERLAP MATRIX. ITS DONE IN
C TREDS AND THEN REBUILD THE SYMMETRIC NB*NB OVERLAP MATRIX
C------------------------------------------------------------------------

c     WRITE(6,*) ' DIAGR BEFORE GIVENS H ', (H(I),I=1,40)
c     call mat_sym(F,nb)
c     CALL GIVENS(F,BB,EIG,NB,NROOT,IBGIV,QQ,NIN)
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
      call diag_para(F,NRR,NB,BB,EIG)
c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "diagr cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

c     WRITE(6,*) ' DIAGR AFTER GIVENS F ',(BB(I),I=1,21)
c     WRITE(6,*) ' DIAGR AFTER GIVENS BB ',(BB(I),I=1,3)
c     WRITE(6,*) ' DIAGR AFTER GIVENS EIG',(EIG(I),I=1,20)
C     THE VECTORS IN BB ARE IN THE ORTHOGONAL BASIS SET
      
C
c     IF ((IDD2.EQ.20.OR.IDD2.EQ.21) .AND. NXT.NE.1 .AND. IT.NE.0) THEN
C      WRITE(6,2222) (EIG(I),I=1,NROOT)
C2222  FORMAT(10F8.6)
c        DO 50 I = NEL/2 + 1, NROOT
c           EIG(I) = EIG(I) - TWO*FLAM
c  50    CONTINUE
C      WRITE(6,2222) (EIG(I),I=1,NROOT)
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
c     WRITE(6,*) ' DIAGR BEFORE ASSIGN',(BB(I),I=1,20)
      CALL ASSIGN(BB,IZERO,NROOT,IONE,IT,NASS,MIM,FOP,ISWANN,ISW3,IPUN
     1 ,NB,EIG,NEL,MULT,W,NCORE)
c     IF (ISW1.EQ.0 .OR. ISW1.EQ.9) CALL EXOUT(RAD,W,BB,H,FA,F,EIG,NB,
c    .    NXT,DRYRUN)
c     IF (SCFTOL.LT.ZERO) CALL BACKTR(BB,H,FA,NB)
c     IF (ISW1.EQ.0.OR.ISW1.EQ.9) CALL BACKTR(BB,H,FA,NB)
c     IF (NXT.EQ.1) GO TO 60
c     IF (ISW3.EQ.0 .AND. IPRINT.GE.1.OR.
c    .     INDEX(CONTRLO,'SCF_ITER').NE.0) THEN
c        WRITE (6,9000) IL
c        IF (NB.GT.999) THEN
c           WRITE (6,9020) (I,EIG(I)+EIGSHF,I=1,NROOT)
c        ELSE
c           WRITE (6,9010) (I,EIG(I)+EIGSHF,I=1,NROOT)
c        END IF
c     END IF
   60 CALL REAR(BB,EIG,W,NB,NROOT,IFLAG)

       if (.not.rtdb_put(rtdb,'semiemp:coefficient',mt_dbl,nb*nb,BB))
     1  call errquit('diagr: failed to write coefficient',0, RTDB_ERR)

       if (.not.rtdb_put(rtdb,'semiemp:occupation',mt_dbl,nb,W))
     1  call errquit('diagr: failed to write ocuupation',0, RTDB_ERR)

      
c     WRITE(6,*) ' DIAGR AFTER REAR ',(Q(I),I=1,20)
   70 CONTINUE
c     IF (NXT.EQ.1)  CALL SYMORB(BB,EIG,NU,NW,IIMO,NB,IONE,IFLAG)
c     write(6,*)"NEED SYMORB"
c     IAB = 1
c     IF (IPUN.GT.8) CALL BACKUP(BB,IT,NR,IAB)
c     WRITE (14) (BB(I),I=1,NR)
c     CALL ZRWRITE (IO_VEC_DEN, BB, 1, NR, 'DIAGR')
c
C      CLOSED SHELL RHF DENSITY
C

C     Linear scaling - calculate B-factors and assign with Fermi
C     function

c     IF (LINFLAG) CALL BFACT(BB,NOWF)
c     NROOT = NB
c     IF (LINFLAG) CALL LINASS(W,IT)

C     JD 1999
c     write(6,*)"NROOT,NB",NROOT,NB
      CALL DENSIT(BB,F,W,NB,NROOT,NRR)
c     CALL DROPDEN(F,RAD,IT)
c     IF (MULT.NE.1) THEN
c        WRITE (6,*)
c    .     ' ***** ERROR IN MULTIPLICITY FOR RESTRICTED RHF ****'
c        CALL ZEXIT(3)
c        RETURN


c     END IF
C
c     write(6,*)"i am here ???"
C     WRITE (14) (F(I),I=1,NRR)
c     CALL ZRWRITE (IO_VEC_DEN, F, 1, NRR, 'DIAGR')
C
C
C     FORM EXTRAPOLATED TOTAL DENSITY
C
C
c  80 CONTINUE
C     REWIND 14
c     CALL ZREWIND (IO_VEC_DEN, 'DIAGR')
C     READ (14)
c     CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'DIAGR')
C     READ (14) (F(J),J=1,NRR)
c     CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'DIAGR')
c     write(6,*)"NXT",NXT
   90 IF (NXT.LE.0) GO TO 120
C
C     VECTOR AND DENSITY OUTPUT SECTION
C     NOTE:  F MUST CONTAIN TOTAL DENSITY WHEN RETURNING.
C     WITH THE NEW MEM SCHEME, IT IS CLEARLY NOT GUARANTEED BY SCFOUT. I
C     HOPE THE ZIO ABOVE TAKES CARE - TT
C
c 100 CONTINUE
c     IF(JSW(12).EQ.1) THEN
c       REWIND(IO_VCORE)
c       READ (IO_VCORE) (BB(I), I = 1,NRR)
c       TRACEV = ZERO
c       DO 105 I = 1, NRR
c       TRACEV = TRACEV + BB(I)*F(I)
c 105   CONTINUE
c       ID = 0
c       TRACEV = TWO*TRACEV
c       DO 108 I = 1,NB
c       ID = ID + I
c       TRACEV = TRACEV - BB(ID)*F(ID)
c 108   CONTINUE
c       WRITE(6,2022) TRACEV
c2022   FORMAT(/,16X,'     VCORE ENERGY =', F12.8,' au',//)
c     ENDIF
C     WRITE(6,*) ' DIAGR BEFORE SCFOUT'
c     CALL SCFOUT(H,IIMO,IT,IXS,DRYRUN)
      IL = 0
      NXT = 0
C     CALL SLITE(IONE)
C      if(nxt.eq.0) go to 5000
      DO 110 I = 1, NB*NB
         BB(I) = ZERO
  110 CONTINUE
C      GO TO 5000
C 4900 WRITE(6,2045)
C      CALL ZEXIT(3)
  120 CONTINUE
      RETURN
C
C     Format Statements
C
 9000 FORMAT (' EIGENVALUES',10X,'DIAG. NO. ',I3)
C 2045 FORMAT(' ***** ERROR - AT LEAST ONE SCF CALC. FOR C.I. CALC IS
C     1 NEEDED *****')
 9010 FORMAT (6 (I4,F9.5))
 9020 FORMAT (6 (I5,F8.4))

      END
C-------------------------------------------------------------------
C
C-------------------------------------------------------------------
      SUBROUTINE SLITET(I,K,JSLIT)
      INTEGER I,K
      INTEGER JSLIT(4)
c     write(6,*)"i am in slitet",JSLIT
      K = JSLIT(I)
      JSLIT(I) = 2
      RETURN

      END
c---------------------------------------------------------------------
C
C--------------------------------------------------------------------
      SUBROUTINE SLITE(I,JSLIT)
      INTEGER I
      INTEGER JSLIT(4)
c     write(6,*)"i am in slite"
      IF (I.EQ.0) GO TO 10
C     TURN SENSE LIGHT I ON
      JSLIT(I) = 1
      GO TO 20
C     TURN ALL SENSE LIGHTS OFF
   10 JSLIT(1) = 2
      JSLIT(2) = 2
      JSLIT(3) = 2
      JSLIT(4) = 2
   20 RETURN

      END


