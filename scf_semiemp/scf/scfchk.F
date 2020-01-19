      SUBROUTINE SCFCHK(F,FA,IT,ITGL,ID,UPDATE_DIFF_ONLY,DE1,NDN,NB,NA
     1 ,ISW1,ISW3,NTYP,NU,JSLIT,NRR,PNOW,ISWANN,PRA,PRB,RA,RB,XBACK1,
     2 XBACK2,PBACK1,slowscf)
C
C     F IS NOW DETERMINED FOR THE MOLECULE FROM THE ATOMIC FACTORS
C     WEIGHTED BY THE NUMBER OF BASIS FUNCTIONS ON THE ATOM
C     ITAPE IS FILE THAT CONTAINS THE LATEST DENSITIES
C     VECTORS AND DENSITIES OBTAINED FROM PREVIOUS CYCLE ARE ON UNIT 14
C     EXTRAPOLATED VALUES TO BE USED ARE ON JTAPE
C     PNOW IS THE CURRENT DENSITY (OUT)
C     PBACK1 IS THE PREVIOUS DENSITY (OUT)
C     XBACK1 IS THE PREVIOUS DENSITY (EXTRAPOLATED)
C     XBACK2 IS THE N-2 DENSITY (EXTRAPOLATED)
C
C            PROCEDURE BASED ON
C                 DYNAMICAL DAMPING...
C                 ZERNER ET HEHENBERGER
C                 CHEM. PHYS. LETT.
C                 62, 3, (1979) P550-554
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
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
     +   ,ITEN = 10,IDD1=0)
      LOGICAL      DIIS, DIISCHK, DDFLAG, UPDATE_DIFF_ONLY
      DIMENSION    F(*), FA(*), DIFF(50),PNOW(NA,2),PBACK1(NA,2)
     + ,FAC(33),XBACK1(NA,2),XBACK2(NA,2),NTYP(*),NU(*),JSLIT(4)
     + ,PRA(*),PRB(*),RA(*),RB(*)
      logical slowscf
      SAVE         DIFF
C     DIMENSION XBACK3(MXA,8)
C     DIMENSION    IWAIT(8)
C
C     DATA         IWAIT/8*0/
C
C   CVT IS THE CONVERGENCE RATIO -> IF DIFF(IT)/DIFF(IT-1) .GT. CVT
C   THEN CONVERGENCE IS ACCELERATED.
C
C     CVT = 0.6D0
C
       DIISCHK=.false.
C   PACK ENERGY CONVERGENCE HISTORY INTO DIFF
       IPRINT = 0
       K = 0
C
             if (.not.MA_Push_Get(MT_Dbl,nb*nb,'TRA',lTRA,iTRA))
     &   call errquit('scfit: cannot allocate RA',0,
     1    MA_ERR)
c     if (ga_nodeid.eq.0)
c     write(6,*)"i am in scfchk" ,ID
      JJ = ID
c     if (IT.eq.0.and.ID.eq.1) JJ=2
c     if (IT.eq.0.and.ID.eq.2) JJ=1
      IF (IT.GT.50) THEN
         DO 10 I = 2, 50
            DIFF(I-1) = DIFF(I)
   10    CONTINUE
         DIFF(50) = DE1
      ELSE
         DIFF(IT) = DE1
      END IF

      IF (UPDATE_DIFF_ONLY) RETURN

c     write(6,*)"scfchk passed here",J
      TEST = 1.0D-8
      IFL = 0
      CALL SLITET(2,J,JSLIT)
c     write(6,*)"J before IFL",J
      IF (J.EQ.2) GO TO 20
              IFL = 1
       CALL SLITE(2,JSLIT)
   20  CONTINUE
C
      DO 30 I = 1, NA
c     write(6,*)"PNOW",PNOW(I,JJ)
         PBACK1(I,JJ) = 0.0D0
         IF (IT.GT.0) PBACK1(I,JJ) = PNOW(I,JJ)
         PNOW(I,JJ) = 0.0D0
   30 CONTINUE
      L = 0
      DO 40 I = 1, NB
         K = NU(I)
         L = L + I
c      write(6,*)"F(L),PNOW",F(L),PNOW(K,JJ)
         PNOW(K,JJ) = PNOW(K,JJ) + F(L)
   40 CONTINUE
      if (.not.slowscf) GO TO 45         
      IF (IT.GT.21) GO TO 80
      IF (IT.GT.20) GO TO 74
      IF (IT.GT.15) GO TO 73
      IF (IT.GT.10) GO TO 72
      IF (IT.GT.5) GO TO 71
      IF (IT.GT.0) GO TO 70

   45 CONTINUE
      IF (IT.GT.1) GO TO 80
      IF (IT.GT.0) GO TO 74      
C
C     CYCLE ZERO
C
      DO 50 I = 1, NA
       XBACK2(I,JJ) = PNOW(I,JJ)
c        write(6,*)"JJ,XBACK2(I,JJ)",JJ,XBACK2(I,JJ)
   50 CONTINUE
c     write(6,*)"NDN, JJ",NDN,JJ
      IF (JJ.LT.NDN) GO TO 220
C     REWIND 12
c     CALL ZREWIND (IO_SCF_DENSITY, 'SCFCHK')
C     REWIND 17
c     CALL ZREWIND (IO_DENSEXT, 'SCFCHK')
C     REWIND 14
c     CALL ZREWIND (IO_VEC_DEN, 'SCFCHK')
c     ITAPE = IO_DENSEXT
c     DO 60 I = 1, NDN
C     READ (14)
c        CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'SCFCHK')
C     READ (14) (FA(J),J=1,NRR)
c        CALL ZRREAD (IO_VEC_DEN, FA, 1, NRR, 'SCFCHK')
C  25 WRITE (17) (FA(J),J=1,NRR)
c        CALL ZRWRITE (IO_DENSEXT, FA, 1, NRR, 'SCFCHK')
c  60 CONTINUE
C     WRITE(6,*) ' SCFCHK', ITAPE
C     CALL PRNT (FA,NB,NRR)
      GO TO 220
C
C     FIRST CYCLE
C
   70 FAC(JJ) = 0.98D0
      IF (JJ.LT.NDN) GO TO 220
      GO TO 100      
   71 FAC(JJ) = 0.95D0
      IF (JJ.LT.NDN) GO TO 220
      GO TO 100      
   72 FAC(JJ) = 0.90D0
      IF (JJ.LT.NDN) GO TO 220
      GO TO 100      
   73 FAC(JJ) = 0.60D0
      IF (JJ.LT.NDN) GO TO 220
      GO TO 100      
   74 FAC(JJ) = 0.40D0   
      IF (ISW1.EQ.0) FAC(JJ) = 0.50D0
      IF (JJ.LT.NDN) GO TO 220
      GO TO 100
C
C     CYCLE 2=>N-1
C
   80 SUM = 0.0D0
c     write(6,*)"here IT",IT
c     write(6,*)"XBACK1",((XBACK1(I,J),I=1,1),J=1,1)
c     write(6,*)"XBACK2",((XBACK2(I,J),I=1,1),J=1,1)
      DO 90 I = 1, NA
c     write(6,*)"XBACK",XBACK1(I,JJ),XBACK2(I,JJ)
         DENOM = XBACK1(I,JJ) - XBACK2(I,JJ)
c     write(6,*)"DENOM",DENOM
C
C   START SAVING CONVERGENCE INFORMATION AFTER SEVENTH CYCLE
C
C     IF (IT .GE. 7) THEN
C        XBACK3(I,JJ) = XBACK2(I,JJ)
C     ENDIF
C
         XBACK2(I,JJ) = XBACK1(I,JJ)
         IF (ABS(DENOM).LT.TEST) GO TO 90
c       write(6,*)"PNOW,PNOW,DENOM",PNOW(I,JJ),PBACK1(I,JJ),DENOM
         SLOPE = (PNOW(I,JJ)-PBACK1(I,JJ))/DENOM
         DIISCHK = .FALSE.
         IF (SLOPE.GT.0.0D0 .OR. SLOPE.LT.-19.0D0) GO TO 90
         DIISCHK = .TRUE.
         L = NTYP(I)
         N = L*L
         IF (L.EQ.4) N = 9
         IF (L.GT.4) N = 16
C
C     STORE ATOMIC FACTOR IN XBACK1 TEMPORARLY
C
c        write(6,*)"SLOPE,N",SLOPE,N
         XBACK1(I,JJ) = SLOPE/ (SLOPE-ONE)
         SUM = SUM + N*XBACK1(I,JJ)
   90 CONTINUE
c     IF(.NOT.DIISCHK) PRINT *,'SLOPE OUT OF RANGE',IT,SLOPE
      FAC(JJ) = SUM/NB
c     write(6,*)"IFL",IFL
      IF (JJ.LT.NDN) GO TO 220
      IF (IFL.EQ.1) GO TO 190
C
C     NOW INTERPOLATE THE DENSITIES
C         ITAPE HOLDS OLD DENSITIES (EXTRAPOLATED)
C         (14)  HOLDS NEW DENSITIES
C
  100 CONTINUE
      IF (ISWANN .EQ. IZERO)  then
      IF (IT.GT.0) THEN
       DO I = 1, NRR
         FA(I) = 2*FA(I)
       END DO
      ELSE
       DO   I = 1, NRR
         FA(I)=F(I)
       END DO
      END IF
      END IF

C     REWIND 12
c     CALL ZREWIND (IO_SCF_DENSITY, 'SCFCHK')
C     REWIND 17
c     CALL ZREWIND (IO_DENSEXT, 'SCFCHK')
C     REWIND 14
c     CALL ZREWIND (IO_VEC_DEN, 'SCFCHK')
c     IF (ITAPE.EQ.IO_SCF_DENSITY) JTAPE = IO_DENSEXT
c     IF (ITAPE.EQ.IO_DENSEXT) JTAPE = IO_SCF_DENSITY
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
      if (IT.ge.2) then
       TEMP = FAC(1)
       FAC(1) = FAC(2)
       FAC(2) = TEMP
      end if
      end if
      DO 180 IOP = 1, NDN
C     READ (ITAPE) (FA(I),I=1,NRR)
c        CALL ZRREAD (ITAPE, FA, 1, NRR, 'SCFCHK')
C     READ(14)
c        CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'SCFCHK')
C     READ (14) (F(I),I=1,NRR)
c        CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'SCFCHK')
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
      if (IOP.eq.1) then
       DO I =1, NRR
        TEMP = FA(I)
        FA(I) = PRB(I)
        PRB(I) = TEMP
       END DO
      else
       DO I =1, NRR
        F(I) = PRB(I)
        FA(I) = PRA(I)
       end do
      end if
      end if
c     write(6,*)"after 1,FA",(FA(I),I=1,nrr)
c     write(6,*)"after 1,F",(F(I),I=1,nrr)
         IJ = 1
         DO 120 I = 1, NB
            DO 110 J = 1, I
c         write(6,*)"F,FAC,FA",F(IJ),FAC(IOP),FA(IJ)
               F(IJ) = FAC(IOP)*FA(IJ) + (1.0D0-FAC(IOP))*F(IJ)
               IJ = IJ + 1
  110       CONTINUE
  120    CONTINUE
C
C     WRITE (JTAPE) (F(I),I=1,NRR)
c        CALL ZRWRITE (JTAPE, F, 1, NRR, 'SCFCHK')
         II = 0
         KOLD = 0
         DO 130 I = 1, NB
            K = NU(I)
            IF (K.NE.KOLD) XBACK1(K,IOP) = 0.0D0
            KOLD = K
            II = II + I
c           write(6,*)"XBACK1",XBACK1(K,IOP),F(II)
            XBACK1(K,IOP) = XBACK1(K,IOP) + F(II)
c           write(6,*)"XBACK1",XBACK1(K,IOP)
  130    CONTINUE
         IF (K.GE.NA) GO TO 150
C
c     point charges are INCLUDEd in na, but have no basis sets
C
         KPLUS = K + 1
         DO 140 I = KPLUS, NA
            XBACK1(I,IOP) = 0.0D0
  140    CONTINUE
  150    CONTINUE
c        IF ((ISW3.EQ.0 .AND. IPRINT.GE.1.).OR.(ISW3.EQ.0.
c    .       AND.INDEX(CONTRLO,'SCF_ITER').NE.0)) THEN
c           WRITE (6,9050)
c           WRITE (6,9030) IOP
c           NSTART = -7
c           DO 160 I = 1, NA, 8
c              NSTART = NSTART + 8
c              NSTOP = NSTART + 7
c              IF (NSTOP.GT.NA) NSTOP = NA
c              WRITE (6,9000) (XBACK2(J,IOP),J=NSTART,NSTOP)
c              WRITE (6,9000) (PNOW(J,IOP),J=NSTART,NSTOP)
c              WRITE (6,9020) FAC(IOP)
c              WRITE (6,9000) (XBACK1(J,IOP),J=NSTART,NSTOP)
c              WRITE (6,9050)
C     IF(NSTOP.EQ.NA) GO TO 109
c 160       CONTINUE
c 170       CONTINUE
c        END IF
c     write(6,*)"F(I)",(F(I),I=1,nrr)
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
      if (IOP.eq.1) then
       DO I =1, NRR
        Dbl_MB(iTRA+I-1) = F(I)
c      write(6,*)"F(I)",F(I)
       END DO
      else
       DO I = 1, NRR
        RA(I) = F(I)
        RB(I) = Dbl_MB(iTRA+I-1)
c      write(6,*)"RA,RB",RA(I),RB(I)
       END DO
      end if 
      end if 
  180 CONTINUE
      ITAPE = JTAPE
      GO TO 220
C
C     LAST CYCLE
C
  190 CONTINUE
C     REWIND 14
C     REWIND 12
c     CALL ZREWIND (IO_SCF_DENSITY, 'SCFCHK')
c     CALL ZREWIND (IO_VEC_DEN, 'SCFCHK')
      DO 210 J = 1, NDN
         IF(ITGL.EQ.0.OR.ITGL.GT.4) THEN
         IF (ISW3.EQ.0 .AND. IPRINT.GE.0) THEN
            IF (ABS(DIFF(IT)).GT.ABS(DIFF(IT-1))) THEN
               DO 200 I = 1, NA
                  IF ((XBACK1(I,J).GT.0.1D0) .AND.
     .                (XBACK1(I,J).LE.1.0D0).AND.ga_nodeid().eq.0) 
     .                WRITE (6,9010) I, J, XBACK1(I,J)
  200          CONTINUE
            END IF
         ENDIF
         END IF
C     READ (14)
C     READ (14) (F(I),I=1,NRR)
c        CALL ZIREAD (IO_VEC_DEN, IDU, 1, 1, 'SCFCHK')
c        CALL ZRREAD (IO_VEC_DEN, F, 1, NRR, 'SCFCHK')
c     WRITE (12) (F(I),I=1,NRR)
c        CALL ZRWRITE (IO_SCF_DENSITY, F, 1, NRR, 'SCFCHK')
  210 CONTINUE
c     ITAPE = IO_SCF_DENSITY
C
C
  220 CONTINUE
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
c     write(6,*)"HERE U"
      if (IT.ge.1.and.JJ.eq.2) then
c     write(6,*)"exchange"
      DO I = 1, NA
       TEMP = XBACK1(I,1)
       XBACK1(I,1)=  XBACK1(I,2)
       XBACK1(I,2)= TEMP
c      TEMP = XBACK2(I,1)
c      XBACK2(I,1)=  XBACK2(I,2)
c      XBACK2(I,2)= TEMP
c      TEMP = PNOW(I,1)
c      PNOW(I,1) = PNOW(I,2)
c      PNOW(I,2) = TEMP
      end do
      end if
c      RA(I) = TRA(I)
c      RB(I) = TRB(I)
c     END DO
c     END IF
c     IF (IT.EQ.0) then
c      DO I = 1, NRR
c       FA(I) = F(I)
c      END DO
c     END IF
      ENDIF
c     write(6,*)"F(I) in scfchk",(F(I),I=1,nrr)
c     write(6,*)"FA(I) in scfchk",(FA(I),I=1,nrr)
C
C     RETURN TO SING
C
c     write(6,*)"J in SCFCHK",J
      if (.not.ma_pop_stack(lTRA))
     &   call errquit('scfit: cannot pop stack lPRA',3, MA_ERR)
      RETURN

 9000 FORMAT (1X,8 (F8.5,1X))
 9010 FORMAT (1X,'***WARNING: DENSITY ON ATOM',I3,'NOT CONVERGED ',
     .       'FOR OPERATOR = ',I2,' F(N-1)= ',F9.5)
 9020 FORMAT (1X,'FACTOR = ',F10.6)
 9030 FORMAT (1X,'ATOMIC DENSITIES FOR OPERATOR = ',I2)
 9040 FORMAT (1X,'*** EXTRAPOLATION FOR CREEP DONE THIS CYCLE....')
 9050 FORMAT ()

      END
                                





