      SUBROUTINE FILLF(BB,FA,F,Q,GAMMA,SUMELEL,NG,NRR,NINTGL,NIN,ITAG,
     1 NB)
C     ZERNER   QTP,   UNIV. OF FL.
C     MODIFIED BY NR   VERSION : 29.10.84
C     CLOSED SHELL F=TOTAL DENSITY, FA= HALF OF F=ALPHA SPIN DENSITY
C     OPEN SHELL F=TOTAL DENSITY, FA=ALPHA SPIN DENSITY
C     Q RETURNS AS THE ONE_ELECTRON MATRIX INCLUDING REACTION FIELD.
C
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
#include "bas.fh"
#include "schwarz.fh"
#include "dftps.fh"
#include "bgj.fh"
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
      dimension  F(*),Q(*),FA(*),GAMMA(*),BB(*),DIP(3)
     +,AD(3),BD(3),BBDD(3),BBD(3),NIN(*),SUMELEL(*),NG(*)
      CHARACTER*8  FLAGM, FLAGU
     
c     write(6,*)"i am in fillf",NINTGL
c     write(6,*)"F",(F(I),I=1,20)
c     write(6,*)"FA",(FA(I),I=1,20)
c     write(6,*)"BB",(BB(I),I=1,20)
c     write(6,*)"Q",(Q(I),I=1,nrr)
          if (.not.MA_Push_Get(MT_INT,NINTGL,'LABI',lLABI,
     1   iLABI))
     &   call errquit('semiemp_energy: cannot allocate LABI',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABJ',lLABJ,
     1   iLABJ))
     &   call errquit('semiemp_energy: cannot allocate LABJ',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABK',lLABK,
     1   iLABK))
     &   call errquit('semiemp_energy: cannot allocate LABK',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABL',lLABL,
     1   iLABL))
     &   call errquit('semiemp_energy: cannot allocate LABL',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,NINTGL,'QQ',lQQ,iQQ))
     &   call errquit('semiemp_energy: cannot allocate QQ',0,
     1    MA_ERR)

        call semiemp_integral_load(NINTGL,INT_MB(iLABI),INT_MB(iLABJ),
     A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ))
           
c      write(6,*)"got the integrals"
      ITDIP = 1
      FLAGU = FLAGU
C     JD 1999
      DO 10 I = 1, NRR
         BB(I) = ZERO
         FA(I) = F(I)/TWO
   10 CONTINUE
c     write(6,*)"BB",(BB(I),I=1,nrr)
c     write(6,*)"F in fillf",(F(I),I=1,nrr)
c     write(6,*)"FA",(FA(I),I=1,nrr)
c     write(6,*)"QQ",(QQ(I),I=1,nb)
      LL = 0
c     NINTGL = 0
C
C
C
c     write(6,*)"LABI",(LABI(I),I=1,20)
c     write(6,*)"LABJ",(LABJ(I),I=1,20)
c     write(6,*)"LABK",(LABK(I),I=1,20)
c      call ga_sync()
c     go to 170
   20 CONTINUE
C     BEGIN TWO _ ELECTRON PART
      LL = LL + 1
c     write(6,*)"LL,NINTGL",LL,NINTGL
      IF (LL.LE.NINTGL) THEN
      GO TO 30
      ELSE
      GO TO 170
      ENDIF
c     IF (NINTGL.LE.0) GO TO 170
c     I = 0
c     J = 0
c     K = 0
c     L = 0
c     LL = 1
   30 CONTINUE
c     write(6,*)"i am here in fillf",LL,NINTGL
      I = int_mb(iLABI+LL-1)
      J = int_mb(iLABJ+LL-1)
      K = int_mb(iLABK+LL-1)
      L = int_mb(iLABL+LL-1)
      V = Dbl_MB(iQQ+LL-1)

c     I = LABI(LL)
c     J = LABJ(LL)
c     K = LABK(LL)
c     L = LABL(LL)
c     V = QQ(LL)
c     WRITE(6,*) 'I J K L V', I,J,K,L,V
C
cyw BB matrix starts with a zero matrix
c      write(6,*)"(I-J)",(I-J)
c      write(6,*)"(K-L)",(K-L)
c      write(6,*)"I-K",(I-K)
c      if ((I-K).eq.0.and.ga_nodeid().eq.0) go to 60
c      if ((I-K).gt.0.and.ga_nodeid().eq.0) go to 180
c      if ((I-K).lt.0.and.ga_nodeid().eq.0) go to 60 
C
      IF ((I-J).LT.0)   GO TO 90
      IF ((I-J).EQ.0)   GO TO 40
      GO TO 180


   40 IF ((K-L).LT.0)  GO TO 80
      IF ((K-L).EQ.0)  GO TO 50
      GOTO 180

   50 IF ((I-K).LT.0)  GO TO 70
      IF ((I-K).EQ.0) GO TO 60
      GOTO 180

c     IF (I-J) 90, 40, 180
c  40 IF (K-L) 80, 50, 180
c  50 IF (I-K) 70, 60, 180
c     (I,I/I,I)
   60 IN = NIN(I) + I
c     write(6,*) "i am here 0",IN,V
c     IF (ITAG.EQ.1) THEN 
      BB(IN) = BB(IN) + (F(IN)-FA(IN))*V
C     FOR NDDO THEORY
c     SUMELEL(I) = SUMELEL(I) + F(IN)*V
C     FOR INDO
c     SUMELEL(I) = SUMELEL(I) + F(IN)*GAMMA(NIN(NG(I))+NG(I))
c     if (ga_nodeid().eq.0) write(6,*) "SUMELEL",SUMELEL(I)
c     GO TO 160
c     ELSE
      GO TO 160
C     (I,I/K,K)
   70 IN = NIN(I) + I
c     write(6,*)"i am here 1",IN
      KN = NIN(K) + K
c     IF (ITAG.EQ.1) THEN
      BB(IN) = BB(IN) + F(KN)*V
      BB(KN) = BB(KN) + F(IN)*V
C     FOR NDDO
C     SUMELEL(I) = SUMELEL(I) + F(KN)*V
C     SUMELEL(K) = SUMELEL(K) + F(IN)*V
C     FOR INDO, CNDO
      IIKK = NIN(NG(K)) + NG(I)
c     SUMELEL(I) = SUMELEL(I) + F(KN)*GAMMA(IIKK)
c     SUMELEL(K) = SUMELEL(K) + F(IN)*GAMMA(IIKK)
      IN = NIN(K) + I
c     write(6,*)"BB,FA,V",BB(IN),FA(IN),V
      BB(IN) = BB(IN) - FA(IN)*V
c     if (ga_nodeid().eq.0) write(6,*)"BB",BB(IN)
c     GO TO 160
c     ELSE
      GO TO 160
C     (I,I/K,L)
   80 IN = NIN(I) + I
c     write(6,*)"i am here 2",IN,KN
      KN = NIN(L) + K
      IF (ITAG.eq.1)
     .BB(IN) = BB(IN) + 2.*F(KN)*V
      BB(KN) = BB(KN) + F(IN)*V
      IN = NIN(K) + I
      KN = NIN(L) + I
      BB(IN) = BB(IN) - FA(KN)*V
      BB(KN) = BB(KN) - FA(IN)*V
      IF (I.NE.K) GO TO 160
      IF (ITAG.eq.1)
     .BB(IN) = BB(IN) - FA(KN)*V
      IF (ITAG.eq.2)
     .BB(IN) = BB(IN) + FA(KN)*V
c     write(6,*)"before leaving"
      GO TO 160

   90 IF (K-L) 110, 100, 180
C     (I,J,/K,K)
  100 IN = NIN(J) + I
c     write(6,*)"i am here 3"
      KN = NIN(K) + K
      BB(IN) = BB(IN) + F(KN)*V
      IF (ITAG.eq.1)
     .BB(KN) = BB(KN) + 2.*F(IN)*V
      IN = NIN(K) + I
      KN = NIN(K) + J
      IF (K.LT.J) KN = NIN(J) + K
      BB(IN) = BB(IN) - FA(KN)*V
      BB(KN) = BB(KN) - FA(IN)*V
      IF (J.NE.K) GO TO 160
      IF (ITAG.eq.1)
     .BB(KN) = BB(KN) - FA(IN)*V
      IF (ITAG.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
      GO TO 160
C     (I,J/K,L)
  110 IN = NIN(J) + I
c     write(6,*)"i am here 4"
      KN = NIN(L) + K
      IF (ITAG.eq.1)
     .BB(IN) = BB(IN) + 2.*F(KN)*V
      IF (IN.EQ.KN) GO TO 120
      IF (ITAG.eq.1)
     .BB(KN) = BB(KN) + 2.*F(IN)*V
  120 IN = NIN(K) + I
c     write(6,*)"i am here 5"
      KN = NIN(L) + J
      IF (J.GT.L) KN = NIN(J) + L
      BB(IN) = BB(IN) - FA(KN)*V
      BB(KN) = BB(KN) - FA(IN)*V
      IF (I.NE.K) GO TO 130
      IF (J.EQ.L) GO TO 140
      IF (ITAG.eq.1)
     .BB(IN) = BB(IN) - FA(KN)*V
      IF (ITAG.eq.2)
     .BB(IN) = BB(IN) + FA(KN)*V
      GO TO 140

  130 IF (J.NE.L) GO TO 140
c     write(6,*)"i am here 6"
      IF (ITAG.eq.1) 
     .BB(KN) = BB(KN) - FA(IN)*V
      IF (ITAG.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
  140 IN = NIN(L) + I
c     write(6,*)"i am here 7"
      KN = NIN(J) + K
      IF (K.GT.J) KN = NIN(K) + J
      BB(IN) = BB(IN) - FA(KN)*V
      IF (IN.EQ.KN) GO TO 150
      BB(KN) = BB(KN) - FA(IN)*V
  150 IF (J.NE.K) GO TO 160
c     write(6,*)"i am here 8"
      IF (ITAG.eq.1)
     .BB(KN) = BB(KN) - FA(IN)*V
      IF (ITAG.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
  160 CONTINUE
      GO TO 20

  170 CONTINUE
c     call ga_sync()
c     write(6,*)"i am here??"
c     write(6,*)"BB",(BB(I),I=1,nrr)
c     write(6,*)"F",(F(I),I=1,nrr)
c     write(6,*)"FA",(FA(I),I=1,nrr)
C     REWIND 16
c     CALL ZREWIND (IO_SCF_MATRICES, 'FILLF')
      ISZ = 1
c     write(6,*)"ERROR"
c     IF(ISW1.EQ.8) THEN
C       ADD EXTRA TERMS FOR XINDO
cyw
c    this way, no diagonal elements are considered.
cyw
c       CALL ZRREAD (IO_SCF_MATRICES, Q, 1, NRR, 'FILLF')
c     ftr = ONE
c     DO 179 I = 2, NB
c        ID = NIN(I)
c        DO 175 J = 1,I-1
c           IF(NU(I).EQ.NU(J)) GO TO 175
c           IJ = ID + J
c           BB(IJ) = BB(IJ) + ftr*Q(IJ)*(SUMELEL(I) + SUMELEL(J))/TWO
c 175    CONTINUE
c 179 CONTINUE
C       CALL XEXCHANGE(Q,FA,BB,GAMMA,dryrun)
c     ELSE
C       READ (16)
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
c     ENDIF
c     READ (16)  SKIP OVERLAP
c     CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
C     READ (16)  SKIP S**(-1/2)
c     CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
C
C     HSCRF modifies Q.
c     if(lscrf) CALL MOMOUT(.false.,DRYRUN)
C
C     READ THE ONE ELECTRON MATRIX
C     READ (16) (Q(I),I=1,NRR)
c     CALL ZRREAD (IO_SCF_MATRICES, Q, 1, NRR, 'FILLF')
C
c     IF(JSW(6).EQ.50) THEN
C     MODIFY THE ONE ELECTRON MATRIX BY AMSOL
c        CALL GMSCRF(NB,NU,KAN,Q)
c     ELSEIF (LSCRF) THEN
C        THIS SECTION ADDS TO THE ONE ELECTRON MATRIX BB THE SELF
C        CONSISTENT REACTION FIELD.
c        CALL HSCRF(Q,DIP,NB,DRYRUN)
c     ENDIF
C        WRITE(6,*) ' FILLF AFTER ELSEIF'
  175 CONTINUE
 
      if (ITAG.eq.2) GO TO 190
      DO 189 I = 1, NRR
c      write(6,*)"BB,Q",BB(I),Q(I)
         BB(I) = BB(I) + Q(I)
  189 CONTINUE
      LL = 0
c     WRITE (6,*)"BB", (BB(M),M=1,NRR)
      GO TO 190

  180 WRITE (6,9000) I, J, K, L
      CALL ZEXIT(3)
CNR
  190 CONTINUE
         if (.not.ma_pop_stack(lQQ))
     & call errquit('semiemp_energy: cannot pop stack QQ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABL))
     & call errquit('semiemp_energy: cannot pop stack LABL',3, MA_ERR)
         if (.not.ma_pop_stack(lLABK))
     & call errquit('semiemp_energy: cannot pop stack LABK',3, MA_ERR)
         if (.not.ma_pop_stack(lLABJ))
     & call errquit('semiemp_energy: cannot pop stack LABJ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABI))
     & call errquit('semiemp_energy: cannot pop stack LABI',3, MA_ERR)


c     call ga_sync()
CNR
      RETURN

 9000 FORMAT (' **** ERROR IN INTEGRAL LIST',4I5,'****')

      END


