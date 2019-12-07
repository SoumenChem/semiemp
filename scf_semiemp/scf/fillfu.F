      SUBROUTINE FILLFU(BB,HB,FA,F,Q,GAMMA,SUMELEL,NG,QQ,NRR,LABI
     1 ,LABJ,LABK,LABL,NINTGL,NIN,ITAG)
C
C     ZERNER   UNIVERSITY OF FLORIDA.  COPYRIGHTED
C     UNRESTRICTED HF.
C     CLOSED SHELL F=TOTAL DENSITY, FA= HALF OF F=ALPHA SPIN DENSITY
C     OPEN SHELL F=TOTAL DENSITY, FA=ALPHA SPIN DENSITY
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
c     INCLUDE '../include/xindo.cmn' 
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/rhf.cmn'
c     INCLUDE '../include/char.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
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
      dimension  F(*),Q(*),FA(*),GAMMA(*),HB(*),BB(*),DIP(3)
     +,AD(3),BD(3),BBDD(3),BBD(3),NIN(*),SUMELEL(*),LABI(*),
     + LABJ(*),LABK(*),LABL(*),NG(*),QQ(*)
      CHARACTER*8  FLAGM, FLAGU
      
c     LSCRF=JSW(6).GT.0.AND.JSW(6).NE.50.AND.JSW(6).NE.40
c     IF (DRYRUN) THEN
c        IF (LSCRF) THEN 
CKA         CALL GRNDIP(DIP,AD,BD,BBBD,BBD,Q,FA,HB,RAD,ITDIP,DRYRUN)
c           IF(JSW(6).EQ.30) THEN
CGMP               CALL COSMO2(BB,DIP,NB,DRYRUN)  
CGMP     THE CALL ARGUMENTS DO NOT MATCH THAT FOR COSMO2 IN GSNDO1 OR COSMO2
CGMP     THERE ALSO DOES NOT APPEAR TO BE ANY REASON TO CALL COSMO2 HERE
CGMP     FOR THE DRYRUN, SINCE IT IS NOT CALLED FROM HERE IN THE REAL RUN
CGMP
c           ELSE
CKA            CALL HSCRF(BB,DIP,NB,DRYRUN)
c           ENDIF
c        END IF
c        RETURN
c     END IF
C     FOR SCRF YOU NEED MOMENTS EACH CYCLE
c     write(6,*)"i am here fillfu"
c     ITDIP = 1
c     IF (LSCRF) 
c    x       CALL GRNDIP(DIP,AD,BD,BBBD,BBD,Q,BB,HB,RAD,ITDIP,DRYRUN)
      
C
C     CONVENTION I.LE.J,K.LE.L,I.LE.K. FORI.EQ.K,J.LE.L
C     V=(I,J/K,L)=(1,1/2,2)
C
      FLAGU = FLAGU
C     REWIND 11
c     CALL ZREWIND (IO_AOINTS, 'FILLFU')
      ISZ = 1
c     CALL FIND(IO_AOINTS,FLAGM,'FILLFU')
C +++++++++
      DO 10 I = 1, NRR
         BB(I) = ZERO
         HB(I) = ZERO
   10 CONTINUE
      LL = 0
c     write(6,*)"F(I)",(F(I),I=1,nrr)
c     write(6,*)"FA(I)",(FA(I),I=1,nrr) 
C
C
c     Go to 170
   20 CONTINUE
C     write(6,*)"LL,NINTGL",LL,NINTGL
      LL = LL + 1
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
      I = LABI(LL)
      J = LABJ(LL)
      K = LABK(LL)
      L = LABL(LL)
      V = QQ(LL)
c     WRITE(6,*) 'I J K L V', I,J,K,L,V
C
C     SEE SUBROUTINE FILLF FOR REMARKS ABOUT ILAB AND EQUIVALENCE
C
C     WRITE(6,*) ' FILLFU ', I,J,K,L,v
      IF (I-J) 90, 40, 200
   40 IF (K-L) 80, 50, 200
   50 IF (I-K) 70, 60, 200
C     (I,I/I,I)
   60 IN = NIN(I) + I
      BB(IN) = BB(IN) + (F(IN)-FA(IN))*V
      HB(IN) = HB(IN) + FA(IN)*V
C     FOR NDDO THEORY
C     SUMELEL(I) = SUMELEL(I) + F(IN)*V
C     FOR INDO
      SUMELEL(I) = SUMELEL(I) + F(IN)*GAMMA(NIN(NG(I))+NG(I))
      GO TO 160
C     (I,I/K,K)
   70 IN = NIN(I) + I
      KN = NIN(K) + K
      BB(IN) = BB(IN) + F(KN)*V
      HB(IN) = HB(IN) + F(KN)*V
      BB(KN) = BB(KN) + F(IN)*V
      HB(KN) = HB(KN) + F(IN)*V
C     FOR NDDO
C     SUMELEL(I) = SUMELEL(I) + F(KN)*V
C     SUMELEL(K) = SUMELEL(K) + F(IN)*V
C     FOR INDO, CNDO
      IIKK = NIN(NG(K)) + NG(I)
      SUMELEL(I) = SUMELEL(I) + F(KN)*GAMMA(IIKK)
      SUMELEL(K) = SUMELEL(K) + F(IN)*GAMMA(IIKK)
      IN = NIN(K) + I
      BB(IN) = BB(IN) - FA(IN)*V
      HB(IN) = HB(IN) - (F(IN)-FA(IN))*V
      GO TO 160
C     (I,I/K,L)
   80 IN = NIN(I) + I
      KN = NIN(L) + K
      if (itag.eq.1)
     .BB(IN) = BB(IN) + 2.*F(KN)*V
      if (itag.eq.1)
     .HB(IN) = HB(IN) + 2.*F(KN)*V
      BB(KN) = BB(KN) + F(IN)*V
      HB(KN) = HB(KN) + F(IN)*V
      IN = NIN(K) + I
      KN = NIN(L) + I
      BB(IN) = BB(IN) - FA(KN)*V
      HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      BB(KN) = BB(KN) - FA(IN)*V
      HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      IF (I.NE.K) GO TO 160
      if (itag.eq.1)
     .BB(IN) = BB(IN) - FA(KN)*V
      if (itag.eq.1)
     .HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      if (itag.eq.2)
     .BB(IN) = BB(IN) + FA(KN)*V
      if (itag.eq.2)
     .HB(IN) = HB(IN) + (F(KN)-FA(KN))*V
      GO TO 160

   90 IF (K-L) 110, 100, 200
C     (I,J,/K,K)
  100 IN = NIN(J) + I
      KN = NIN(K) + K
      BB(IN) = BB(IN) + F(KN)*V
      HB(IN) = HB(IN) + F(KN)*V
      if (itag.eq.1)
     .BB(KN) = BB(KN) + 2.*F(IN)*V
      if (itag.eq.1)
     .HB(KN) = HB(KN) + 2.*F(IN)*V
      IN = NIN(K) + I
      KN = NIN(K) + J
      IF (K.LT.J) KN = NIN(J) + K
      BB(IN) = BB(IN) - FA(KN)*V
      HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      BB(KN) = BB(KN) - FA(IN)*V
      HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      IF (J.NE.K) GO TO 160
      if (itag.eq.1)
     .BB(KN) = BB(KN) - FA(IN)*V
      if (itag.eq.1)
     .HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      if (itag.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
      if (itag.eq.2)
     .HB(KN) = HB(KN) + (F(IN)-FA(IN))*V
      GO TO 160
C     (I,J/K,L)
  110 IN = NIN(J) + I
      KN = NIN(L) + K
      if (itag.eq.1)
     .BB(IN) = BB(IN) + 2.*F(KN)*V
      if (itag.eq.1)
     .HB(IN) = HB(IN) + 2.*F(KN)*V
      IF (IN.EQ.KN) GO TO 120
      if (itag.eq.1)
     .BB(KN) = BB(KN) + 2.*F(IN)*V
      if (itag.eq.1)
     .HB(KN) = HB(KN) + 2.*F(IN)*V
  120 IN = NIN(K) + I
      KN = NIN(L) + J
      IF (J.GT.L) KN = NIN(J) + L
      BB(IN) = BB(IN) - FA(KN)*V
      HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      BB(KN) = BB(KN) - FA(IN)*V
      HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      IF (I.NE.K) GO TO 130
      IF (J.EQ.L) GO TO 140
      if (itag.eq.1)
     .BB(IN) = BB(IN) - FA(KN)*V
      if (itag.eq.1)
     .HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      if (itag.eq.2)
     .BB(IN) = BB(IN) + FA(KN)*V
      if (itag.eq.2)
     .HB(IN) = HB(IN) + (F(KN)-FA(KN))*V
      GO TO 140

  130 IF (J.NE.L) GO TO 140
      if (itag.eq.1)
     .BB(KN) = BB(KN) - FA(IN)*V
      if (itag.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
      if (itag.eq.1)
     .HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      if (itag.eq.2)
     .HB(KN) = HB(KN) + (F(IN)-FA(IN))*V
CSG   BB(KN) = BB(KN) - FA(IN)*V
CSG   HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
  140 IN = NIN(L) + I
      KN = NIN(J) + K
      IF (K.GT.J) KN = NIN(K) + J
      BB(IN) = BB(IN) - FA(KN)*V
      HB(IN) = HB(IN) - (F(KN)-FA(KN))*V
      IF (IN.EQ.KN) GO TO 150
      BB(KN) = BB(KN) - FA(IN)*V
      HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
  150 IF (J.NE.K) GO TO 160
      if (itag.eq.1)
     .BB(KN) = BB(KN) - FA(IN)*V
      if (itag.eq.1)
     .HB(KN) = HB(KN) - (F(IN)-FA(IN))*V
      if (itag.eq.2)
     .BB(KN) = BB(KN) + FA(IN)*V
      if (itag.eq.2)
     .HB(KN) = HB(KN) + (F(IN)-FA(IN))*V
  160 CONTINUE
      GO TO 20
C
C     NOW ADD THE ONE ELECTRON MATRIX
  170 CONTINUE
C     REWIND 16
c     CALL ZREWIND (IO_SCF_MATRICES, 'FILLFU')
c     IF(ISW1.EQ.8) THEN
C       ADD EXTRA TERMS FOR XINDO
C       READ (16)  WEIGHTED OVERLAP
c       CALL ZRREAD (IO_SCF_MATRICES, Q, 1, NRR, 'FILLF')
c       DO 175 I = 2, NB
c       ID = NIN(I)
c       DO 172 J = 1,I-1
c       IF(NU(I).EQ.NU(J)) GO TO 175
c       IJ = ID + J
c       BB(IJ) = BB(IJ) + Q(IJ)*(SUMELEL(I) + SUMELEL(J))/TWO
c       HB(IJ) = HB(IJ) + Q(IJ)*(SUMELEL(I) + SUMELEL(J))/TWO
  172   CONTINUE
  175   CONTINUE
C       CALL XEXCHANGE(Q,FA,BB,GAMMA)
c     ELSE
C       READ (16)
c       CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
c     ENDIF


C     READ 16 - SKIP OVERLAP
c     CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
C     READ (16)  SKIP S**(-1/2)
c     CALL ZIREAD (IO_SCF_MATRICES, IDU, 1, ISZ, 'FILLF')
C
C     THIS SECTION ADDS TO THE ONE ELECTRON MATRIX BB THE SELF
C     CONSISTENT REACTION FIELD.
C     IF (LSCRF) 
C    x    CALL GRNDIP(DIP,AD,BD,BBBD,BBD,Q,FA,HB,RAD,ITDIP,DRYRUN)
C     READ THE ONE ELECTRON MATRIX
C      READ (16) (Q(I),I=1,NRR)
c     CALL ZRREAD (IO_SCF_MATRICES, Q, 1, NRR, 'FILLFU')
C     WRITE(6,*) ' ONE EL IN FILLFU'
C     CALL PRNT(Q,NB,NRR)
C     WRITE(6,*) ' TOTAL DENSITY IN FULLFU'
C     CALL PRNT(F,NB,NRR)
C     WRITE(6,*) ' ALPHA DENSITY IN FULLFU'
C     CALL PRNT(FA,NB,NRR)
CJiabo
c     IF(JSW(6).EQ.50) CALL GMSCRF(NB,NU,KAN,Q)
CLi
c     IF (LSCRF) CALL HSCRF(Q,DIP,NB,DRYRUN)
      if (itag.eq.2) GO TO 210
      DO I = 1, NRR
c     write(6,*)"BB,HB,Q",BB(I),HB(I),Q(I)
         BB(I) = BB(I) + Q(I)
c     write(6,*)"HB",HB(I)
         HB(I) = HB(I) + Q(I)
      END DO
C-----------------------------------------------------------------------
C
C     SPIN CONSTRAINED UHF
C
C     SEE [N. HANDY, CHEM. PHYS. LETT.  VOL. 183, NO. 5, P. 423 (1991)]
C
C-----------------------------------------------------------------------
c     IF (SCFTYP.EQ.'SUHF  ') THEN
C
c     NOPEN = MULT - 1
c     NCLOSED = ( NEL - NOPEN ) / 2

c        CALL GETIT('VECSCFB',SVEC ,NR,1)

c        DO 300 I=1,NRR
c300        FZ(I)=0.0D0
           
c     DO 310 I=1,NOPEN
c       IMO1=NCLOSED+I
c310    CALL CISD(FZ,SVEC ,1,IMO1,IMO1,1.0D0,1.0D0,1.0D0)

c        DO 180 I = 1, NRR
c           HB(I) = HB(I) - AMBDA*( 2.0D0 * FA(I) - F(I) - FZ(I) )
c180        CONTINUE

c        CALL GETIT('VECSCFA',SVEC ,NR,1)

c        DO 320 I=1,NRR
c320       FZ(I)=0.0D0

c        DO 330 I=1,NOPEN
c          IMO1=NCLOSED+I
c330       CALL CISD(FZ,SVEC ,1,IMO1,IMO1,1.0D0,1.0D0,1.0D0)

c        DO 190 I = 1, NRR
c           BB(I) = BB(I) - AMBDA*( F(I) - 2.0D0*FA(I) + FZ(I) )
c190        CONTINUE

c     END IF

C
C      WRITE (16) (HB(I),I=1,NRR)
c     CALL ZRWRITE (IO_SCF_MATRICES, HB, 1, NRR, 'FILLFU')
C      WRITE (16) (BB(I),I=1,NRR)
c     CALL ZRWRITE (IO_SCF_MATRICES, BB, 1, NRR, 'FILLFU')
C     ******UNIT 16  *******
C       (OVERLAP),(WEIGHTED OVERLAP),(OVERAP)**(-1/2),
C       (ONE ELECTRON MATRIX),
C       (HB=BETA FOCK MATRIX),(BB=ALPHA FOCK MATRIX)
C       WRITE(6,*) ' FOCK MATRICES IN FILLFU'
C       CALL PRNT(HB,NB,NRR)
C       CALL PRNT(BB,NB,NRR)
      GO TO 210

  200 WRITE (6,9000) I, J, K, L
      CALL ZEXIT(3)
  210 CONTINUE
c     write(6,*)"after fock F",(BB(I),I=1,nrr)
      call ga_sync()
      RETURN

 9000 FORMAT (' **** ERROR IN INTEGRAL LIST',4I5,'****')

      END
