      SUBROUTINE DIPOLE(rtdb,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,NB,
     1 NRR,ISW1,Z12,D12,P,NEL,D,FX,FY,FZ,ITAG,PME,IHOMO,ILUMO,IORB)

C
C     (17-JUN-1983) /MCZ/
C
C     THIS SUBROUTINE CALC. THE DIPOLE MATRIX ELEMENTS AND STORES
C     THEM ON UNIT 10, AND CALC. THE GRD. STATE DIPOLE.
C
C     CODE IS FOR ZDO BASIS. TO GENERALIZE ACTIVATE ALL 'CT' CARDS
C
C     RAD ARE THE INTERATOMIC DISTANCES IN ANGSTROMS, FG ARE THE
C     EULERIAN TRANSFORMATION COEF. NEEDED FOR SUB MOMENT.
C
C     JDIP=1   CALCULATE THE INTEGRALS
C     JDIP=2   THE INTEGRALS HAVE BEEN STORED, CALC THE DIPOLE
C              MOMENT OF A GIVEN STATE PRINTING ALL THE COMPONENTS
C              THIS TAKES CONSIDERABLE STORAGE.
C     JDIP=3   THE INTEGRALS HAVE BEEN STORED. CALC ONLY THE TOTAL
C              DIPOLE U(X), U(Y) AND U(Z).
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
c     INCLUDE '../include/const.cmn'
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
      CHARACTER*8  LABEL
C
C****************************************************************
c     INCLUDE '../include/switch.cmn'
C****************************************************************
C***      COMMON/SWITCH/
C***     .  ISWANN,III,IPRINT,IPUN,ISAVE,ISPEC,ISW1,ISW2,ISW3,ISW4,
C***     .  JSW(20),IBETA,IDD1,IDD2
C****************************************************************
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     .  ATOM(100),NIN(MXB),
C***     .  NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     .  NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     .  Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
c     INCLUDE '../include/dipgrn.cmn'
c     INCLUDE '../include/summary.cmn'
C****************************************************************
C     COMMON/DIPGRN/
C    . DIPGR(3),DPM(3),POLL,DIPNUC(3)
C****************************************************************
C
c     INCLUDE '../include/debug.cmn'
c     INCLUDE '../include/outcntrl.cmn'
      integer rtdb
      integer orbspace(2), aorb(3) 
       character*8 ATOM(105)
      DIMENSION    H(*),  RAD(NA,NA), A(3),
     .             BB(3), E(3), BBB(3),NLL(16),NIN(*),NU(*),
     .             B(3),NW(*),CO(3,NA),NP(*),ALPHA(NB,6),KAN(*),DIPP(6)
     .             ,DPM(3),D12(*),Z12(*),P(*),D(3),FX(*),FY(*)
     .             ,FZ(*), PME(3*iorb)
      DIMENSION    LABEL(3), MPG(3)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
      PARAMETER (DEBYE    = 2.541748047D0,BOHRS    = 5.29177249D-1)
      LOGICAL      DRYRUN
      SAVE         SQRT2, SQRT3, SQRT5, SQRT10, SRT2D5, SRT6D5, SRT3D2,
     .             TEST, LABEL, NPROP
C---
      DATA         SQRT2/1.414213562D0/, SQRT3/1.732050808D0/,
     .             SQRT5/2.236067977D0/, SQRT10/3.16227766D0/,
     .             SRT2D5/0.632455532D0/, SRT6D5/1.095445115D0/,
     .             SRT3D2/1.224744871D0/
      DATA         NLL/0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3/
      DATA         TEST/1.0D-4/
      double precision FACT(35)
      integer ehinfo(5)
      logical got_orbspace, got_elh
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
      DATA ATOM/' H','Ct','Li','Be',' B',' C',' N',' O',' F','Ne','Na',
     1 'Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',' V','Cr',
     2 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',
     3 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
     4 'Sn','Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm',
     5 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re',
     6 'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     7 'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     8 'Md','No','Lr','XX','XX'/
C---
C DEFINITIONS FOR WRTPRP (WRITE PROPERTIES)
C 22-JUN-1983 /HJAAJ
C
      DATA         LABEL/'       X', '       Y', '       Z'/
C---
C     IDIP=1 FOR CNDO APPROX. OF POPLE TO BE USED, IDIP=2 FOR THE
C        DIPOLE TO BE CALC. OVER STO'S AND TRANSFORMED TO A SYM.
C        ORTHOG. BASIS, IDIP=3 FOR ONLY CHARGE TERM TO BE USED.
C        IDIP = 4 FOR EHT OR AB-INITIO
C---
C
      DX = ZERO
      DY = ZERO
      DZ = ZERO 
             if (.not.MA_Push_Get(MT_Dbl,3*na*na,'FG',lFG,iFG))
     &   call errquit('dipole: cannot allocate FG',0,
     1    MA_ERR)

      got_orbspace = .false.
      got_elh = .false.

      if (itag.eq.2) then
      if (.not. rtdb_get (rtdb, "rt_semiemp:got_orbspace", mt_log, 1,
     1 got_orbspace))        
     1 call errquit ('dipole:failed to get got_orbspace target from 
     1 rtdb',0,RTDB_ERR)

      if (.not. rtdb_get (rtdb, "rt_semiemp:got_elh", mt_log, 1,
     1 got_elh))
     1 call errquit ('dipole:failed to get got_elh target from
     1 rtdb',0,RTDB_ERR)

       end if

       if (got_elh) then
        if (.not.rtdb_get(rtdb,'rtsemiemp:ehinfo',mt_int,5, ehinfo))
     1 call errquit('eh_state: failed to get ehinfo',0,RTDB_ERR)
       end if

       if (got_orbspace) then
c       write(6,*)"here"
        if (.not.rtdb_get(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1  call errquit('dipole: failed to read orbspace',0,
     2   RTDB_ERR)
c       if (.not.rtdb_get(rtdb,'rt_semiemp:actsize',mt_int,3,aorb))
c    1  call errquit('dipole: failed to read aorb',0,
c    2   RTDB_ERR)
c       if (.not.MA_Push_Get(MT_Dbl,3*aorb(1)*aorb(2),'PME',lPME,iPME))
c    &   call errquit('dipole: cannot allocate FG',0,
c    1    MA_ERR)
c       if (.not.MA_Push_Get(MT_Int,aorb(1),'IHOMO',lhomo,ihomo))
c    &   call errquit('dipole: cannot allocate IHOMO',0,
c    1    MA_ERR)
c       if (.not.MA_Push_Get(MT_Int,aorb(2),'ILUMO',llumo,ilumo))
c    &   call errquit('dipole: cannot allocate ILUMO',0,
c    1    MA_ERR)

c       write(6,*)"in dipole",orbspace(1),orbspace(2)
      end if 

c     write(6,*)"i am in dipole"!,(P(I),I=1,nrr)
      JDIP = IONE
C     NPROP=IO_MOMENTS
      IDIP = IONE
c     IDIP = JSW(10)
   10 IF (IDIP.LT.10) GO TO 20
      IDIP = IDIP - 10
      GO TO 10

   20 CONTINUE
c     IF (III.NE.0 .AND. IDIP.EQ.2) THEN
c        ISCRF = JSW(6)
c        IF (ISCRF.GT.0) THEN
c           WRITE (6,9000)
c           IDIP = 1
c        END IF
c     END IF
      IF(ISW1.GT.50) THEN
       IDIP = IDIP
      ELSE
        IF (ISW1.EQ.0 .OR. ISW1.GT.7) IDIP = 4
      ENDIF
      DIPFAC = ONE
      IF (IDIP.EQ.3) DIPFAC = ZERO
C---
      CALL CMS(C1,C2,C3,CO,KAN,NA)
C     DO 22 I = 1,NA
C     CO(1,I) = CO(1,I)-C1
C     CO(2,I) = CO(2,I)-C2
C     CO(3,I) = CO(3,I)-C3
C  22 CONTINUE
      FF1 = (DEBYE/BOHRS)
      NA1 = NA*NA
      NA2 = 2*NA1
      IF (ISW1.EQ.0.OR.ISW1.EQ.8.OR.ISW1.EQ.9) GO TO 40
      IF (JDIP.EQ.2) THEN
C     FIRST THE POPULATION IN THE ZDO BASIS
         WRITE (6,9100)
         WRITE (6,9110)
         FCHG = ZERO
         DO 30 I = 1, NA
            FF = Z12(I) - D12(I)
            KK = KAN(I)
            IF (KK.EQ.0) KK = 2
            CX = CO(1,I)*BOHRS
            CY = CO(2,I)*BOHRS
            CZ = CO(3,I)*BOHRS
            WRITE (6,9120) I, ATOM(KK), FF, D12(I), CX, CY, CZ
            FCHG = FCHG + FF
   30    CONTINUE
         GO TO 660
      END IF
      IF (JDIP.GT.1) GO TO 660
C      IF(III.EQ.0) GO TO 150
   40 CONTINUE
C      REGENERATE GEOMETRIC MATRICES IF NOT ALREADY THERE.
C
C
      NA1 = NA*NA
      THRESH = 1.0D-7
      NA2 = 2*NA1
      DO 100 I = 1, NA
         DO 90 J = 1, NA
            IF (I.GT.J) IJ = J + (I* (I-1))/2
            IF (I.LE.J) IJ = I + (J* (J-1))/2
            ITJ = I + (J-1)*NA
            IF (ABS(RAD(I,J)).LT.1.0D-4) GO TO 60
            IF (I-J) 50, 90, 50
   50       Dbl_MB(iFG+ITJ-1) = (CO(3,J)-CO(3,I))/RAD(I,J)
            G = ONE - Dbl_MB(iFG+ITJ-1)*Dbl_MB(iFG+ITJ-1)
            IF (THRESH-G) 80, 70, 70
   60       Dbl_MB(iFG+ITJ-1) = ONE
   70       Dbl_MB(iFG+ITJ+NA1-1) = ONE
C
C  CHANGE ADDED BY JDHD JAN 18 83
C
            IF (Dbl_MB(iFG+ITJ-1).LT.ZERO) Dbl_MB(iFG+ITJ+NA1-1) = -ONE
            Dbl_MB(iFG+ITJ+NA2-1) = ZERO
            GO TO 90

   80       G = ONE/SQRT(G)
            Dbl_MB(iFG+ITJ+NA1-1) = (CO(1,J)-CO(1,I))*G/RAD(I,J)
            Dbl_MB(iFG+ITJ+NA2-1) = (CO(2,J)-CO(2,I))*G/RAD(I,J)
   90    CONTINUE
  100 CONTINUE
  110 CONTINUE
      if (itag.eq.1) then
CTDIP    2
C     RETRIEVE THE ORBITAL OVERLAPS FOR SUB. MOMENT
CTDIP    2
C ***** FLAG
C      REWIND 16
      DO 620 I = 1, NB
         K = NW(I)
         KA = NU(I)
         KX = KAN(KA)
         NT1 = NP(I)
         II = NIN(I) + I
C     ONE CENTRE DIAGONAL TERMS
        FX(II) = -(CO(1,KA)*BOHRS-C1)*FF1
        FY(II) = -(CO(2,KA)*BOHRS-C2)*FF1
        FZ(II) = -(CO(3,KA)*BOHRS-C3)*FF1
c      write(6,*)"FX,FY,FZ",FX(II),FY(II),FZ(II)
         NL = NLL(K+1)
  120    CONTINUE
         LL = I + 1
         IF (NB.LT.LL) GO TO 620
         CONTINUE
         DO 610 J = LL, NB
            IJ = I + NIN(J)
            M = NW(J)
            MA = NU(J)
            MX = KAN(MA)
            FX(IJ) = ZERO
            FY(IJ) = ZERO
            FZ(IJ) = ZERO
            NT2 = NP(J)
            SS = ZERO
            MM = NLL(M+1) + 1
C *****FLAG
            IF (MA.NE.KA) GO TO 560
C     ONE CENTRE OFF DIAGONAL TERMS
            IF (MM.NE. (NL+2)) GO TO 600
            MMM = NIN(MM) + NL + 1
            DO 180 ISS = 1, 2
               AMU = ALPHA(I,ISS+1)
               IF (AMU.LT.TEST) GO TO 180
               ACU = ALPHA(I,ISS+4)
               DO 170 JSS = 1, 2
                  BMU = ALPHA(J,JSS+1)
                  IF (BMU.LT.TEST) GO TO 170
                  BCU = ALPHA(J,JSS+4)
                  ABC = ACU*BCU
                  GO TO (600,130,600,600,140,600,600,
     .                   600,150,600), MMM
C     SP POLARIZATION TERM
  130             J2 = 2*NT2 + 1
                  K2 = 2*NT1 + 1
                  T = SQRT((TWO*AMU)**K2* (TWO*BMU)**J2/
     .                (FACT(K2)*FACT(J2)*THREE))
                  K3 = NT1 + NT2 + 2
                  SSS = T*FACT(K3)/ (AMU+BMU)**K3
                  GO TO 160
C     PD POLARIZATION TERM
  140             J3 = 2*NT2 + 1
                  K3 = NT1 + NT2 + 2
                  J2 = 2*NT1 + 1
C                 T=SQRT((TWO*AMU)**J2*(TWO*BMU)**J3/(FACT(J2)*FACT(J3)
C    X               *FIVE))
                  T = SQRT((TWO*AMU)**J2* (TWO*BMU)**J3/
     .                (FACT(J2)*FACT(J3)*FIVE))
                  SSS = T*FACT(K3)/ (AMU+BMU)**K3
                  GO TO 160
C     DF POLARIZATION
  150             CONTINUE
                  J4 = 2*NT2 + 1
                  J2 = 2*NT1 + 1
                  K3 = NT1 + NT2 + 2
                  T = SQRT((TWO*AMU)**J2* (TWO*BMU)**J4/
     .                (FACT(J2)*FACT(J4)*SEVEN))
                  SSS = T*FACT(K3)/ (AMU+BMU)**K3
  160             CONTINUE
                  SS = SS + SSS*ABC
  170          CONTINUE
  180       CONTINUE
            IF (MMM.NE.2) GO TO 220
C      SP SECTION
            GO TO (190,200,210), M

  190       FX(IJ) = -DEBYE*SS
            GO TO 600

  200       FY(IJ) = -DEBYE*SS
c           write(6,*)"here 2"
            GO TO 600

  210       FZ(IJ) = -DEBYE*SS
c           write(6,*)"here 3"
            GO TO 600

  220       IF (MMM.EQ.9) GO TO 500
C     PD SECTION
            MMMM = K + 3* (M-4)
            GO TO (230,240,250,260,270,600,280,260,
     .             600,290,600,260,600,290,280), MMMM

  230       FX(IJ) = DEBYE*SS/SQRT3
            GO TO 600

  240       FY(IJ) = DEBYE*SS/SQRT3
            GO TO 600

  250       FZ(IJ) = -DEBYE*SS*TWO/SQRT3
            GO TO 600

  260       FX(IJ) = -DEBYE*SS
            GO TO 600

  270       FY(IJ) = DEBYE*SS
            GO TO 600

  280       FY(IJ) = -DEBYE*SS
            GO TO 600

  290       FZ(IJ) = -DEBYE*SS
            GO TO 600

C 300       FZ(IJ) = DEBYE*SS*THREE/SQRT5
C           GO TO 600

  310       FX(IJ) = -DEBYE*SS*SRT6D5
            GO TO 600

  320       FX(IJ) = DEBYE*SS*SRT6D5
            GO TO 600

  330       FY(IJ) = -DEBYE*SS*SRT6D5
            GO TO 600

  340       FY(IJ) = DEBYE*SS*SRT6D5
            GO TO 600

C 350       FZ(IJ) = DEBYE*SS*TWO*SRT2D5
C           GO TO 600

  360       FX(IJ) = DEBYE*SS
            GO TO 600

  370       FY(IJ) = -DEBYE*SS/SQRT10
            GO TO 600

C 380       FX(IJ) = -DEBYE*SS/SQRT10
C           GO TO 600

  390       FZ(IJ) = DEBYE*SS
            GO TO 600

  400       FY(IJ) = -DEBYE*SS*SRT3D2
            GO TO 600

  410       FX(IJ) = DEBYE*SS*SRT3D2
            GO TO 600

C 420       FY(IJ) = DEBYE*SS/SQRT10
C           GO TO 600

  430       FY(IJ) = DEBYE*SS*SRT3D2
            GO TO 600

C 440       FX(IJ) = -DEBYE*SS*SRT2D5
C           GO TO 600

C 450       FX(IJ) = DEBYE*SS*SRT2D5/SQRT3
C           GO TO 600

  460       FZ(IJ) = DEBYE*SS*FOUR/SQRT5
            GO TO 600

  470       FX(IJ) = -DEBYE*SS/SQRT10
            GO TO 600

  480       FY(IJ) = DEBYE*SS/SQRT10
            GO TO 600

  490       FZ(IJ) = DEBYE*SS*TWO*SQRT2
            GO TO 600

  500       CONTINUE
C      DF SECTION
            MMMMD = K - 3
            MMMMF = M - 8
            GO TO (510,520,530,540,550), MMMMD

  510       GO TO (460,320,340,600,600,600,600), MMMMF

  520       GO TO (600,470,480,390,600,410,430), MMMMF

  530       GO TO (600,370,470,600,390,400,410), MMMMF

  540       GO TO (310,490,600,360,270,600,600), MMMMF

  550       GO TO (330,600,490,280,360,600,600), MMMMF

  560       CONTINUE
            GO TO (600,570,600,570), IDIP

  570       CONTINUE
C     THIS IS THE TWO CENTRE BOND TERM
C
CTDIP   26
C ***** FLAG
            KAMA = KA + (MA* (MA-1))/2
            R = RAD(KA,MA)
            NT2 = NP(J)
            DDZ = ZERO
            DDY = ZERO
            DDX = ZERO
            S = H(IJ)
            DO 590 ISS = 1, 2
               AMU = ALPHA(I,ISS+1)
               IF (AMU.LT.TEST) GO TO 590
               ACU = ALPHA(I,ISS+4)
               DO 580 JSS = 1, 2
                  BMU = ALPHA(J,JSS+1)
                  IF (BMU.LT.TEST) GO TO 580
                  BCU = ALPHA(J,JSS+4)
                  MMM1 = MM - 1
c                 CALL MOMENT(I,NT1,NL,AMU,J,NT2,MMM1,BMU,R,S,DX,DY,DZ,
c    .                        FACT,IDIP,NU,NW,CO,FG,C1,C2,C3,NA)
                  ABC = ACU*BCU
                  DDX = DDX + DX*ABC
                  DDY = DDY + DY*ABC
                  DDZ = DDZ + DZ*ABC
  580          CONTINUE
  590       CONTINUE
             FX(IJ) = -DEBYE*DDX*DIPFAC
             FY(IJ) = -DEBYE*DDY*DIPFAC
             FZ(IJ) = -DEBYE*DDZ*DIPFAC
C ***** FLAG
  600       CONTINUE
  610    CONTINUE
  620 CONTINUE
      end if
C     NOW STORE THE DIPOLE INTEGRALS ON UNIT 10
C---
C WRITE DIPOLE INTEGRALS ON UNIT NPROP IN THE FORMAT OF THE
C POLYATOM PROPERTIES PROGRAM.
C
      GO TO (640,630,640,640), IDIP

  630 CONTINUE
C     READ S**(-1/2)
C      READ(16) (H(I),I=1,NRR)
c     CALL ZRREAD (IO_SCF_MATRICES, H, 1, NRR, 'DIPOLE')
C     FMCK=1.0D0
C     DO 1750 J=2,NB
C     JL=J-1
C     DO 1748 K=1,JL
C     JK=NIN(J)+K
C     H(JK)=H(JK)*FMCK
C1748 CONTINUE
C1750 CONTINUE
c     IF(INDEX(CONTRLO,'DIPOLE').NE.0) THEN
c       WRITE(6,*) ' DIPOLE MOMENT COMPONENTS '
c       WRITE(6,*) '      X'
c       CALL PRNT(FX,NB,NRR)
c       WRITE(6,*) '      Y'
c       CALL PRNT(FY,NB,NRR)
c       WRITE(6,*) '      Z'
c       CALL PRNT(FZ,NB,NRR)
c     ENDIF
c     CALL FOCKTR(FX,H,NB,DRYRUN)
c     CALL FOCKTR(FY,H,NB,DRYRUN)
c     CALL FOCKTR(FZ,H,NB,DRYRUN)
  640 CONTINUE
C
C     INTEGRALS ARE STORED IN DEBYES,
C
C     REWIND NPROP
c     CALL ZREWIND (IO_MOMENTS, 'DIPOLE')
      IF (IDIP.EQ.3) THEN
         IDD = 1
         IDDD = 1
         DO 650 II = 1, NRR
            IF (II.EQ.IDD) THEN
               IDDD = IDDD + 1
               IDD = NIN(IDDD+1)
            ELSE
                FX(IJ) = 0.0D0
                FY(IJ) = 0.0D0
                FZ(IJ) = 0.0D0
            END IF
  650    CONTINUE
      END IF
c     CALL WRTPRP(NPROP,LABEL(1),FX,NRR,.FALSE.,'DIPOLE')
c     CALL WRTPRP(NPROP,LABEL(2),FY,NRR,.FALSE.,'DIPOLE')
c     CALL WRTPRP(NPROP,LABEL(3),FZ,NRR,.FALSE.,'DIPOLE')
c     IF(INDEX(CONTRLO,'DIPOLE').NE.0) THEN
c       WRITE(6,*) ' DIPOLE MOMENT COMPONENTS '
c       WRITE(6,*) '      X'
c       CALL PRNT(FX,NB,NRR)
c       WRITE(6,*) '      Y'
c       CALL PRNT(FY,NB,NRR)
c       WRITE(6,*) '      Z'
c       CALL PRNT(FZ,NB,NRR)
c     ENDIF
C---
C     NOW CALCULATE THE GROUND STATE DIPOLE
C
  660 CONTINUE
CTEST
      ITDIP = 0
      JDIP = 2
      IDIP = 1
      IF (JDIP.EQ.3) ITDIP = 1
      IF (JDIP.EQ.1) GO TO 710
c     IHOMO = NEL/2
c     ILUMO = NEL/2+1
c     write(6,*)"HOMO,LUMO",IHOMO,ILUMO
C     IF(ABS(FCHG).GT.1.0D-7.AND.JSW(6).EQ.0) GO TO 3300
C     FIRST TIME THROUGH THE DIPOLE INTEGRALS ARE STORED, BUT THERE
C     IS NOT YET A DENSITY MATRIX FOR CALCULATING A DIPOLE MOMENT.
C
C      - TT May 1995 -
C     If it later seems, that GRNDIP is better off locally allocating
C     FX, FY, FZ, we may as well free them up right here, not up in the
c     wrapper, but the wrapper needs to call the free's in the dryrun
c     case still.
C  
      CALL GRNDIP(P,E,A,B,BBB,BB,FX,FY,FZ,RAD,ITDIP,NA,NU,NRR,ISW1,NIN,
     1 CO,Z12,KAN,NB,NW,PME,orbspace,got_orbspace,got_elh,IHOMO,ILUMO,
     1 IORB,NEL,ehinfo)

c      call ga_sync()
C
c     write(6,*)"after grndip"
      IF (JDIP.NE.2) GO TO 710
C

      FNACH = ZERO
      DO 670 I = 1, NA
         FNACH = FNACH + Z12(I)
  670 CONTINUE
      FNECH = FNACH - DFLOAT(NEL)
      IF ((ABS(FNECH).GT.1.0D-6).AND.(ga_nodeid().eq.0)) THEN
         WRITE (6,9040)
      END IF
      IDIP = IONE
      IF (ga_nodeid().eq.0.and.itag.eq.1) WRITE (6,9050) IONE, IONE
      IF (IDIP.EQ.1.and.ga_nodeid().eq.0.and.itag.eq.1) THEN
         WRITE (6,9010) (E(J),J=1,3)
         WRITE (6,9020) (A(J),J=1,3)
         WRITE (6,9030) (B(J),J=1,3)
         WRITE (6,9140) (BBB(J),J=1,3)
         WRITE (6,9130) (BB(J),J=1,3)
      END IF
      DO 690 J = 1, 3
c     write(6,*)"E,A,B,BB,BBB"!,E(J),A(J),B(J),BB(J),BBB(J)
         E(J) = E(J) + A(J) + B(J) + BB(J) + BBB(J)
         DIPP(J) = E(J)
  690 CONTINUE
      NECH = FNECH + 0.00001
      if (ga_nodeid().eq.0.and.itag.eq.1)
     & WRITE (6,9060) (E(J),J=1,3)
  700 CONTINUE
      AAA = E(1)**2 + E(2)**2 + E(3)**2
      DPM(1) = E(1)/DEBYE
      DPM(2) = E(2)/DEBYE
      DPM(3) = E(3)/DEBYE
c     write(6,*)"E",E(1),E(2),E(3)
       D(1) = E(1)*0.393430307
       D(2) = E(2)*0.393430307
       D(3) = E(3)*0.393430307
      DTOT = SQRT(AAA)
      IF (NECH.EQ.0.and.ga_nodeid().eq.0.and.itag.eq.1) THEN
         WRITE (6,9090) DTOT
      END IF
c     for summary...
      dipole_moment=dtot
  710 CONTINUE
C     RESTORE THE ABSOLUTE COORDINATES
C     DO 4422 I = 1,NA
C     CO(1,I) = CO(1,I)+C1
C     CO(2,I) = CO(2,I)+C2
C     CO(3,I) = CO(3,I)+C3
C4422 CONTINUE
 9000 FORMAT (/,'  ***** FOR GEOMETRY OPTIMISATION AND REACTION',
     .       ' FIELD, ONLY ONE CENTER DIPOLE APPROXIMATION IS USED ****'
     .       ,/)
 9010 FORMAT (' CHARGE DENSITY CONTRIBUTION ',3F12.4)
 9020 FORMAT (' SP LOCAL ATOMIC POLARIZATION',3F12.4)
 9030 FORMAT (' PD LOCAL ATOMIC POLARIZATION',3F12.4)
 9040 FORMAT (11X,'  > CHARGED SPECIES < ')
 9050 FORMAT (1X,2I3,' STATE DIPOLE MOMENT',10X,'X',11X,'Y',11X,'Z')
 9060 FORMAT (1X,'TOTAL DIPOLE MOMENT (DEBYES)',3F12.4)
 9090 FORMAT (' TOTAL DIPOLE=',F14.6,' DEBYES'//)
 9100 FORMAT (9X,'SUMMARY OF POPULATION STUDY IN THE ZDO BASIS',//)
 9110 FORMAT (9X,'ATOM',3X,'TYPE',4X,'FORMAL',5X,'TOTAL',14X,
     .       'COORDINAT','ES'/24X,'CHARGE',5X,'CHARGE',8X,'X',10X,'Y',
     .       10X,'Z'/)
 9120 FORMAT (10X,I3,5X,A2,5F11.3)
 9130 FORMAT (1X,'TWO CENTRE BOND CONTRIBUTION',3F12.4)
 9140 FORMAT (1X,'DF LOCAL ATOMIC POLARIZATION',3F12.4)

c      WRITE (6,*) 'DIPOLE out'
c      write(6,*)"ILUMO",(ILUMO(I),I=1,IORB)
c      write(6,*)"IHOMO",(IHOMO(I),I=1,IORB)
c      if (got_orbspace) then
c       if (.not.ma_pop_stack(lLUMO))
c    &   call errquit('dipole: cannot pop stack lPME',3, MA_ERR)
c       if (.not.ma_pop_stack(lHOMO))
c    &   call errquit('dipole: cannot pop stack lPME',3, MA_ERR)
c       if (.not.ma_pop_stack(lPME))
c    &   call errquit('dipole: cannot pop stack lPME',3, MA_ERR)
c      end if
      if (.not.ma_pop_stack(lFG))
     &   call errquit('dipole: cannot pop stack lFG',3, MA_ERR)

      RETURN

      END
C------------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------------
      SUBROUTINE CMS(X,Y,Z,CO,KAN,NA)
C
C     calculates coordinates of mass centre
C
c     IMPLICIT REAL*8(A-H,O-Z)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/z.cmn'
c     INCLUDE '../include/mass.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION KAN(*),CO(3,NA), FMASS(105)
      PARAMETER (ZERO = 0.0d0,BOHRS    = 5.29177249D-1)
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

c     write(6,*)"i am in cms"
      X = ZERO
      Y = ZERO
      Z = ZERO
      TM = ZERO
      DO 10 I = 1, NA
         IF (KAN(I).NE.0) THEN
            B = FMASS(KAN(I))
            TM = TM + B
            X = X + CO(1,I)*B*BOHRS
            Y = Y + CO(2,I)*B*BOHRS
            Z = Z + CO(3,I)*B*BOHRS
         END IF
   10 CONTINUE
      X = X/TM
      Y = Y/TM
      Z = Z/TM
      RETURN

      END

