      SUBROUTINE PURVIS(F,NB,RAD,NA,KAN,NTYP,CO,NCORE,ISW1)

C
C                     (C) COPYRIGHT 1986
C
C                     QUANTUM THEORY PROJECT
C                      UNIVERSITY OF FLORIDA
C                      GAINESVILLE, FL 32611
C
C                        ALL RIGHTS RESERVED
C
C    ****************************************************************
C    *                                                              *
C    *                      N  O  T  I  C  E                        *
C    *                                                              *
C    *       THE SET OF PROGRAMS BUILT BY THIS PROCEDURE            *
C    *       CONTAINS PROPRIETARY INFORMATION OF                    *
C    *       THE AUTHORS AND THE QUANTUM THEORY PROJECT.            *
C    *       THE INFORMATION HEREIN                                 *
C    *       SHALL NOT BE DISCLOSED TO OTHERS, AND SHALL NOT BE     *
C    *       REPRODUCED IN WHOLE OR IN PART, WITHOUT WRITTEN        *
C    *       PERMISSION.                                            *
C    *                                                              *
C    *        QTP AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANT*
C    *        WHATSOEVER IN RESPECT TO THIS PROGRAM, INCLUDING THE  *
C    *        IMPLIED WARRANTY OF MERCHANTABILITY OR FITNESS FOR    *
C    *        A PARTICULAR PURPOSE; AND THE USER ASSUMES ALL RISKS, *
C    *        INCLUDING CONSEQUENTIAL LOSS OR DAMAGE,IN RESPECT TO  *
C    *        USE OF THIS PROGRAM.  IN ADDITION, QTP AND THE AUTHORS*
C    *        SHALL NOT BE OBLIGATED TO CORRECT OR MAINTAIN THE     *
C    *        PROGRAM, OR NOTIFY USER OF MODIFICATIONS OR UPDATES   *
C    *        THAT IT MAY MAKE FROM TIME TO TIME.                   *
C    *                                                              *
C    ****************************************************************
C
C
C     CALCULATES INNER SHELL M.O.'S NEEDED FOR CERTAIN PROPERTIES.
C     ASSUMES SLATER RULES FOR EXPONENTS, AND DIAGONALIZES CORE
C     MATRIX ASSUMED NON INTERACTING WITH RESPECT TO THE VALENCE
C     SHELL

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
c     Integer rtdb
C     WARNING   RPA PGM MUST HAVE THESE DIMENSIONS
      integer iNPC
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
     +   ,ITEN = 10, IPRINT  = 0)
c     DIMENSION  BB(*), F(*), H(*), FZ(*), 
      DIMENSION F(*), RAD(NA,NA), NTYP(*), KAN(*), CO(3,NA) 
C
C
      DATA         THRESH/1.0D-7/
      double precision FACT(35)
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

C
C     THE EXPONENTS ARE NOW BEING FILLED IN. ALPC(I)
C
C     EXCEPT WHERE NOTED THESE ARE SLATER EXPONENTS.
C     DO NOT CALL DEFCOR IF HCORE HAS ALREADY GENERATED THE CORE
C     BASIS
C     SAVE NFIRST 
c     write(6,*)"i am in purvis",ncore
      if (NB.GE.NCORE) NCORE = NB
            if (.not.MA_Push_Get(MT_int,ncore,'NFIRST',lnfirst,infirst))
     &   call errquit('purvis: cannot allocate nfirst',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_int,ncore,'NUC',lnuc,inuc))
     &   call errquit('purvis: cannot allocate nuc',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_int,ncore,'NWC',lnwc,inwc))
     &   call errquit('purvis: cannot allocate nwc',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_int,ncore,'NPC',lnpc,inpc))
     &   call errquit('purvis: cannot allocate nuc',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_dbl,ncore,'ALPC',lalpc,ialpc))
     &   call errquit('purvis: cannot allocate alpc',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_int,ncore*ncore,'NIN',lNIN,iNIN))
     &   call errquit('purvis: cannot allocate NIN',0,
     1    MA_ERR)
        L=0
       do I=0,ncore*ncore-1
        INT_MB(iNIN+I)= L
c       write(6,*)"NIN",INT_MB(iNIN+I)
        L=L+I+1
       end do
c     DO 120 I = 1, NA
c        LCFRST(I) = NFIRST(I)
c 120 CONTINUE
      CALL DEFCOR(NA,KAN,NTYP,LC,dbl_mb(iALPC),int_mb(iNWC),int_mb(iNUC)
     1 ,int_mb(iNPC),int_mb(iNFIRST),ISW1)
c     LC = NCORE
c      write(6,*)"LC =",LC 
      iF (LC.EQ.0) THEN
         RETURN
      END IF
      NTEST = NB
      if (.not.MA_Push_Get(MT_dbl,LC,'CEIG',lceig,iceig))
     &   call errquit('purvis: cannot allocate ceig',0,
     1    MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,LC*LC,'H',lH,iH))
     &   call errquit('purvis: cannot allocate H',0,
     1    MA_ERR)
      if (.not.MA_Push_Get(MT_dbl,LC*LC,'BB',lBB,iBB))
     &   call errquit('purvis: cannot allocate BB',0,
     1    MA_ERR)    
c 110 IF (LC.GT.NTEST) THEN
c        WRITE (6,9000) LC, NTEST, LC

c9000    FORMAT (' ***** THE NUMBER OF CORE ORBITALS',I5,
c    .          ' IS GREATER THAN THE NUMBER OF',/,
c    .          ' ***** VALENCE ORBITALS',I5,//,
c    .          ' *****THE CORE M.O.s WILL NOT BE CALCULATED UNLESS',
c    .          ' THE DYNAMIC ALLOCATION '/
c    .          ' ***** (DYNAL IN $CONTRL) IS INCREASED TO TOTAL NBASIS'
c    .         ,'=',I5)

c        RETURN
c     END IF
C     CALCULATE THE FIRST CORE ORBITAL LCFRST(I) FOR ATOM I
      NUCF = 0
      DO 130 I = 0, LC-1
         IA = int_mb(iNUC+I)
c        write(6,*)"IA",IA
         IF (IA.EQ.NUCF) GO TO 130
         int_mb(iNFIRST+IA-1) = I
         NUCF = IA
  130 CONTINUE
C     IF (IPRINT.GE.0) THEN
C        WRITE (6,9020)
C        WRITE (6,9030)
C        DO 150 I = 1, LC
C           J = NUC(I)
C           L = NWC(I)
C           IF (ISW1.NE.7) GO TO 140
C           WRITE (6,9040) I, J, AOS(L+4), L, NPC(I), NTYP(J), ALPC(I),
C    .        ALPC(I), ZERO, ONE, ONE, ZERO
C           GO TO 150

C 140       WRITE (6,9040) I, J, AOS(L+1), L, NPC(I), NTYP(J), ALPC(I),
C    .        ALPC(I), ZERO, ONE, ONE, ZERO
C 150    CONTINUE
C     END IF
C
C     THIS PORTION OF THE PROGRAMME IS CONCERNED WITH SETTING UP VARIOUS
C     GEOMETRICAL MATRICES
C
      NA1 = NA*NA
      NA2 = 2*NA1
      ITJ = 0
      DO 210 J = 1, NA
         DO 200 I = 1, NA
            IF (I.GT.J) IJ = J + (I* (I-1))/2
            IF (I.LE.J) IJ = I + (J* (J-1))/2
            ITJ = ITJ + 1
            IF (ABS(RAD(I,J)).LT.1.0D-4) GO TO 170
            IF (I-J) 160, 200, 160
  160       F(ITJ) = (CO(3,J)-CO(3,I))/RAD(I,J)
            G = ONE - F(ITJ)*F(ITJ)
            IF (THRESH-G) 190, 180, 180
  170       F(ITJ) = ONE
  180       F(ITJ+NA1) = ONE
            IF (F(ITJ).LT.ZERO) F(ITJ+NA1) = -ONE
            F(ITJ+NA2) = ZERO
            GO TO 200

  190       G = ONE/SQRT(G)
            F(ITJ+NA1) = (CO(1,J)-CO(1,I))*G/RAD(I,J)
            F(ITJ+NA2) = (CO(2,J)-CO(2,I))*G/RAD(I,J)
  200    CONTINUE
  210 CONTINUE
C     CALL CORONE TO FORM THE OVERLAP MATRIX OVER CORE ORBITALS
C
      CALL CORONE(Dbl_MB(iH),RAD,F,dbl_mb(iALPC),FACT,LC,int_mb(iNUC),NA
     1 ,int_mb(iNWC),KAN,int_mb(iNPC),Int_MB(iNIN),ISW1)
C
      LCLC = Int_MB(iNIN+LC-1) + LC
c     WRITE(6,*) ' CORE OVERLAP MATRIX IN PURVIS'
C     CALL PRNT(H,LC,LCLC)
C
C 3500 IF (ISAVE.EQ.5) CALL MOLTRN(H,SH,FA,REP)
C
C  SINCE SH IS A SUB BLOCK OF BB MOVE SH INTO F
C
C
C
  220 CONTINUE
C
C     THIS PORTION OF THE PROGRAMME SETS UP THE HAMILTONIAN H(IJ).
C
C     THE DIAGONAL OF H ARE TREATED AS PARAMETERS AND THESE
C     PARAMETERS ARE NOW BEING INSERTED. (A.U.)
C
      LL = 0
      DO 230 I = 1, LC
         KA = int_mb(iNUC+I-1)
c      write(6,*)"KA",KA
         CALL PARAMC(HHHHH,I,int_mb(iNPC+I-1),KA,KAN(KA),
     1   int_mb(iNWC+I-1))
         LL = LL + I
c      write(6,*)"HHHHH",HHHHH
         Dbl_MB(iH+LL-1) = HHHHH
  230 CONTINUE
C
C     WE NOW FILL IN THE DIAGONAL ELEMENTS OF THE HAMILTONIAN AND SCALE
C     THE OFF-DIAGONAL ELEMENTS.
C
      HUKFAC = 2.2D0
C
      DO 250 I = 2, LC
         JS = I - 1
         II = Int_MB(iNIN+I-1) + I
         DO 240 J = 1, JS
            JJ = Int_MB(iNIN+J-1) + J
            IJ = II - I + J
      Dbl_MB(iH+IJ-1) = HUKFAC* (Dbl_MB(iH+II-1)+Dbl_MB(iH+JJ-1))
     1 *Dbl_MB(iH+IJ-1)/TWO
  240    CONTINUE
  250 CONTINUE
C
c     WRITE(6,*) ' HAMILTONIAN OVER CORE'
C     CALL PRNT(H,LC,LCLC)
C
C     ZERO OUT SMALL OFF-DIAGONALIZATION ELEMENTS
C     DO 260 I = 1, LCLC
C        IF (ABS(H(I)).LE.1.0D-8) H(I) = 0.0D0
C 260 CONTINUE
c     IBGIV = .TRUE.
      call diag_para(Dbl_MB(iH),LC*(LC+1)/2,LC,Dbl_MB(iBB),
     1 dbl_mb(iCEIG))
      if (ga_nodeid().eq.0) write(6,*)"here 1"
c     CALL GIVENS(H,BB,CEIG,LC,LC,IBGIV,QQ,NIN)
C     SYMORB WOULD NEED NUC, NWC, NFIRST AND OTHERS

c     CALL SYMORB(BB,CEIG,NUC,NWC,JMO,LC,3,IFLAG)

C     WRITE(6,2222) (JMO(IM), IM=1, LC)
C2222 FORMAT(10I5)
c     IF (IPRINT.GE.0.and.ga_nodeid().eq.0) THEN
c        WRITE (6,9010)
c9010    FORMAT (//)

C        DO 270 I = 1, LC
C           JMO(I) = 0
C 270    CONTINUE
C        IFLAG = -1 SAYS NO SYMMETRY
C        IFLAG = -1
c        CALL MOOUT(JMO,CEIG,LC,IFLAG)
C        WRITE(6,*) ' ISYM IN PURVIS 3', ISYM, IFLAG
c        CALL MOPRNT(JMO,CEIG,BB,LC,LC,NUC,NWC,NPC,KAN,ATOM,IFLAG)
C        WRITE(6,*) ' ISYM IN PURVIS 4', ISYM, IFLAG
c     END IF
c     IF (IELEC.EQ.1) CALL INB(BB,ALPC,NPC,NWC,NUC,LC)
C     RESTORE NFIRST
c     DO 280 I = 1, NA
c        NFIRST(I) = LCFRST(I)
c 280 CONTINUE
C
C9020 FORMAT (//,18X,' CORE ORBITALS - FROZEN FOR THE CALCULATION '//)
 9030 FORMAT ('  AO.  ATM  SYM. TYPE    N  TYPE  ',
     .       '3 EXPONENTS',14X,'3 COEFFICIENTS'//)
 9040 FORMAT (1X,I3,2X,I3,2X,A8,I3,1X,I3,1X,I3,2X,6F8.3)

      if (.not.ma_pop_stack(lBB))
     &   call errquit('purvis: cannot pop stack lBB',3, MA_ERR)
      if (.not.ma_pop_stack(lH))
     &   call errquit('purvis: cannot pop stack lH',3, MA_ERR)
      if (.not.ma_pop_stack(lceig))
     &   call errquit('purvis: cannot pop stack lceig',3, MA_ERR)
      if (.not.ma_pop_stack(lNIN))
     &   call errquit('purvis: cannot pop stack lNIN',3, MA_ERR)
      if (.not.ma_pop_stack(lalpc))
     &   call errquit('purvis: cannot pop stack lalpc',3, MA_ERR)
      if (.not.ma_pop_stack(lnpc))
     &   call errquit('purvis: cannot pop stack lnpc',3, MA_ERR)
      if (.not.ma_pop_stack(lnwc))
     &   call errquit('purvis: cannot pop stack lnwc',3, MA_ERR)
      if (.not.ma_pop_stack(lnuc))
     &   call errquit('purvis: cannot pop stack lnuc',3, MA_ERR)
      if (.not.ma_pop_stack(lnfirst))
     &   call errquit('purvis: cannot pop stack lnfirst',3, MA_ERR)
c     IF (USEUNIT0) WRITE (0,*) 'PURVIS out'

      RETURN

      END
C-----------------------------------------------------------------------------------
C
C-----------------------------------------------------------------------------------
      SUBROUTINE  DEFCOR(NA,KAN,NTYP,LC,ALPC,NWC,NUC,NPC,NFIRST,ISW1)
C
C     *********************************************************
C     *        QUANTUM THEORY PROJECT, P.O. BOX 1184-35       *
C     *      UNIVERSITY OF FLORIDA, GAINESVILLE, FL 32611     *
C     *                                                       *
C     *  This subroutine was part of the old subroutine       *
C     *  PURVIS.F. It has been partially rewritten by J.D.    *
C     *  da Motta to calculate the inner shell MO's required  *
C     *  evaluation of the effective core potential WB in     *
C     *  the  ZINDO framework, in order to implement a        *
C     *   S-NDDO/1 model.                                     *
C     *                                                       *
C     *      Last revision:  Joaquim D. da Motta Neto,        *
C     *      QTP-UF, Gainesville, FL, October 5th, 1994.      *
C     *                                                       *
C     *********************************************************
C
C     This subroutine is called by GSNDO1.F, just before the
C     section that forms the one electron matrix. Apparently,
C     it could be called afterwards (below line labeled 4500) 
C     without difference, since the sections are pretty 
C     independent from each other. Must check this later.
      IMPLICIT DOUBLE PRECISION( A-H, O-Z )
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      CHARACTER*8 AOS
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0, IPRINT = 0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
      DIMENSION  NTYP(*),ALPC(*),NUC(*),NPC(*),NWC(*),
     .           KAN(*),NFIRST(*),AOS(16)
      DATA         AOS/' S      ', ' P(X)   ', ' P(Y)   ', ' P(Z)   ',
     X   ' D(Z2)  ', 'D(X2-Y2)', ' D(XY)  ', ' D(XZ)  ',
     X   ' D(YZ)  ', ' F(Z3)  ', ' F(XZ2) ', ' F(YZ2) ',
     X   'Z(X2-Y2)', ' F(XYZ) ', ' X3-3XY2', ' 3YX2-Y3'/

C     THE EXPONENTS ARE NOW BEING FILLED IN. ALPC(I)
C
C     EXCEPT WHERE NOTED THESE ARE SLATER EXPONENTS.
      write(6,*)"i am in defcor"
      LC = 0
C
      DO 1100 I = 1, NA
      KA = KAN(I)
      NT = NTYP(I)
      IF (KA.LE.2) GO TO 1100
C                                            
C     ++++++++++++++++++++++  1S  ++++++++++++++++++++++++++++
C                                                                
C                                            
C     ++++++++++++++++++++++  1S  ++++++++++++++++++++++++++++
C                                                                
      LC       = LC + 1
      ALPC(LC) = REAL(KA) - 0.30D0
      NWC(LC)  = 0
      NUC(LC)  = I
      NPC(LC)  = 1
      IF (KA.LE.10) GO TO 1100
C                                           
C     +++++++++++++++++++   2S 2P   ++++++++++++++++++++++++++
C
      ALP = (REAL(KA) - 1.7200D0 - 2.45D0) / 2.0D0
C
      DO 110 J = 1, 4
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = J - 1
      NUC(LC)  = I
  110 NPC(LC)  = 2
      IF (KA.LE.18) GO TO 1100
C
C     +++++++++++++++++++   3S 3P 3D   ++++++++++++++++++++++++
C
      ALP = ( REAL(KA) - 2.0D0 - 6.80D0 - 2.45D0 ) / 3.0D0
C
C     K = NU(I)   ( Used only in PURVIS.F   JDMN94 )
C
      DO 120 J = 1, 4
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = J - 1
      NUC(LC)  = I
  120 NPC(LC)  = 3
C        WRITE(6,*) ' DEFCOR KA, NT', KA, NT
         IF (KA.LE.30.AND.NT.GE.3) GO TO 1100

C       THESE ARE 3D
        ALP = (DFLOAT(KA)-13.589D0-2.43D0)/3.0D0
C       THIS IS THE WAHLGREN-ODD 3D INNER CORE.
c       IF(KA.LE.30) ALP = AD(KA,1)
        DO 130 J = 1, 5
         LC = LC + 1
         ALPC(LC) = ALP
         NWC(LC) = 3 + J
         NUC(LC) = I
         NPC(LC) = 3
  130   CONTINUE

C                                                         
C     ++++++++++++++++++++  4S 4P 4D   ++++++++++++++++++++++++
C
      IF (KA.LE.54) GO TO 1100
      ALP = ( REAL(KA) - 10.0 - 15.3 - 2.45 ) / 3.7D0
      DO 140 J = 1, 4
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = J - 1
      NUC(LC)  = I
  140 NPC(LC)  = 4
C
      IF (KA.LT.48) GO TO 1100
      IF (KA.EQ.48.AND.NT.GE.3) GO TO 1100
      ALP = ( REAL(KA) - 36.00 - 3.15D0 ) / 3.70D0
      DO 160 J = 1, 5
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = 3 + J
      NUC(LC)  = I
  160 NPC(LC)  = 4
      IF (KA.LE.54) GO TO 1100
C
C
C     +++++++++++++++++++++  5S 5P   +++++++++++++++++++++++++++
C
      ALP = ( REAL(KA) - 28.0D0 - 15.30D0 - 2.45D0 ) / 4.0D0
      DO 170 J = 1, 4
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = J - 1
      NUC(LC)  = I
  170 NPC(LC)  = 5
      IF (KA.LE.72) GO TO 1100
C
C     +++++++++++++++++++++++  4F  +++++++++++++++++++++++++++++
C
      ALP = ( REAL(KA) - 42.4 ) / 4.0D0
      DO 180 J = 1, 7
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = J + 8
      NUC(LC)  = I
  180 NPC(LC)  = 4
      IF (KA.LT.80) GO TO 1100
      IF (KA.EQ.80.AND.NT.GT.2) GO TO 1100
C
C     ++++++++++++++++++++++++  5D   ++++++++++++++++++++++++++++
C
C
      ALP = ( REAL(KA) - 48.00 ) / 5.0D0
      DO 190 J = 1, 5
      LC       = LC + 1
      ALPC(LC) = ALP
      NWC(LC)  = 3 + J
      NUC(LC)  = I
  190 NPC(LC)  = 5
      IF (KA.LE.86) GO TO 1100
C
C     +++++++++++++++++++++++  6S 6P  +++++++++++++++++++++++++++
C
      ALP = ( REAL(KA) - 48.0D0 ) / 6.0D0
      DO 200 J = 1, 4
      LC = LC + 1
C                                                           
      ALPC(LC) = ALP
      NWC(LC)  = J - 1
      NUC(LC)  = I
  200 NPC(LC)  = 6
C
 1100 CONTINUE
      NCORE = LC
      IF (LC.EQ.0) THEN
c        IF (USEUNIT0) WRITE (0,*) 'DEFCOR return one'
         RETURN
      ENDIF
      K0  = IZERO
      K1  = IZERO
      K4  = IZERO
      K9  = IZERO
      K16 = IZERO
c     WRITE(6,*) ' IN GSNDO 1', NA, LC

      DO  I=1,NA
       IF(NTYP(I) .LE. IZERO) THEN
        K0 = K0 + IONE
       ELSEIF(NTYP(I) .LE. 1) THEN
        K1 = K1 + IONE
       ELSEIF(NTYP(I) .LE. 2) THEN
        K4 = K4 + IONE
       ELSEIF(NTYP(I) .LE. 4) THEN
        K9 = K9 + IONE
       ELSE
        K16 = K16 + IONE
       ENDIF
      END DO
c     write(6,*)"K1,K4,K9,K16",K1,K4,K9,K16
      NTEST = K1 + 4 * K4 + 9 * K9 + 16 * K16
C
C     Calculate the first core orbital LCFRST(I) for atom I
      NUCF = 0
C
      DO 1150 I = 1, LC 
      IA = NUC(I) 
c     write(6,*)"IA",IA      
      IF(IA.EQ.NUCF) GO TO 1150
c     write(6,*)"IA =", IA
      NFIRST(IA) = I                                                       
      NUCF      = IA
 1150 CONTINUE
C                                                          
      IF (IPRINT.GE.0.AND.ga_nodeid().eq.0) THEN                                              
        WRITE  (6, 14000)
14000   FORMAT (//,18X,'CORE ORBITALS - FROZEN FOR THE CALCULATION '//)
        WRITE  (6, 16000)
16000   FORMAT ('  AO.   ATM  SYM. TYPE     N ',
     x  ' TYPE   3 EXPONENTS', 14X, '3 COEFFICIENTS'//)
C                                                 
        NCORE = LC

        DO 2300 I = 1, LC                                                  
        J = NUC(I)                                                        
        L = NWC(I)                                                        
        IF (ISW1.NE.7) GO TO 2200                                        
        WRITE (6,52000)  I, J, AOS(L+4), L, NPC(I), NTYP(J), ALPC(I),
     1    ALPC(I), ZERO, ONE, ONE, ZERO
52000   FORMAT(1X,I3,2X,I3,2X,A8,I3,1X,I3,1X,I3,2X,6F8.3)
        GO TO 2300
 
 2200   WRITE (6,52000)  I, J, AOS(L+1), L, NPC(I), NTYP(J), ALPC(I),
     1    ALPC(I), ZERO, ONE, ONE, ZERO
 2300   CONTINUE
      ENDIF
C
      RETURN
      END

