      SUBROUTINE GFORMAM(NB,NA,GAMMA,NU,KAN,LL,LABI,LABJ,LABK,LABL,
     1           QQ,GSS,GPP,GSP,GP2,HSP,HPP,PO,DAD,NFIRST,CO)
C     From voityuk@theochem.tu-muenchen.de Mon Jul 28 10:49:30 1997
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/z.cmn'
C     *
C     W- ALL 2-electron 2-centren INTEGRALS
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MXCI = 6000,
     + ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     + BOHRS    = 5.29177249D-1)
      DIMENSION  SP(3,3),PP(6,3,3)
      DIMENSION MAT(NA),NLAST(NA),IND2(45,45),INDX(9,9),INDEX(9,9)
      DIMENSION WW(2025),LABI(*),LABJ(*),LABK(*),LABL(*),QQ(*)
      DIMENSION NFIRST(NA),GSS(86),GPP(105),GSP(105),GP2(86),HSP(86),
     1         HPP(86),KAN(*),NU(*),GAMMA(*), PO(*), DAD(*), CO(*)
C
c     write(6,*)"i am in gformam"
      CALL INDXXX(IND2,INDX,INDEX)
C
      NA1=0
      NUMAT=0
      DO I=1, NB
        NA2=NU(I)
        IF(NA2.NE.NA1)THEN
          NA1=NA2
          NUMAT=NUMAT+1
          MAT(NUMAT)=KAN(I)
          NLAST(NUMAT)=I
         ELSE
          NLAST(NUMAT)=I
        ENDIF
      ENDDO
c     write(6,*)"NUMAT",NUMAT
      LL = 0
C     LL IS THE INTEGRAL COUNTER

      CALL GFORMA1(NUMAT,MAT,NFIRST,LABI,LABJ,LABK,LABL,QQ,GSS,GPP,
     1     GSP,GP2,HSP,HPP,LL)

c     write(6,*)"LL in gforma1",LL
C *** TWO-ELECTRON REPULSIONS 
      IJPAIR = (NUMAT*(NUMAT+1))/2
      KINT=LL
c     write(6,*)"QQ",(QQ(I),I=1,nb*nb)
C     SUM OVER PAIRS OF ATOMS
C     WRITE(6,*) ' GFORMA IJPAIR, NUMAT', IJPAIR, NUMAT

c     write(6,*)"NFIRST",(NFIRST(I),I=1,NA)
      DO 90 IJ=2,IJPAIR
C     I IS FIRST ATOM, J IS SECOND
      I      = NINT(SQRT(TWO*FLOAT(IJ)))
      J      = IJ-(I*(I-1))/2
c     write(6,*)"IJ,I",IJ,I
      IF(I.EQ.J) GO TO 90
      NI     = MAT(I)
      NJ     = MAT(J)
      IA     = NFIRST(I)
      JA     = NFIRST(J)
      IB     = NLAST(I)
      JB     = NLAST(J)
C     WRITE(6,*) ' GFORMAM', NI, NJ, IA,IB,JA, JB
C
C *** DISTANCE R (AU)  ROTATION MATRIX
c     write(6,*)"I,J",I,J
      CALL ROTAT (J,I,NJ,NI,CO,R,SP,PP)
C     CONVERT DISTANCE TO ATOMIC UNITS.
c     RAU=R/BOHRS
      RAU = R
c     write(6,*)"RAU",RAU
C     WRITE(*,'(1X,A,6I5,F8.4)')
C    .   'IA,IB,JA,JB,NI,NJ,RAU',IA,IB,JA,JB,NI,NJ,RAU
      CALL ROTINT (IA,IB,JA,JB,NI,NJ,RAU,WW,PO,DAD,INDX,INDEX,IND2,SP,PP
     . )
C
c     write(6,*)"LL before",LL
      IY = 0
      DO  I=IA,IB
      DO  J=IA,I
      DO  K=JA,JB
      DO  L=JA,K
       IY = IY + 1
       LL = LL + 1
c      write(6,*)"WW",WW(IY)
       QQ(LL)=WW(IY)
c      WRITE(6,'(1x,I4,''<'',I2,'','',I2,'','',I2,'',''I2,
c    .      ''> ='',F12.5)')LL,L,K,J,I,QQ(LL)
       LABI(LL) = L
       LABJ(LL) = K
       LABK(LL) = J
       LABL(LL) = I
      ENDDO
      ENDDO
      ENDDO
      ENDDO
c     WRITE(6,*)"LL",LL,nb*(nb+1)/2
c     write(6,*)"QQ",(QQ(I),I=1,nb*(nb+1)/2)
C     ENDDO
C
C     DEBUG
c     IF (LL.LT.10000) GO TO  80
c           NINTGL = 10000
C 
c           CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, NINTGL, 1, 1, 'GFORM')
c           WRITE (6,*)  ' NO. OF AO. INTEGRALS', NINTGL
c           LL = LL - 10000
c           IF (LL.EQ.0) GO TO  80
c
c           DO  70 K = 1, LL
c              QQ(K) = QQ(K+10000)
c              LABI(K) = LABI(K+10000)
c              LABJ(K) = LABJ(K+10000)
c              LABK(K) = LABK(K+10000)
c              LABL(K) = LABL(K+10000)
c  70       CONTINUE
c  80    CONTINUE
   90 CONTINUE

C 
c        IF (LL.EQ.0) GO TO 160
c     CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LL, 1, 1, 'GFORM')
C
  160 CONTINUE
c     LL = -1
C
c     CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LL, 1, 1, 'GFORM')


      RETURN
      END

C-----------------------------------------------------------------------------
C
C-----------------------------------------------------------------------------
      SUBROUTINE INDXXX(IND2,INDX,INDEX)
C     *
C     DEFINITION OF INDICES
C     *
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION IND2(45,45)
      DIMENSION INDEX(9,9),INDX(9,9)
C
c      write(6,*)"i am in indxxx"
         DO I=1,9
         DO J=1,I
          INDEX(I,J) = -(J*(J-1))/2 + I + 9*(J-1)
          INDX(I,J) = (I*(I-1))/2 + J
          INDEX(J,I)=INDEX(I,J)
          INDX(J,I)=INDX(I,J)
         ENDDO
         ENDDO
C
         DO I=1,25
         DO J=1,25
           IND2( I, J) =  0
         ENDDO
         ENDDO
C
      IND2( 1, 1) =  1
      IND2( 1, 2) =  2
      IND2( 1,10) =  3
      IND2( 1,18) =  4
      IND2( 1,25) =  5
      IND2( 2, 1) =  6
      IND2( 2, 2) =  7
      IND2( 2,10) =  8
      IND2( 2,18) =  9
      IND2( 2,25) = 10
      IND2(10, 1) = 11
      IND2(10, 2) = 12
      IND2(10,10) = 13
      IND2(10,18) = 14
      IND2(10,25) = 15
      IND2( 3, 3) = 16
      IND2( 3,11) = 17
      IND2(11, 3) = 18
      IND2(11,11) = 19
      IND2(18, 1) = 20
      IND2(18, 2) = 21
      IND2(18,10) = 22
      IND2(18,18) = 23
      IND2(18,25) = 24
      IND2( 4, 4) = 25
      IND2( 4,12) = 26
      IND2(12, 4) = 27
      IND2(12,12) = 28
      IND2(19,19) = 29
      IND2(25, 1) = 30
      IND2(25, 2) = 31
      IND2(25,10) = 32
      IND2(25,18) = 33
      IND2(25,25) = 34
      RETURN
      END
C-------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
      SUBROUTINE GFORMA1(NUMAT,MAT,NFIRST,LABI,LABJ,LABK,LABL,QQ,GSS,
     1 GPP,GSP,GP2,HSP,HPP,LL)
C
C      One-center two-electron integrals for all atoms
C      From voityuk@theochem.tu-muenchen.de Wed Jul 30 08:53:39 1997
C      transferred to W. Indexes LI,LJ,LK,LL defined.
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/am1.cmn'
      DIMENSION MAT(*),NFIRST(*),LABI(*),LABK(*),LABJ(*)
      DIMENSION LP2(4,3),LHP(4,3),LABL(*),QQ(*)
      DIMENSION GSS(86),GPP(105),GSP(105),GP2(86),HSP(86),
     1         HPP(86)
      PARAMETER  (EVOLT    = 2.72113961D1)
      DATA LP2/1,1,2,2, 1,1,3,3, 2,2,3,3/
      DATA LHP/1,2,1,2, 1,3,1,3, 2,3,2,3/
c     SAVE LP2, LHP
C
      write(6,*)"i am in gforma1"
      DO 100 I=1,NUMAT
         NI=MAT(I)
         IA=NFIRST(I)
         LL=LL+1
C    GSS
         QQ(LL)=GSS(NI)/EVOLT
         LABI(LL)=IA
         LABJ(LL)=IA
         LABK(LL)=IA
         LABL(LL)=IA
         IF(NI.LT.3) GOTO 100
C    GSP
          DO J=1,3
           LL=LL+1
           QQ(LL)=GSP(NI)/EVOLT
           LABI(LL)=IA
           LABJ(LL)=IA
           LABK(LL)=IA+J
           LABL(LL)=IA+J
          ENDDO
C    HSP
          DO J=1,3
           LL=LL+1
           QQ(LL)=HSP(NI)/EVOLT
           LABI(LL)=IA
           LABJ(LL)=IA+J
           LABK(LL)=IA
           LABL(LL)=IA+J
          ENDDO
C    GPP
          DO J=1,3
           LL=LL+1
           QQ(LL)=GPP(NI)/EVOLT
           LABI(LL)=IA+J
           LABJ(LL)=IA+J
           LABK(LL)=IA+J
           LABL(LL)=IA+J
          ENDDO
C    GP2
          DO J=1,3
           LL=LL+1
           QQ(LL)=GP2(NI)/EVOLT
           LABI(LL)=IA+LP2(1,J)
           LABJ(LL)=IA+LP2(2,J)
           LABK(LL)=IA+LP2(3,J)
           LABL(LL)=IA+LP2(4,J)
          ENDDO
C    HPP
          DO J=1,3
           LL=LL+1
           QQ(LL)=HPP(NI)/EVOLT
           LABI(LL)=IA+LHP(1,J)
           LABJ(LL)=IA+LHP(2,J)
           LABK(LL)=IA+LHP(3,J)
           LABL(LL)=IA+LHP(4,J)
          ENDDO
c           IF (LL.LT.10000) GO TO 90
c           NINTGL = 10000
C
c           CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, NINTGL, 1, 1, 'GFORM')
c           WRITE (6,*)  ' NO. OF AO. INTEGRALS', NINTGL
C

c           LL = LL - 10000
c           IF (LL.EQ.0) GO TO 90
C
c           DO 80  K = 1, LL
c              QQ(K) = QQ(K+10000)
c              LABI(K) = LABI(K+10000)
c              LABJ(K) = LABJ(K+10000)
c              LABK(K) = LABK(K+10000)
c              LABL(K) = LABL(K+10000)
c  80       CONTINUE

C
   90    CONTINUE

C
 100  CONTINUE
C
C      DEBUG PRINT
C     WRITE(6,'(1x,''One-center  integrals for all atoms'')')
C      DO J=1,LL
C        WRITE(6,'(1x,I5,''<'',I2,'','',I2,'','',I2,'',''I2,
C    .      ''> ='',F12.5)')J,LABI(J),LABJ(J),LABK(J),LABL(J),QQ(J)
C      ENDDO
      RETURN
      END
C-----------------------------------------------------------------------------
C
C----------------------------------------------------------------------------

      SUBROUTINE ROTAT (J,I,NJ,NI,COORD,R,SP,PP)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/const.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      PARAMETER (MXCI = 6000,
     + ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0)

CSG
C     *
C     ROTATION MATRIX FOR A GIVEN ATOM PAIR I-J (I.GT.J).
C     *
      PARAMETER (SMALL=1.0D-07, BOHRS    = 5.29177249D-1)
      DIMENSION  SP(3,3),PP(6,3,3)
      DIMENSION P(3,3)
      DIMENSION COORD(3,*)

C     CA  = COS(PHI)    , SA  = SIN(PHI)
C     CB  = COS(THETA)  , SB  = SIN(THETA)
C     C2A = COS(2*PHI)  , S2A = SIN(2*PHI)
C     C2B = COS(2*THETA), S2B = SIN(2*PHI)
      X11    = (COORD(1,J)-COORD(1,I))!*BOHRS
      X22    = (COORD(2,J)-COORD(2,I))!*BOHRS
      X33    = (COORD(3,J)-COORD(3,I))!*BOHRS
C     WRITE(6,'(1X,''COORDINATES '')')
C     WRITE(6,'(1X,I5,3F8.3 )')J,(COORD(K,J),K=1,3)
C     WRITE(6,'(1X,I5,3F8.3 )')I,(COORD(K,I),K=1,3)
      B2     = X11*X11+X22*X22
      B      = SQRT(B2)
      B3     = B*B2
      R2     =(B2+X33*X33)
c     write(6,*)"B2,X33",B2,X33
      R      = SQRT(R2)
      IF(NI.LT.3.AND.NJ.LT.3) GOTO 200
C
c     write(6,*)"B,X11,X22,X33",b,X11,X22,X33
      R3     = R*R2
      SB     = B/R
C     CHECK FOR SPECIAL CASE (BOTH ATOMS ON Z AXIS).
      IF(SB.GT.SMALL) THEN
         CA  = X11/B
         SA  = X22/B
         CB  = X33/R
      ELSE
         SA  = ZERO
         SB  = ZERO
         IF(X33.LT.ZERO) THEN
            CA  =-ONE
            CB  =-ONE
         ELSE IF(X33.GT.ZERO) THEN
            CA  = ONE
            CB  = ONE
         ELSE
            CA  = ZERO
            CB  = ZERO
         ENDIF
      ENDIF
C *** CALCULATE ROTATION MATRIX ELEMENTS.
c     write(6,*)"CA,SB,CB,SA",CA,SB,CB,SA
      P(1,1) = CA*SB
      P(2,1) = CA*CB
      P(3,1) =-SA
      P(1,2) = SA*SB
      P(2,2) = SA*CB
      P(3,2) = CA
      P(1,3) = CB
      P(2,3) =-SB
      P(3,3) = ZERO
C  S-P
      DO  K=1,3
      DO  L=1,3
      SP(K,L) =   P(K,L)
      ENDDO
      ENDDO
C  P-P
      DO  K=1,3
      PP(1,K,K) =       P(K,1)*P(K,1)
      PP(2,K,K) =       P(K,2)*P(K,2)
      PP(3,K,K) =       P(K,3)*P(K,3)
      PP(4,K,K) =       P(K,1)*P(K,2)
      PP(5,K,K) =       P(K,1)*P(K,3)
      PP(6,K,K) =       P(K,2)*P(K,3)
      IF(K.GT.1) THEN
      DO  L=1,K-1
      PP(1,K,L)=        TWO*  P(K,1)*P(L,1)
      PP(2,K,L)=        TWO*  P(K,2)*P(L,2)
      PP(3,K,L)=        TWO*  P(K,3)*P(L,3)
      PP(4,K,L)=        P(K,1)*P(L,2)+ P(K,2)*P(L,1)
      PP(5,K,L)=        P(K,1)*P(L,3)+ P(K,3)*P(L,1)
      PP(6,K,L)=        P(K,2)*P(L,3)+ P(K,3)*P(L,2)
      ENDDO
      ENDIF
      ENDDO
C
200   CONTINUE
      RETURN
      END
C-----------------------------------------------------------------------------
C
C----------------------------------------------------------------------------
      SUBROUTINE ROTINT(IA,IB,JA,JB,NI,NJ,R,WW,PO,DAD,INDX,INDEX,IND2
     .  ,SP,PP)
C     *
C     CALCULATION OF TWO-CENTER TWO-ELECTRON INTEGRALS
C     IN THE MOLECULAR COODINATE SYSTEM BY TWO-STEP PROCEDURE
C     *
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      LOGICAL LOGV
      DIMENSION  SP(3,3),PP(6,3,3)
      DIMENSION INDEX(9,9),INDX(9,9)
      DIMENSION  V(45,45),LOGV(45,45),LORB(9)
      DIMENSION  RI(22),REP(491)
C
      DIMENSION WW(*),PO(*),DAD(*),IND2(*)
      DIMENSION MET(45)
      DIMENSION IPOS(34)
      DATA IPOS/ 1, 5, 11, 12, 12, 2, 6, 13, 14, 14, 3, 8, 16, 18, 18,
     -           7,15, 10, 20,  4, 9,17, 19, 21,  7,15,10, 20, 22,  4,
     -           9,17, 21, 19 /
      DATA LORB  /0,3*1,5*2/
      DATA MET/ 1,2,3,2,3,3,2,3,3,3,4,5,5,5,6,4,5,5,5,6,6,
     1          4,5,5,5,6,6,6,4,5,5,5,6,6,6,6,4,5,5,5,5*6/
C
      INDW(I,J)=(INDX(I,J)-1)*LIMKL+KL
C
c     write(6,*)"INDW",INDW(i,J)
      CALL RESP(NI,NJ,R,RI,PO,DAD)
      DO  I=1,34
      INT=IPOS(I)
c     write(6,*)"RI",RI(INT)
      REP(I)=RI(INT)
      ENDDO
C
      II=IB-IA+1
      KK=JB-JA+1
      LIMIJ=INDX(II,II)
      LIMKL=INDX(KK,KK)
      ISTEP=LIMKL*LIMIJ
      DO I=1,ISTEP
      WW(I)=0.0D0
      ENDDO
C
      CALL RI1(II,KK,REP,LOGV,V,INDX,INDEX,IND2,SP,PP,LL)
c
      DO 100 I1=1,II
      DO  80 J1=1,I1
C
      IJ=INDEX(I1,J1)
      JJ=INDX(I1,J1)
      MM=MET(JJ)
C
      DO  50 K=1,KK
      DO  40 L=1,K
      KL=INDX(K,L)
      IF(.NOT.LOGV(IJ,KL))GOTO 40
c     write(6,*)"V in rotint",V(IJ,KL)
      WREPP=V(IJ,KL)
C
C
      IF(MM.EQ.1)THEN
      IW=INDW(1,1)
      WW(IW)= WREPP
      ELSEIF(MM.EQ.2)THEN
      DO I=1,3
      IW  =INDW(I+1,1)
      WW(IW)=WW(IW)+SP(I1-1,i) * WREPP
      ENDDO
      ELSEIF(MM.EQ.3)THEN
      DO I=1,3
      CC=PP(I,I1-1,J1-1)
      IW=INDW(I+1,I+1)
      WW(IW)=WW(IW)+CC * WREPP
      IMINUS=I-1
      IF(IMINUS.NE.0)THEN
      DO J=1,IMINUS
      CC=PP(1+I+J,I1-1,J1-1)
      IW=INDW(I+1,J+1)
      WW(IW)=WW(IW)+CC * WREPP
      ENDDO
      ENDIF
      ENDDO
      ENDIF
C
 40   CONTINUE
 50   CONTINUE
 80   CONTINUE
 90   CONTINUE
 100  CONTINUE
C
      RETURN
      END

C--------------------------------------------------------------------
C
C--------------------------------------------------------------------

      SUBROUTINE RI1(II,KK,REP,LOGV,V,INDX,INDEX,IND2,SP,PP,LL)
C     *
C     ROTATION OF TWO-ELECTRON TWO-CENTER INTEGRALS
C     FIRST STEP
C     *
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)

      LOGICAL LOGV
      DIMENSION SP(3,3),PP(6,3,3)
      DIMENSION IND2(45,45)
      DIMENSION INDEX(9,9),INDX(9,9)
      DIMENSION V(45,45),LOGV(45,45)
      DIMENSION REP(491)
      DIMENSION MET(45)
      DATA MET / 1,2,3,2,3,3,2,3,3,3,4,5,5,5,6,4,5,5,5,6,6,
     -           4,5,5,5,6,6,6,4,5,5,5,6,6,6,6,4,5,5,5,5*6/
C
      LIMKL=INDX(KK,KK)
      DO   I=1,45
      DO   K=1,LIMKL
      LOGV(I,K)=.FALSE.
      V(I,K)=0.D0
      ENDDO
      ENDDO
C
      DO 100 I1=1,II
      DO  80 J1=1,I1
      IJ=INDEX(I1,J1)
c
      DO  60  K1=1,KK
C
      DO  50 L1=1,K1
      KL=INDEX(K1,L1)
      ND=IND2(IJ,KL)
      IF(ND.NE.0) THEN
C
      WREPP=REP(ND)
      LL=INDX(K1,L1)
      MM=MET(LL)
C
c     write(6,*)"WREPP,MM",WREPP,MM
      IF(MM.EQ.1)THEN
      V(IJ,1 )= WREPP
      ELSEIF(MM.EQ.2)THEN
      K=K1-1
      V(IJ,2)=V(IJ,2)+SP(K,1) * WREPP
      V(IJ,4)=V(IJ,4)+SP(K,2) * WREPP
      V(IJ,7)=V(IJ,7)+SP(K,3) * WREPP
      ELSEIF(MM.EQ.3)THEN
      K=K1-1
      L=L1-1
      V(IJ,3)= V(IJ,3)  + PP(1,K,L)* WREPP
      V(IJ,6)= V(IJ,6)  + PP(2,K,L)* WREPP
      V(IJ,10)=V(IJ,10) + PP(3,K,L)* WREPP
      V(IJ,5)= V(IJ,5)  + PP(4,K,L)* WREPP
      V(IJ,8)= V(IJ,8)  + PP(5,K,L)* WREPP
      V(IJ,9)= V(IJ,9)  + PP(6,K,L)* WREPP
      ENDIF
      ENDIF
  50  CONTINUE
  60  CONTINUE
      DO KL=1,LIMKL
c     write(6,*)"V",V(IJ,KL)
      IF(V(IJ,KL).NE.0.00D00)LOGV(IJ,KL)=.TRUE.
      ENDDO
  80  CONTINUE
 100  CONTINUE
C
      RETURN
      END
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
      SUBROUTINE RESP (NI,NJ,R,A,PO,DAD)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/am1.cmn'
C     *
C     THIS ROUTINE COMPUTES THE TWO-CENTRE REPULSION INTEGRALS AND THE
C     NUCLEAR ATTRACTION INTEGRALS IN LOCAL COORDINATES.
C     *
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
     
      PARAMETER (MXCI = 6000,
     + ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0)
 
      DIMENSION A(22), PO(9,105), DAD(9,105)
      DIMENSION X(69),PXX(33),PXY(69)
      DATA PXX/ 1.0D0   ,-0.5D0   ,-0.5D0   , 0.5D0   , 0.25D0  ,
     1          0.25D0  , 0.5D0   , 0.25D0  ,-0.25D0  , 0.25D0  ,
     2         -0.25D0  ,-0.125D0 ,-0.125D0 , 0.50D0  , 0.125D0 ,
     3          0.125D0 ,-0.5D0   ,-0.25D0  ,-0.25D0  ,-0.25D0  ,
     4          0.25D0  ,-0.125D0 , 0.125D0 ,-0.125D0 , 0.125D0 ,
     5          0.125D0 , 0.25D0  , 0.0625D0, 0.0625D0,-0.25D0  ,
     6          0.25D0  , 0.25D0  ,-0.25D0  /
      DATA PXY/ 1.0D0   ,-0.5D0   ,-0.5D0   , 0.5D0   , 0.25D0  ,
     1          0.25D0  , 0.5D0   , 0.25D0  ,-0.25D0  , 0.25D0  ,
     2         -0.25D0  ,-0.125D0 ,-0.125D0 ,-0.50D0  ,-0.50D0  ,
     3          0.5D0   , 0.25D0  , 0.25D0  , 0.5D0   , 0.25D0  ,
     4         -0.25D0  ,-0.25D0  ,-0.125D0 ,-0.125D0 , 0.5D0   ,
     5         -0.5D0   , 0.25D0  , 0.25D0  ,-0.25D0  ,-0.25D0  ,
     6         -0.25D0  , 0.25D0  ,-0.25D0  , 0.25D0  ,-0.125D0 ,
     7          0.125D0 ,-0.125D0 , 0.125D0 ,-0.125D0 , 0.125D0 ,
     8         -0.125D0 , 0.125D0 , 0.125D0 , 0.125D0 , 0.25D0  ,
     9          0.125D0 , 0.125D0 , 0.125D0 , 0.125D0 , 0.0625D0,
     A          0.0625D0, 0.0625D0, 0.0625D0,-0.25D0  , 0.25D0  ,
     B          0.25D0  ,-0.25D0  ,-0.25D0  , 0.25D0  , 0.25D0  ,
     C         -0.25D0  , 0.125D0 ,-0.125D0 ,-0.125D0 , 0.125D0 ,
     D         -0.125D0 , 0.125D0 , 0.125D0 ,-0.125D0 /
c     write(6,*)"i am in resp",NI,NJ
C *** INITIALIZATION.
      R2     = R*R
      AEE    = (PO(1,NI)+PO(1,NJ))**2
C *** HYDROGEN - HYDROGEN.
      IF(NI.LT.3 .AND. NJ.LT.3) THEN
         EE     = ONE/SQRT(R2+AEE)
         A(1)   = EE
C *** HEAVY ATOM - HYDROGEN
      ELSE IF(NI.GT.2 .AND. NJ.LT.3) THEN
         DA     = DAD(2,NI)
         QA     = DAD(3,NI)/SQRT(TWO)
         TWOQA  = QA+QA
         ADE    = (PO(2,NI)+PO(1,NJ))**2
         AQE    = (PO(3,NI)+PO(1,NJ))**2
         X(1)   = R2+AEE
         X(2)   = R2+AQE
         X(3)   = (R+DA)**2+ADE
         X(4)   = (R-DA)**2+ADE
         X(5)   = (R-TWOQA)**2+AQE
         X(6)   = (R+TWOQA)**2+AQE
         X(7)   = R2+TWOQA*TWOQA+AQE
         DO  I=1,7
         X(I)   = PXY(I)/SQRT(X(I))
         ENDDO
         A(1)   =  X(1)
         A(2)   = (X(3)+X(4))
         A(3)   = (X(1)+X(2)+X(5)+X(6))
         A(4)   = (X(1)+X(2)+X(7))
C *** HYDROGEN - HEAVY ATOM
      ELSE IF(NI.LT.3 .AND. NJ.GT.2) THEN
         DB     = DAD(2,NJ)
         QB     = DAD(3,NJ)/SQRT(TWO)
         TWOQB  = QB+QB
         AED    = (PO(1,NI)+PO(2,NJ))**2
         AEQ    = (PO(1,NI)+PO(3,NJ))**2
c        write(6,*)"DB,QB,AED,AEQ",DB,QB,AED,AEQ
         X(1)   = R2+AEE
         X(2)   = R2+AEQ
         X(3)   = (R-DB)**2+AED
         X(4)   = (R+DB)**2+AED
         X(5)   = (R-TWOQB)**2+AEQ
         X(6)   = (R+TWOQB)**2+AEQ
         X(7)   = R2+TWOQB*TWOQB+AEQ
         DO  I=1,7
         X(I)   = PXY(I)/SQRT(X(I))
         ENDDO
         A(1)   =  X(1)
         A(5)   = (X(3)+X(4))
         A(11)  = (X(1)+X(2)+X(5)+X(6))
         A(12)  = (X(1)+X(2)+X(7))
C *** HEAVY ATOM - HEAVY ATOM
      ELSE
         DA     = DAD(2,NI)
         QA     = DAD(3,NI)/SQRT(TWO)
         TWOQA  = QA+QA
         TWOQA2 = TWOQA*TWOQA
         RPDA2  = (R+DA)**2
         RMDA2  = (R-DA)**2
         RP2QA2 = (R+TWOQA)**2
         RM2QA2 = (R-TWOQA)**2
         ADE    = (PO(2,NI)+PO(1,NJ))**2
         AQE    = (PO(3,NI)+PO(1,NJ))**2
         ADD    = (PO(2,NI)+PO(2,NJ))**2
         ADQ    = (PO(2,NI)+PO(3,NJ))**2
         AQQ    = (PO(3,NI)+PO(3,NJ))**2
         TWOQAQ = TWOQA2+AQQ
         X(1)   = R2+AEE
         X(2)   = R2+AQE
         X(3)   = RPDA2+ADE
         X(4)   = RMDA2+ADE
         X(5)   = RM2QA2+AQE
         X(6)   = RP2QA2+AQE
         X(7)   = R2+TWOQA2+AQE
         X(8)   = RPDA2+ADQ
         X(9)   = RMDA2+ADQ
         X(10)  = R2+AQQ
         X(11)  = R2+TWOQAQ
         X(12)  = RP2QA2+AQQ
         X(13)  = RM2QA2+AQQ
         IF(NI.NE.NJ) THEN
            DB     = DAD(2,NJ)
            QB     = DAD(3,NJ)/SQRT(TWO)
            TWOQB  = QB+QB
            TWOQB2 = TWOQB*TWOQB
            TWOQBQ = TWOQB2+AQQ
            RPDB2  = (R+DB)**2
            RMDB2  = (R-DB)**2
            RP2QB2 = (R+TWOQB)**2
            RM2QB2 = (R-TWOQB)**2
            AED    = (PO(1,NI)+PO(2,NJ))**2
            AEQ    = (PO(1,NI)+PO(3,NJ))**2
            AQD    = (PO(3,NI)+PO(2,NJ))**2
            X(14)  = R2+AEQ
            X(15)  = RMDB2+AED
            X(16)  = RPDB2+AED
            X(17)  = RM2QB2+AEQ
            X(18)  = RP2QB2+AEQ
            X(19)  = R2+TWOQB2+AEQ
            X(20)  = RMDB2+AQD
            X(21)  = RPDB2+AQD
            X(22)  = R2+TWOQBQ
            X(23)  = RP2QB2+AQQ
            X(24)  = RM2QB2+AQQ
            X(25)  = R2+(DA-DB)**2+ADD
            X(26)  = R2+(DA+DB)**2+ADD
            X(27)  = (R+DA-DB)**2+ADD
            X(28)  = (R-DA+DB)**2+ADD
            X(29)  = (R-DA-DB)**2+ADD
            X(30)  = (R+DA+DB)**2+ADD
            X(31)  = RPDA2+TWOQB2+ADQ
            X(32)  = RMDA2+TWOQB2+ADQ
            X(33)  = RMDB2+TWOQA2+AQD
            X(34)  = RPDB2+TWOQA2+AQD
            X(35)  = (R+DA-TWOQB)**2+ADQ
            X(36)  = (R-DA-TWOQB)**2+ADQ
            X(37)  = (R+DA+TWOQB)**2+ADQ
            X(38)  = (R-DA+TWOQB)**2+ADQ
            X(39)  = (R+TWOQA-DB)**2+AQD
            X(40)  = (R+TWOQA+DB)**2+AQD
            X(41)  = (R-TWOQA-DB)**2+AQD
            X(42)  = (R-TWOQA+DB)**2+AQD
            X(43)  = R2+FOUR*(QA-QB)**2+AQQ
            X(44)  = R2+FOUR*(QA+QB)**2+AQQ
            X(45)  = R2+TWOQA2+TWOQBQ
            X(46)  = RM2QB2+TWOQAQ
            X(47)  = RP2QB2+TWOQAQ
            X(48)  = RP2QA2+TWOQBQ
            X(49)  = RM2QA2+TWOQBQ
            X(50)  = (R+TWOQA-TWOQB)**2+AQQ
            X(51)  = (R+TWOQA+TWOQB)**2+AQQ
            X(52)  = (R-TWOQA-TWOQB)**2+AQQ
            X(53)  = (R-TWOQA+TWOQB)**2+AQQ
            X(54)  = (R-QB)**2+(DA-QB)**2+ADQ
            X(55)  = (R+QB)**2+(DA-QB)**2+ADQ
            X(56)  = (R-QB)**2+(DA+QB)**2+ADQ
            X(57)  = (R+QB)**2+(DA+QB)**2+ADQ
            X(58)  = (R+QA)**2+(QA-DB)**2+AQD
            X(59)  = (R-QA)**2+(QA-DB)**2+AQD
            X(60)  = (R+QA)**2+(QA+DB)**2+AQD
            X(61)  = (R-QA)**2+(QA+DB)**2+AQD
            QMADD  = (QA-QB)**2+AQQ
            QPADD  = (QA+QB)**2+AQQ
            X(62)  = (R+QA-QB)**2+QMADD
            X(63)  = (R+QA+QB)**2+QMADD
            X(64)  = (R-QA-QB)**2+QMADD
            X(65)  = (R-QA+QB)**2+QMADD
            X(66)  = (R+QA-QB)**2+QPADD
            X(67)  = (R+QA+QB)**2+QPADD
            X(68)  = (R-QA-QB)**2+QPADD
            X(69)  = (R-QA+QB)**2+QPADD
            DO  I=1,69
            X(I)   = PXY(I)/SQRT(X(I))
            ENDDO
            EE     = X(1)
            DZE    = X(3) +X(4)
            QZZE   = X(2) +X(5) +X(6)
            QXXE   = X(2) +X(7)
            EDZ    = X(15)+X(16)
            EQZZ   = X(14)+X(17)+X(18)
            EQXX   = X(14)+X(19)
            DXDX   = X(25)+X(26)
            DZDZ   = X(27)+X(28)+X(29)+X(30)
            X89    = X(8) +X(9)
            X2021  = X(20)+X(21)
            DZQXX  = X89  +X(31)+X(32)
            QXXDZ  = X2021+X(33)+X(34)
            DZQZZ  = X89  +X(35)+X(36)+X(37)+X(38)
            QZZDZ  = X2021+X(39)+X(40)+X(41)+X(42)
            X1011  = X(10)+X(11)
            X1213  = X(12)+X(13)
            X2324  = X(23)+X(24)
            QXXQXX = X1011+X(22)+X(43)+X(44)
            QXXQYY = X1011+X(22)+X(45)
            QXXQZZ = X1011+X2324+X(46)+X(47)
            QZZQXX = X(10)+X1213+X(22)+X(48)+X(49)
            QZZQZZ = X(10)+X1213+X2324+X(50)+X(51)+X(52)+X(53)
            DXQXZ  = X(54)+X(55)+X(56)+X(57)
            QXZDX  = X(58)+X(59)+X(60)+X(61)
            QXZQXZ = 0.0D0
            DO  I=62,69
            QXZQXZ = QXZQXZ+X(I)
            ENDDO
         ELSE
            TWODA  = DA+DA
            X(14)  = R2+ADD
            X(15)  = RP2QA2+TWOQAQ
            X(16)  = RM2QA2+TWOQAQ
            X(17)  = R2+TWODA**2+ADD
            X(18)  = (R-TWODA)**2+ADD
            X(19)  = (R+TWODA)**2+ADD
            X(20)  = RPDA2+TWOQA2+ADQ
            X(21)  = RMDA2+TWOQA2+ADQ
            X(22)  = (R+DA-TWOQA)**2+ADQ
            X(23)  = (R-DA-TWOQA)**2+ADQ
            X(24)  = (R+DA+TWOQA)**2+ADQ
            X(25)  = (R-DA+TWOQA)**2+ADQ
            X(26)  = R2+FOUR*TWOQA2+AQQ
            X(27)  = R2+TWOQA2+TWOQAQ
            X(28)  = (R+TWOQA+TWOQA)**2+AQQ
            X(29)  = (R-TWOQA-TWOQA)**2+AQQ
            RMQA2  = (R-QA)**2
            RPQA2  = (R+QA)**2
            DMADD  = (DA-QA)**2+ADQ
            DPADD  = (DA+QA)**2+ADQ
            X(30)  = RMQA2+DMADD
            X(31)  = RPQA2+DMADD
            X(32)  = RMQA2+DPADD
            X(33)  = RPQA2+DPADD
            DO  I=1,33
            X(I)   = PXX(I)/SQRT(X(I))
            ENDDO
            EE     = X(1)
            DZE    = X(3) +X(4)
            QZZE   = X(2) +X(5) +X(6)
            QXXE   = X(2) +X(7)
            EDZ    =-DZE
            EQZZ   = QZZE
            EQXX   = QXXE
            DXDX   = X(14)+X(17)
            DZDZ   = X(14)+X(18)+X(19)
            X89    = X(8) +X(9)
            DZQXX  = X89  +X(20)+X(21)
            QXXDZ  =-DZQXX
            DZQZZ  = X89  +X(22)+X(23)+X(24)+X(25)
            QZZDZ  =-DZQZZ
            X1010  = X(10)+X(10)*HALF
            X1111  = X(11)+X(11)
            X1213  = X(12)+X(13)
            QXXQXX = X1010+X1111+X(26)
            QXXQYY = X1111+X(10)+X(27)
            QXXQZZ = X(10)+X(11)+X1213+X(15)+X(16)
            QZZQXX = QXXQZZ
            QZZQZZ = X1010+X1213+X1213+X(28)+X(29)
            DXQXZ  = X(30)+X(31)+X(32)+X(33)
            QXZDX  =-DXQXZ
            QXZQXZ = QXXQZZ
         ENDIF
         A(1)  = EE
         A(2)  = DZE
         A(3)  = EE + QZZE
         A(4)  = EE + QXXE
         A(5)  = EDZ
         A(6)  = DZDZ
         A(7)  = DXDX
         A(8)  = EDZ + QZZDZ
         A(9)  = EDZ + QXXDZ
         A(10) = QXZDX
         A(11) = EE  + EQZZ
         A(12) = EE  + EQXX
         A(13) = DZE + DZQZZ
         A(14) = DZE + DZQXX
         A(15) = DXQXZ
         A(16) = EE + EQZZ + QZZE + QZZQZZ
         A(17) = EE + EQZZ + QXXE + QXXQZZ
         A(18) = EE + EQXX + QZZE + QZZQXX
         A(19) = EE + EQXX + QXXE + QXXQXX
         A(20) = QXZQXZ
         A(21) = EE + EQXX + QXXE + QXXQYY
         A(22) = HALF*(A(19)-A(21))
      ENDIF
      RETURN
      END











