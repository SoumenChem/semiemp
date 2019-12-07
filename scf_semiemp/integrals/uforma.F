
      SUBROUTINE UFORMA(NB,NA,UCORE,F,FA,RAD,LQNO,NIN,NU,NW,KAN,PO,DAD
     1     ,Z12)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/z.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/const.cmn'
c     INCLUDE '../include/datain.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0,
     +  THREE =3.0D0,FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0,
     +  SEVEN = 7.0D0,EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0)
      DIMENSION F(*),FA(*), RAD(NA,NA),H1EC(10), PO(9,105), DAD(9,105)
      DIMENSION E(8),D(8),UCORE(4,86),LQNO(*),NIN(*),NU(*),KAN(*),NW(*)
     1 ,Z12(*)
C
      TF = 0.0d0
      TD = 0.0d0
C     ADD THE CORE INTEGFALS
c     write(6,*)"i am in uforma",DAD(2,8),DAD(3,8)
      DO I=1,NB
        L=LQNO(I)
        ID=NIN(I)+I
        FA(ID) = ZERO
c       write(6,*)"FA,ID,KAN",FA(ID),ID,KAN(NU(I))
        FA(ID)=FA(ID)+UCORE(L+1,KAN(NU(I)))
c       WRITE(6,*) ' UFORMA I ,L, UCORE' !, I, LQNO(I),
C    x          UCORE(L+1,KAN(NU(I)))
      ENDDO  
C
c     CALL GEOMAT(F,RAD)
C
      NA1 = NA*NA
      NA2 = NA1 + NA1
      KIF=0
      DO 200 K=1,NA
      DO 100 I=1,NB
        LI=NW(I)
        KI=NU(I)
        IF(K.EQ.KI) GO TO 100
        IF(KI.GT.K)THEN
          KIKK=NIN(KI)+K
        ELSE
          KIKK=NIN(K)+KI
        ENDIF
C
        IF(KIF.NE.KI)THEN
            KIF=KI
c           RAU=RAD(KIKK)/BOHRS
            RAU=RAD(KI,K)
c          write(6,*)"RAU",RAU
            CALL IJCORE(KAN(KI),KAN(K),RAU,H1EC,PO,DAD)
        ENDIF
        DO 90 J=I,NB
           IJ = NIN(J) + I
           KJ=NU(J)
           IF(KI.NE.KJ)GO TO 95
           LJ=NW(J)
           KORM=KI+(K-1)*NA
c       write(6,*)"KORM",KORM,KORM+NA1,KORM+NA2
           U1=F(KORM)
           U2=F(KORM+NA1)
           U3=F(KORM+NA2)
c       write(6,*)"U1,U2,U3,LJ",U1,U2,U3,LJ
           CALL TRANSMAT(U1,U2,U3,LJ,D)
           CALL TRANSMAT(U1,U2,U3,LI,E)
c       write(6,*)"D,E",D(1),E(1)
C          MM=LANGL(NW(J)+1)+1
           MM = LQNO(J)+ 1
           GO TO (27,28,29,30), MM
C
  30       EPHI= ONE
           EDEL=-ONE
C
  28       EPI= ONE
           ESIG=-ONE
           GOTO 31
C
  29       EDEL= ONE
           EPI =-ONE
C
  27       ESIG= ONE
  31       CONTINUE
C..........................................
C     THE WAY THE ONE CENTER IS ARRANGED THE LOWER L VALUE
C     IS ALWAYS FIRST
 
      MINKM= LQNO(I)  + 1
C          SIGMA PART
           IND=LQNO(I)+LQNO(J)+1
c     write(6,*)"E,D,H1EC,Z12,IND",E(1),D(1),H1EC(IND),Z12(k),IND
           FA(IJ)=FA(IJ)+E(1)*D(1)*H1EC(IND)*Z12(K)
c     write(6,*)"here 1"
c     write(6,*)"FA",FA(IJ)   
C            WRITE (6,9012)
C9012       FORMAT(' IJ  IND, E,           D,       LI      LJ')
c           WRITE(6,9011) IJ ,IND, E(1), D(1),LQNO(I),LQNO(J)
c9011       FORMAT(2I5,2F10.6,2I5)
C
C  CHECK HOW MANY TERMS IN OVERLAP
C
      FACTOR = ONE
C     IF(I.EQ.J) FACTOR = HALF
      GO TO (50,47,48,49),MINKM
C
C  PHI   
C
   49 CONTINUE   
c     write(6,*) "here 2"
      FA(IJ)=(E(6)*D(6)+E(7)*D(7))*TF*FACTOR +FA(IJ)
c     write(6,*)"FA",FA(IJ)
C
C  DELTA  
C
   48 CONTINUE 
c     write(6,*)"here 3"
      FA(IJ)=(E(4)*D(4)+E(5)*D(5))*TD*FACTOR +FA(IJ)
c     write(6,*)"FA",FA(IJ)
      WRITE(6,*) ' ***** ERROR IN UFORMA *****'
      CALL ZEXIT(3) 
C
C  PI (X**2 AND Y**2)
C
   47 CONTINUE
      IND = 4
c     write(6,*)"here 4"
      FA(IJ)=(E(2)*D(2)+E(3)*D(3))*H1EC(IND)*Z12(K)*FACTOR + FA(IJ)
c     write(6,*)"FA",FA(IJ)
C           WRITE(6,9014) IJ,IND,E(2), D(2),E(3),D(3)
C9014       FORMAT(2I5,4F10.6)
C
  50       CONTINUE
  90       CONTINUE
  95       CONTINUE
 100       CONTINUE
 200       CONTINUE
c     write(6,*)"F in uforma",(F(I),I=1,nb*(nb+1)/2)
           
      RETURN
      END
C------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
      SUBROUTINE IJCORE(NI,NJ,R,H1EC,PO,DAD)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/am1.cmn'
C     *
C     CALCULATE THE NUCLEAR ATTRACTION INTEGRALS IN LOCAL COORDINATES.
C     CALCUALTION OF (u(A)v(A)|1/R(B)) FOR AM1.  NI IS THE ATOMIC
C     NUMBER FOR ATOM A, NJ IS THE ATOMIC NUMBER OF ATOM B
C     R in AU
C     H1EC IN AU
C     H1EC (*): FIRST INDEX 1,2,3,4 corresponds to ss, sp,pp-sig, pp-del
C     of NI
C                           5,6,7,8  sd-sig, pd-sig, dd-sig,pd-pi
C                           9,10, dd-pi, dd-del 
C     *  

C     *
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0,
     +  THREE =3.0D0,FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0,
     +  SEVEN = 7.0D0,EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0)
      DIMENSION PXY(7),X(7),H1EC(10), PO(9,105), DAD(9,105)
      DATA PXY/ 1.0 D0,-0.5 D0,-0.5 D0, 0.5 D0, 0.25D0,
     .          0.25D0, 0.5D0 /
C *** INITIALIZATION.
      SAVE PXY
c     write(6,*)"i am in ijcore,NI,NJ",NI,NJ,DAD(2,8),DAD(3,8)
      DO   I=1,10
      H1EC(I)=ZERO
      ENDDO
      IF (ABS(R).LT.1.0D-7) RETURN
      R2     = R*R
      POCORE=PO(9,NJ)
C     WRITE(6,9010) ' IJCORE', NI, PO(1,NI), PO(2,NI), PO(3,NI),
C    x       PO(7,NI), PO(9,NI)
C     WRITE(6,9010) ' IJCORE', NI, DAD(1,NI), DAD(2,NI), 
C    x        DAD(3,NI), DAD(4,NI)
C     WRITE(6,9010) ' IJCORE', NJ, PO(1,NJ), PO(2,NJ), PO(3,NJ), 
C    x          PO(7,NJ), PO(9,NJ)
C     WRITE(6,9010) ' IJCORE', NJ, DAD(1,NJ), DAD(2,NJ), 
C    x    DAD(3,NJ), DAD(4,NJ)
C9010 FORMAT(1X,A,I5,5F10.6)

C *** SS -CORE INTERACTION
         SS    = (POCORE +PO(1,NI))**2
         H1EC(1) =-ONE/SQRT(R2+SS)
c      write(6,*)"H1EC",H1EC(1)
      IF(NI.LT.3) RETURN
C *** NI -  HEAVY ATOM
      IF(NI.GE.3 ) THEN
         PP     = (POCORE +PO(7,NI))**2
         DA     = DAD(2,NI)
         QA     = DAD(3,NI)/SQRT(TWO)
c      write(6,*)"DA,QA",DA,QA,NI
         TWOQA  = QA+QA
         AD    = (PO(2,NI)+ POCORE)**2
         AQ    = (PO(3,NI)+ POCORE)**2
         X(1)   = R2+PP
         X(2)   = R2+AQ
         X(3)   = (R+DA)**2+AD
         X(4)   = (R-DA)**2+AD
         X(5)   = (R-TWOQA)**2+AQ
         X(6)   = (R+TWOQA)**2+AQ
         X(7)   = R2+TWOQA*TWOQA+AQ
         DO   I=1,7
            X(I)   = PXY(I)/SQRT(X(I))
         ENDDO
         H1EC(2)= -(X(3)+X(4))
         H1EC(3)= -(X(1)+X(2)+X(5)+X(6))
         H1EC(4)= -(X(1)+X(2)+X(7))
c     write(6,*)"H1EC",H1EC(1),H1EC(2),H1EC(3)
      ENDIF
      RETURN
      END

