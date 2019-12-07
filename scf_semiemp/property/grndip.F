      SUBROUTINE GRNDIP(P,E,A,B,BBB,BB,FX,FY,FZ,RAD,ITDIP,NA,NU,NRR,ISW1
     .     ,NIN,CO,Z12,KAN,NB,NW,PME,orbspace,got_orbspace,got_elh,IHOMO
     .     ,ILUMO,IORB,NEL,ehinfo)
C
C     THIS SUBROUTINE CALC. THE DIPOLE MATRIX ELEMENTS FOR A GIVEN
C     STATE.
C
C     ITDIP=0, THEN ALL COMPONENTS OF THE DIPOLE ARE CALCULATED
C              SEPARATELY
C          =1, THEN THE THREE TOTAL COMPONENTS ARE CALCULATED
C              USING ONLY THE LARGE MATRICES FZ AND P. THE DIPOLES
C              COMPONENTS RETURN AS E(3) IN A.U.
C
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
c     INCLUDE '../include/const.cmn'

c     INCLUDE '../include/rhf.cmn'
C****************************************************************
C***      COMMON/RHF/
C***     .  AR(21),BR(21),FOP(8),NOP,NDN,MIM(14),ITAPE
C****************************************************************
c     INCLUDE '../include/switch.cmn'
C****************************************************************
C***      COMMON/SWITCH/
C***     .  ISWANN,III,IPRINT,IPUN,ISAVE,ISPEC,ISW1,ISW2,ISW3,ISW4,
C***     .  JSW(10),IBETA,IDD1,IDD2
C****************************************************************
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     .  ATOM(100),NIN(MXB),
C***     .  NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     .  NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     .  Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      integer orbspace(2)
      double precision PME(3*IORB)
      DIMENSION   EX(NA),EY(NA),EZ(NA)

      DIMENSION    E(3), A(3), B(3), BBB(3), BB(3), FX(*), FY(*), FZ(*),
     .             P(*), RAD(NA,NA), NIN(*), CO(3,NA), Z12(*),NU(*)
      DIMENSION    LABEL(3),KAN(*),NW(*), MPG(3),IHOMO(*)
      DIMENSION    ILUMO(*)
ctt
      DIMENSION    TOTDIP(3),COM(3)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER (DEBYE    = 2.541748047D0,BOHRS    = 5.29177249D-1)
C
      CHARACTER*8  LABEL
      DATA         LABEL/'       X', '       Y', '       Z'/
      integer ehinfo(5)
      LOGICAL      FFLG, got_orbspace, got_elh
C
c     write(6,*)"i am in grndip"!,E!,FELEC,FFFF,ITDIP
      FELEC = ONE  
      FNUCL = ONE  
c     NPROP = IO_MOMENTS
C     RELATIVE TO THE CENTER OF MASS!
      CALL CMS(COM(1),COM(2),COM(3),CO,KAN,NA)
c     iF (JSW(6).EQ.0) GO TO 40
C     FOR SELF CONSISTENT REACTION FIELD
c     NELN = 0
c     DO 10 I = 1, NA
c        NELN = NELN + Z12(I)
c  10 CONTINUE
c     ICHARG = NELN - NEL
c     IF (ICHARG) 20, 40, 30
c  20 FELEC = DFLOAT(NEL+ICHARG)/DFLOAT(NEL)
c     GO TO 40

c  30 FNUCL = DFLOAT(NEL)/DFLOAT(NEL+ICHARG)
   40 CONTINUE
      FF1 = (DEBYE/BOHRS)
      DO 50 J = 1, 3
         E(J) = ZERO
         A(J) = ZERO
         B(J) = ZERO
         BBB(J) = ZERO
         BB(J) = ZERO
c        PME(J) =ZERO
   50 CONTINUE
      IF (got_orbspace) THEN
      Do I = 1, IORB
       DO J=1, 3
       PME(3*(I-1)+J) = ZERO
       END DO
      END DO
      END IF
C
C     NOW FIND THE FIRST ORDER DENSITY
C
c     IF(ISW1.NE.8) THEN
c       CALL FINDP(P,NRR,ISWANN,NDN,1)
c     ELSE
C       FOR XINDO THIS IS THE TOTAL DENSITY IN THE STO BASIS
c       REWIND (IO_XINDO)
c       READ (IO_XINDO) (P(I) , I = 1,NRR)
c     ENDIF

c     WRITE(6,*) ' GRNDIP - DENSITY '
C     CALL PRNT(P,NB,NRR)
C     REWIND NPROP
c     CALL ZREWIND (IO_MOMENTS, 'GRNDIP')
      IF (ITDIP.NE.0) GO TO 140
C--   DECOMPOSES THE DIPOLE MOMENT- UNITS ARE IN DEBYE
CMCZ  FFLG=.TRUE.
c     FFLG = .FALSE.
c     CALL RDPRP(NPROP,LABEL(1),FX,NRR,FFLG,'GRNDIP')
c     CALL RDPRP(NPROP,LABEL(2),FY,NRR,FFLG,'GRNDIP')
c     CALL RDPRP(NPROP,LABEL(3),FZ,NRR,FFLG,'GRNDIP')
      DO 58 I = 1,NA
         EX(I) = ZERO
         EY(I) = ZERO
         EZ(I) = ZERO
   58 CONTINUE
      DO 60 I = 1, NB
         NUIII= NU(I)
C
C        II = NIN(I) + I has been moved up here now.			!CK1105
         II = NIN(I) + I                                                !CK1105
         EX(NUIII) = EX(NUIII) + P(II)*FX(II)*FELEC
         EY(NUIII) = EY(NUIII) + P(II)*FY(II)*FELEC
         EZ(NUIII) = EZ(NUIII) + P(II)*FZ(II)*FELEC
C
C     E(J) IS THE CHARGE CONTRIBUTION TO THE GROUND STATE DIPOLE.
C        II = NIN(I) + I						!CK1105
         E(1) = E(1) + P(II)*FX(II)*FELEC
         E(2) = E(2) + P(II)*FY(II)*FELEC
         E(3) = E(3) + P(II)*FZ(II)*FELEC
   60 CONTINUE
      DO 70 I = 1, NA
         FFFF = Z12(I)*FF1*FNUCL
         E(1) = E(1) + FFFF*(CO(1,I)*BOHRS-COM(1))
         E(2) = E(2) + FFFF*(CO(2,I)*BOHRS-COM(2))
         E(3) = E(3) + FFFF*(CO(3,I)*BOHRS-COM(3))
         EX(I) = EX(I) + FFFF*(CO(1,I)*BOHRS-COM(1))
         EY(I) = EY(I) + FFFF*(CO(2,I)*BOHRS-COM(2))
         EZ(I) = EZ(I) + FFFF*(CO(3,I)*BOHRS-COM(3))
c      IF (got_orbspace) THEN
c       DO J = 1, 105
c        PME(1,J) = PME(1,J) + FFFF*(CO(1,I)*BOHRS-COM(1))
c        PME(2,J) = PME(2,J) + FFFF*(CO(2,I)*BOHRS-COM(2))
c        PME(3,J) = PME(3,J) + FFFF*(CO(3,I)*BOHRS-COM(3))
c       END DO
c      END IF
   70 CONTINUE
       L = 0
      DO 120 J = 2, NB
         JA = NU(J)
         JB = NW(J)
         JD = NIN(J)
         JS = J - 1
         DO 115 K = 1, JS
            KA = NU(K)
            KJ = JD + K
            KAJA = NIN(JA) + KA
            IF (ABS(RAD(KA,JA)).GT.0.001) GO TO 100
            IF (JB.GT.8 .OR. NW(K).GT.8) GO TO 90
            IF (JB.GT.3 .OR. NW(K).GT.3) GO TO 80
            A(1) = A(1) + P(KJ)*FX(KJ)
            A(2) = A(2) + P(KJ)*FY(KJ)
            A(3) = A(3) + P(KJ)*FZ(KJ)
            GO TO 110

   80       B(1) = B(1) + P(KJ)*FX(KJ)
            B(2) = B(2) + P(KJ)*FY(KJ)
            B(3) = B(3) + P(KJ)*FZ(KJ)
            GO TO 110

   90       BBB(1) = BBB(1) + P(KJ)*FX(KJ)
            BBB(2) = BBB(2) + P(KJ)*FY(KJ)
            BBB(3) = BBB(3) + P(KJ)*FZ(KJ)
            GO TO 110

  100       BB(1) = BB(1) + P(KJ)*FX(KJ)
            BB(2) = BB(2) + P(KJ)*FY(KJ)
            BB(3) = BB(3) + P(KJ)*FZ(KJ)
  110       CONTINUE
c           write(6,*)"P,FX,FY,FZ",P(KJ),FX(KJ),FY(KJ),FZ(KJ)
            EX(KA) = EX(KA) + P(KJ)*FX(KJ)
            EY(KA) = EY(KA) + P(KJ)*FY(KJ)
            EZ(KA) = EZ(KA) + P(KJ)*FZ(KJ)
            IF (got_orbspace) THEN
             IF (got_elh) THEN
              IF ((K.GE.orbspace(1).AND.K.LE.ehinfo(4)).AND.
     1   (J.GE.ehinfo(1).AND.J.LE.orbspace(2)).AND.(K.NE.J)) THEN
             L = L + 1
             IHOMO(L) = K
             ILUMO(L) = J
c            write(6,*)"L,K,J",L,K,J,IHOMO(L),ILUMO(L)
c            write(6,*)"HERE,KJ",KJ,P(KJ),FZ(KJ)
             PME((L-1)*3+1) = PME((L-1)*3+1) + P(KJ)*FX(KJ)
             PME((L-1)*3+2) = PME((L-1)*3+2) + P(KJ)*FY(KJ)
             PME((L-1)*3+3) = PME((L-1)*3+3) + P(KJ)*FZ(KJ)
            ELSE IF ((K.GE.orbspace(1).AND.K.LE.NEL/2).AND.
     1   (J.LE.orbspace(2).AND.J.GT.NEL/2)) THEN
             L = L + 1
             IHOMO(L) = K
             ILUMO(L) = J
c            write(6,*)"L,K,J",L,K,J,IHOMO(L),ILUMO(L)
c            write(6,*)"HERE,KJ",KJ,P(KJ),FZ(KJ)
             PME((L-1)*3+1) = PME((L-1)*3+1) + P(KJ)*FX(KJ)
             PME((L-1)*3+2) = PME((L-1)*3+2) + P(KJ)*FY(KJ)
             PME((L-1)*3+3) = PME((L-1)*3+3) + P(KJ)*FZ(KJ)
              END IF
            ELSE IF ((K.GE.orbspace(1).AND.K.LE.NEL/2).AND.
     1   (J.LE.orbspace(2).AND.J.GT.NEL/2)) THEN
             L = L + 1
             IHOMO(L) = K
             ILUMO(L) = J
c            write(6,*)"L,K,J",L,K,J,IHOMO(L),ILUMO(L)
c            write(6,*)"HERE,KJ",KJ,P(KJ),FZ(KJ)
             PME((L-1)*3+1) = PME((L-1)*3+1) + P(KJ)*FX(KJ)
             PME((L-1)*3+2) = PME((L-1)*3+2) + P(KJ)*FY(KJ)
             PME((L-1)*3+3) = PME((L-1)*3+3) + P(KJ)*FZ(KJ)
             END IF
            END IF
  115 CONTINUE
  120 CONTINUE
c     write(6,*)"L",L

      DO 130 J = 1, 3
         A(J) = A(J)*TWO
         BB(J) = BB(J)*TWO
         BBB(J) = BBB(J)*TWO
         B(J) = B(J)*TWO
c        write(6,*)"A,BB,BBB,B",A(j),BB(J),BBB(J),B(J)
c        PME(J) = PME(J)*TWO
  130 CONTINUE
       IF (got_orbspace) THEN
        DO J= 1, 3
         DO I= 1, IORB
         PME(3*(I-1)+J)=TWO*PME(3*(I-1)+J)
         END DO
        END DO
       END IF
c     write(6,*)"HERE"
      RETURN

  140 CONTINUE
C--   ONLY THE THREE COMPONENTS OF THE TOTAL DIPOLE - IN A.U.(ELECTRON-B
c     ITER = 1
c     FFLG = .TRUE.
c 150 CONTINUE
c     CALL RDPRP(NPROP,LABEL(ITER),FZ,NRR,FFLG,'GRNDIP')
c     FFLG = .FALSE.
c     DO 160 I = 1, NB
C
C     E(J) IS THE CHARGE CONTRIBUTION TO THE GROUND STATE DIPOLE.
c        II = NIN(I) + I
c        E(ITER) = E(ITER) + P(II)*FZ(II)*FELEC
c 160 CONTINUE
c     DO 170 I = 1, NA
c        FFFF = Z12(I)*FF1*FNUCL
c        E(ITER) = E(ITER) + FFFF*(CO(ITER,I)-COM(ITER))*BOHRS
c 170 CONTINUE
c     DO 190 J = 2, NB
c        JD = NIN(J)
c        JS = J - 1
c        DO 180 K = 1, JS
c           KJ = JD + K
c           BB(ITER) = BB(ITER) + P(KJ)*FZ(KJ)
c 180    CONTINUE
c 190 CONTINUE
c     BB(ITER) = BB(ITER)*TWO
c     ITER = ITER + 1
c     IF (ITER.LE.3) GO TO 150
c     DO 200 I = 1, 3
c        E(I) = (E(I)+BB(I))/DEBYE
c 200 CONTINUE
ctt
c     IF (JSW(6).EQ.9 .OR. JSW(6).EQ.10.OR.JSW(6).EQ.19.OR.JSW(6).EQ.20)
c    .     THEN
c        FFLG = .TRUE.
c        CALL RDPRP(NPROP,LABEL(1),FX,NRR,FFLG,'GRNDIP')
c        CALL RDPRP(NPROP,LABEL(2),FY,NRR,FFLG,'GRNDIP')
c        CALL RDPRP(NPROP,LABEL(3),FZ,NRR,FFLG,'GRNDIP')
c        CALL MOLDIP(NA,KAN,CO,P,FX,FY,FZ,TOTDIP,ICHARG)
C           DIPOLES AND TOTDIP ARE IN DEBYE FROM MOLDIP
c     END IF
CGP
C      WRITE(6,2000) (EX(I),I=1,NA)
C      WRITE(6,2000) (EY(I),I=1,NA)
C      WRITE(6,2000) (EZ(I),I=1,NA)
CGP
2000  FORMAT(3F20.10)
      RETURN

      END
