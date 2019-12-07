      SUBROUTINE  core_valence_charge(NAT,KAN,VAL,CORE,NTYP,DK)
C     THIS SUBROUTINE GENERATES THE MVS = MINIMUM VALENCE SET BASIS.
C     DK = NUMBER OF VALENCE ELECTRONS, NTYP = 1, 2, 3, 4, 5,
C     FOR S/ S,P / NS,NP,ND / NS, NP, (N-1)D/ NS,NP,(N-1)D,(N-2)F/
C     RESPECTIVELY
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer NAT
      integer NTYP(NAT)
      integer NTYPE
      double precision VAL(NAT)
      double precision CORE(NAT)
      integer KAN(*)
      double precision DK(NAT),X(NAT)
      integer I
      double precision  KA
C     WRITE(6,2222) (DK(I), I=1,NAT), (Z12(I),I=1,NAT)
C     WRITE(6,2244) (KAN(I),I=1,NAT)
C     WRITE(6,2244) (NTYP(I),I=1,NAT)
C2244 FORMAT (' DEFMVB',/,7I10)
C2222 FORMAT( ' DEFMVB',/,7F10.6)
       NTYPE = 0
       DO I=1, NAT
           DK(I)=0.0d0
           VAL(I)=0.0d0
           NTYP(I)= 0
           CORE(I)=0.0d0
           X(I)=0.0d0
        ENDDO
c       write(6,*) "i am here for valence electron"        
        DO 500 I = 1, NAT
c       write(6,*) NAT,DK(I)
        KA = KAN(I)
        IF (KA.EQ.0) GO TO 500
        IF (KA.GT.2) GO TO 10
        DK(I) = KA
        X(I)=0.0
        NTYPE = 1
        GO TO 140
   10   IF (KA.GT.10) GO TO 20
c       write(6,*) KA
        DK(I) = KA-2
        X(I)= 2
        NTYPE = 2
        GO TO 140
   20   IF (KA.GT.18) GO TO 30
        X(I)= 10
        DK(I) = KA-10
        NTYPE = 2
        GO TO 140
   30   IF (KA.GT.30) GO TO 40
C       FIRST TRANSITION SERIES
        DK(I) = KA-18
        NTYPE = 4
        X(I)= 18
        IF (KA.EQ.19) NTYPE = 2
C       CHECK FOR ZN SP OR SPD BASIS
c       IF (NTYP(I).LT.4.AND.NTYP(I).NE.0) DK(I) = VAL(I)
        GO TO 140
   40   IF (KA.GT.36) GO TO 50
        DK(I) = KA-28
        NTYPE = 2
        X(I)= 28
        GO TO 140
  
   50   IF (KA.GT.48) GO TO 60
C       SECOND TRANSITION SERIES
        DK(I) = KA-36
        NTYPE = 4
        X(I)= 36
         IF (KA.EQ.37) NTYPE = 2
c        IF(NTYP(I).LT.4.AND.NTYP(I).NE.0) DK(I) = VAL(I) 
        GO TO 140

   60   IF (KA.GT.54) GO TO 70
        X(I)= 46
        DK(I) = KA-46
        IF (KA.EQ.48 .AND. NTYP(I).GT.2) DK(I) = 12.0
        NTYPE = 2
C       IF (KA.EQ.53) GO TO 2000
        GO TO 140

   70   CONTINUE
        IF (KA.GT.56) GO TO 80
C       CS AND BA
        DK(I) = KA-54
        X(I)= 54
        IF(KA.EQ.55) NTYPE = 2
        IF(KA.EQ.56) NTYPE = 4
        GO TO 140
C
   80   IF (KA.GT.71) GO TO 90
C       LANTHANIDES
        X(I)= 54
        DK(I) = KA-54
        NTYPE = 5
        GO TO 140

   90   IF (KA.GT.80) GO TO 100
C       THIRD TRANSITION SERIES
        DK(I) = KA-68
        x(I)= 68
        NTYPE = 4
c        IF(NTYP(I).LT.4.AND.NTYP(I).NE.0) DK(I) = VAL(I) 
        GO TO 140

  100   IF (KA.GT.85) GO TO 110
        DK(I) = KA-78
        X(I)= 78
        NTYPE = 2
        GO TO 140

  110   IF (KA.GT.88) GO TO 120
        DK(I) = KA-86
        X(I)= 86
        NTYPE = 2
        GO TO 140

  120   IF (KA.GT.105) then
         CALL  errquit('semiemp_energy: input error',74, GEOM_ERR)
        DK(I) = KA-86
        X(I)= 86
        NTYPE = 5
        GO TO 140
        endif
C 
  140   CONTINUE
        IF(NTYP(I).EQ.0) NTYP(I) = NTYPE
  500   CONTINUE
C     WRITE(6,2222) (DK(I), I=1,NAT), (Z12(I),I=1,NAT)
C     WRITE(6,2244) (NTYP(I),I=1,NAT)
      DO I = 1, NAT
c     write(6,*) "NTYP(I)",NTYP(I)
      VAL(I) = DK(I)
      CORE(I)= X(I)
      ENDDO
      RETURN
      END
