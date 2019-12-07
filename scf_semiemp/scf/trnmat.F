      SUBROUTINE TRNMAT(F,EIG,NB,NRR,BB)
C
C     ALLOCATE A TEMPORARY MATRIX
C      
C     INCLUDE '../include/param.cmn'
C     INCLUDE '../include/mem.cmn'
      IMPLICIT DOUBLE PRECISION (A-F, O-Z)
      IMPLICIT INTEGER (I-N)
      parameter (ZERO = 0.0d0, ONE = 1.0d0)
C-----------------------------------------------------------------------
C     THIS SUBROUTINE FORMS THE FOCK TRANSFORMATION MATRIX, S**-1/2, AND
C     STORES IT IN F.
C-----------------------------------------------------------------------

C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL      IBGIV
      DIMENSION    F(*), BB(*), EIG(*)
C     SAVE D, EK
C-----------------------------------------------------------------------
c     write(6,*)"i am in trnmat"
      IBGIV = .TRUE.
c     BB= ZERO
c     CALL GIVENS(F,BB,EIG,NB,NB,IBGIV,QQ,NIN)
c     write(6,*)"BB",(BB(I),I=1,20)      
      DO 10 I = 1, NB
         IF (EIG(I).LT.1.0D-10) THEN
            WRITE (6,9000) EIG(I)
            EIG(I) = 1.0D-10
         END IF
         EIG(I) = ONE/SQRT(EIG(I))
   10 CONTINUE

      DO 20 IJ = 1, NRR
         F(IJ) = ZERO
   20 CONTINUE
c     write(6,*)"here"
      KK = 0
      DO 50 K = 1, NB
         IJ = 0
         DO 40 I = 1, NB
            XX = BB(KK+I)*EIG(K)
            DO 30 J = 1, I
               IJ = IJ + 1
c     write(6,*)"F(IJ),BB,XX",F(IJ),BB(KK+J),XX
               F(IJ) = F(IJ) + BB(KK+J)*XX
c     write(6,*)"F(IJ)",F(IJ)
   30       CONTINUE
   40    CONTINUE
         KK = KK + NB
   50 CONTINUE
c     write(6,*)"F",(F(I),I=1,nrr)
      RETURN
 9000 FORMAT ('POSSIBLE LINEAR DEPENDENCY IN BASIS',' S EIGENVALUE =',
     .       D16.10)
      END
c-------------------------------------------------------------------------------
c-
c------------------------------------------------------------------------------
      SUBROUTINE GIVENS(A,VEC,W,N,NEIG,IBGIV,QQ,NIN)

c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      parameter (MXCI = 6000)
      LOGICAL      IBGIV
      DIMENSION    A(*), VEC(N,NEIG), W(*), QQ(*)
      DIMENSION    E2(N),IND(N),RV1(N),RV2(N),RV3(N)
      DIMENSION    RV4(N),RV6(N),D(N), EK(N), IQ(N)
      DIMENSION    NIN(*)
c     write(6,*)"i am in givens"
      IF(N.LE.2) THEN
         CALL HDIAG(N,N,-1,VEC,A,QQ,NIN,IQ)
         II = 0
         DO 10 I = 1, N
            II = II + I
            W(I) = A(II)
   10    CONTINUE
         RETURN

      END IF
      NV = (N*N+N)/2
      ND = N
      WRITE(6,*)ND
      write(6,*)N
      write(6,*)NEIG
c     write(6,*)"N befor tred3",N
c     write(6,*)"D befor tred3",(D(I),I=1,10)
c     write(6,*)"E befor tred3",(EK(I),I=1,10)
      CALL TRED3_SG(N,NV,A,D,EK,E2)
c     write(6,*)"N befor imtqlv",N
c     write(6,*)"D befor imtqlv",(D(I),I=1,10)
c     write(6,*)"E befor imtqlv",(EK(I),I=1,10)
      CALL IMTQLV(N,D,EK,E2,W,IND,IERR,RV1)
      
      IF(IERR.NE.0) THEN
         WRITE (6,9000)

 9000    FORMAT (/////,5X,'DIAGONALIZATION FAILED, IMTQLV CHECK MACHEP',
     .          /////)

         CALL ZEXIT(3)
      END IF
c     write(6,*)"E",E
      CALL TINVIT(ND,N,D,EK,E2,NEIG,W,IND,VEC,IERR,RV1,RV2,RV3,RV4,RV6)
      IF(IERR.NE.0) THEN
         WRITE (6,9010)

 9010    FORMAT (/////,5X,'DIAGONALIZATION FAILED, TOO MANY ITERATIONS',
     .          /////)

c        CALL ZEXIT(3)
      END IF
      CALL TRBAK3(ND,N,NV,A,NEIG,VEC)
c     write(6,*)"VEC",VEC
      RETURN

 9020 FORMAT ('   EISPTR  ND :',I5,'  N:',I5,'  NEIG :',I5)

      END
c--------------------------------------------------------------
c
c--------------------------------------------------------------
      SUBROUTINE HDIAG(N,IEGEN,ISW,U,H,QQ,NIN,IQ)
C     IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MXCI =  6000,
     + ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
   
      dimension  QQ(*),NIN(*)
CAT     INTEGER*2 LABI,LABJ,LABK,LABL
C     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/z.cmn'
C
      DIMENSION    IQ(N),U(N,N), H(1)
C
C     WE USE A SECOND DIAGONALIZATION ROUTINE HERE SO AS TO LEAVE THE
C     ORDER OF THE ATOMIC ORBITALS UNPERMUTED WHEN THE IMPROVED
C     NATURAL ORBITALS ARE COMPUTED.  THIS ROUTINE IS SLOW BUT THE
C     INCOMING MATRIX IS ALMOST DIAGONAL ANYWAY, WHILE GVNS SERIOUSLY
C     PERMUTES THE ORDER OF THE ATOMIC ORBITALS EVEN THOUGH IT IS MUCH
C     FASTER FOR THE SCF AND CI ENERGY DIAGONALIZATIONS. MATRIX TO BE
C     DIAGONALIZED IS H (VALUES ON DIAG) VECTORS IN COL. IN U.
C
c9000 FORMAT (10F12.6)
C
      L=0
c     write(6,*) "i am in hdiag",L
c     write(6,*)"N",N
c     write(6,*)"U",U
c     write(6,*)"H",(H(I),I=1,10)
c     DO  I = 1, MXCI
c      NIN(I) = L
c      L = L + I
c     END DO
      RAP = 1.0D-12
   10 IF (IEGEN) 20, 70, 20
   20 DO 60 I = 1, N
         DO 50 J = 1, N
            IF (I-J) 40, 30, 40
   30       U(I,J) = ONE
            GO TO 50

   40       U(I,J) = ZERO
   50    CONTINUE
   60 CONTINUE
c     write(6,*)"i am here 1"
   70 NNR = 0
      IF (N-1) 450, 450, 80
   80 NMI1 = N - 1
      DO 110 I = 1, NMI1
         QQ(I) = ZERO
         IQ(I) = N
         IPL1 = I + 1
         DO 100 J = IPL1, N
            IJ = NIN(J) + I
            ABHIJ = ABS(H(IJ))
            IF (QQ(I)-ABHIJ) 90, 100, 100
   90       QQ(I) = ABHIJ
            IQ(I) = J
  100    CONTINUE
  110 CONTINUE
  120 XMAX = QQ(1)
      IPIV = 1
      JPIV = IQ(1)
      DO 140 I = 1, NMI1
         IF (XMAX-QQ(I)) 130, 140, 140
  130    XMAX = QQ(I)
         IPIV = I
         JPIV = IQ(I)
  140 CONTINUE
c     write(6,*)"i am here2"
      IF (XMAX-RAP) 450, 450, 150
  150 NNR = NNR + 1
      IPIP = NIN(IPIV) + IPIV
c     write(6,*)"JPIV",JPIV
c     write(6,*)"NIN(JPIV)",NIN(JPIV)
      HII = H(IPIP)
      IPJP = NIN(JPIV) + IPIV
      JPJP = NIN(JPIV) + JPIV
c     write(6,*)"JPJP",JPJP
      HJJ = H(JPJP)
      HIJ = H(IPJP)
      HIMHJ = HII - HJJ
c     write(6,*)"HIMHJ,HIJ,FOUR",HIMHJ,HIJ,FOUR
      DENOM = ABS(HIMHJ) + SQRT((HIMHJ**2)+ (FOUR*HIJ*HIJ))
c     write(6,*)"i am here3",DENOM
      IF (DENOM) 170, 160, 170
  160 COSIN = ONE
      SINE = ZERO
      GO TO 220

  170 TANG = SIGN(TWO,HIMHJ)*HIJ/DENOM
      COSIN = ONE/SQRT(ONE+TANG**2)
      SINE = TANG*COSIN
      H(IPIP) = COSIN**2* (HII+TANG* (TWO*HIJ+TANG*HJJ))
      H(JPJP) = COSIN**2* (HJJ-TANG* (TWO*HIJ-TANG*HII))
      H(IPJP) = ZERO
C
C     IF ISW=0 WE DO NOT WISH TO ORDER THE EIGENVALUES AND VECTORS.
C           > 0 GREATEST FIRST
C           < 0 SMALLEST FIRST
C
      IF (ISW) 200, 220, 180
  180 IF (H(IPIP)-H(JPJP)) 190, 220, 220
  190 HTEMP = H(IPIP)
      H(IPIP) = H(JPJP)
      H(JPJP) = HTEMP
      HTEMP = SIGN(ONE,-SINE)*COSIN
      COSIN = ABS(SINE)
      SINE = HTEMP
      GO TO 220

  200 IF (H(IPIP)-H(JPJP)) 220, 220, 210
  210 HTEMP = H(IPIP)
      H(IPIP) = H(JPJP)
      H(JPJP) = HTEMP
      HTEMP = SIGN(ONE,-SINE)*COSIN
      COSIN = ABS(SINE)
      SINE = HTEMP
  220 CONTINUE
      DO 290 I = 1, NMI1
         IF (I-IPIV) 240, 290, 230
  230    IF (I-JPIV) 240, 290, 240
  240    IF (IQ(I)-IPIV) 250, 260, 250
  250    IF (IQ(I)-JPIV) 290, 260, 290
  260    K = IQ(I)
         IK = NIN(K) + I
         HTEMP = H(IK)
         H(IK) = ZERO
         IPL1 = I + 1
         QQ(I) = ZERO
         DO 280 J = IPL1, N
            IJ = NIN(J) + I
            IF (QQ(I)-ABS(H(IJ))) 270, 280, 280
  270       QQ(I) = ABS(H(IJ))
            IQ(I) = J
  280    CONTINUE
         H(IK) = HTEMP
  290 CONTINUE
      QQ(IPIV) = ZERO
      QQ(JPIV) = ZERO
      DO 420 I = 1, N
         IF (I-IPIV) 300, 420, 340
  300    IIP = NIN(IPIV) + I
         IJP = NIN(JPIV) + I
         HTEMP = H(IIP)
         H(IIP) = COSIN*HTEMP + SINE*H(IJP)
         IF (QQ(I)-ABS(H(IIP))) 310, 320, 320
  310    QQ(I) = ABS(H(IIP))
         IQ(I) = IPIV
  320    H(IJP) = -SINE*HTEMP + COSIN*H(IJP)
         IF (QQ(I)-ABS(H(IJP))) 330, 420, 420
  330    QQ(I) = ABS(H(IJP))
         IQ(I) = JPIV
         GO TO 420

  340    CONTINUE
         JPI = NIN(I) + JPIV
         IPI = NIN(I) + IPIV
         IJP = NIN(JPIV) + I
         IF (I-JPIV) 350, 420, 380
  350    CONTINUE
         HTEMP = H(IPI)
         H(IPI) = COSIN*HTEMP + SINE*H(IJP)
         IF (QQ(IPIV)-ABS(H(IPI))) 360, 370, 370
  360    QQ(IPIV) = ABS(H(IPI))
         IQ(IPIV) = I
  370    H(IJP) = -SINE*HTEMP + COSIN*H(IJP)
         IF (QQ(I)-ABS(H(IJP))) 330, 420, 420
  380    HTEMP = H(IPI)
         H(IPI) = COSIN*HTEMP + SINE*H(JPI)
         IF (QQ(IPIV)-ABS(H(IPI))) 390, 400, 400
  390    QQ(IPIV) = ABS(H(IPI))
         IQ(IPIV) = I
  400    H(JPI) = -SINE*HTEMP + COSIN*H(JPI)
         IF (QQ(JPIV)-ABS(H(JPI))) 410, 420, 420
  410    QQ(JPIV) = ABS(H(JPI))
         IQ(JPIV) = I
  420 CONTINUE
      IF (IEGEN) 430, 120, 430
  430 DO 440 I = 1, N
         HTEMP = U(I,IPIV)
c       write(6,*)"COSIN",COSIN
c       write(6,*)"HTEMP",HTEMP
c       write(6,*)"SINE",SINE
c       write(6,*)"U(I,JPIV)",U(I,JPIV)
         U(I,IPIV) = COSIN*HTEMP + SINE*U(I,JPIV)
         U(I,JPIV) = -SINE*HTEMP + COSIN*U(I,JPIV)
c       write(6,*)"U(I,I)",U(I,IPIV)
c       write(6,*)"U(I,J)",U(I,JPIV)
c     write(6,*)"i am here11"
  440 CONTINUE
      GO TO 120

  450 CONTINUE
c     write(6,*)"U",U
  460 RETURN

      END
c-----------------------------------------------------------------------------------
c
c----------------------------------------------------------------------------------
      SUBROUTINE TRED3_SG(N,NV,A,D,E,E2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C     INCLUDE '../include/param.cmn'
      Double precision   A(NV), D(N), E2(N),E(N)
C
C     THIS ROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS ROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING ROUTINE DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ---------- FOR I=N STEP -1 UNTIL 1 DO -- ----------
c     write(6,*)"i am in tred3"
c     write(6,*)"A",(A(I),I=1,NV )
      DO 120 II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = (I*L)/2
         H = 0.0D+00
         SCALE = 0.0D+00
         IF (L.LT.1) GO TO 20
C     ---------- SCALE ROW (ALGOL TOL THEN NOT NEEDED) ----------
         DO 10 K = 1, L
            IZ = IZ + 1
c        write(6,*)"A(IZ),IZ",A(IZ),IZ
            D(K) = A(IZ)
            SCALE = SCALE + ABS(D(K))
   10    CONTINUE
C
         IF (SCALE.NE.0.0D+00) GO TO 30
   20    E(I) = 0.0D+00
         E2(I) = 0.0D+00
         GO TO 110
C
   30    DO 40 K = 1, L
            D(K) = D(K)/SCALE
            H = H + D(K)*D(K)
   40    CONTINUE
C
         E2(I) = SCALE*SCALE*H
         F = D(L)
         G = -SIGN(SQRT(H),F)
c        write(6,*)"E(I)",E(I)
         E(I) = SCALE*G
c        write(6,*)"SCALE,G",SCALE,G
         H = H - F*G
         D(L) = F - G
         A(IZ) = SCALE*D(L)
         IF (L.EQ.1) GO TO 110
         F = 0.0D+00
C
         JK = 1
         DO 70 J = 1, L
            JM1 = J - 1
            DT = D(J)
            G = 0.0D+00
C     ---------- FORM ELEMENT OF A*U ----------
            IF (JM1.EQ.0) GO TO 60
            DO 50 K = 1, JM1
               E(K) = E(K) + DT*A(JK)
               G = G + D(K)*A(JK)
               JK = JK + 1
   50       CONTINUE
   60       E(J) = G + A(JK)*DT
            JK = JK + 1
C     ---------- FORM ELEMENT OF P ----------
   70    CONTINUE
         F = 0.0D+00
         DO 80 J = 1, L
            E(J) = E(J)/H
            F = F + E(J)*D(J)
   80    CONTINUE
C
         HH = F/ (H+H)
         JK = 0
C     ---------- FORM REDUCED A ----------
         DO 100 J = 1, L
            F = D(J)
            G = E(J) - HH*F
            E(J) = G
C
            DO 90 K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F*E(K) - G*D(K)
   90       CONTINUE
  100    CONTINUE
C
  110    D(I) = A(IZ+1)
         A(IZ+1) = SCALE*SQRT(H)
  120 CONTINUE
c     write(6,*)"E",E
c     write(6,*)"D",D
c     write(6,*)"D",(D(I),I=1,N )
C
      RETURN
C     ---------- LAST CARD OF TRED3 ----------
      END
c----------------------------------------------------------------------
C
c----------------------------------------------------------------------
      SUBROUTINE IMTQLV(N,D,E,E2,W,IND,IERR,RV1)
C
c     INCLUDE '../include/param.cmn'
      INTEGER      I, J, K, L, M, N, II, MML, TAG, IERR
      DOUBLE PRECISION D(*), E(*), E2(*), W(*), RV1(N)
      DOUBLE PRECISION B, C, F, G, P, R, S, MACHEP
C     DOUBLE PRECISION SQRT,ABS,DSIGN
      INTEGER      IND(N)
C
C     THIS SUBROUTINE IS A VARIANT OF  IMTQL1  WHICH IS A TRANSLATION OF
C     ALGOL PROCEDURE IMTQL1, NUM. MATH. 12, 377-383(1968) BY MARTIN AND
C     WILKINSON, AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC TRIDIAGONAL
C     MATRIX BY THE IMPLICIT QL METHOD AND ASSOCIATES WITH THEM
C     THEIR CORRESPONDING SUBMATRIX INDICES.
C
C     ON INPUT:
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2(1) IS ARBITRARY.
C
C     ON OUTPUT:
C
C        D AND E ARE UNALTERED;
C
C        ELEMENTS OF E2, CORRESPONDING TO ELEMENTS OF E REGARDED
C          AS NEGLIGIBLE, HAVE BEEN REPLACED BY ZERO CAUSING THE
C          MATRIX TO SPLIT INTO A DIRECT SUM OF SUBMATRICES.
C          E2(1) IS ALSO SET TO ZERO;
C
C        W CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES;
C
C        IND CONTAINS THE SUBMATRIX INDICES ASSOCIATED WITH THE
C          CORRESPONDING EIGENVALUES IN W -- 1 FOR EIGENVALUES
C          BELONGING TO THE FIRST SUBMATRIX FROM THE TOP,
C          2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS;
C
C        RV1 IS A TEMPORARY STORAGE ARRAY.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
      DATA         MACHEP/16.0D-13/
C*SUN  DATA MACHEP/1.0D-14/
C     DATA MACHEP/1.0D-13/
c     DATA MACHEP/-1.0D-14/ 
      write(6,*)"i am in imtqlv"
c      write(6,*)"E(I)", (E(I),I=1,10)
      IERR = 0
      K = 0
      TAG = 0
C
      DO 10 I = 1, N
         W(I) = D(I)
         IF (I.NE.1) RV1(I-1) = E(I)
   10 CONTINUE
C
      E2(1) = 0.0D0
      RV1(N) = 0.0D0
C
      DO 140 L = 1, N
         J = 0
C     :::::::::: LOOK FOR SMALL SUB-DIAGONAL ELEMENT ::::::::::
   20    DO 30 M = L, N
            IF (M.EQ.N) GO TO 40
            IF (ABS(RV1(M)).LE.MACHEP* (ABS(W(M))+ABS(W(M+
     .          1)))) GO TO 40
C     :::::::::: GUARD AGAINST UNDERFLOWED ELEMENT OF E2 ::::::::::
            IF (E2(M+1).EQ.0.0D0) GO TO 50
   30    CONTINUE
C
   40    IF (M.LE.K) GO TO 60
         IF (M.NE.N) E2(M+1) = 0.0D0
   50    K = M
         TAG = TAG + 1
   60    P = W(L)
         IF (M.EQ.L) GO TO 100
         IF (J.EQ.30) GO TO 150
         J = J + 1
C     :::::::::: FORM SHIFT ::::::::::
         G = (W(L+1)-P)/ (2.0D0*RV1(L))
         R = SQRT(G*G+1.0D0)
         G = W(M) - P + RV1(L)/ (G+SIGN(R,G))
         S = 1.0D0
         C = 1.0D0
         P = 0.0D0
         MML = M - L
C     :::::::::: FOR I=M-1 STEP -1 UNTIL L DO -- ::::::::::
         DO 90 II = 1, MML
            I = M - II
            F = S*RV1(I)
            B = C*RV1(I)
            IF (ABS(F).LT.ABS(G)) GO TO 70
            C = G/F
            R = SQRT(C*C+1.0D0)
            RV1(I+1) = F*R
            S = 1.0D0/R
            C = C*S
            GO TO 80

   70       S = F/G
            R = SQRT(S*S+1.0D0)
            RV1(I+1) = G*R
            C = 1.0D0/R
            S = S*C
   80       G = W(I+1) - P
            R = (W(I)-G)*S + 2.0D0*C*B
            P = S*R
            W(I+1) = G + P
            G = C*R - B
   90    CONTINUE
C
         W(L) = W(L) - P
         RV1(L) = G
         RV1(M) = 0.0D0
         GO TO 20
C     :::::::::: ORDER EIGENVALUES ::::::::::
  100    IF (L.EQ.1) GO TO 120
C     :::::::::: FOR I=L STEP -1 UNTIL 2 DO -- ::::::::::
         DO 110 II = 2, L
            I = L + 2 - II
            IF (P.GE.W(I-1)) GO TO 130
            W(I) = W(I-1)
            IND(I) = IND(I-1)
  110    CONTINUE
C
  120    I = 1
  130    W(I) = P
         IND(I) = TAG
  140 CONTINUE
C
      GO TO 160
C     :::::::::: SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS ::::::::::
  150 IERR = L
  160 RETURN
C     :::::::::: LAST CARD OF IMTQLV ::::::::::
      END
c----------------------------------------------------------------------------------
c
c----------------------------------------------------------------------------------
      SUBROUTINE TINVIT(NM,N,D,E,E2,M,W,IND,Z,IERR,RV1,RV2,RV3,RV4,RV6)
C
c     INCLUDE '../include/param.cmn'
C     WARNING   RPA PGM MUST HAVE THESE DIMENSIONS
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 60000)
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
      INTEGER      I, J, M, N, P, Q, R, S, II, IP, JJ, NM, ITS, TAG,
     .             IERR, GROUP
      DOUBLE PRECISION D(N), E(N), E2(N), W(M), Z(NM,M), RV1(N), RV2(N)
      DOUBLE PRECISION RV3(N), RV4(N), RV6(N)
      DOUBLE PRECISION U, V, UK, XU, X0, X1, EPS2, EPS3, EPS4, NORM
      DOUBLE PRECISION ORDER, MACHEP
C     DOUBLE PRECISION SQRT,ABS,DFLOAT
      DOUBLE PRECISION BIG, SCALE
      INTEGER      IND(M)

      U = ZERO
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE INVERSE ITERATION TECH-
C     NIQUE IN THE ALGOL PROCEDURE TRISTURM BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
C
C     THIS SUBROUTINE FINDS THOSE EIGENVECTORS OF A TRIDIAGONAL
C     SYMMETRIC MATRIX CORRESPONDING TO SPECIFIED EIGENVALUES,
C     USING INVERSE ITERATION.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX;
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY;
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E,
C          WITH ZEROS CORRESPONDING TO NEGLIGIBLE ELEMENTS OF E.
C          E(I) IS CONSIDERED NEGLIGIBLE IF IT IS NOT LARGER THAN
C          THE PRODUCT OF THE RELATIVE MACHINE PRECISION AND THE SUM
C          OF THE MAGNITUDES OF D(I) AND D(I-1).  E2(1) MUST CONTAIN
C          0.0D0 IF THE EIGENVALUES ARE IN ASCENDING ORDER, OR 2.0D0
C          IF THE EIGENVALUES ARE IN DESCENDING ORDER.  IF  BISECT,
C          TRIDIB, OR  IMTQLV  HAS BEEN USED TO FIND THE EIGENVALUES,
C          THEIR OUTPUT E2 ARRAY IS EXACTLY WHAT IS EXPECTED HERE;
C
C        M IS THE NUMBER OF SPECIFIED EIGENVALUES;
C
C        W CONTAINS THE M EIGENVALUES IN ASCENDING OR DESCENDING ORDER;
C
C        IND CONTAINS IN ITS FIRST M POSITIONS THE SUBMATRIX INDICES
C          ASSOCIATED WITH THE CORRESPONDING EIGENVALUES IN W --
C          1 FOR EIGENVALUES BELONGING TO THE FIRST SUBMATRIX FROM
C          THE TOP, 2 FOR THOSE BELONGING TO THE SECOND SUBMATRIX, ETC.
C
C     ON OUTPUT:
C
C        ALL INPUT ARRAYS ARE UNALTERED;
C
C        Z CONTAINS THE ASSOCIATED SET OF ORTHONORMAL EIGENVECTORS.
C          ANY VECTOR WHICH FAILS TO CONVERGE IS SET TO ZERO;
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          -R         IF THE EIGENVECTOR CORRESPONDING TO THE R-TH
C                     EIGENVALUE FAILS TO CONVERGE IN 5 ITERATIONS;
C
C        RV1, RV2, RV3, RV4, AND RV6 ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     :::::::::: MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C                MACHEP = 16.0D0**(-13) FOR LONG FORM ARITHMETIC
C                ON S360 ::::::::::
c      DATA MACHEP/16.0D-13/
c      DATA MACHEP/-1.0D-16/
c     DATA MACHEP/1.D35/
      DATA         MACHEP/16.0D-13/
C     DATA MACHEP/1.0D-13/
C
      write(6,*)"i am in tinvit"
      IERR = 0
      IF (M.EQ.0) GO TO 350
      TAG = 0
      ORDER = 1.0D0 - E2(1)
      Q = 0
C     :::::::::: ESTABLISH AND PROCESS NEXT SUBMATRIX ::::::::::
   10 P = Q + 1
C
      DO 20 Q = P, N
         IF (Q.EQ.N) GO TO 30
         IF (E2(Q+1).EQ.0.0D0) GO TO 30
   20 CONTINUE
C     :::::::::: FIND VECTORS BY INVERSE ITERATION ::::::::::
   30 TAG = TAG + 1
      S = 0
C
      DO 340 R = 1, M
         IF (IND(R).NE.TAG) GO TO 340
         ITS = 1
         X1 = W(R)
         IF (S.NE.0) GO TO 70
C     :::::::::: CHECK FOR ISOLATED ROOT ::::::::::
         XU = 1.0D0
         IF (P.NE.Q) GO TO 40
         RV6(P) = 1.0D0
         GO TO 310

   40    NORM = ABS(D(P))
         IP = P + 1
C
         DO 50 I = IP, Q
            NORM = NORM + ABS(D(I)) + ABS(E(I))
   50    CONTINUE
C     :::::::::: EPS2 IS THE CRITERION FOR GROUPING,
C                EPS3 REPLACES ZERO PIVOTS AND EQUAL
C                ROOTS ARE MODIFIED BY EPS3,
C                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ::::::::::
         EPS2 = 1.0D-3*NORM
         EPS3 = MACHEP*NORM
         UK = REAL(Q-P+1)
         EPS4 = UK*EPS3
         UK = EPS4/SQRT(UK)
         S = P
   60    GROUP = 0
         GO TO 80
C     :::::::::: LOOK FOR CLOSE OR COINCIDENT ROOTS ::::::::::
   70    IF (ABS(X1-X0).GE.EPS2) GO TO 60
         GROUP = GROUP + 1
         IF (ORDER* (X1-X0).LE.0.0D0) X1 = X0 + ORDER*EPS3
C     :::::::::: ELIMINATION WITH INTERCHANGES AND
C                INITIALIZATION OF VECTOR ::::::::::
   80    V = 0.0D0
C
         DO 110 I = P, Q
            RV6(I) = UK
            IF (I.EQ.P) GO TO 100
            IF (ABS(E(I)).LT.ABS(U)) GO TO 90
C     :::::::::: WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
C                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ::::::::::
            XU = U/E(I)
            RV4(I) = XU
            RV1(I-1) = E(I)
            RV2(I-1) = D(I) - X1
            RV3(I-1) = 0.0D0
            IF (I.NE.Q) RV3(I-1) = E(I+1)
            U = V - XU*RV2(I-1)
            V = -XU*RV3(I-1)
            GO TO 110

   90       XU = E(I)/U
            RV4(I) = XU
            RV1(I-1) = U
            RV2(I-1) = V
            RV3(I-1) = 0.0D0
  100       U = D(I) - X1 - XU*V
            IF (I.NE.Q) V = E(I+1)
  110    CONTINUE
C
         IF (U.EQ.0.0D0) U = EPS3
         RV1(Q) = U
         RV2(Q) = 0.0D0
         RV3(Q) = 0.0D0
C     :::::::::: BACK SUBSTITUTION
C                FOR I=Q STEP -1 UNTIL P DO -- ::::::::::
  120    DO 130 II = P, Q
            I = P + Q - II
c           write(6,*)"I",I
c           write(6,*)"U",U
c           write(6,*)"V",V
c           write(6,*)"RV2",RV2(I)
c           write(6,*)"RV3",RV3(I)
c           write(6,*)"RV1",RV1(I)
            RV6(I) = (RV6(I)-U*RV2(I)-V*RV3(I))/RV1(I)
            V = U
            U = RV6(I)
c           write(6,*)"U",U
c           write(6,*)"V",V
  130    CONTINUE
C     :::::::::: ORTHOGONALIZE WITH RESPECT TO PREVIOUS
C                MEMBERS OF GROUP ::::::::::
         IF (GROUP.EQ.0) GO TO 180
         J = R
C
         DO 170 JJ = 1, GROUP
  140       J = J - 1
            IF (IND(J).NE.TAG) GO TO 140
            XU = 0.0D0
C
            DO 150 I = P, Q
               XU = XU + RV6(I)*Z(I,J)
  150       CONTINUE
C
            DO 160 I = P, Q
               RV6(I) = RV6(I) - XU*Z(I,J)
  160       CONTINUE
C
  170    CONTINUE
C
  180    NORM = 0.0D0
C
         DO 190 I = P, Q
            NORM = NORM + ABS(RV6(I))
  190    CONTINUE
C
         IF (NORM.GE.1.0D0) GO TO 260
C     :::::::::: FORWARD SUBSTITUTION ::::::::::
         IF (ITS.EQ.5) GO TO 250
CMCZ     IF (ITS.EQ.6) GO TO 250
         IF (NORM.NE.0.0D0) GO TO 200
         RV6(S) = EPS4
         S = S + 1
         IF (S.GT.Q) S = P
         GO TO 220

  200    XU = EPS4/NORM
C
         DO 210 I = P, Q
            RV6(I) = RV6(I)*XU
  210    CONTINUE
C     :::::::::: ELIMINATION OPERATIONS ON NEXT VECTOR
C                ITERATE ::::::::::
  220    DO 240 I = IP, Q
            U = RV6(I)
C     :::::::::: IF RV1(I-1) .EQ. E(I), A ROW INTERCHANGE
C                WAS PERFORMED EARLIER IN THE
C                TRIANGULARIZATION PROCESS ::::::::::
            IF (RV1(I-1).NE.E(I)) GO TO 230
            U = RV6(I-1)
            RV6(I-1) = RV6(I)
  230       RV6(I) = U - RV4(I)*RV6(I-1)
  240    CONTINUE
C
         ITS = ITS + 1
         GO TO 120
C     :::::::::: SET ERROR -- NON-CONVERGED EIGENVECTOR ::::::::::
  250    IERR = -R
         XU = 0.0D0
         GO TO 310
C     :::::::::: NORMALIZE SO THAT SUM OF SQUARES IS
C                1 AND EXPAND TO FULL ORDER ::::::::::
  260    U = 0.0D0
C
CMARK**** 1/26/87 ****************************************************
C     THE NEXT FEW LINES ARE TO SCALE DOWN THE MAGNITUDE
C     OF THE RV6(I) VALUES TO PREVENT AN OVERFLOW WHICH SOMETIMES
C     OCCURS WITH LARGER MOLECULES.
C*********************************************************************
         BIG = ABS(RV6(P))
C******************************************
C     FIND THE LARGEST ELEMENT OF RV6(I)  *
C******************************************
         DO 270 I = P, Q
            IF (ABS(RV6(I)).GT.BIG) BIG = ABS(RV6(I))
  270    CONTINUE
C***********************************************************************
C     IF BIG TOO SMALL...NO PROBLEM...GO TO NORMAL CODE WITH NO SCALE DO
C***********************************************************************
         IF (BIG.LE.1.0D10) GO TO 280
C*****************************************************************
C     SCALE DOWN THE RV6(I) ELEMENTS AND FORM THE SUM OF SQUARES *
C*****************************************************************
C        DO 850 I = P,Q
CMT           IF(ABS(RV6(I)).LE.1.0D0) RV6(I) = 0.0D0
C          RV6(I) = RV6(I)/BIG
C          U = U + RV6(I)**2
C 850    CONTINUE
C        SCALE = 1.0D0/BIG
C        GOTO 1100
C**********************************************************************
C     IF NO ELEMENT OF RV6(I) IS TOO BIG WE JUST FORM THIS SIMPLE SUM *
C**********************************************************************
  280    CONTINUE
         DO 290 I = P, Q
            U = U + RV6(I)**2
  290    CONTINUE
         SCALE = 1.0D0
C
  300    CONTINUE
         XU = (1.0D0/SQRT(U))*SCALE
C
  310    DO 320 I = 1, N
            Z(I,R) = 0.0D0
  320    CONTINUE
C
         DO 330 I = P, Q
c        write(6,*)"I",I,"R",R
            Z(I,R) = RV6(I)*XU
c        write(6,*)"i",Z(I,R)
  330    CONTINUE
C
         X0 = X1
  340 CONTINUE
C
      IF (Q.LT.N) GO TO 10
  350 RETURN
C     :::::::::: LAST CARD OF TINVIT ::::::::::
      END
c---------------------------------------------------------------------------------------
c
c--------------------------------------------------------------------------------------
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
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
      INTEGER      I, J, K, L, M, N, IK, IZ, NM, NV
      DOUBLE PRECISION A(*), Z(NM,M)
      DOUBLE PRECISION H, S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT:
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT;
C
C        N IS THE ORDER OF THE MATRIX;
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT;
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS;
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED;
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT:
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
c     write(6,*)"i am in trbak3"
c     write(6,*)"A",A
c     write(6,*)"Z",Z
      IF (M.EQ.0) GO TO 50
      IF (N.EQ.1) GO TO 50
C
      DO 40 I = 2, N
         L = I - 1
         IZ = (I*L)/2
         IK = IZ + I
         H = A(IK)
         IF (H.EQ.0.0D0) GO TO 40
C
         DO 30 J = 1, M
            S = 0.0D0
            IK = IZ
C
            DO 10 K = 1, L
               IK = IK + 1
               S = S + A(IK)*Z(K,J)
   10       CONTINUE
C     :::::::::: DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ::::::::::
            S = (S/H)/H
            IK = IZ
C
            DO 20 K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S*A(IK)
   20       CONTINUE
C
   30    CONTINUE
C
   40 CONTINUE
C
   50 RETURN
C     :::::::::: LAST CARD OF TRBAK3 ::::::::::
      END
                                                                                






