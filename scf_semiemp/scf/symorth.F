      SUBROUTINE SYMORTH(N,BB,T,QQ,NIN)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION BB(*),T(*),QQ(*),NIN(*)
      parameter (ZERO = 0.0d0, ONE = 1.0d0)

      CALL HDIAG(N,N,0,BB,T,QQ,NIN)
         L = 0
         DO 10 I = 1, N
            L = L + I
          write(6,*)"T",T(L)
            IF (T(L).LE.ZERO) GO TO 170
            QQ(I) = ONE/SQRT(T(L))
   10    CONTINUE
C
C
         DO 40 I = 1, N
            ID = NIN(I)
            DO 30 J = 1, I
               IJ = ID + J
               T(IJ) = 0.0d0
               DO 20 K = 1, N
                  IK = (K-1)*N + I
                  JK = (K-1)*N + J
                  T(IJ) = T(IJ) + BB(IK)*QQ(K)*BB(JK)
   20          CONTINUE
   30       CONTINUE
   40    CONTINUE


      RETURN
  170 WRITE (6,9000) I, I, T(L)
      CALL ZEXIT(3)
      RETURN

 9000 FORMAT (25X,'IMPOSSIBLE S EIGENVALUE S(',I3,',',I3,') = ',F8.3,
     .       25X,'CHECK DYNAMIC ARRAY SIZES AGAINST CI BASIS')

      END

