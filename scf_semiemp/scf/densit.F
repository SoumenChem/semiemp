      SUBROUTINE DENSIT(BB,F,W,NB,NORBS,NRR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      DIMENSION    W(*), BB(*), F(*)

c     write(6,*)"i am in densit",(W(I),I=1,2)
      NB1 = NB
      NRR = (NB1* (NB1+1))/2
      DO 10 I = 1, NRR
         F(I) = 0.0D0
   10 CONTINUE
      DO 40 K = 1, NORBS
         IF (abs(W(K)).LT.0.001D0) GO TO 40
         KNB = (K-1)*NB1
         IJ = 1
         DO 30 I = 1, NB1
            XX = BB(KNB+I)*W(K)
            DO 20 J = 1, I
               KJ = KNB + J
c        write(6,*)"F,BB,XX",F(IJ),BB(KJ),XX
               F(IJ) = F(IJ) + BB(KJ)*XX
c       write(6,*)"F(IJ),IJ",F(IJ),IJ
               IJ = IJ + 1
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
c     write(6,*)"here",(F(I),I=1,nrr)
      RETURN

      END


