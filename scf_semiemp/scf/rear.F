      SUBROUTINE REAR(BB,EIG,W,NB,NORBS,IFLAG)
C
C     SUBROUTINE ORDERS VECTORS AND EIGENVALUES
C     OCCUPIED ORBITALS FIRST, THEN EMPTY ORBITALS
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      DIMENSION    W(*), Q(NB), BB(*), EIG(*)

C
C
c     write(6,*)"i am in rear",IFLAG
c     write(6,*)"BB",(BB(I),I=1,21)
c     write(6,*)"W",(W(I),I=1,21)
c     write(6,*)"EIG",(EIG(I),I=1,6)
      IFLAG = 0
      DO 40 I = 1, NORBS
         IF (W(I).GT.1.0D-3) GO TO 10
         IF (IFLAG.EQ.0) IL = I
         IFLAG = 1
         GO TO 40
   10    IF (IFLAG.EQ.0) GO TO 40
         IA = (I-1)*NB
         IB = (IL-1)*NB
         X = EIG(IL)
         EIG(IL) = EIG(I)
         EIG(I) = X
         Y = W(IL)
         W(IL) = W(I)
         W(I) = Y
         DO 20 J = 1, NB
            Q(J) = BB(IB+J)
   20    CONTINUE
         DO 30 J = 1, NB
            BB(IB+J) = BB(IA+J)
            BB(IA+J) = Q(J)
   30    CONTINUE
         IL = IL + 1
   40 CONTINUE
c     write(6,*)"Q",(Q(I),I=1,20)
      RETURN

      END
                          
 
