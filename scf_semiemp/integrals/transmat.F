      subroutine  transmat(CT,CP,SP,M,H)
C
C     THIS SUBROUTINE GIVES THE TRANSFORMATION MATRIX NECESSARY FOR THE
C     EXPRESSION OF P, D, AND F TENSORS IN THE SIMPLEST CANONICAL TERMS.
C
C     M. C. Zerner
C     For ZINDO
C
C     M=0,1,2,3,4,5,6,7,8 FOR S,X,Y,Z,3Z2-R2,X2-Y2,XY,XZ,YZ RESPECTIVELY
C     M=9,10,11,12,13,14,15 FOR Z3,XZ2,YZ2,Z(X2-Y2),XYZ,X(X2-3Y2),AND
C        Y(3X2-Y2)
C     M=51,52,53 FOR X2,Y2,Z2
C
C     H(I) RETURNS AS SIGMA,PI-X,PI-Y,DELTA-(X**2-Y**2),DELTA-XY,PHI-
C          X(X2-3Y2),PHI-Y(3X2-Y2), AND S COMPONENT IF ANY.
C
C     THE PHASES OF THE REAL FNS. ARE THOSE OF FRANK HARRIS IN THE UPPS
C     ALA MONOGRAPH.
C     SPHERICAL COORDINATES. CT=COS(THETA), CP=COS(PHI), SP=SIN(PHI).
C     H RETURNS AS TRANS ELEMENTS OF 3-DIM. SPACE FOR P ORB. AND OF 5-DI
C     M SPACE FOR THE D ORB., ETC.
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (S,C,R,Z)
#include "global.fh"
      Double precision    H(8)
      integer M
      DATA         SQRT3/1.732050807D0/
      parameter   (ZERO=0.0d0, ONE=1.0D0, TWO= 2.0d0, THREE=3.0d0)
c     
c     
c     if (ga_nodeid().eq.0)
c    .write(6,*) "i am  in tansmat"  
      H(1) = ONE
      H(2) = ZERO
      H(3) = ZERO
      H(4) = ZERO
      H(5) = ZERO
      H(6) = ZERO
      H(7) = ZERO
      H(8) = ZERO
      IF (M) 10, 240, 10
   10 CONTINUE
      R = TWO*SP*CP
      ST = ABS(ONE-CT**2)
      IF (1.0D-7-ST) 30, 30, 20
   20 ST = ZERO
      GO TO 40

   30 CONTINUE
      ST = SQRT(ST)
   40 CONTINUE
      Z = CP**2 - SP**2
      IF (M.GT.15) GO TO 200
      GO TO (50,60,70,80,90,100,110,120,130,
     .       140,150,160,170,180,190), M
C     P TYPE FN.
   50 H(1) = CP*ST
      H(2) = CP*CT
      H(3) = -SP
      GO TO 240

   60 H(1) = SP*ST
      H(2) = SP*CT
      H(3) = CP
      GO TO 240

   70 H(1) = CT
      H(2) = -ST
      H(3) = ZERO
      GO TO 240
C     D TYPE FN .
   80 H(1) = (THREE*CT*CT-ONE)/TWO
      H(2) = -SQRT3*CT*ST
      H(3) = ZERO
      H(4) = SQRT3* (ONE-CT**2)/TWO
      H(5) = ZERO
      GO TO 240

   90 H(1) = SQRT3*Z* (ONE-CT*CT)/TWO
      H(2) = CT*Z*ST
      H(3) = -ST*R
      H(4) = Z* (ONE+CT*CT)/TWO
      H(5) = -R*CT
      GO TO 240

  100 H(1) = SQRT3*CP*SP* (ONE-CT*CT)
      H(2) = CT*ST*R
      H(3) = Z*ST
      H(4) = CP*SP* (ONE+CT*CT)
      H(5) = CT*Z
      GO TO 240

  110 H(1) = SQRT3*CP*CT*ST
      H(2) = CP* (TWO*CT*CT-ONE)
      H(3) = -SP*CT
      H(4) = -H(1)/SQRT3
      H(5) = SP*ST
      GO TO 240

  120 H(1) = SQRT3*SP*CT*ST
      H(2) = SP* (TWO*CT*CT-ONE)
      H(3) = CP*CT
      H(4) = -H(1)/SQRT3
      H(5) = -CP*ST
      GO TO 240
C     F TYPE FN.
  130 H(1) = (5.0D0*CT*CT-3.0D0)*CT/TWO
      H(2) = -0.612372436D0* (5.0D0*CT*CT-ONE)*ST
      H(3) = ZERO
      H(4) = 1.936491673D0*CT*ST*ST
      H(5) = ZERO
      H(6) = -0.790569415D0*ST**3
      H(7) = ZERO
      GO TO 240

  140 H(1) = 0.612372436D0* (5.0D0*CT*CT-ONE)*ST*CP
      H(2) = CT* (15.0D0*CT*CT-11.0D0)*CP/4.0D0
      H(3) = - (5.0D0*CT*CT-ONE)*SP/4.0D0
      H(4) = -0.790569415D0*ST* (3.0D0*CT*CT-ONE)*CP
      H(5) = 1.581138830D0*ST*CT*SP
      H(6) = 0.968245837D0*ST*ST*CT*CP
      H(7) = -0.968245837D0*ST*ST*SP
      GO TO 240

  150 H(1) = 0.612372436D0* (5.0D0*CT*CT-ONE)*ST*SP
      H(2) = CT* (15.0D0*CT*CT-11.0D0)*SP/4.0D0
      H(3) = (5.0D0*CT*CT-ONE)*CP/4.0D0
      H(4) = -0.790569415D0*ST* (THREE*CT*CT-ONE)*SP
      H(5) = -1.581138830D0*ST*CT*CP
      H(6) = 0.968245837D0*ST*ST*CT*SP
      H(7) = 0.968245837D0*ST*ST*CP
      GO TO 240

  160 H(1) = 1.936491673D0*CT*ST*ST*Z
      H(2) = -0.790569415D0*ST* (ONE-THREE*CT*CT)*Z
      H(3) = -1.581138830D0*ST*CT*R
      H(4) = CT* (THREE*CT*CT-ONE)*Z/TWO
      H(5) = - (TWO*CT*CT-ONE)*R
      H(6) = -0.612372436D0*ST* (ONE+CT*CT)*Z
C     H(7)=1.224744871D0*ST*CT*R
      H(7) = -1.224744871D0*ST*CT*R
      GO TO 240

  170 H(1) = 1.936491673D0*CT*ST*ST*R
      H(2) = -0.790569415D0*ST* (ONE-THREE*CT*CT)*R
      H(3) = 1.581138830D0*ST*CT*Z
      H(4) = CT* (THREE*CT*CT-ONE)*R/TWO
      H(5) = (TWO*CT*CT-ONE)*Z
      H(6) = -0.612372436D0*ST* (ONE+CT*CT)*R
      H(7) = -1.224744871D0*ST*CT*Z
C     H(6)= 0.612372436D0*ST*(ONE+CT*CT)*R
C     H(7)= 1.224744871D0*ST*CT*Z
      GO TO 240

  180 CONTINUE
      S = R*CP + SP*Z
      C = Z*CP - R*SP
      H(1) = 0.790569415D0*ST**3*C
      H(2) = 0.968245837D0*ST*ST*CT*C
      H(3) = -0.968245837D0*ST*ST*S
      H(4) = 0.612372436D0*ST* (ONE+CT*CT)*C
      H(5) = -1.224744871D0*ST*CT*S
      H(6) = (THREE+CT*CT)*C*CT/4.0D0
      H(7) = - (ONE+THREE*CT*CT)*S/4.0D0
      GO TO 240

  190 CONTINUE
      S = R*CP + SP*Z
      C = Z*CP - R*SP
      H(1) = 0.790569415D0*ST**3*S
      H(2) = 0.968245837D0*ST*ST*CT*S
      H(3) = 0.968245837D0*ST*ST*C
      H(4) = 0.612372436D0*ST* (ONE+CT*CT)*S
      H(5) = 1.224744871D0*ST*CT*C
      H(6) = (THREE+CT*CT)*S*CT/4.0D0
      H(7) = (ONE+THREE*CT*CT)*C/4.0D0
      GO TO 240

  200 CONTINUE
      MM = M - 50
      GO TO (210,220,230), MM
C     QUADRATIC TRANSFORM, X2,Y2,Z2
  210 H(1) = Z*SQRT3* (ONE-CT*CT)/TWO - (THREE*CT*CT-ONE)/ (TWO*SQRT3)
      H(2) = CT*ST* (ONE+Z)
      H(3) = -ST*R
      H(4) = (Z* (ONE+CT*CT)-ONE+CT*CT)/TWO
      H(5) = -CT*R
      H(6) = ZERO
      H(7) = ZERO
      H(8) = 1.29099445D0
      GO TO 240

  220 H(1) = - (THREE*CT*CT-ONE)/ (TWO*SQRT3) - Z* (ONE-CT*CT)*SQRT3/TWO
      H(2) = CT*ST* (ONE-Z)
      H(3) = ST*R
      H(4) = - (ONE-CT*CT+Z* (ONE+CT*CT))/TWO
      H(5) = CT*R
      H(6) = ZERO
      H(7) = ZERO
      H(8) = 1.29099445D0
      GO TO 240

  230 H(1) = (THREE*CT*CT-ONE)/SQRT3
      H(4) = (ONE-CT*CT)
      H(5) = ZERO
      H(2) = -TWO*ST*CT
      H(3) = ZERO
      H(6) = ZERO
      H(7) = ZERO
      H(8) = 1.29099445D0
C     ABOVE NUMBER IS SQRT(5/3) AND IS S COMPONENT OF X2 OR Y2 OR Z2
  240 CONTINUE
      RETURN

      END
