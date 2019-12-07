      SUBROUTINE DGEOM(CT,CP,SP,M,HT,HP)
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
C
C  THIS SUBROUTINE IS ANALAGOUS TO ROUTINE GEOM
C  EXCEPYT NOW WE DETERMINE THE DIFFERENTIAL OF THE TRANSFORMATION
C  PARTS       FEB 27/1981
C
C  HT = DH / DTHETA
C  HP = DH / DPHI
C
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)

      DIMENSION    HT(10), HP(10)
      DATA         DSQ3/1.7320508075688773D0/
      DATA         THRESH/1.D-10/

c     write(6,*)"i am in dgeom"
      S2P = ZERO
      S2T = ZERO
      C2P = ZERO
      C2T = ZERO
      IF (M.EQ.0) THEN
         HT(1) = ZERO
         HP(1) = ZERO
         GO TO 200

      END IF
      ST = ONE - CT*CT
      IF (THRESH-ST) 20, 20, 10
   10 ST = ZERO
      GO TO 30

   20 ST = SQRT(ST)
   30 IF (M.LT.4) GO TO 40
      C2P = CP*CP - SP*SP
      S2P = TWO*SP*CP
      C2T = CT*CT - ST*ST
      S2T = TWO*ST*CT
   40 GO TO (50,60,70,80,90,100,110,120,130,
     .       140,150,160,170,180,190), M
C  PX
   50 HT(1) = CT*CP
      HT(2) = -ST*CP
      HT(3) = ZERO
      HP(1) = -ST*SP
      HP(2) = -CT*SP
      HP(3) = -CP
      GO TO 200
C  PY
   60 HT(1) = CT*SP
      HT(2) = -ST*SP
      HT(3) = ZERO
      HP(1) = ST*CP
      HP(2) = CT*CP
      HP(3) = -SP
      GO TO 200
C  PZ
   70 HT(1) = -ST
      HT(2) = -CT
      HT(3) = ZERO
      HP(1) = ZERO
      HP(2) = ZERO
      HP(3) = ZERO
      GO TO 200
C  DZ2
   80 HT(1) = -1.5D0*S2T
      HT(2) = -DSQ3*C2T
      HT(3) = ZERO
      HT(4) = -HT(1)/DSQ3
      HT(5) = ZERO
      HP(1) = ZERO
      HP(2) = ZERO
      HP(3) = ZERO
      HP(4) = ZERO
      HP(5) = ZERO
      GO TO 200
C  DX2-Y2
   90 HT(1) = DSQ3*S2T*C2P/TWO
      HT(2) = C2P*C2T
      HT(3) = -CT*S2P
      HT(4) = -C2P*S2T/TWO
      HT(5) = ST*S2P
      HP(1) = -DSQ3*S2P*ST*ST
      HP(2) = -S2T*S2P
      HP(3) = -TWO*C2P*ST
      HP(4) = -S2P* (ONE+CT*CT)
      HP(5) = -TWO*CT*C2P
      GO TO 200
C  DXY
  100 HT(1) = DSQ3*S2T*S2P/TWO
      HT(2) = C2T*S2P
      HT(3) = CT*C2P
      HT(4) = -S2T*S2P/TWO
      HT(5) = -ST*C2P
      HP(1) = DSQ3*C2P*ST*ST
      HP(2) = S2T*C2P
      HP(3) = -TWO*S2P*ST
      HP(4) = C2P* (ONE+CT*CT)
      HP(5) = -TWO*CT*S2P
      GO TO 200
C  DXZ
  110 HT(1) = DSQ3*C2T*CP
      HT(2) = -TWO*CP*S2T
      HT(3) = ST*SP
      HT(4) = -C2T*CP
      HT(5) = CT*SP
      HP(1) = -DSQ3*S2T*SP/TWO
      HP(2) = -C2T*SP
      HP(3) = -CT*CP
      HP(4) = S2T*SP/TWO
      HP(5) = ST*CP
      GO TO 200
C  DYZ
  120 HT(1) = DSQ3*C2T*SP
      HT(2) = -TWO*S2T*SP
      HT(3) = -ST*CP
      HT(4) = -C2T*SP
      HT(5) = -CT*CP
      HP(1) = DSQ3*S2T*CP/TWO
      HP(2) = C2T*CP
      HP(3) = -CT*SP
      HP(4) = -S2T*CP/TWO
      HP(5) = ST*SP
      GO TO 200
C    ***G9
  130 HT(1) = -7.500*ST*CT**2 + 1.500*ST
      HT(2) = 6.1237243600*ST**2*CT - 3.0618621800*CT**3 +
     .        0.6123724360*CT
      HT(3) = ZERO
      HT(4) = -1.9364916730*ST**3 + 3.8729833460*ST*CT**2
      HT(5) = ZERO
      HT(6) = -2.3717082450*ST**2*CT
      HT(7) = ZERO
C    ***G9
      HP(1) = ZERO
      HP(2) = ZERO
      HP(3) = ZERO
      HP(4) = ZERO
      HP(5) = ZERO
      HP(6) = ZERO
      HP(7) = ZERO
      GO TO 200
C    ***G10
  140 HT(1) = -6.1237243600*ST**2*CT*CP + 3.0618621800*CT**3*CP -
     .        0.6123724360*CT*CP
      HT(2) = -11.250*ST*CT**2*CP + 2.750*ST*CP
      HT(3) = 2.500*ST*CT*SP
      HT(4) = 4.7434164900*ST**2*CT*CP - 2.3717082450*CT**3*CP +
     .        0.7905694150*CT*CP
      HT(5) = -1.5811388300*ST**2*SP + 1.5811388300*CT**2*SP
      HT(6) = -0.9682458370*ST**3*CP + 1.9364916740*ST*CT**2*CP
      HT(7) = -1.9364916740*ST*CT*SP
C    ***G10
      HP(1) = -3.0618621800*ST*CT**2*SP + 0.6123724360*ST*SP
      HP(2) = -3.750*CT**3*SP + 2.750*CT*SP
      HP(3) = -1.250*CT**2*CP + 0.250*CP
      HP(4) = 2.3717082450*ST*CT**2*SP - 0.7905694150*ST*SP
      HP(5) = 1.5811388300*ST*CT*CP
      HP(6) = -0.9682458370*ST**2*CT*SP
      HP(7) = -0.9682458370*ST**2*CP
      GO TO 200
C    ***G11
  150 HT(1) = -6.1237243600*ST**2*CT*SP + 3.0618621800*CT**3*SP -
     .        0.6123724360*CT*SP
      HT(2) = -11.250*ST*CT**2*SP + 2.750*ST*SP
      HT(3) = -2.500*ST*CT*CP
      HT(4) = 4.7434164900*ST**2*CT*SP - 2.3717082450*CT**3*SP +
     .        0.7905694150*CT*SP
      HT(5) = 1.5811388300*ST**2*CP - 1.5811388300*CT**2*CP
      HT(6) = -0.9682458370*ST**3*SP + 1.9364916740*ST*CT**2*SP
      HT(7) = 1.9364916740*ST*CT*CP
C    ***G11
      HP(1) = 3.0618621800*ST*CT**2*CP - 0.6123724360*ST*CP
      HP(2) = 3.750*CT**3*CP - 2.750*CT*CP
      HP(3) = -1.250*CT**2*SP + 0.250*SP
      HP(4) = -2.3717082450*ST*CT**2*CP + 0.7905694150*ST*CP
      HP(5) = 1.5811388300*ST*CT*SP
      HP(6) = 0.9682458370*ST**2*CT*CP
      HP(7) = -0.9682458370*ST**2*SP
      GO TO 200
C    ***G12
  160 HT(1) = 1.9364916730*ST**3*SP**2 - 1.9364916730*ST**3*CP**2 -
     .        3.8729833460*ST*CT**2*SP**2 + 3.8729833460*ST*CT**2*CP**2
      HT(2) = 4.7434164900*ST**2*CT*SP**2 -
     .        4.7434164900*ST**2*CT*CP**2 - 2.3717082450*CT**3*SP**2 +
     .        2.3717082450*CT**3*CP**2 + 0.7905694150*CT*SP**2 -
     .        0.7905694150*CT*CP**2
      HT(3) = 3.1622776600*ST**2*SP*CP - 3.1622776600*CT**2*SP*CP
      HT(4) = 4.500*ST*CT**2*SP**2 - 4.500*ST*CT**2*CP**2 -
     .        0.500*ST*SP**2 + 0.500*ST*CP**2
      HT(5) = 8.000*ST*CT*SP*CP
      HT(6) = -1.2247448720*ST**2*CT*SP**2 +
     .        1.2247448720*ST**2*CT*CP**2 + 0.6123724360*CT**3*SP**2 -
     .        0.6123724360*CT**3*CP**2 + 0.6123724360*CT*SP**2 -
     .        0.6123724360*CT*CP**2
      HT(7) = -2.4494897420*ST**2*SP*CP + 2.4494897420*CT**2*SP*CP
C    ***G12
      HP(1) = -7.7459666920*ST**2*CT*SP*CP
      HP(2) = -9.4868329800*ST*CT**2*SP*CP + 3.1622776600*ST*SP*CP
      HP(3) = 3.1622776600*ST*CT*SP**2 - 3.1622776600*ST*CT*CP**2
      HP(4) = -6.000*CT**3*SP*CP + 2.000*CT*SP*CP
      HP(5) = 4.000*CT**2*SP**2 - 4.000*CT**2*CP**2 - 2.000*SP**2 +
     .        2.000*CP**2
      HP(6) = 2.4494897440*ST*CT**2*SP*CP + 2.4494897440*ST*SP*CP
      HP(7) = -2.4494897420*ST*CT*SP**2 + 2.4494897420*ST*CT*CP**2
      GO TO 200
C    ***G13
  170 HT(1) = -3.8729833460*ST**3*SP*CP + 7.7459666920*ST*CT**2*SP*CP
      HT(2) = -9.4868329800*ST**2*CT*SP*CP + 4.7434164900*CT**3*SP*CP -
     .        1.5811388300*CT*SP*CP
      HT(3) = 1.5811388300*ST**2*SP**2 - 1.5811388300*ST**2*CP**2 -
     .        1.5811388300*CT**2*SP**2 + 1.5811388300*CT**2*CP**2
      HT(4) = -9.000*ST*CT**2*SP*CP + 1.000*ST*SP*CP
      HT(5) = 4.000*ST*CT*SP**2 - 4.000*ST*CT*CP**2
      HT(6) = -2.4494897420*ST**2*CT*SP*CP + 1.2247448710*CT**3*SP*CP +
     .        1.2247448710*CT*SP*CP
      HT(7) = 1.2247448710*ST**2*SP**2 - 1.2247448710*ST**2*CP**2 -
     .        1.2247448710*CT**2*SP**2 + 1.2247448710*CT**2*CP**2
C    ***G13
      HP(1) = -3.8729833460*ST**2*CT*SP**2 + 3.8729833460*ST**2*CT*CP**2
      HP(2) = -4.7434164900*ST*CT**2*SP**2 +
     .        4.7434164900*ST*CT**2*CP**2 + 1.5811388300*ST*SP**2 -
     .        1.5811388300*ST*CP**2
      HP(3) = -6.3245553200*ST*CT*SP*CP
      HP(4) = -3.000*CT**3*SP**2 + 3.000*CT**3*CP**2 + 1.000*CT*SP**2 -
     .        1.000*CT*CP**2
      HP(5) = -8.000*CT**2*SP*CP + 4.000*SP*CP
      HP(6) = -1.2247448710*ST*CT**2*SP**2 +
     .        1.2247448710*ST*CT**2*CP**2 - 1.2247448710*ST*SP**2 +
     .        1.2247448710*ST*CP**2
      HP(7) = -4.8989794840*ST*CT*SP*CP
      GO TO 200
C    ***G14
  180 HT(1) = -7.1151247350*ST**2*CT*SP**2*CP +
     .        2.3717082450*ST**2*CT*CP**3
      HT(2) = 2.9047375110*ST**3*SP**2*CP - 0.9682458370*ST**3*CP**3 -
     .        5.8094750220*ST*CT**2*SP**2*CP +
     .        1.9364916740*ST*CT**2*CP**3
      HT(3) = 1.9364916740*ST*CT*SP**3 - 5.8094750220*ST*CT*SP*CP**2
      HT(4) = 3.6742346160*ST**2*CT*SP**2*CP -
     .        1.2247448720*ST**2*CT*CP**3 -
     .        1.8371173080*CT**3*SP**2*CP + 0.6123724360*CT**3*CP**3 -
     .        1.8371173080*CT*SP**2*CP + 0.6123724360*CT*CP**3
      HT(5) = -1.2247448710*ST**2*SP**3 + 3.6742346130*ST**2*SP*CP**2 +
     .        1.2247448710*CT**2*SP**3 - 3.6742346130*CT**2*SP*CP**2
      HT(6) = 2.250*ST*CT**2*SP**2*CP - 0.750*ST*CT**2*CP**3 +
     .        2.250*ST*SP**2*CP - 0.750*ST*CP**3
      HT(7) = -1.500*ST*CT*SP**3 + 4.500*ST*CT*SP*CP**2
C    ***G14
      HP(1) = 2.3717082450*ST**3*SP**3 - 7.1151247350*ST**3*SP*CP**2
      HP(2) = 2.9047375110*ST**2*CT*SP**3 -
     .        8.7142125330*ST**2*CT*SP*CP**2
      HP(3) = 8.7142125330*ST**2*SP**2*CP - 2.9047375110*ST**2*CP**3
      HP(4) = 1.8371173080*ST*CT**2*SP**3 -
     .        5.5113519240*ST*CT**2*SP*CP**2 + 1.8371173080*ST*SP**3 -
     .        5.5113519240*ST*SP*CP**2
      HP(5) = 11.0227038390*ST*CT*SP**2*CP - 3.6742346130*ST*CT*CP**3
      HP(6) = 0.750*CT**3*SP**3 - 2.250*CT**3*SP*CP**2 +
     .        2.250*CT*SP**3 - 6.750*CT*SP*CP**2
      HP(7) = 6.750*CT**2*SP**2*CP - 2.250*CT**2*CP**3 +
     .        2.250*SP**2*CP - 0.750*CP**3
      GO TO 200
C    ***G15
  190 HT(1) = -2.3717082450*ST**2*CT*SP**3 +
     .        7.1151247350*ST**2*CT*SP*CP**2
      HT(2) = 0.9682458370*ST**3*SP**3 - 2.9047375110*ST**3*SP*CP**2 -
     .        1.9364916740*ST*CT**2*SP**3 +
     .        5.8094750220*ST*CT**2*SP*CP**2
      HT(3) = -5.8094750220*ST*CT*SP**2*CP + 1.9364916740*ST*CT*CP**3
      HT(4) = 1.2247448720*ST**2*CT*SP**3 -
     .        3.6742346160*ST**2*CT*SP*CP**2 -
     .        0.6123724360*CT**3*SP**3 + 1.8371173080*CT**3*SP*CP**2 -
     .        0.6123724360*CT*SP**3 + 1.8371173080*CT*SP*CP**2
      HT(5) = 3.6742346130*ST**2*SP**2*CP - 1.2247448710*ST**2*CP**3 -
     .        3.6742346130*CT**2*SP**2*CP + 1.2247448710*CT**2*CP**3
      HT(6) = 0.750*ST*CT**2*SP**3 - 2.250*ST*CT**2*SP*CP**2 +
     .        0.750*ST*SP**3 - 2.250*ST*SP*CP**2
      HT(7) = 4.500*ST*CT*SP**2*CP - 1.500*ST*CT*CP**3
C    ***G15
      HP(1) = -7.1151247350*ST**3*SP**2*CP + 2.3717082450*ST**3*CP**3
      HP(2) = -8.7142125330*ST**2*CT*SP**2*CP +
     .        2.9047375110*ST**2*CT*CP**3
      HP(3) = 2.9047375110*ST**2*SP**3 - 8.7142125330*ST**2*SP*CP**2
      HP(4) = -5.5113519240*ST*CT**2*SP**2*CP +
     .        1.8371173080*ST*CT**2*CP**3 - 5.5113519240*ST*SP**2*CP +
     .        1.8371173080*ST*CP**3
      HP(5) = 3.6742346130*ST*CT*SP**3 - 11.0227038390*ST*CT*SP*CP**2
      HP(6) = -2.250*CT**3*SP**2*CP + 0.750*CT**3*CP**3 -
     .        6.750*CT*SP**2*CP + 2.250*CT*CP**3
      HP(7) = 2.250*CT**2*SP**3 - 6.750*CT**2*SP*CP**2 + 0.750*SP**3 -
     .        2.250*SP*CP**2
      GO TO 200

  200 RETURN

      END
