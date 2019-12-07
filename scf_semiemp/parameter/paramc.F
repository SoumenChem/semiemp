      SUBROUTINE PARAMC(HHHHH,IORB,NPQ,KA,KZ,NWL)
C
C     THIS SUBROUTINE INSERTS PARAMETERS FOR THE CORE (INACTIVE)
C     ORBITALS
C
C     NPQ=NP(I)
C     KA=NU(I)
C     KZ=KAN(KA)
C     NWL=NW(I)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      DIMENSION    SCORE(92,6), PCORE(82,5), DCORE(72,4), FCORE(62,3)
C     THESE ARE EXPERIMENTAL CORE IONIZATION POTENTIALS IN EV.
C      1S                          Li   Be     B     C    N     O
      DATA         SCORE/14., 25., 55., 111., 188., 284., 399., 532.,
C      1S/2S       F     Ne                  H    He    Li    Be
     .             686., 867., 82*999.9999, 0.00, 0.00, 0.00, 0.00,
C      2S           B     C      N   O    F    Ne   Na   Mg   Al
     .             0.00, 0.00, 0.00, 24., 31., 45., 63., 89., 118.,
C      2S           Si    P    S     Cl    Ar     K    Ca
     .             149., 189., 229., 270., 320., 377., 438.,
C      3S                                Na  Mg  Al Si   P    S    Cl
     .             72*999.9999, 10*0.00, 1., 2., 1., 8., 16., 16., 18.,
     .             25., 34., 44., 54., 59., 66., 74., 84., 95., 101.,
     .             112., 120., 137., 62*999.9999, 30*0.00, 0.00, 0.00,
     .             0.00, 0.00, 27., 24., 30., 38., 46., 52.,
     .             52*999.9999, 50*0.00, 7., 12., 14., 18., 23., 40.,
     .             33., 38., 38., 38., 38., 39., 32., 36., 40., 63.,
     .             51., 60., 53., 53., 57., 65., 71., 77., 83., 84.,
     .             96., 102., 108., 120., 12*999.9999, 80*0.00, 0.00,
     .             3., 8., 12., 18., 26., 34., 44., 52., 60.,
     .             2*999.9999/
C      2P                               Na   Mg   Al    Si   P     S
      DATA         (PCORE(I,1),I=1,82)/38.1, 54., 74., 100., 136., 165.,
     .             202., 247., 297., 350., 407., 461., 520., 584., 652.,
     .             723., 794., 872., 951., 1044., 62*999.9999/
      DATA         (PCORE(I,2),I=1,82)/0.00, 0.00, 0.00, 3.00, 10.0,
     .             8.0, 7.00, 12.0, 18.0, 26.0, 32.0, 34.0, 38.0, 43.0,
     .             49.0, 56.0, 60.0, 68.0, 74.0, 87.0, 107., 129., 147.,
     .             168., 189., 223., 248., 280., 313., 345.,
     .             52*999.9999/
      DATA         (PCORE(I,3),I=1,82)/20*0.00, 1., 3., 3., 6., 5., 11.,
     .             15., 20., 26., 29., 34., 35., 39., 43., 48., 51.,
     .             59., 67., 77., 89., 99., 110., 123., 147., 177.,
     .             192., 206., 224., 237., 244., 255., 267., 284., 289.,
     .             311., 332., 343., 366., 386., 396., 22*999.9999/
      DATA         (PCORE(I,4),I=1,82)/40*0.000, 2., 2., 3., 7., 12.,
     .             16., 15., 20., 23., 22., 22., 22., 22., 21., 26.,
     .             26., 20., 29., 32., 23., 28., 34., 41., 42., 41.,
     .             52., 57., 57, 63., 70., 88., 96., 105., 114., 132.,
     .             145., 161., 177., 191., 206., 223., 227./
      DATA         (PCORE(I,5),I=1,82)/70*0.00, 1., 1., 3., 5., 8., 11.,
     .             15., 19., 33., 46., 38., 38.0/
      DATA         DCORE/7., 3., 2., 2., 4., 6., 3., 4., 2., 9., 18.,
     .             29., 41., 57., 70., 89., 111., 134., 159., 182.,
     .             206., 229., 255., 282., 309., 337., 370., 408., 447.,
     .             489., 533., 577., 625., 678., 733., 788., 840., 893.,
     .             941., 989., 32*999.9999, 20*0.00, 4., 2., 2., 2., 3.,
     .             1., 3., 9., 16., 24., 32., 40., 50., 63., 78., 92.,
     .             99., 111., 114., 118., 121., 130., 134., 141., 148.,
     .             154., 161., 172., 180., 191., 200., 219., 236., 252.,
     .             267., 281., 303., 322., 343., 370., 396., 424., 452.,
     .             487., 520., 554., 590., 620., 657., 696., 726., 759.,
     .             72*0.00, 72*0.00/
      DATA         FCORE/62*0.00, 62*0.00, 62*0.00/

c     if (ga_nodeid().eq.0)
c    .write(6,*)"i am in paramc"
      FKZ = DFLOAT(KZ)
      NL = (NWL+2)/3
      IF (NWL.GT.3) NL = 2
      IF (NWL.GT.8) NL = 3
      NLP1 = NL + 1
      GO TO (10,80,140,200), NLP1
C     NS ELECTRONS
   10 GO TO (20,30,40,50,60,70), NPQ

   20 CONTINUE
      IF (KZ.LT.11) THEN
         HHHHH = SCORE(KZ,1)
      ELSE
         HHHHH = 955.7462 - 133.85218*FKZ + 14.03841*FKZ*FKZ
      END IF
      GO TO 190

   30 CONTINUE
      IF (KZ.LT.21) THEN
         HHHHH = SCORE(KZ,2)
      ELSE
         HHHHH = 897.20 - 85.24430*FKZ + 3.14643*FKZ*FKZ
      END IF
      GO TO 190

   40 CONTINUE
      IF (KZ.LT.31) THEN
         HHHHH = SCORE(KZ,3)
      ELSE
         HHHHH = 589.708 - 45.85252*FKZ + 1.03380*FKZ*FKZ
      END IF
      GO TO 190

   50 CONTINUE
      IF (KZ.LT.41) THEN
         HHHHH = SCORE(KZ,4)
      ELSE
C     THIS IS A 1% FIT
         HHHHH = 335.7508 - 19.85173*FKZ + 0.319371*FKZ*FKZ
      END IF
      GO TO 190

   60 CONTINUE
      IF (KZ.LT.81) THEN
         HHHHH = SCORE(KZ,5)
      ELSE
         HHHHH = -2195.847 + 39.212613*FKZ - 0.128557*FKZ*FKZ
      END IF
      GO TO 190

   70 CONTINUE
      IF (KZ.LT.86) THEN
         HHHHH = SCORE(KZ,6)
      ELSE
         HHHHH = -4867.74456 + 102.7500*FKZ - 0.5333333*FKZ*FKZ
      END IF
      GO TO 190
C     NP ELECTRONS
   80 CONTINUE
      GO TO (200,90,100,110,120,130), NPQ

   90 CONTINUE
      IF (KZ.LT.31) THEN
         HHHHH = PCORE(KZ-10,1)
      ELSE
         HHHHH = 1249.15 - 104.345*FKZ + 3.2575*FKZ*FKZ
      END IF
      GO TO 190

  100 CONTINUE
      IF (KZ.LT.41) THEN
         HHHHH = PCORE(KZ-10,2)
      ELSE
         HHHHH = 1596.686 - 83.805811*FKZ + 1.323555*FKZ*FKZ
      END IF
      GO TO 190

  110 CONTINUE
      IF (KZ.LT.71) THEN
         HHHHH = PCORE(KZ-10,3)
      ELSE
         HHHHH = 2089.72 - 73.191313*FKZ + 0.696723*FKZ*FKZ
      END IF
      GO TO 190

  120 HHHHH = PCORE(KZ-10,4)
      GO TO 190

  130 HHHHH = PCORE(KZ-10,5)
      GO TO 190
C     ND ELECTRONS
  140 GO TO (200,200,150,160,170,180), NPQ

  150 CONTINUE
      IF (KZ.LT.60) THEN
         HHHHH = DCORE(KZ-20,1)
      ELSE
         HHHHH = 2420.0 - 93.2100*FKZ + 1.15500*FKZ*FKZ
      END IF
      GO TO 190

  160 HHHHH = DCORE(KZ-20,2)
      GO TO 190

  170 HHHHH = 0.00
      GO TO 190

  180 CONTINUE
      GO TO 200
C  40 GO TO (5000,5000,5000,44,5000,5000), NPQ
C  44 HHHHH=0.00
  190 CONTINUE
C
      IF(ABS(HHHHH).LT.1.0) GO TO 200
      HHHHH = -HHHHH/27.2114D0
      GO TO 210

  200 WRITE (6,9000) IORB, NPQ, NL, KA

 9000 FORMAT (/' ***** NO CORE ORBITAL PARAMETERS FOR',3I5,' ON ATOM',
     .       I5,' ERROR SUB PARAMC*****'/)

  210 RETURN

      END
