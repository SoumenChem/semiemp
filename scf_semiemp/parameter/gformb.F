       SUBROUTINE GFORMB(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
     1  I,J,LL,N1,N2,NT,N1N2
     1 ,SCG1SP,SCF2PP,SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,
     2 SCG3SF,SCG2PF,SCG4PF,SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,
     3 SCF2FF,SCF4FF,SCF6FF,R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,
     4 R2PFFF,R2PPPF,R2SDFF,R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF)
CPK** -----------------------------------------------------------
CPK   ONE PART OF THE FORMER SUBROUTINE GFORM
CPK   NEW SUBROUTINES ARE NAMED GFORMA TO GFORME.
CPK   NEW COMMON 'ALLGF' INSERTED FOR COMMUNICATION WITH THE
CPK   NEW SUBROUTINES.
CPK** -----------------------------------------------------------
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MXB = 4000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
c     INCLUDE '../include/allgf.cmn'
c     INCLUDE '../include/ftch.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/z.cmn'
      DIMENSION    GAMMA(*)
      dimension  QQ(*),LABI(*),LABJ(*)
      dimension   LABK(*),LABL(*), NIN(*)
C****************************************************************
      Dimension SCG1SP(*),SCF2PP(*),SCG2SD(*),SCG1PD(*),
     . SCF2PD(*),SCG3PD(*),SCF2DD(*),SCF4DD(*),SCG3SF(*),
     . SCG2PF(*),SCG4PF(*),SCG1DF(*),SCG3DF(*),SCG5DF(*),
     . SCF2PF(*),SCF2DF(*),SCF4DF(*),SCF2FF(*),SCF4FF(*),
     . SCF6FF(*),
C***C D type R-integrals
     .   R1SPPD(*),R2SDDD(*),R2SDPP(*),
C***C F type R-integrals
     .   R1PDDF(*),R1SPDF(*),R2PFDD(*),R2PFFF(*),
     .   R2PPPF(*),R2SDFF(*),R2SDPF(*),R3PDDF(*),
     .   R3SFDF(*),R3SFPD(*),R4PFDD(*),R4PFFF(*)
C
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
      PARAMETER    (ILEVEN=11,ITWELV=12,ITHRTN=13,IFRTEN=14,IFVTEN=15)
      DATA         FP258/0.2581989D0/FP218/0.2182179D0/
      DATA         SQRT3/1.7320508D0/F3P0/3.0D0/
      DATA         F2P0/2.0D0/F4P0/4.0D0/
C
C     ...................................................
C
C     THIS BIT OF NONSENSE IS NECESSARY SINCE UNIX FORTRAN
C     GENERATES ONLY "SHORT" JUMPS FOR COMPUTED GOTOS
C     GOOD FOR ONLY A FEW HUNDRED LINES.
C     THE STANDARD GOTOS GENERATE A PROPER JUMP.
C
C     ...................................................
C
c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in gformb"
      ICS1 = ICS - 14
      THRE = THREE
      GO TO (10,20,30,40,50,60,70,80,90,
     .       100,110,120,130,140,150,160,170,180,
     .       190,200,210,220), ICS1
      
      WRITE (6,9000) ICS1

 9000 FORMAT (' ERROR SUBROUTINE GFORMB.  ICS1= ',I5)

      CALL ZEXIT(3)
   10 GO TO 230

   20 GO TO 240

   30 GO TO 250

   40 GO TO 260

   50 GO TO 270

   60 GO TO 280

   70 GO TO 290

   80 GO TO 300

   90 GO TO 310

  100 GO TO 320

  110 GO TO 330

  120 GO TO 340

  130 GO TO 350

  140 GO TO 360

  150 GO TO 370

  160 GO TO 380

  170 GO TO 390

  180 GO TO 400

  190 GO TO 410

  200 GO TO 420

  210 GO TO 430

  220 GO TO 440

  230 X = GAMMA(N1N2)
      QQ(LL) = X + FOUR*SCF2DD(IP) + 36.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X - FOUR*SCF2DD(IP) + SIX*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X - FOUR*SCF2DD(IP) + SIX*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = X + TWO*SCF2DD(IP) - 24.00D0*SCF4DD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = X + TWO*SCF2DD(IP) - 24.00D0*SCF4DD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + IFOUR
      LL = LL + 4
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X + EIGHT*SCF2DF(IP) + 36.00D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = X + SIX*SCF2DF(IP) + SIX*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + ISIX
      QQ(LL+3) = X + SIX*SCF2DF(IP) + SIX*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = X - 42.00D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IEIGHT
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = X - 42.00D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ININE
      LABL(LL+5) = J + ININE
      QQ(LL+6) = X - TEN*SCF2DF(IP) + 18.00D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITEN
      LABL(LL+6) = J + ITEN
      QQ(LL+7) = X - TEN*SCF2DF(IP) + 18.00D0*SCF4DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ILEVEN
      LABL(LL+7) = J + ILEVEN
      LL = LL + 7
      GO TO 450

  240 X = GAMMA(N1N2)
      QQ(LL) = 0.1549193D0*R2SDPP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - IFOUR
      QQ(LL+1) = -0.1549193D0*R2SDPP(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - ITHREE
      QQ(LL+2) = ONE*SCG2SD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J
      QQ(LL+3) = -0.1277753D0*R2SDDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J
      QQ(LL+4) = 0.1106567D0*R2SDDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + ITWO
      QQ(LL+5) = -0.1106567D0*R2SDDD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IEIGHT
      LABL(LL+5) = J + ITHREE
      LL = LL + 5
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0414039D0*R2SDPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = 0.1603567D0*R2SDPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.0414039D0*R2SDPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = 0.1603567D0*R2SDPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + ITEN
      QQ(LL+5) = 0.1309307D0*R2SDPF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = -0.1333333D0*R2SDFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ININE
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = 0.1032796D0*R2SDFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ITEN
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = -0.0666667D0*R2SDFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITEN
      LABL(LL+8) = J + ININE
      QQ(LL+9) = -0.1032796D0*R2SDFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ILEVEN
      LABL(LL+9) = J + ISIX
      QQ(LL+10) = -0.0666667D0*R2SDFF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ILEVEN
      LABL(LL+10) = J + ITEN
      LL = LL + 10
      GO TO 450

  250 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = THRE*SCG1PD(IP) - 21.00D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = THRE*SCG1PD(IP) - SIX*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      LL = LL + 2
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1851640D0*R1PDDF(IP) - 0.0132260D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = -0.0853735D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = -0.0534522D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ININE
      QQ(LL+5) = -0.0534522D0*R1PDDF(IP) + 0.0229081D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ISIX
      QQ(LL+6) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFIVE
      LABL(LL+6) = J + ITEN
      QQ(LL+7) = -0.1309307D0*R1PDDF(IP) - 0.0062348D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISIX
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + ISEVEN
      QQ(LL+9) = 0.1690309D0*R1PDDF(IP) - 0.0724418D0*R3PDDF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISEVEN
      LABL(LL+9) = J + IEIGHT
      LL = LL + 9
      GO TO 450

  260 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -THRE*SCG1PD(IP) + SIX*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITHREE
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.1851640D0*R1PDDF(IP) + 0.0132260D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = -0.0853735D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITEN
      QQ(LL+3) = -0.0534522D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = -0.2070197D0*R1PDDF(IP) - 0.0098581D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ITEN
      QQ(LL+5) = 0.0534522D0*R1PDDF(IP) - 0.0229081D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFOUR
      LABL(LL+5) = J + IFIVE
      QQ(LL+6) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFOUR
      LABL(LL+6) = J + ININE
      QQ(LL+7) = -0.1690309D0*R1PDDF(IP) + 0.0724418D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFIVE
      LABL(LL+7) = J + IEIGHT
      QQ(LL+8) = 0.1309307D0*R1PDDF(IP) + 0.0062348D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + IFOUR
      QQ(LL+9) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISIX
      LABL(LL+9) = J + ISEVEN
      LL = LL + 9
      GO TO 450

  270 X = GAMMA(N1N2)
      QQ(LL) = 15.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0623480D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + ININE
      QQ(LL+4) = -0.0381802D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ITEN
      LL = LL + 5
      GO TO 450

  280 X = GAMMA(N1N2)
      QQ(LL) = FOUR*SCF2DD(IP) + 15.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -3.4641016D0*SCF2DD(IP) + 17.3205081D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = 3.4641016D0*SCF2DD(IP) - 17.3205081D0*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + ITHREE
      LL = LL + 2
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 8.9442719D0*SCF2DF(IP) - 6.7082039D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + ISEVEN
      QQ(LL+2) = -6.9282032D0*SCF2DF(IP) + 17.3205081D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 4.4721360D0*SCF2DF(IP) + 20.1246118D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISIX
      LABL(LL+3) = J + ININE
      QQ(LL+4) = 6.9282032D0*SCF2DF(IP) - 17.3205081D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 4.4721360D0*SCF2DF(IP) + 20.1246118D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ITEN
      LL = LL + 5
      GO TO 450

  290 X = GAMMA(N1N2)
      QQ(LL) = X + FOUR*SCF2DD(IP) + 36.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X + FOUR*SCF2DD(IP) - 34.00D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X - TWO*SCF2DD(IP) - FOUR*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = X - TWO*SCF2DD(IP) - FOUR*SCF4DD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      LL = LL + 3
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X - EIGHT*SCF2DF(IP) + SIX*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = X - SIX*SCF2DF(IP) + ONE*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = -27.1108834D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + ININE
      QQ(LL+4) = X - SIX*SCF2DF(IP) + ONE*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 27.1108834D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ITEN
      QQ(LL+6) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = X - 42.00D0*SCF4DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IEIGHT
      QQ(LL+8) = X + TEN*SCF2DF(IP) + THRE*SCF4DF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ININE
      LABL(LL+8) = J + ININE
      QQ(LL+9) = X + TEN*SCF2DF(IP) + THRE*SCF4DF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITEN
      LABL(LL+9) = J + ITEN
      LL = LL + 9
      GO TO 450

  300 X = GAMMA(N1N2)
      QQ(LL) = 0.1549193D0*R2SDPP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - IFOUR
      QQ(LL+1) = ONE*SCG2SD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J
      QQ(LL+2) = -0.1277753D0*R2SDDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J
      QQ(LL+3) = 0.1106567D0*R2SDDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + ITWO
      LL = LL + 3
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0414039D0*R2SDPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = 0.1603567D0*R2SDPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = -0.0414039D0*R2SDPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = -0.1603567D0*R2SDPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = 0.1309307D0*R2SDPF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = -0.1333333D0*R2SDFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ININE
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = 0.1032796D0*R2SDFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ITEN
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = -0.0666667D0*R2SDFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITEN
      LABL(LL+8) = J + ININE
      QQ(LL+9) = 0.0666667D0*R2SDFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ILEVEN
      LABL(LL+9) = J + IEIGHT
      LL = LL + 9
      GO TO 450

  310 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -1.7320508D0*SCG1PD(IP) - 5.1961524D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J - ITWO
      QQ(LL+2) = -THRE*SCG1PD(IP) + 21.00D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J - IONE
      QQ(LL+3) = THRE*SCG1PD(IP) - SIX*SCG3PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ITWO
      LL = LL + 3
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1851640D0*R1PDDF(IP) - 0.0132260D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = -0.0853735D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.0534522D0*R1PDDF(IP) - 0.0229081D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ININE
      QQ(LL+5) = -0.0534522D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IFOUR
      QQ(LL+6) = -0.2070197D0*R1PDDF(IP) - 0.0098581D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFIVE
      LABL(LL+6) = J + IEIGHT
      QQ(LL+7) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISIX
      LABL(LL+7) = J + ISEVEN
      QQ(LL+8) = -0.1309307D0*R1PDDF(IP) - 0.0062348D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISEVEN
      LABL(LL+8) = J + ITHREE
      QQ(LL+9) = -0.1690309D0*R1PDDF(IP) + 0.0724418D0*R3PDDF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISEVEN
      LABL(LL+9) = J + ISIX
      LL = LL + 9
      GO TO 450

  320 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = THRE*SCG1PD(IP) - SIX*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1851640D0*R1PDDF(IP) - 0.0132260D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = 0.0853735D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IEIGHT
      QQ(LL+3) = -0.0534522D0*R1PDDF(IP) + 0.0229081D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = -0.0534522D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFOUR
      LABL(LL+5) = J + IFIVE
      QQ(LL+6) = 0.2070197D0*R1PDDF(IP) + 0.0098581D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFOUR
      LABL(LL+6) = J + ININE
      QQ(LL+7) = -0.1309307D0*R1PDDF(IP) - 0.0062348D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFIVE
      LABL(LL+7) = J + ITHREE
      QQ(LL+8) = 0.1690309D0*R1PDDF(IP) - 0.0724418D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IFIVE
      LABL(LL+8) = J + ISIX
      QQ(LL+9) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISIX
      LABL(LL+9) = J + ISEVEN
      LL = LL + 9
      GO TO 450

  330 X = GAMMA(N1N2)
      QQ(LL) = 15.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0623480D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + ININE
      QQ(LL+4) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = -0.0492904D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IEIGHT
      LL = LL + 5
      GO TO 450

  340 X = GAMMA(N1N2)
      QQ(LL) = FOUR*SCF2DD(IP) + 15.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -3.4641016D0*SCF2DD(IP) + 17.3205081D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITWO
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 8.9442719D0*SCF2DF(IP) - 6.7082039D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + ISEVEN
      QQ(LL+2) = -6.9282032D0*SCF2DF(IP) + 17.3205081D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 4.4721360D0*SCF2DF(IP) + 20.1246118D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISIX
      LABL(LL+3) = J + ININE
      QQ(LL+4) = -4.4721360D0*SCF2DF(IP) - 20.1246118D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + IEIGHT
      LL = LL + 4
      GO TO 450

  350 X = GAMMA(N1N2)
      QQ(LL) = 35.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -27.1108834D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + ININE
      QQ(LL+2) = -27.1108834D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + IEIGHT
      QQ(LL+3) = 35.00D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + ISEVEN
      LL = LL + 3
      GO TO 450

  360 X = GAMMA(N1N2)
      QQ(LL) = X + FOUR*SCF2DD(IP) + 36.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X - TWO*SCF2DD(IP) - FOUR*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X - TWO*SCF2DD(IP) - FOUR*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      LL = LL + 2
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X - EIGHT*SCF2DF(IP) + SIX*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = X - SIX*SCF2DF(IP) + ONE*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = 27.1108834D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IEIGHT
      QQ(LL+4) = X - SIX*SCF2DF(IP) + ONE*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = -27.1108834D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ININE
      QQ(LL+6) = X - 42.00D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + ISIX
      QQ(LL+7) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + ISEVEN
      QQ(LL+8) = X + TEN*SCF2DF(IP) + THRE*SCF4DF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + IEIGHT
      QQ(LL+9) = X + TEN*SCF2DF(IP) + THRE*SCF4DF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ININE
      LABL(LL+9) = J + ININE
      LL = LL + 9
      GO TO 450

  370 X = GAMMA(N1N2)
      QQ(LL) = 0.1549193D0*R2SDPP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - IFOUR
      QQ(LL+1) = ONE*SCG2SD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J
      QQ(LL+2) = 0.0638877D0*R2SDDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J
      QQ(LL+3) = 0.1106567D0*R2SDDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J
      QQ(LL+4) = 0.1106567D0*R2SDDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + IONE
      LL = LL + 4
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.1014185D0*R2SDPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = 0.1309307D0*R2SDPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 0.1309307D0*R2SDPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = 0.1656157D0*R2SDPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ITHREE
      QQ(LL+5) = 0.0421637D0*R2SDFF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ININE
      LABL(LL+5) = J + ITHREE
      QQ(LL+6) = 0.0816497D0*R2SDFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITEN
      LABL(LL+6) = J + IFIVE
      QQ(LL+7) = 0.0816497D0*R2SDFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ILEVEN
      LABL(LL+7) = J + ISIX
      QQ(LL+8) = 0.1054093D0*R2SDFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITWELV
      LABL(LL+8) = J + ISEVEN
      QQ(LL+9) = 0.1054093D0*R2SDFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITHRTN
      LABL(LL+9) = J + IEIGHT
      LL = LL + 9
      GO TO 450

  380 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = THRE*SCG1PD(IP) - SIX*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = 3.4641016D0*SCG1PD(IP) - 15.5884573D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J - ITHREE
      QQ(LL+3) = 15.00D0*SCG3PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J - ITWO
      LL = LL + 3
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2267787D0*R1PDDF(IP) - 0.0431959D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = -0.0623480D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 0.2138090D0*R1PDDF(IP) + 0.0229081D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ITHREE
      QQ(LL+6) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = 0.2138090D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISEVEN
      LABL(LL+8) = J + IEIGHT
      LL = LL + 8
      GO TO 450

  390 X = GAMMA(N1N2)
      QQ(LL) = 15.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 15.00D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J - IONE
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0623480D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + IEIGHT
      QQ(LL+4) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + ITHREE
      QQ(LL+5) = -0.0492904D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ISEVEN
      LL = LL + 5
      GO TO 450

  400 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1851640D0*R1PDDF(IP) + 0.0529040D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = -0.0534522D0*R1PDDF(IP) + 0.0610883D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITHREE
      QQ(LL+3) = 0.2070197D0*R1PDDF(IP) - 0.0394323D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = -0.0534522D0*R1PDDF(IP) + 0.0610883D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = 0.2070197D0*R1PDDF(IP) - 0.0394323D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J + IEIGHT
      QQ(LL+6) = -0.1309307D0*R1PDDF(IP) + 0.0249392D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFOUR
      LABL(LL+6) = J + ITWO
      QQ(LL+7) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFOUR
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IFIVE
      LABL(LL+8) = J + ISIX
      LL = LL + 8
      GO TO 450

  410 X = GAMMA(N1N2)
      QQ(LL) = ONE*SCF2DD(IP) + 30.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 1.7320508D0*SCF2DD(IP) - 8.6602540D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J
      QQ(LL+2) = 1.7320508D0*SCF2DD(IP) - 8.6602540D0*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IONE
      LL = LL + 2
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 1.4142136D0*SCF2DF(IP) + 21.2132034D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = 2.7386128D0*SCF2DF(IP) + 21.9089023D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 2.7386128D0*SCF2DF(IP) + 21.9089023D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = 3.5355339D0*SCF2DF(IP) - 21.2132034D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IEIGHT
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = 3.5355339D0*SCF2DF(IP) - 21.2132034D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ININE
      LABL(LL+5) = J + IEIGHT
      LL = LL + 5
      GO TO 450

  420 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCF2DD(IP) + 20.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = THRE*SCF2DD(IP) - 15.00D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 2.4494897D0*SCF2DF(IP) - 6.1237244D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = -33.2039154D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + ISEVEN
      QQ(LL+3) = 4.7434165D0*SCF2DF(IP) + 4.7434165D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 4.7434165D0*SCF2DF(IP) - 17.3925271D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IEIGHT
      LABL(LL+6) = J + IEIGHT
      LL = LL + 6
      GO TO 450

  430 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCF2DD(IP) + 20.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 2.4494897D0*SCF2DF(IP) - 6.1237244D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = -33.2039154D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + IEIGHT
      QQ(LL+3) = 4.7434165D0*SCF2DF(IP) + 4.7434165D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = -4.7434165D0*SCF2DF(IP) + 17.3925271D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + IEIGHT
      QQ(LL+6) = -6.1237244D0*SCF2DF(IP) - 6.1237244D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ISEVEN
      LL = LL + 6
      GO TO 450

  440 X = GAMMA(N1N2)
      QQ(LL) = X + FOUR*SCF2DD(IP) + 36.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X - TWO*SCF2DD(IP) - FOUR*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      LL = LL + 1
      IF (NT.LT.5) GO TO 450
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X + FOUR*SCF2DF(IP) - 24.00D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = -7.7459667D0*SCF2DF(IP) - 7.7459667D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = X + FNINE*SCF2DF(IP) + 16.00D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = -3.8729833D0*SCF2DF(IP) + 23.2379001D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = X - THRE*SCF2DF(IP) - 24.00D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFOUR
      LABL(LL+5) = J + IFOUR
      QQ(LL+6) = -3.8729833D0*SCF2DF(IP) + 23.2379001D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFOUR
      LABL(LL+6) = J + IEIGHT
      QQ(LL+7) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFIVE
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + ISIX
      QQ(LL+9) = X - FIVE*SCF2DF(IP) - 12.00D0*SCF4DF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISEVEN
      LABL(LL+9) = J + ISEVEN
      QQ(LL+10) = X - FIVE*SCF2DF(IP) - 12.00D0*SCF4DF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + IEIGHT
      LABL(LL+10) = J + IEIGHT
      LL = LL + 10
      GO TO 450

  450 CONTINUE
      RETURN

      END
