      SUBROUTINE GFORMC(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
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
CPK** -----------------------------------------------------------
CPK   ONE PART OF THE FORMER SUBROUTINE GFORM
CPK   NEW SUBROUTINES ARE NAMED GFORMA TO GFORME.
CPK   NEW COMMON 'ALLGF' INSERTED FOR COMMUNICATION WITH THE
CPK   NEW SUBROUTINES.
CPK** -----------------------------------------------------------
c     INCLUDE '../include/param.cmn'
c     INTEGER    ILEVEN, ITWELV, ITHRTN, IFRTEN,IFVTEN
c     PARAMETER  (ILEVEN=11,ITWELV=12,ITHRTN=13,IFRTEN=14,IFVTEN=15)
c     INCLUDE '../include/allgf.cmn'
c     INCLUDE '../include/ftch.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/z.cmn'
c     DIMENSION    GAMMA(1)
C
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
c     write(6,*)"i am in gformc"
      ICS1 = ICS - 36
      THRE = THREE
      SEVN = SEVEN
      GO TO (10,20,30,40,50,60,70,80,90,
     .       100,110,120,130,140,150,160,170,180,
     .       190,200,210,220,230,240,250,260), ICS1

      GO TO 530

   10 GO TO 270

   20 GO TO 280

   30 GO TO 290

   40 GO TO 300

   50 GO TO 310

   60 GO TO 320

   70 GO TO 330

   80 GO TO 340

   90 GO TO 350

  100 GO TO 360

  110 GO TO 370

  120 GO TO 380

  130 GO TO 390

  140 GO TO 400

  150 GO TO 410

  160 GO TO 420

  170 GO TO 430

  180 GO TO 440

  190 GO TO 450

  200 GO TO 460

  210 GO TO 470

  220 GO TO 480

  230 GO TO 490

  240 GO TO 500

  250 GO TO 510

  260 GO TO 520

  270 X = GAMMA(N1N2)
      QQ(LL) = 0.1549193D0*R2SDPP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + ITWO
      LABL(LL) = J - IFIVE
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
      QQ(LL+3) = -0.1106567D0*R2SDDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J
      QQ(LL+4) = 0.1106567D0*R2SDDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J - IONE
      LL = LL + 4
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1309307D0*R2SDPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = -0.1014185D0*R2SDPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IONE
      QQ(LL+3) = -0.1309307D0*R2SDPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + IFOUR
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
      QQ(LL+7) = -0.0816497D0*R2SDFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ILEVEN
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.1054093D0*R2SDFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITWELV
      LABL(LL+8) = J + ISEVEN
      QQ(LL+9) = -0.1054093D0*R2SDFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITHRTN
      LABL(LL+9) = J + ISIX
      LL = LL + 9
      GO TO 530

  280 X = GAMMA(N1N2)
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
      QQ(LL+2) = 15.00D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J - ITWO
      LL = LL + 2
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0623480D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + ITHREE
      QQ(LL+3) = 0.0492904D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISIX
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = 0.0381802D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + ITWO
      QQ(LL+5) = -0.0492904D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ISIX
      LL = LL + 5
      GO TO 530

  290 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 3.4641016D0*SCG1PD(IP) - 15.5884573D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J - IFOUR
      QQ(LL+2) = -15.00D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J - ITHREE
      LL = LL + 2
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2267787D0*R1PDDF(IP) - 0.0431959D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IONE
      QQ(LL+2) = 0.0623480D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + IONE
      QQ(LL+3) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = 0.2138090D0*R1PDDF(IP) - 0.0534522D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ITWO
      QQ(LL+6) = -0.0492904D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFIVE
      LABL(LL+6) = J + ISIX
      QQ(LL+7) = 0.2138090D0*R1PDDF(IP) + 0.0229081D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISIX
      LABL(LL+7) = J + ITHREE
      QQ(LL+8) = -0.0492904D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + ISEVEN
      LL = LL + 8
      GO TO 530

  300 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1PD(IP) + 24.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1851640D0*R1PDDF(IP) + 0.0529040D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = 0.0534522D0*R1PDDF(IP) - 0.0610883D0*R3PDDF(IP)
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
      LABL(LL+4) = J + ITWO
      QQ(LL+5) = -0.2070197D0*R1PDDF(IP) + 0.0394323D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J + ISIX
      QQ(LL+6) = 0.1690309D0*R1PDDF(IP) + 0.0482945D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFOUR
      LABL(LL+6) = J + IFIVE
      QQ(LL+7) = -0.1309307D0*R1PDDF(IP) + 0.0249392D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFIVE
      LABL(LL+7) = J + IONE
      QQ(LL+8) = -0.1690309D0*R1PDDF(IP) - 0.0482945D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IFIVE
      LABL(LL+8) = J + IFOUR
      LL = LL + 8
      GO TO 530

  310 X = GAMMA(N1N2)
      QQ(LL) = ONE*SCF2DD(IP) + 30.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -1.7320508D0*SCF2DD(IP) + 8.6602540D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J
      QQ(LL+2) = 1.7320508D0*SCF2DD(IP) - 8.6602540D0*SCF4DD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J - IONE
      LL = LL + 2
      IF (NT.LT.5) GO TO 530
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
      QQ(LL+3) = -2.7386128D0*SCF2DF(IP) - 21.9089023D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = 3.5355339D0*SCF2DF(IP) - 21.2132034D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IEIGHT
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = -3.5355339D0*SCF2DF(IP) + 21.2132034D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ININE
      LABL(LL+5) = J + ISIX
      LL = LL + 5
      GO TO 530

  320 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCF2DD(IP) + 20.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -THRE*SCF2DD(IP) + 15.00D0*SCF4DD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J - IONE
      LL = LL + 1
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -2.4494897D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = -33.2039154D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + ISEVEN
      QQ(LL+3) = -4.7434165D0*SCF2DF(IP) + 17.3925271D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 4.7434165D0*SCF2DF(IP) + 4.7434165D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = -6.1237244D0*SCF2DF(IP) - 6.1237244D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IEIGHT
      LABL(LL+6) = J + ISIX
      LL = LL + 6
      GO TO 530

  330 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCF2DD(IP) + 20.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 2.4494897D0*SCF2DF(IP) - 6.1237244D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = 33.2039154D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + ISIX
      QQ(LL+3) = 4.7434165D0*SCF2DF(IP) - 17.3925271D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = 4.7434165D0*SCF2DF(IP) + 4.7434165D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ISIX
      QQ(LL+6) = 6.1237244D0*SCF2DF(IP) + 6.1237244D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ISEVEN
      LL = LL + 6
      GO TO 530

  340 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCF2DD(IP) + 20.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -7.7459667D0*SCF2DF(IP) - 7.7459667D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = SIX*SCF2DF(IP) + 20.00D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + ITHREE
      QQ(LL+3) = -3.8729833D0*SCF2DF(IP) + 23.2379001D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = 3.8729833D0*SCF2DF(IP) - 23.2379001D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ISIX
      LL = LL + 4
      GO TO 530

  350 X = GAMMA(N1N2)
      QQ(LL) = X + FOUR*SCF2DD(IP) + 36.00D0*SCF4DD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 530
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X + FOUR*SCF2DF(IP) - 24.00D0*SCF4DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = 7.7459667D0*SCF2DF(IP) + 7.7459667D0*SCF4DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = X - THRE*SCF2DF(IP) - 24.00D0*SCF4DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ITWO
      QQ(LL+4) = 3.8729833D0*SCF2DF(IP) - 23.2379001D0*SCF4DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = X + FNINE*SCF2DF(IP) + 16.00D0*SCF4DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J + ITHREE
      QQ(LL+6) = 3.8729833D0*SCF2DF(IP) - 23.2379001D0*SCF4DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITHREE
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFOUR
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = X + 28.00D0*SCF4DF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IFIVE
      LABL(LL+8) = J + IFIVE
      QQ(LL+9) = X - FIVE*SCF2DF(IP) - 12.00D0*SCF4DF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISIX
      LABL(LL+9) = J + ISIX
      QQ(LL+10) = X - FIVE*SCF2DF(IP) - 12.00D0*SCF4DF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ISEVEN
      LABL(LL+10) = J + ISEVEN
      LL = LL + 10
      GO TO 530

  360 X = GAMMA(N1N2)
      QQ(LL) = -0.0724418D0*R3SFPD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - ITWO
      QQ(LL+1) = -0.0724418D0*R3SFPD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - IONE
      QQ(LL+2) = 0.1254729D0*R3SFPD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J - IFIVE
      QQ(LL+3) = ONE*SCG3SF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J
      QQ(LL+4) = 0.0851835D0*R3SFDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J
      QQ(LL+5) = -0.0952381D0*R3SFDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ITHREE
      QQ(LL+6) = -0.0952381D0*R3SFDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + IFOUR
      QQ(LL+7) = 0.0301169D0*R3SFDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + IONE
      QQ(LL+8) = 0.0301169D0*R3SFDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ITWO
      LL = LL + 8
      GO TO 530

  370 X = GAMMA(N1N2)
      QQ(LL) = -0.0323970D0*R2PFDD(IP) + 0.0599944D0*R4PFDD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + ITHREE
      LABL(LL) = J - ITWO
      QQ(LL+1) = -0.0561132D0*R2PFDD(IP) - 0.0173189D0*R4PFDD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J - ITWO
      QQ(LL+2) = -0.0561132D0*R2PFDD(IP) - 0.0173189D0*R4PFDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J - IONE
      QQ(LL+3) = FNINE*SCG2PF(IP) + TEN*SCG4PF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J
      QQ(LL+4) = -11.6189500D0*SCG2PF(IP) - 3.8729833D0*SCG4PF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I
      LABL(LL+4) = J + ITHREE
      QQ(LL+5) = -11.6189500D0*SCG2PF(IP) - 3.8729833D0*SCG4PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IONE
      LABL(LL+5) = J + IFOUR
      QQ(LL+6) = -14.6969385D0*SCG2PF(IP) + 12.2474487D0*SCG4PF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITWO
      LABL(LL+6) = J + IONE
      QQ(LL+7) = -0.0213809D0*R2PFFF(IP) + 0.0269961D0*R4PFFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IONE
      QQ(LL+8) = -0.0414039D0*R2PFFF(IP) + 0.0278814D0*R4PFFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ININE
      LABL(LL+8) = J + ITHREE
      QQ(LL+9) = -0.0414039D0*R2PFFF(IP) + 0.0278814D0*R4PFFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITEN
      LABL(LL+9) = J + IFOUR
      QQ(LL+10) = -0.0534522D0*R2PFFF(IP) - 0.0269961D0*R4PFFF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ILEVEN
      LABL(LL+10) = J + IFIVE
      QQ(LL+11) = -0.0534522D0*R2PFFF(IP) - 0.0269961D0*R4PFFF(IP)
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + ITWELV
      LABL(LL+11) = J + ISIX
      LL = LL + 11
      GO TO 530

  380 X = GAMMA(N1N2)
      QQ(LL) = -0.0323970D0*R2PFDD(IP) + 0.0599944D0*R4PFDD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + ITWO
      LABL(LL) = J - IONE
      QQ(LL+1) = 0.0561132D0*R2PFDD(IP) + 0.0173189D0*R4PFDD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J - IONE
      QQ(LL+2) = -0.0561132D0*R2PFDD(IP) - 0.0173189D0*R4PFDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J - ITWO
      QQ(LL+3) = FNINE*SCG2PF(IP) + TEN*SCG4PF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J
      QQ(LL+4) = 11.6189500D0*SCG2PF(IP) + 3.8729833D0*SCG4PF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I
      LABL(LL+4) = J + ITHREE
      QQ(LL+5) = -14.6969385D0*SCG2PF(IP) + 12.2474487D0*SCG4PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IONE
      LABL(LL+5) = J + ITWO
      QQ(LL+6) = -0.0213809D0*R2PFFF(IP) + 0.0269961D0*R4PFFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ITWO
      QQ(LL+7) = -0.0414039D0*R2PFFF(IP) + 0.0278814D0*R4PFFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.0414039D0*R2PFFF(IP) - 0.0278814D0*R4PFFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ININE
      LABL(LL+8) = J + ITHREE
      QQ(LL+9) = -0.0534522D0*R2PFFF(IP) - 0.0269961D0*R4PFFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITEN
      LABL(LL+9) = J + ISIX
      QQ(LL+10) = 0.0534522D0*R2PFFF(IP) + 0.0269961D0*R4PFFF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ILEVEN
      LABL(LL+10) = J + IFIVE
      LL = LL + 10
      GO TO 530

  390 X = GAMMA(N1N2)
      QQ(LL) = 0.1122263D0*R2PFDD(IP) + 0.0831306D0*R4PFDD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - IFIVE
      QQ(LL+1) = -0.1122263D0*R2PFDD(IP) + 0.0138551D0*R4PFDD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - IFOUR
      QQ(LL+2) = -0.1122263D0*R2PFDD(IP) + 0.0138551D0*R4PFDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J - ITHREE
      QQ(LL+3) = 0.0561132D0*R2PFDD(IP) - 0.0554204D0*R4PFDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J - ITWO
      QQ(LL+4) = 0.0561132D0*R2PFDD(IP) - 0.0554204D0*R4PFDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J - IONE
      QQ(LL+5) = 27.00D0*SCG2PF(IP) + 16.00D0*SCG4PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I
      LABL(LL+5) = J
      QQ(LL+6) = 0.1047446D0*R2PFFF(IP) + 0.0529013D0*R4PFFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J
      QQ(LL+7) = 0.0785584D0*R2PFFF(IP) + 0.0088169D0*R4PFFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + IONE
      QQ(LL+8) = 0.0785584D0*R2PFFF(IP) + 0.0088169D0*R4PFFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ITWO
      QQ(LL+9) = -0.0617182D0*R4PFFF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ININE
      LABL(LL+9) = J + ITHREE
      QQ(LL+10) = -0.0617182D0*R4PFFF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ITEN
      LABL(LL+10) = J + IFOUR
      QQ(LL+11) = -0.1309307D0*R2PFFF(IP) + 0.0264507D0*R4PFFF(IP)
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + ILEVEN
      LABL(LL+11) = J + IFIVE
      QQ(LL+12) = -0.1309307D0*R2PFFF(IP) + 0.0264507D0*R4PFFF(IP)
      LABI(LL+12) = I
      LABJ(LL+12) = J
      LABK(LL+12) = I + ITWELV
      LABL(LL+12) = J + ISIX
      LL = LL + 12
      GO TO 530

  400 X = GAMMA(N1N2)
      QQ(LL) = FNINE*SCG1DF(IP) + 16.00D0*SCG3DF(IP) +
     .         1000.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 6.7082039D0*SCG1DF(IP) - 17.8885438D0*SCG3DF(IP) +
     .           223.6067977D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = 6.7082039D0*SCG1DF(IP) - 17.8885438D0*SCG3DF(IP) +
     .           223.6067977D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = 8.4852814D0*SCG1DF(IP) + 5.6568542D0*SCG3DF(IP) -
     .           707.1067812D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + IONE
      QQ(LL+4) = 8.4852814D0*SCG1DF(IP) + 5.6568542D0*SCG3DF(IP) -
     .           707.1067812D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ITWO
      LL = LL + 4
      GO TO 530

  410 X = GAMMA(N1N2)
      QQ(LL) = 20.00D0*SCG3DF(IP) + 350.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -12.2474487D0*SCG3DF(IP) + 428.6607050D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IONE
      QQ(LL+2) = -15.8113883D0*SCG3DF(IP) - 110.6797181D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 12.2474487D0*SCG3DF(IP) - 428.6607050D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITWO
      QQ(LL+4) = -15.8113883D0*SCG3DF(IP) - 110.6797181D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ISIX
      LL = LL + 4
      GO TO 530

  420 X = GAMMA(N1N2)
      QQ(LL) = 20.00D0*SCG3DF(IP) + 350.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -12.2474487D0*SCG3DF(IP) + 428.6607050D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = -15.8113883D0*SCG3DF(IP) - 110.6797181D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + ISIX
      QQ(LL+3) = -12.2474487D0*SCG3DF(IP) + 428.6607050D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + IONE
      QQ(LL+4) = 15.8113883D0*SCG3DF(IP) + 110.6797181D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + IFIVE
      LL = LL + 4
      GO TO 530

  430 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1DF(IP) + TWO*SCG3DF(IP) + 800.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -3.8729833D0*SCG1DF(IP) + 3.8729833D0*SCG3DF(IP) -
     .           309.8386677D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = -3.8729833D0*SCG1DF(IP) + 3.8729833D0*SCG3DF(IP) -
     .           309.8386677D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFOUR
      LL = LL + 2
      GO TO 530

  440 X = GAMMA(N1N2)
      QQ(LL) = THRE*SCG1DF(IP) + TWO*SCG3DF(IP) + 800.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 3.8729833D0*SCG1DF(IP) - 3.8729833D0*SCG3DF(IP) +
     .           309.8386677D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ITHREE
      LL = LL + 1
      GO TO 530

  450 X = GAMMA(N1N2)
      QQ(LL) = X + 16.00D0*SCF2FF(IP) + 36.00D0*SCF4FF(IP) +
     .         40000.00D0*SCF6FF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X + 12.00D0*SCF2FF(IP) + SIX*SCF4FF(IP) -
     .           30000.00D0*SCF6FF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X + 12.00D0*SCF2FF(IP) + SIX*SCF4FF(IP) -
     .           30000.00D0*SCF6FF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = X - 42.00D0*SCF4FF(IP) + 12000.00D0*SCF6FF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = X - 42.00D0*SCF4FF(IP) + 12000.00D0*SCF6FF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = X - 20.00D0*SCF2FF(IP) + 18.00D0*SCF4FF(IP) -
     .           2000.00D0*SCF6FF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IFIVE
      QQ(LL+6) = X - 20.00D0*SCF2FF(IP) + 18.00D0*SCF4FF(IP) -
     .           2000.00D0*SCF6FF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + ISIX
      LL = LL + 6
      GO TO 530

  460 X = GAMMA(N1N2)
      QQ(LL) = 0.1024482D0*R3SFPD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - ISIX
      QQ(LL+1) = -0.0295742D0*R3SFPD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J - IFIVE
      QQ(LL+2) = -0.0295742D0*R3SFPD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J - IFOUR
      QQ(LL+3) = 0.1182970D0*R3SFPD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J - ITHREE
      QQ(LL+4) = ONE*SCG3SF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I
      LABL(LL+4) = J
      QQ(LL+5) = 0.0638877D0*R3SFDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFOUR
      LABL(LL+5) = J
      QQ(LL+6) = 0.0737711D0*R3SFDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFIVE
      LABL(LL+6) = J
      QQ(LL+7) = -0.0476190D0*R3SFDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFIVE
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.0737711D0*R3SFDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + IONE
      QQ(LL+9) = -0.0476190D0*R3SFDF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ISIX
      LABL(LL+9) = J + IFIVE
      QQ(LL+10) = 0.0301169D0*R3SFDF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ISEVEN
      LABL(LL+10) = J - IONE
      QQ(LL+11) = 0.0583212D0*R3SFDF(IP)
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + ISEVEN
      LABL(LL+11) = J + ITWO
      QQ(LL+12) = 0.0583212D0*R3SFDF(IP)
      LABI(LL+12) = I
      LABJ(LL+12) = J
      LABK(LL+12) = I + IEIGHT
      LABL(LL+12) = J + ITHREE
      LL = LL + 12
      GO TO 530

  470 X = GAMMA(N1N2)
      QQ(LL) = -0.0320713D0*R2PPPF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - IEIGHT
      QQ(LL+1) = 0.1282854D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - ISEVEN
      QQ(LL+2) = 0.0916324D0*R2PFDD(IP) - 0.0509069D0*R4PFDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J - ISIX
      QQ(LL+3) = 0.0264520D0*R2PFDD(IP) + 0.0367389D0*R4PFDD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J - IFIVE
      QQ(LL+4) = -0.0916324D0*R2PFDD(IP) - 0.0084845D0*R4PFDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J - IFIVE
      QQ(LL+5) = -0.0916324D0*R2PFDD(IP) - 0.0084845D0*R4PFDD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J - IFOUR
      QQ(LL+6) = 0.0229081D0*R2PFDD(IP) + 0.0763604D0*R4PFDD(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J - ITHREE
      QQ(LL+7) = 0.0687243D0*R2PFDD(IP) - 0.0084845D0*R4PFDD(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J - ITWO
      QQ(LL+8) = 19.5000000D0*SCG2PF(IP) + 13.5000000D0*SCG4PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I
      LABL(LL+8) = J
      QQ(LL+9) = -5.8094750D0*SCG2PF(IP) - 1.9364917D0*SCG4PF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I
      LABL(LL+9) = J + IFOUR
      QQ(LL+10) = 16.5000000D0*SCG2PF(IP) - 1.5000000D0*SCG4PF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + IONE
      LABL(LL+10) = J + IONE
      QQ(LL+11) = -5.8094750D0*SCG2PF(IP) - 1.9364917D0*SCG4PF(IP)
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + IONE
      LABL(LL+11) = J + IFIVE
      QQ(LL+12) = 22.0454077D0*SCG2PF(IP) - 9.7979590D0*SCG4PF(IP)
      LABI(LL+12) = I
      LABJ(LL+12) = J
      LABK(LL+12) = I + ITWO
      LABL(LL+12) = J - IONE
      QQ(LL+13) = -4.7434165D0*SCG2PF(IP) + 9.4868330D0*SCG4PF(IP)
      LABI(LL+13) = I
      LABJ(LL+13) = J
      LABK(LL+13) = I + ITWO
      LABL(LL+13) = J + ITWO
      QQ(LL+14) = 0.0855236D0*R2PFFF(IP) - 0.0323953D0*R4PFFF(IP)
      LABI(LL+14) = I
      LABJ(LL+14) = J
      LABK(LL+14) = I + IEIGHT
      LABL(LL+14) = J - IONE
      QQ(LL+15) = 0.0276026D0*R2PFFF(IP) - 0.0104555D0*R4PFFF(IP)
      LABI(LL+15) = I
      LABJ(LL+15) = J
      LABK(LL+15) = I + IEIGHT
      LABL(LL+15) = J + ITWO
      QQ(LL+16) = 0.0427618D0*R2PFFF(IP) + 0.0215969D0*R4PFFF(IP)
      LABI(LL+16) = I
      LABJ(LL+16) = J
      LABK(LL+16) = I + ININE
      LABL(LL+16) = J
      QQ(LL+17) = 0.0138013D0*R2PFFF(IP) + 0.0313666D0*R4PFFF(IP)
      LABI(LL+17) = I
      LABJ(LL+17) = J
      LABK(LL+17) = I + ININE
      LABL(LL+17) = J + IFOUR
      QQ(LL+18) = 0.0855236D0*R2PFFF(IP) - 0.0323953D0*R4PFFF(IP)
      LABI(LL+18) = I
      LABJ(LL+18) = J
      LABK(LL+18) = I + ITEN
      LABL(LL+18) = J + IONE
      QQ(LL+19) = 0.0138013D0*R2PFFF(IP) + 0.0313666D0*R4PFFF(IP)
      LABI(LL+19) = I
      LABJ(LL+19) = J
      LABK(LL+19) = I + ITEN
      LABL(LL+19) = J + IFIVE
      QQ(LL+20) = 0.0377945D0*R4PFFF(IP)
      LABI(LL+20) = I
      LABJ(LL+20) = J
      LABK(LL+20) = I + ILEVEN
      LABL(LL+20) = J + ITWO
      QQ(LL+21) = 0.0377945D0*R4PFFF(IP)
      LABI(LL+21) = I
      LABJ(LL+21) = J
      LABK(LL+21) = I + ITWELV
      LABL(LL+21) = J + ITHREE
      QQ(LL+22) = -0.1069045D0*R2PFFF(IP) - 0.0161977D0*R4PFFF(IP)
      LABI(LL+22) = I
      LABJ(LL+22) = J
      LABK(LL+22) = I + ITHRTN
      LABL(LL+22) = J + IFOUR
      QQ(LL+23) = -0.1069045D0*R2PFFF(IP) - 0.0161977D0*R4PFFF(IP)
      LABI(LL+23) = I
      LABJ(LL+23) = J
      LABK(LL+23) = I + IFRTEN
      LABL(LL+23) = J + IFIVE
      LL = LL + 23
      GO TO 530

  480 X = GAMMA(N1N2)
      QQ(LL) = 0.0264520D0*R2PFDD(IP) + 0.0367389D0*R4PFDD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + ITWO
      LABL(LL) = J - IFOUR
      QQ(LL+1) = -0.0229081D0*R2PFDD(IP) + 0.0424224D0*R4PFDD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFIVE
      LABL(LL+1) = J - ITWO
      QQ(LL+2) = 1.5000000D0*SCG2PF(IP) + 7.5000000D0*SCG4PF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J
      QQ(LL+3) = 5.8094750D0*SCG2PF(IP) + 1.9364917D0*SCG4PF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = -4.7434165D0*SCG2PF(IP) + 9.4868330D0*SCG4PF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IONE
      LABL(LL+4) = J + ITHREE
      QQ(LL+5) = 0.0276026D0*R2PFFF(IP) - 0.0104555D0*R4PFFF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ITHREE
      QQ(LL+6) = -0.0213809D0*R2PFFF(IP) + 0.0269961D0*R4PFFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IEIGHT
      LABL(LL+6) = J + IONE
      QQ(LL+7) = 0.0138013D0*R2PFFF(IP) + 0.0313666D0*R4PFFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = -0.0138013D0*R2PFFF(IP) - 0.0313666D0*R4PFFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ININE
      LABL(LL+8) = J + IFOUR
      LL = LL + 8
      GO TO 530

  490 X = GAMMA(N1N2)
      QQ(LL) = 0.0529040D0*R2PFDD(IP) + 0.0734778D0*R4PFDD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - ITHREE
      QQ(LL+1) = 0.0916324D0*R2PFDD(IP) - 0.0212112D0*R4PFDD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - ITHREE
      QQ(LL+2) = 0.0916324D0*R2PFDD(IP) - 0.0212112D0*R4PFDD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J - ITWO
      QQ(LL+3) = 24.00D0*SCG2PF(IP) + 15.00D0*SCG4PF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J
      QQ(LL+4) = 0.0349149D0*R2PFFF(IP) + 0.0330633D0*R4PFFF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J
      QQ(LL+5) = 0.0676123D0*R2PFFF(IP) + 0.0341476D0*R4PFFF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ITWO
      QQ(LL+6) = 0.0676123D0*R2PFFF(IP) + 0.0341476D0*R4PFFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IEIGHT
      LABL(LL+6) = J + ITHREE
      QQ(LL+7) = 0.0872872D0*R2PFFF(IP) - 0.0330633D0*R4PFFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ININE
      LABL(LL+7) = J + IFOUR
      QQ(LL+8) = 0.0872872D0*R2PFFF(IP) - 0.0330633D0*R4PFFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITEN
      LABL(LL+8) = J + IFIVE
      LL = LL + 8
      GO TO 530

  500 X = GAMMA(N1N2)
      QQ(LL) = SIX*SCG1DF(IP) + FNINE*SCG3DF(IP) + 900.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -1.7320508D0*SCG1DF(IP) + 10.3923048D0*SCG3DF(IP) -
     .           259.8076211D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J
      QQ(LL+2) = 6.7082039D0*SCG1DF(IP) - 6.7082039D0*SCG3DF(IP) +
     .           67.0820393D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = -1.7320508D0*SCG1DF(IP) + 10.3923048D0*SCG3DF(IP) -
     .           259.8076211D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + IONE
      QQ(LL+4) = 6.7082039D0*SCG1DF(IP) - 6.7082039D0*SCG3DF(IP) +
     .           67.0820393D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = -4.2426407D0*SCG1DF(IP) + 4.2426407D0*SCG3DF(IP) +
     .           848.5281374D0*SCG5DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHREE
      LABL(LL+5) = J - IONE
      QQ(LL+6) = 5.4772256D0*SCG1DF(IP) + 8.2158384D0*SCG3DF(IP) -
     .           328.6335345D0*SCG5DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITHREE
      LABL(LL+6) = J + ITWO
      QQ(LL+7) = 5.4772256D0*SCG1DF(IP) + 8.2158384D0*SCG3DF(IP) -
     .           328.6335345D0*SCG5DF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFOUR
      LABL(LL+7) = J + ITHREE
      LL = LL + 7
      GO TO 530

  510 X = GAMMA(N1N2)
      QQ(LL) = 0.5000000D0*SCG1DF(IP) + 17.00D0*SCG3DF(IP) +
     .         425.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -1.9364917D0*SCG1DF(IP) - 7.7459667D0*SCG3DF(IP) -
     .           19.3649167D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = 0.5000000D0*SCG1DF(IP) + SEVN*SCG3DF(IP) -
     .           275.00D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IONE
      QQ(LL+3) = -1.9364917D0*SCG1DF(IP) - 7.7459667D0*SCG3DF(IP) -
     .           19.3649167D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 1.2247449D0*SCG1DF(IP) + 4.8989795D0*SCG3DF(IP) -
     .           244.9489743D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J - IONE
      QQ(LL+5) = -1.5811388D0*SCG1DF(IP) + 1.5811388D0*SCG3DF(IP) +
     .           537.5872022D0*SCG5DF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITWO
      LABL(LL+5) = J + ITWO
      QQ(LL+6) = -1.5811388D0*SCG1DF(IP) + 17.3925271D0*SCG3DF(IP) -
     .           347.8505426D0*SCG5DF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITHREE
      LABL(LL+6) = J + ITHREE
      LL = LL + 6
      GO TO 530

  520 X = GAMMA(N1N2)
      QQ(LL) = 0.5000000D0*SCG1DF(IP) + 17.00D0*SCG3DF(IP) +
     .         425.00D0*SCG5DF(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = 1.9364917D0*SCG1DF(IP) + 7.7459667D0*SCG3DF(IP) +
     .           19.3649167D0*SCG5DF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = -1.5811388D0*SCG1DF(IP) + 1.5811388D0*SCG3DF(IP) +
     .           537.5872022D0*SCG5DF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + ITHREE
      QQ(LL+3) = 1.2247449D0*SCG1DF(IP) + 4.8989795D0*SCG3DF(IP) -
     .           244.9489743D0*SCG5DF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J - IONE
      QQ(LL+4) = 1.5811388D0*SCG1DF(IP) - 17.3925271D0*SCG3DF(IP) +
     .           347.8505426D0*SCG5DF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + ITWO
      LL = LL + 4
      GO TO 530

  530 CONTINUE
c     write(6,*)"NIN in gformc",(NIN(I),I=1,NB)
      RETURN

      END
