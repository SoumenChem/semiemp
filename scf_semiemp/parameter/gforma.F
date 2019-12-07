      SUBROUTINE GFORMA(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
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
C     ...................................................
C
C     THIS BIT OF NONSENSE IS NECESSARY SINCE UNIX FORTRAN
C     GENERATES ONLY "SHORT" JUMPS FOR COMPUTED GOTOS
C     GOOD FOR ONLY A FEW HUNDRED LINES.
C     THE STANDARD GOTOS GENERATE A PROPER JUMP.
C
c      if (ga_nodeid().eq.0)
c    .write(6,*)"i am in gforma",I,J,LL,N1,N2,NT,N1N2
      GO TO (10,20,30,40,50,60,70,80,90,
     .       100,110,120,130,140), ICS

      WRITE (6,9000) ICS

 9000 FORMAT (' ERROR SUBROUTINE GFORMA.  ICS= ',I5)

      CALL ZEXIT(3)
   10 GO TO 150

   20 GO TO 160

   30 GO TO 170

   40 GO TO 180

   50 GO TO 190

   60 GO TO 200

   70 GO TO 210

   80 GO TO 220

   90 GO TO 230

  100 GO TO 240

  110 GO TO 250

  120 GO TO 260

  130 GO TO 270

  140 GO TO 280

  150 X = GAMMA(N1N2)
c     write(6,*)"i am here 1",X,LL,I,J
      QQ(LL) = X
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.1) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = X
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      LL = LL + 3
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = X
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = X
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISIX
      LABL(LL+3) = J + ISIX
      QQ(LL+4) = X
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = X
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IEIGHT
      LABL(LL+5) = J + IEIGHT
      LL = LL + 5
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 3
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ININE
      LABL(LL+1) = J + ININE
      QQ(LL+2) = X
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITEN
      LABL(LL+2) = J + ITEN
      QQ(LL+3) = X
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ILEVEN
      LABL(LL+3) = J + ILEVEN
      QQ(LL+4) = X
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWELV
      LABL(LL+4) = J + ITWELV
      QQ(LL+5) = X
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITHRTN
      LABL(LL+5) = J + ITHRTN
      QQ(LL+6) = X
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFRTEN
      LABL(LL+6) = J + IFRTEN
      QQ(LL+7) = X
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFVTEN
      LABL(LL+7) = J + IFVTEN
      LL = LL + 7
      GO TO 290

  160 X = GAMMA(N1N2)
c     write(6,*)"i am here 2"      
      QQ(LL) = SCG1SP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.1490712D0*R1SPPD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = FP258*R1SPPD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = FP258*R1SPPD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = FP258*R1SPPD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ISIX
      LL = LL + 4
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2390457D0*R1SPDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ININE
      QQ(LL+2) = -0.0690066D0*R1SPDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.2672612D0*R1SPDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + ITHRTN
      QQ(LL+4) = -0.0690066D0*R1SPDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + ITEN
      QQ(LL+5) = 0.2672612D0*R1SPDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + IFRTEN
      QQ(LL+6) = -0.1690309D0*R1SPDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + IEIGHT
      QQ(LL+7) = FP218*R1SPDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + ILEVEN
      QQ(LL+8) = FP218*R1SPDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ITWELV
      LL = LL + 8
      GO TO 290

  170 X = GAMMA(N1N2)
c     write(6,*)"i am here 3"
      QQ(LL) = X + F4P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X - F2P0*SCF2PP(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X - F2P0*SCF2PP(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      LL = LL + 2
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X - F2P0*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ITHREE
      QQ(LL+2) = -3.4641016D0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = X + F2P0*SCF2PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFOUR
      QQ(LL+4) = X + F2P0*SCF2PD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = X + F2P0*SCF2PD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ISIX
      QQ(LL+6) = X - F4P0*SCF2PD(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ISEVEN
      LL = LL + 6
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0962140D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ININE
      QQ(LL+2) = 0.1242118D0*R2PPPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J + ITHRTN
      QQ(LL+3) = -0.0320713D0*R2PPPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + ITEN
      QQ(LL+4) = 0.1242118D0*R2PPPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IONE
      LABL(LL+4) = J + IFRTEN
      QQ(LL+5) = -0.0785584D0*R2PPPF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITWO
      LABL(LL+5) = J + IEIGHT
      QQ(LL+6) = 0.1014185D0*R2PPPF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITWO
      LABL(LL+6) = J + ILEVEN
      QQ(LL+7) = X - F4P0*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IEIGHT
      QQ(LL+8) = -7.7459667D0*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ILEVEN
      QQ(LL+9) = X + F3P0*SCF2PF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ININE
      LABL(LL+9) = J + ININE
      QQ(LL+10) = -3.8729833D0*SCF2PF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ININE
      LABL(LL+10) = J + ITHRTN
      QQ(LL+11) = X - FNINE*SCF2PF(IP)
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + ITEN
      LABL(LL+11) = J + ITEN
      QQ(LL+12) = -3.8729833D0*SCF2PF(IP)
      LABI(LL+12) = I
      LABJ(LL+12) = J
      LABK(LL+12) = I + ITEN
      LABL(LL+12) = J + IFRTEN
      QQ(LL+13) = X
      LABI(LL+13) = I
      LABJ(LL+13) = J
      LABK(LL+13) = I + ILEVEN
      LABL(LL+13) = J + ILEVEN
      QQ(LL+14) = X
      LABI(LL+14) = I
      LABJ(LL+14) = J
      LABK(LL+14) = I + ITWELV
      LABL(LL+14) = J + ITWELV
      QQ(LL+15) = X + FIVE*SCF2PF(IP)
      LABI(LL+15) = I
      LABJ(LL+15) = J
      LABK(LL+15) = I + ITHRTN
      LABL(LL+15) = J + ITHRTN
      QQ(LL+16) = X + FIVE*SCF2PF(IP)
      LABI(LL+16) = I
      LABJ(LL+16) = J
      LABK(LL+16) = I + IFRTEN
      LABL(LL+16) = J + IFRTEN
      LL = LL + 16
      GO TO 290

  180 X = GAMMA(N1N2)
c     write(6,*)"i am here 4"
      QQ(LL) = ONE*SCG1SP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = FP258*R1SPPD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = -0.1490712D0*R1SPPD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = -FP258*R1SPPD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = FP258*R1SPPD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITHREE
      LABL(LL+4) = J + ISIX
      LL = LL + 4
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2390457D0*R1SPDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ININE
      QQ(LL+2) = 0.0690066D0*R1SPDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.2672612D0*R1SPDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + ITHRTN
      QQ(LL+4) = -0.0690066D0*R1SPDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISIX
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = -0.2672612D0*R1SPDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISIX
      LABL(LL+5) = J + ITWELV
      QQ(LL+6) = FP218*R1SPDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ILEVEN
      QQ(LL+7) = -0.1690309D0*R1SPDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + ISEVEN
      QQ(LL+8) = -FP218*R1SPDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ITEN
      LL = LL + 8
      GO TO 290

  190 X = GAMMA(N1N2)
c     write(6,*)"i am here 5"
      QQ(LL) = F3P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -3.4641016D0*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = F3P0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + ISIX
      LL = LL + 2
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0320713D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ININE
      QQ(LL+2) = 0.1242118D0*R2PPPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J + ITHRTN
      QQ(LL+3) = -0.0320713D0*R2PPPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + IEIGHT
      QQ(LL+4) = -0.1242118D0*R2PPPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IONE
      LABL(LL+4) = J + ITWELV
      QQ(LL+5) = 0.1014185D0*R2PPPF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITWO
      LABL(LL+5) = J + ILEVEN
      QQ(LL+6) = -7.7459667D0*SCF2PF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IEIGHT
      LABL(LL+6) = J + ILEVEN
      QQ(LL+7) = SIX*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ININE
      LABL(LL+7) = J + ININE
      QQ(LL+8) = -3.8729833D0*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ININE
      LABL(LL+8) = J + ITHRTN
      QQ(LL+9) = 3.8729833D0*SCF2PF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITEN
      LABL(LL+9) = J + ITWELV
      LL = LL + 9
      GO TO 290

  200 X = GAMMA(N1N2)
c     write(6,*)"i am here 6"
      QQ(LL) = X + F4P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = X - F2P0*SCF2PP(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      LL = LL + 1
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X - F2P0*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + ITWO
      QQ(LL+2) = 3.4641016D0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITHREE
      QQ(LL+3) = X + F2P0*SCF2PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = X + F2P0*SCF2PD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = X - F4P0*SCF2PD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IFIVE
      QQ(LL+6) = X + F2P0*SCF2PD(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + ISIX
      LL = LL + 6
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0962140D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ININE
      QQ(LL+2) = -0.1242118D0*R2PPPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J + ITHRTN
      QQ(LL+3) = -0.0785584D0*R2PPPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = -0.1014185D0*R2PPPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IONE
      LABL(LL+4) = J + ITEN
      QQ(LL+5) = X - F4P0*SCF2PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ISEVEN
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = 7.7459667D0*SCF2PF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISEVEN
      LABL(LL+6) = J + ITEN
      QQ(LL+7) = X - FNINE*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IEIGHT
      LABL(LL+7) = J + IEIGHT
      QQ(LL+8) = 3.8729833D0*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + ITWELV
      QQ(LL+9) = X + F3P0*SCF2PF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ININE
      LABL(LL+9) = J + ININE
      QQ(LL+10) = 3.8729833D0*SCF2PF(IP)
      LABI(LL+10) = I
      LABJ(LL+10) = J
      LABK(LL+10) = I + ININE
      LABL(LL+10) = J + ITHRTN
      QQ(LL+11) = X
      LABI(LL+11) = I
      LABJ(LL+11) = J
      LABK(LL+11) = I + ITEN
      LABL(LL+11) = J + ITEN
      QQ(LL+12) = X
      LABI(LL+12) = I
      LABJ(LL+12) = J
      LABK(LL+12) = I + ILEVEN
      LABL(LL+12) = J + ILEVEN
      QQ(LL+13) = X + FIVE*SCF2PF(IP)
      LABI(LL+13) = I
      LABJ(LL+13) = J
      LABK(LL+13) = I + ITWELV
      LABL(LL+13) = J + ITWELV
      QQ(LL+14) = X + FIVE*SCF2PF(IP)
      LABI(LL+14) = I
      LABJ(LL+14) = J
      LABK(LL+14) = I + ITHRTN
      LABL(LL+14) = J + ITHRTN
      LL = LL + 14
      GO TO 290

  210 X = GAMMA(N1N2)
c     write(6,*)"i am here 7"
      QQ(LL) = ONE*SCG1SP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = FP258*R1SPPD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = FP258*R1SPPD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = 0.2981424D0*R1SPPD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + IONE
      LL = LL + 3
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2927700D0*R1SPDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IFOUR
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = FP218*R1SPDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFIVE
      LABL(LL+2) = J + ININE
      QQ(LL+3) = FP218*R1SPDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISIX
      LABL(LL+3) = J + ITEN
      QQ(LL+4) = 0.2760262D0*R1SPDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = 0.2760262D0*R1SPDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IEIGHT
      LABL(LL+5) = J + IEIGHT
      LL = LL + 5
      GO TO 290

  220 X = GAMMA(N1N2)
c     write(6,*)"i am here 8"
      QQ(LL) = F3P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = SQRT3*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + IFOUR
      QQ(LL+2) = F3P0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + IFOUR
      QQ(LL+3) = F3P0*SCF2PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFIVE
      LABL(LL+3) = J + IFIVE
      LL = LL + 3
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0785584D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = 0.1014185D0*R2PPPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.1014185D0*R2PPPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + ITEN
      QQ(LL+4) = 0.1282854D0*R2PPPF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ITWO
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = 2.4494897D0*SCF2PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IEIGHT
      LABL(LL+5) = J + ISEVEN
      QQ(LL+6) = 4.7434165D0*SCF2PF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ININE
      LABL(LL+6) = J + ININE
      QQ(LL+7) = 4.7434165D0*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ITEN
      LABL(LL+7) = J + ITEN
      QQ(LL+8) = 6.1237244D0*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ILEVEN
      LABL(LL+8) = J + ILEVEN
      QQ(LL+9) = 6.1237244D0*SCF2PF(IP)
      LABI(LL+9) = I
      LABJ(LL+9) = J
      LABK(LL+9) = I + ITWELV
      LABL(LL+9) = J + ITWELV
      LL = LL + 9
      GO TO 290

  230 X = GAMMA(N1N2)
c     write(6,*)"i am here 9"
      QQ(LL) = F3P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = SQRT3*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = -F3P0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + IFIVE
      QQ(LL+3) = F3P0*SCF2PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + IFOUR
      LL = LL + 3
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.0785584D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = -0.1014185D0*R2PPPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I
      LABL(LL+2) = J + ININE
      QQ(LL+3) = 0.1282854D0*R2PPPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IONE
      LABL(LL+3) = J + IEIGHT
      QQ(LL+4) = 2.4494897D0*SCF2PF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ISEVEN
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = 4.7434165D0*SCF2PF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IEIGHT
      LABL(LL+5) = J + ITEN
      QQ(LL+6) = -4.7434165D0*SCF2PF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ININE
      LABL(LL+6) = J + ININE
      QQ(LL+7) = 6.1237244D0*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ITEN
      LABL(LL+7) = J + ITWELV
      QQ(LL+8) = -6.1237244D0*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ILEVEN
      LABL(LL+8) = J + ILEVEN
      LL = LL + 8
      GO TO 290

  240 X = GAMMA(N1N2)
c     write(6,*)"i am here 10"
      QQ(LL) = X + F4P0*SCF2PP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.EQ.2) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = X + F4P0*SCF2PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IONE
      QQ(LL+2) = X - F4P0*SCF2PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = X - F4P0*SCF2PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ITHREE
      QQ(LL+4) = X + F2P0*SCF2PD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + IFOUR
      QQ(LL+5) = X + F2P0*SCF2PD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IFIVE
      LL = LL + 5
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 2
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1571169D0*R2PPPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = X + EIGHT*SCF2PF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ISIX
      LABL(LL+2) = J + ISIX
      QQ(LL+3) = X + SIX*SCF2PF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ISEVEN
      LABL(LL+3) = J + ISEVEN
      QQ(LL+4) = X + SIX*SCF2PF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IEIGHT
      LABL(LL+4) = J + IEIGHT
      QQ(LL+5) = X
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ININE
      LABL(LL+5) = J + ININE
      QQ(LL+6) = X
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ITEN
      LABL(LL+6) = J + ITEN
      QQ(LL+7) = X - TEN*SCF2PF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ILEVEN
      LABL(LL+7) = J + ILEVEN
      QQ(LL+8) = X - TEN*SCF2PF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ITWELV
      LABL(LL+8) = J + ITWELV
      LL = LL + 8
      GO TO 290

  250 X = GAMMA(N1N2)
c     write(6,*)"i am here 11"
      QQ(LL) = -0.0894427D0*R2SDPP(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I + IONE
      LABL(LL) = J - ITHREE
      QQ(LL+1) = -0.0894427D0*R2SDPP(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J - ITWO
      QQ(LL+2) = 0.1788854D0*R2SDPP(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J - IONE
      QQ(LL+3) = ONE*SCG2SD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I
      LABL(LL+3) = J
      QQ(LL+4) = 0.1277753D0*R2SDDD(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J
      QQ(LL+5) = -0.1277753D0*R2SDDD(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + IONE
      QQ(LL+6) = -0.1277753D0*R2SDDD(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + ITWO
      QQ(LL+7) = 0.0638877D0*R2SDDD(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISEVEN
      LABL(LL+7) = J + ITHREE
      QQ(LL+8) = 0.0638877D0*R2SDDD(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IEIGHT
      LABL(LL+8) = J + IFOUR
      LL = LL + 8
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.1434274D0*R2SDPF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = 0.1434274D0*R2SDPF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + ISEVEN
      QQ(LL+3) = 0.1756620D0*R2SDPF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + IFIVE
      QQ(LL+4) = 0.1192570D0*R2SDFF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + ININE
      LABL(LL+4) = J + IFIVE
      QQ(LL+5) = 0.0894427D0*R2SDFF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + ITEN
      LABL(LL+5) = J + ISIX
      QQ(LL+6) = 0.0894427D0*R2SDFF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ILEVEN
      LABL(LL+6) = J + ISEVEN
      QQ(LL+7) = -0.1490712D0*R2SDFF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + IFRTEN
      LABL(LL+7) = J + ITEN
      QQ(LL+8) = -0.1490712D0*R2SDFF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + IFVTEN
      LABL(LL+8) = J + ILEVEN
      LL = LL + 8
      GO TO 290

  260 X = GAMMA(N1N2)
c     write(6,*)"i am here 12"
      QQ(LL) = ONE*SCG1PD(IP) + 18.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = -SQRT3*SCG1PD(IP) - 5.1961524D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + IONE
      QQ(LL+2) = -SQRT3*SCG1PD(IP) - 5.1961524D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + ITWO
      QQ(LL+3) = -SQRT3*SCG1PD(IP) + 20.7846097D0*SCG3PD(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITWO
      LABL(LL+3) = J + ITHREE
      LL = LL + 3
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.1069045D0*R1PDDF(IP) + 0.0458162D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITHREE
      LABL(LL+1) = J + ISIX
      QQ(LL+2) = 0.0308607D0*R1PDDF(IP) + 0.0529040D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IFOUR
      LABL(LL+2) = J + ISIX
      QQ(LL+3) = -0.1195229D0*R1PDDF(IP) - 0.0341494D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + IFOUR
      LABL(LL+3) = J + ITEN
      QQ(LL+4) = 0.0308607D0*R1PDDF(IP) + 0.0529040D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFIVE
      LABL(LL+4) = J + ISEVEN
      QQ(LL+5) = -0.1195229D0*R1PDDF(IP) - 0.0341494D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ILEVEN
      QQ(LL+6) = 0.0755929D0*R1PDDF(IP) + 0.0215980D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + ISIX
      LABL(LL+6) = J + IFIVE
      QQ(LL+7) = -0.0975900D0*R1PDDF(IP) + 0.0418243D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISIX
      LABL(LL+7) = J + IEIGHT
      QQ(LL+8) = -0.0975900D0*R1PDDF(IP) + 0.0418243D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISEVEN
      LABL(LL+8) = J + ININE
      LL = LL + 8
      GO TO 290

  270 X = GAMMA(N1N2)
c     write(6,*)"i am here 13"
      QQ(LL) = ONE*SCG1PD(IP) + 18.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      QQ(LL+1) = SQRT3*SCG1PD(IP) + 5.1961524D0*SCG3PD(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I
      LABL(LL+1) = J + IONE
      QQ(LL+2) = -SQRT3*SCG1PD(IP) + 20.7846097D0*SCG3PD(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + IONE
      LABL(LL+2) = J + IFOUR
      LL = LL + 2
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = -0.1069045D0*R1PDDF(IP) + 0.0458162D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + ITWO
      LABL(LL+1) = J + ISEVEN
      QQ(LL+2) = -0.0308607D0*R1PDDF(IP) - 0.0529040D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITHREE
      LABL(LL+2) = J + ISEVEN
      QQ(LL+3) = -0.1195229D0*R1PDDF(IP) - 0.0341494D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ILEVEN
      QQ(LL+4) = 0.0308607D0*R1PDDF(IP) + 0.0529040D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 0.1195229D0*R1PDDF(IP) + 0.0341494D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFOUR
      LABL(LL+5) = J + ITEN
      QQ(LL+6) = -0.0975900D0*R1PDDF(IP) + 0.0418243D0*R3PDDF(IP)
      LABI(LL+6) = I
      LABJ(LL+6) = J
      LABK(LL+6) = I + IFIVE
      LABL(LL+6) = J + ININE
      QQ(LL+7) = 0.0755929D0*R1PDDF(IP) + 0.0215980D0*R3PDDF(IP)
      LABI(LL+7) = I
      LABJ(LL+7) = J
      LABK(LL+7) = I + ISIX
      LABL(LL+7) = J + IFIVE
      QQ(LL+8) = 0.0975900D0*R1PDDF(IP) - 0.0418243D0*R3PDDF(IP)
      LABI(LL+8) = I
      LABJ(LL+8) = J
      LABK(LL+8) = I + ISIX
      LABL(LL+8) = J + IEIGHT
      LL = LL + 8
      GO TO 290

  280 X = GAMMA(N1N2)
c     write(6,*)"i am here 14"
      QQ(LL) = F4P0*SCG1PD(IP) + 27.00D0*SCG3PD(IP)
      LABI(LL) = I
      LABJ(LL) = J
      LABK(LL) = I
      LABL(LL) = J
      IF (NT.LT.5) GO TO 290
      N3 = N2 + 1
      N1N2 = NIN(N3) + N1
      X = GAMMA(N1N2)
      QQ(LL+1) = 0.2618615D0*R1PDDF(IP) + 0.0748176D0*R3PDDF(IP)
      LABI(LL+1) = I
      LABJ(LL+1) = J
      LABK(LL+1) = I + IONE
      LABL(LL+1) = J + IFIVE
      QQ(LL+2) = 0.1951800D0*R1PDDF(IP) - 0.0836486D0*R3PDDF(IP)
      LABI(LL+2) = I
      LABJ(LL+2) = J
      LABK(LL+2) = I + ITWO
      LABL(LL+2) = J + IEIGHT
      QQ(LL+3) = 0.1951800D0*R1PDDF(IP) - 0.0836486D0*R3PDDF(IP)
      LABI(LL+3) = I
      LABJ(LL+3) = J
      LABK(LL+3) = I + ITHREE
      LABL(LL+3) = J + ININE
      QQ(LL+4) = 0.2468854D0*R1PDDF(IP) + 0.0264520D0*R3PDDF(IP)
      LABI(LL+4) = I
      LABJ(LL+4) = J
      LABK(LL+4) = I + IFOUR
      LABL(LL+4) = J + ISIX
      QQ(LL+5) = 0.2468854D0*R1PDDF(IP) + 0.0264520D0*R3PDDF(IP)
      LABI(LL+5) = I
      LABJ(LL+5) = J
      LABK(LL+5) = I + IFIVE
      LABL(LL+5) = J + ISEVEN
      LL = LL + 5
  290 CONTINUE
c     write(6,*)"i am i gforma"
c     write(6,*)"LABI",(LABI(I),I=1,30)
c     write(6,*)"LABJ",(LABJ(I),I=1,30)
c     write(6,*)"LABK",(LABK(I),I=1,30)
c     write(6,*)"LABL",(LABL(I),I=1,30)
      RETURN

      END
