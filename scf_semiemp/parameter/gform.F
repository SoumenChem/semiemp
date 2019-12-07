       SUBROUTINE GFORM(NB,NG,GAMMA,H,NU,NW,KAN,NTYP,ISW1,SCG1SP,SCF2PP,
     1 SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,SCG3SF,SCG2PF,SCG4PF,
     2 SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,SCF2FF,SCF4FF,SCF6FF,
     3 R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,R2PFFF,R2PPPF,R2SDFF,
     4 R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF,QQ,LABI,LABJ,LABK,
     5 LABL,LL,NIN) 
C
C     THIS SUBROUTINE ASSEMBLES ALL INTEGRALS  AND STORES THEM AND THEIR
C     CODE 
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
      PARAMETER (MXCI = 6000,
     + ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0, MXA=1200,MXB=4000,
     +  MPE=105,ISPEC = 0)
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
      dimension  QQ(*),LABI(*),LABJ(*)
      dimension   LABK(*),LABL(*), NIN(*)
      DIMENSION NTYP(*),NU(*),NW(*),NG(*),KAN(*)
      DIMENSION  H(*),GAMMA(*)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
c     INCLUDE '../include/allgf.cmn'
c     INCLUDE '../include/nind.cmn'
c     INCLUDE '../include/switch.cmn'
c     INCLUDE '../include/z.cmn'
C
C
c     L=0
c     if (ga_nodeid.eq.0)
c    .write(6,*) "i am in gform"
      DO  I = 1, NB*NB
            QQ(I) = ZERO
            LABI(I) = 0
            LABJ(I) = 0
            LABK(I) = 0
            LABL(I) = 0
      END DO

      RT3 = SQRT(THREE)
      LL = 0
      I = 0
      J = 0
      IF (ISW1.GT.2 .AND. ISW1.NE.7) GO TO 30
C
C     **** CNDO and PPP integrals only ****
C
      DO 20 I = 1, NB
         N1 = NG(I)
c        write(6,*)"N1",N1
C
         DO 10 J = I, NB
            N2 = NG(J)
c      write(6,*)"N2",N2
c           IF (NEXTR(I).LT.0 .OR. NEXTR(J).LT.0) THEN
C
C     IF(I.NE.J.AND.NU(I).EQ.NU(J)) GO TO 760
C
c              IF (NU(I).EQ.NU(J)) GO TO 10
            LL = LL + 1
            N1N2 = NIN(N2) + N1
c           write(6,*)"N1N2",N1N2
c           write(6,*)"GAMMA",GAMMA(N1N2)
            QQ(LL) = GAMMA(N1N2)
            LABI(LL) = I
            LABJ(LL) = I
            LABK(LL) = J
            LABL(LL) = J
c           write(6,*)"LL,I,J",LL,I,J
c           IF (LL.LT.10000) GO TO 10
c           NINTGL = 10000
C
C      WRITE (11) (QQ(K),K=1,10000),
C    1           (LABI(K),LABJ(K),LABK(K),LABL(K),K=1,10000),NINTGL
C           WRITE(6,*) ' GFORM ',(QQ(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABI(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABJ(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABK(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABL(K),K=1,100)
c           CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, NINTGL, 1, 1, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, NINTGL, 1, 1, 'GFORM')
C
c           IF(IPRINT.EQ.1) WRITE (6,*)' NO. OF AO. INTEGRALS', NINTGL
C           WRITE (6,*) LABI, LABJ, LABK, LABL
C
c           LL = 0
   10    CONTINUE
   20 CONTINUE
c     write(6,*)"LABI",LABI
c     write(6,*)"LABJ",LABJ
c     write(6,*)"LABK",LABK
c     write(6,*)"LABL",LABL
      GO TO 150

   30 CONTINUE
C
C     +++++++++++++++  INDO INTEGRALS  +++++++++++++++++++++++++++
C
c     write(6,*)"NIN",(NIN(I),I=1,NB)
      DO 140 I = 1, NB
         IP = NU(I)
         I1 = NW(I) + 1
         N1 = NG(I)
         NZ = KAN(I)
         NT = NTYP(IP)
c       write(6,*)"before J cycle"
C
         DO 130 J = I, NB
            N2 = NG(J)
c     write(6,*)"NIN",(NIN(IJ),IJ=1,NB)
c         write(6,*)"NIN(N2),N2",NIN(9),NIN(N2),N2
            N1N2 = NIN(N2) + N1
            LL = LL + 1
c         write(6,*)"before 50",NU(J),IP
            IF (NU(J).EQ.IP) GO TO 50
C
C           TWO CENTER INTEGRALS
c           write(6,*)"here",N1N2!,GAMMA(N1N2)
c           write(6,*)"QQ",QQ(LL)
   40       QQ(LL) = GAMMA(N1N2)
c           write(6,*)"here 1"
            LABI(LL) = I
            LABJ(LL) = I
            LABK(LL) = J
            LABL(LL) = J
c           write(6,*)"LL,I,J",LL,I,J
            GO TO 110
C
C     One center integrals
C
   50       CONTINUE
c           IF (NEXTR(I).LT.0 .OR. NEXTR(J).LT.0) THEN
C       IF(I.NE.J) THEN
c              LL = LL - 1
c              GO TO 110
C       ENDIF
C
C       N1N2     = NIN(N2) + N1
C       QQ(LL)   = GAMMA(N1N2)
C       LABI(LL) = I
C       LABJ(LL) = I
C       LABK(LL) = J
C       LABL(LL) = J
C       GO TO 185
C
            J1 = NW(J) + 1
            ICS = NIN(J1) + I1
c        write(6,*)"ICS",ICS
C
C   Beginning of a new integral code
C
   60       IF (ICS.LT.15) THEN
         CALL GFORMA(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
     1  I,J,LL,N1,N2,NT,N1N2
     1 ,SCG1SP,SCF2PP,SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,
     2 SCG3SF,SCG2PF,SCG4PF,SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,
     3 SCF2FF,SCF4FF,SCF6FF,R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,
     4 R2PFFF,R2PPPF,R2SDFF,R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF)
 9000          FORMAT (5I5,F10.6)
c       write(6,*)"LABL",LABL(LL)
               GO TO 110

            END IF
   70       IF (ICS.LT.37) THEN
            CALL GFORMB(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
     1  I,J,LL,N1,N2,NT,N1N2
     1 ,SCG1SP,SCF2PP,SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,
     2 SCG3SF,SCG2PF,SCG4PF,SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,
     3 SCF2FF,SCF4FF,SCF6FF,R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,
     4 R2PFFF,R2PPPF,R2SDFF,R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF)

               GO TO 110

            END IF
   80       IF (ICS.LT.63) THEN
           CALL GFORMC(GAMMA,ICS,NIN,QQ,LABI,LABJ,LABK,LABL,IP,
     1  I,J,LL,N1,N2,NT,N1N2
     1 ,SCG1SP,SCF2PP,SCG2SD,SCG1PD,SCF2PD,SCG3PD,SCF2DD,SCF4DD,
     2 SCG3SF,SCG2PF,SCG4PF,SCG1DF,SCG3DF,SCG5DF,SCF2PF,SCF2DF,SCF4DF,
     3 SCF2FF,SCF4FF,SCF6FF,R1SPPD,R2SDDD,R2SDPP,R1PDDF,R1SPDF,R2PFDD,
     4 R2PFFF,R2PPPF,R2SDFF,R2SDPF,R3PDDF,R3SFDF,R3SFPD,R4PFDD,R4PFFF)
c      write(6,*)"HERE"
               GO TO 110

            END IF
   90       IF (ICS.LT.95) THEN
            write(6,*)"gforme"
              CALL GFORME(GAMMA,ICS)
               GO TO 110

            END IF
  100       IF (ICS.LT.137) THEN
            write(6,*)"gforme"
               CALL GFORME(GAMMA,ICS)
               GO TO 110

            END IF
C
C   End of the integral code
C
  110       CONTINUE
c           write(6,*)"LL here",LL
c           IF (LL.LT.10000) GO TO 130
C
C     ---------  Storage of INDO integrals in unit 11 ----------
C
C     NINTGL = 10000 is the currently used number. It is already
C     too small, and should be changed to  100,000 or similar.
C     But this change  must be done over  the entire program !
C
c           NINTGL = 10000
C
C     WRITE (11) (QQ(K), K = 1, 10000 ),
C    1           (LABI(K),LABJ(K),LABK(K),LABL(K),K=1,10000),NINTGL
C
c           CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
C           WRITE(6,*) ' GFORM ',(QQ(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABI(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABJ(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABK(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABL(K),K=1,100)
c           CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c           CALL ZIWRITE (IO_AOINTS, NINTGL, 1, 1, 'GFORM')
C
C     WRITE (6,*)' GFORM '
c           IF(IPRINT.EQ.1) WRITE (6,*)' NO. OF AO. INTEGRALS', NINTGL
C     WRITE (6,*) LABI, LABJ, LABK, LABL
c           LL = LL - 10000
c           IF (LL.EQ.0) GO TO 130
C
c           DO 120 K = 1, LL
c              QQ(K) = QQ(K+10000)
c              LABI(K) = LABI(K+10000)
c              LABJ(K) = LABJ(K+10000)
c              LABK(K) = LABK(K+10000)
c              LABL(K) = LABL(K+10000)
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
C
  150 CONTINUE

C
C     +++++++++++++++  NDDO  INTEGRALS   ++++++++++++++++++++
C
C     For ISW1 = 5 (NDDO/1) or ISW1 = 6 (NDDO/2),
C     call NDDO two electron routine
C
c     IF(ISW1.EQ.5.OR.ISW1.EQ.6) !CALL NDDOIN( H, GAMMA, LL)
C
c     IF (LL.EQ.0) GO TO 160
C      WRITE (11) (QQ(K),K=1,10000),
C    1           (LABI(K),LABJ(K),LABK(K),LABL(K),K=1,10000),LL
c     CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LL, 1, 1, 'GFORM')
C
C     WRITE (6,*)' GFORM '
c     IF(IPRINT.EQ.1) WRITE (6,*)' NO. OF AO. INTEGRALS', LL
C     WRITE (6,*) LABI, LABJ, LABK, LABL
C
c 160 CONTINUE
c     LL = -1
C
C      WRITE (11) (QQ(K),K=1,10000),
C    1           (LABI(K),LABJ(K),LABK(K),LABL(K),K=1,10000),LL
C
c     CALL ZRWRITE (IO_AOINTS, QQ, 1, 10000, 'GFORM')
c           WRITE(6,*) ' GFORM ',(QQ(K),K=1,100)
c           WRITE(6,*) ' GFORM ',(LABI(K),K=1,100)
c           WRITE(6,*) ' GFORM ',(LABJ(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABK(K),K=1,100)
C           WRITE(6,*) ' GFORM ',(LABL(K),K=1,100)
c     CALL ZIWRITE (IO_AOINTS, LABI, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABJ, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABK, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LABL, 1, 10000, 'GFORM')
c     CALL ZIWRITE (IO_AOINTS, LL, 1, 1, 'GFORM')
C     WRITE (6,*)' GFORM '
c     IF(IPRINT.EQ.1) WRITE (6,*)' NO. OF AO. INTEGRALS', LL
C     WRITE (6,*) LABI, LABJ, LABK, LABL
c      WRITE (6,*)' NO. OF AO. INTEGRALS', LL
C
c     write(6,*)"LABI",(LABI(I),I=1,20)
c     write(6,*)"LABJ",(LABJ(I),I=1,20)
c     write(6,*)"LABK",(LABK(I),I=1,20)
c     write(6,*)"LABL",(LABL(I),I=1,20)
  170 RETURN

      END




