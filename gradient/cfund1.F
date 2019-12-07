      SUBROUTINE CFUND1(IA,IB,IC,ID,IE,RHOA,RHOB,SNAG,DSNAG,RR,NIN,FACT
     1 , A,B)
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER(I-N)
C
C     THIS SUBROUTINE CALCULATES THE C-FUNCTIONS OF ROOTHAAN FIR OVERLAP
C     INTEGRALS. ONLY C FNS WITH POSITIVE INDICES ARE NEEDED AND THEY
C     WILL BE CALC. BY THE BINOMIAL THEOREM.
C     ORIGINAL BY DICK STEVENS, ADJUSTED BY M.ZERNER.
C
C    UPDATED BY JDHD   OCT 26 84
C
c     INCLUDE '../include/z.cmn'
C****************************************************************
C***      COMMON/Z/
C***     . ATOM(100),NIN(MXB),
C***     . NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     X NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     X Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
c     INCLUDE '../include/ev.cmn'
C****************************************************************
C***      COMMON/EV/
C***     . A(35),B(35),DD(35)
C****************************************************************
c     INCLUDE '../include/binc.cmn'
C****************************************************************
C***      COMMON/BINC/BINCOE(465)
C****************************************************************
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)

       dimension BINCOE(465)
       dimension A(35),B(35)
       dimension NIN(*)
       dimension FACT(35)

c     write(6,*)"i am in cfund1",IA,IB,IC,ID,IE
C
C     FORM BINOMIAL COEFFECIENT
C
C     ( N )
C     ( M )  INDEX = N*(N+1)/2 + M+1
C
      BINCOE(1) = ONE
      IND = 2
      DO 70 I = 1, 29
       I1 = I + 1
       DO 60 J = 1, I1
        JJ = I1 - J + 1
        IF(J.LT.26) GO TO 40
        BIN = 1.D0
        IF(J.GE.JJ) THEN
         DO 20 K = J, I
          BIN = BIN*DFLOAT(K)
   20    CONTINUE
         BINCOE(IND) = BIN/FACT(JJ)
c       write(6,*)"BINCOE",BINCOE(IND)
        ELSE
         DO 30 K = JJ, I
          BIN = BIN*DFLOAT(K)
   30    CONTINUE
         BINCOE(IND) = BIN/FACT(J)
c       write(6,*)"BINCOE",BINCOE(IND)
        END IF
        GO TO 50

   40   BINCOE(IND) = FACT(I1)/FACT(J)/FACT(JJ)
   50   IND = IND + 1
   60  CONTINUE
   70 CONTINUE
      COUNT = 0.0D0
      DCOUNT = 0.D0
c     write(6,*)"HERE"
      TAU = (RHOA-RHOB)/ (RHOA+RHOB)
      IAB = IA + 1
      IBB = IB + 1
      ICB = IC + 1
      IDB = ID + 1
      IEB = IE + 1
c     write(6,*)"IEB,IE",IEB,IE
c     write(6,*)"NIN",(NIN(I),I=1,10)
      INDE = NIN(IEB)
      INDD = NIN(IDB)
      INDC = NIN(ICB)
      INDB = NIN(IBB)
      INDA = NIN(IAB)
      DO 600 I6 = 1, IEB
c     write(6,*)"HERE 1",(INDE+I6)
         B6 = BINCOE(INDE+I6)
c      write(6,*)"B6",B6
         DO 500 I5 = 1, IEB
            B5 = BINCOE(INDE+I5)
c      write(6,*)"B5",B5
            DO 400 I4 = 1, IDB
               B4 = BINCOE(INDD+I4)
c      write(6,*)"B4",B4
               DO 300 I3 = 1, ICB
                  B3 = BINCOE(INDC+I3)
c      write(6,*)"B3",B3
                  DO 200 I2 = 1, IBB
                     B2 = BINCOE(INDB+I2)
c      write(6,*)"B2",B2
                     DO 100 I1 = 1, IAB
                        B1 = BINCOE(INDA+I1)
                        TERM = B1*B2*B3*B4*B5*B6*
     .                         (-1)** (I2+I5+I6+I4+IE+ID)
                        IR = I1 + I2 - I3 - I4 + IE + IE - I6 - I6 +
     .                       IC + ID + 3
                        IP = IA - I1 + IB - I2 + IE + IE - I5 - I5 +
     .                       IC - I3 + ID - I4 + 7
c       write(6,*)"A,B,TERM",A(IP),B(IR),TERM
                        COUNT = COUNT + A(IP)*B(IR)*TERM
                        DCOUNT = DCOUNT - (A(IP+1)*B(IR)+
     .                           TAU*A(IP)*B(IR+1))*TERM
c        write(6,*)"COUNT,DCOUNT",COUNT,DCOUNT
  100                CONTINUE
  200             CONTINUE
  300          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
      N = IA + IB + IC + ID + IE + IE + 1
C      RHOBN=(RHOB/2.0D0)**N
C      RHOBN=1.D0
C      SNAG=COUNT*RHOBN
C      DSNAG=DFLOAT(N)*SNAG/RR+(RHOA+RHOB)/RR/2.D0*RHOBN*DCOUNT
c     write(6,*)"COUNT",COUNT
      SNAG = COUNT
      DSNAG = DFLOAT(N)*SNAG/RR + 0.5D0* (RHOA+RHOB)/RR*DCOUNT
  700  RETURN

      END
