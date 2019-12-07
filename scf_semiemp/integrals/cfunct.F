      SUBROUTINE CFUNCT(IA,IB,IC,ID,IE,RHOA,RHOB,SNAG,FACT,A,B,NIN)
C
C     THIS SUBROUTINE CALCULATES THE C-FUNCTIONS OF ROOTHAAN FOR OVERLAP
C     INTEGRALS. ONLY C FNS WITH POSITIVE INDICES ARE NEEDED AND THEY
C     WILL BE CALC. BY THE BINOMIAL THEOREM.
C     ORIGINAL BY DICK STEVENS, ADJUSTED BY M.ZERNER.
C
C    UPDATED BY JDHD   OCT 26 84
C
c     INCLUDE '../include/param.cmn'
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "global.fh"
      parameter (MXCI = 6000, MXB = 4000)
      parameter (ONE=1.0d0, ZERO=0.0d0)
C
c     INCLUDE '../include/binc.cmn'
C****************************************************************
     
       dimension    BINCOE(465)
C****************************************************************
c     INCLUDE '../include/ev.cmn'
C****************************************************************
c         COMMON/EV/
        dimension    A(35),B(35)
C****************************************************************
c     INCLUDE '../include/z.cmn'
        dimension NIN(*)
        dimension FACT(35)
C****************************************************************
C***      COMMON/Z/
C***     . ATOM(100),NIN(MXB),
C***     . NA,NB,NR,NRR,NCONF,NREE,NEL,NTYP(MXA),KAN(MXA),NU(MXB),
C***     X NW(MXB),NG(MXB),MULT,NP(MXB),CO(3,MXA),ALPHA(MXB,6),
C***     X Z12(MXA),EIG(MXB),DIPP(6),FACT(30),EIGA(MXB),D12(MXA)
C****************************************************************
c     write(6,*)" i am in bin cof"
      L = 0
      EIGSHF = ZERO
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
        ELSE
         DO 30 K = JJ, I
          BIN = BIN*DFLOAT(K)
   30    CONTINUE
         BINCOE(IND) = BIN/FACT(J)
        END IF
        GO TO 50

   40   BINCOE(IND) = FACT(I1)/FACT(J)/FACT(JJ)
   50   IND = IND + 1
   60  CONTINUE
   70 CONTINUE
c      write(6,*)"MXCI",MXCI
c     DO 80 I = 1, MXCI
c      NIN(I) = L
c      L = L + I
   80 CONTINUE
c      if (ga_nodeid().eq.0)
c    .write(6,*) "i am in cfunct"
   90 COUNT = 0.0D0
      IAB = IA + 1
      IBB = IB + 1
      ICB = IC + 1
      IDB = ID + 1
      IEB = IE + 1
c     write(6,*)"IEB",IEB,NIN(1)
c     write(6,*) NIN(IEB)
      INDE = NIN(IEB)
c     write(6,*) NIN(IEB)
      INDD = NIN(IDB)
c     write(6,*) NIN(IDB)
      INDC = NIN(ICB)
      INDB = NIN(IBB)
      INDA = NIN(IAB)
c     write(6,*)"IEB","IDB","ICB","IBB","IAB"
c     write(6,*)IEB, IDB, ICB, IAB
      DO 600 I6 = 1, IEB
         B6 = BINCOE(INDE+I6)
         DO 500 I5 = 1, IEB
            B5 = BINCOE(INDE+I5)
            DO 400 I4 = 1, IDB
               B4 = BINCOE(INDD+I4)
               DO 300 I3 = 1, ICB
                  B3 = BINCOE(INDC+I3)
                  DO 200 I2 = 1, IBB
                     B2 = BINCOE(INDB+I2)
                     DO 100 I1 = 1, IAB
                        B1 = BINCOE(INDA+I1)
                        TERM = B1*B2*B3*B4*B5*B6*
     .                         (-1)** (I2+I5+I6+I4+IE+ID)
                        IR = I1 + I2 - I3 - I4 + IE + IE - I6 - I6 +
     .                       IC + ID + 3
                        IP = IA - I1 + IB - I2 + IE + IE - I5 - I5 +
     .                       IC - I3 + ID - I4 + 7
c          write(6,*) "here",COUNT,A(IP),B(IR),TERM
                        COUNT = COUNT + A(IP)*B(IR)*TERM
  100                CONTINUE
  200             CONTINUE
  300          CONTINUE
  400       CONTINUE
  500    CONTINUE
  600 CONTINUE
C      SNAG=COUNT*(RHOB/2.0D0)**(IA+IB+IC+ID+IE+IE+1)
c     write(6,*)"COUNT",COUNT
      SNAG = COUNT
  700 RETURN

      END


