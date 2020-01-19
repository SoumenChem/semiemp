      SUBROUTINE SCFIT(rtdb,F,H,Q,SH,FT,RAD,GAMMA,ISW1,NB,NA,ISWANN,
     1 scftol,itmax,MIM,ISW3,EIG,FOP,NEL,MULT,NRR,BB,NDN,NTYP,NU
     2,NG,JSLIT,SUMELEL,QQ,LABI,LABJ,LABK,LABL,LL,REP,FZ,NIN,E1,D12,AA,
     3 RA,RB,NCORE,FA)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT INTEGER (I-N)
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
#include "inp.fh"
      integer rtdb
      double precision mult,FOP(*)
c     Integer rtdb
C     WARNING   RPA PGM MUST HAVE THESE DIMENSIONS
      PARAMETER (MXB = 40000000, MXA = 1200,MXF = 100, MPE = 105)
      PARAMETER (MXCI = 6000)
      PARAMETER (MATCL = 8000, MAXTYP = 50)
      PARAMETER (IADIM=272402)
      PARAMETER
     + (ZERO =  0.0D0, ONE  =  1.0D0, TWO  =  2.0D0, THREE = 3.0D0,
     +  FOUR =  4.0D0, FIVE =  5.0D0, SIX  =  6.0D0, SEVEN = 7.0D0,
     +  EIGHT = 8.0D0,FNINE =  9.0D0, TEN  = 10.0D0, HALF =  0.5D0,
     +  THRD =  ONE  / THREE, QURT = 0.25D0)
      PARAMETER
     +   (IZERO = 0, IONE = 1, ITWO = 2, ITHREE = 3, IFOUR = 4,
     +   IFIVE = 5, ISIX = 6, ISEVEN = 7, IEIGHT = 8, ININE = 9
     +   ,ITEN = 10)
      dimension  F(*),H(*),FA(*),GAMMA(*),SH(*),RAD(*),NU(*),NG(*)
     .    ,D13(NA),NASS(7),MIM(31),EIG(*),BB(*),FT(*),!SC(NB,33)
     . JSLIT(4),SUMELEL(*),QQ(*),LABI(*),LABJ(*),LABK(*),!W(2*NB)
     . LABL(*), D12(NA),FZ(*),NIN(*),!RA(NRR),RB(NRR)
     . AA(*),RA(*),RB(*),Q(*)
      CHARACTER*8  BLOCK
      CHARACTER*6  ASSINP
      CHARACTER*20 INTPOL
C     CHARACTER*20 ONAME
      CHARACTER*80 LINE,BLANK2
      logical slowscf

C     WARNING: SC EXISTS ONLY IN ITERATIVE EXTENDED HUCKEL (ISW1=0) !

c     LOGICAL      NEWROT_CALLED

c     DATA         ASSINP/'ASSINP'/
       E1 = ZERO
       IQMMM = IZERO
       IPRINT = IZERO
C     write(6,*)"i am in scfit"
C
       slowscf = .false.
       if (.not.rtdb_cget(rtdb,'semiemp:intpol',1,intpol)) then
       call errquit('scfit: failed to read intpol',0, RTDB_ERR)
       endif
       if (inp_compare(.false.,intpol,'slow'))  then
        slowscf = .true.
       end if
       
C     H IS EXTENDED HUCKEL STARTING MATRIX
             if (.not.MA_Push_Get(MT_Dbl,2*nb,'W',lW,iW))
     &   call errquit('scfit: cannot allocate W',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nb*nb,'PRA',lPRA,iPRA))
     &   call errquit('scfit: cannot allocate PRA',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nb*nb,'PRB',lPRB,iPRB))
     &   call errquit('scfit: cannot allocate PRB',0,
     1    MA_ERR)
c            if (.not.MA_Push_Get(MT_Dbl,nrr,'Q',lQ,iQ))
c    &   call errquit('scfit: cannot allocate RB',0,
c    1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,2*NA,'PNOW',lPNOW,iPNOW))
     &   call errquit('scfit: cannot allocate PNOW',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,2*NA,'XBACK1',lXBACK1,iXBACK1
     &   )) call errquit('scfit: cannot allocate XBACK1',0, MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,2*NA,'XBACK2',lXBACK2,iXBACK2
     &   )) call errquit('scfit: cannot allocate XBACK2',0, MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,2*NA,'PBACK1',lPBACK1,iPBACK1
     &   )) call errquit('scfit: cannot allocate PBACK1',0, MA_ERR)

      if (ga_nodeid().eq.0)
     & WRITE (6,*) 'SCFIT IS STARTING SCF ITERATION'
c     write(6,*)"FA"!,(FA(I),I=1,40)
C     NEWROT_CALLED = .FALSE.
      KKPP = 0
c     IF (ISW1 .EQ. IZERO) THEN
c      DO 5 I=1,33
c       DO 4 J=1,NB
c        SC(J,I) = ZERO
c   4   CONTINUE
c   5  CONTINUE
c     END IF
      DO 10 I=1,NA
       D13(I) = ZERO
       DO J = 1, 2
        Dbl_MB(iPNOW+i*j-1)= ZERO
       END DO
   10 CONTINUE
   20 PREVS = ZERO
      IT    = IZERO
c     IF ((ISW3 .NE. IZERO) .OR. (JSW(6) .EQ. 4)) GO TO 200
c     DO 30 IJ = 1, 80
c      BLANK2(IJ:IJ) = ' '
   30 CONTINUE
   50 DO 60 MGC=1,7
       NASS(MGC) = 0
   60 CONTINUE
      GO TO 170
  170 CONTINUE
c     if (.not.rtdb_get(rtdb,'semiemp:scftol',mt_dbl,1,scftol))
c    &  call errquit('scfit failed read to scftol',0,RTDB_ERR)
c     end if
C     
C     SCF LOOP STARTS HERE (LOOPS TO LABEL 200).
C
C
C --- LASTIT Needs to be initialized further up ---
C !JT1201
C --- Optimizations dont work otherwise        ---
C !JT1201
C
      LASTIT=0
C     !JT1201
  200 ITT = IT + 1
      IF (IT.EQ.0) THEN
       DO I= 1, NRR
c        Dbl_MB(iQ+I-1) = FA(I)
c        FA(I) = Q(I)
         BB(I) = SH(I)
c      write(6,*)"BB",BB(I)
c      write(6,*)"Q(I)",Dbl_MB(iQ+I-1)
       END DO
c      DO I= 1, NB*NB
c        BB(I)=SH(I)
c      END DO
      END IF
c     write(6,*)"F",(F(I),I=1,NRR)
c     write(6,*)"H",(H(I),I=1,NRR)
c     write(6,*)"FA",(FA(I),I=1,NRR)
c     write(6,*)"BB",(BB(I),I=1,NRR)
c     write(6,*)"SH",(SH(I),I=1,NRR)
      N_scf_cycles=n_scf_cycles+1
CSG   IF(IT.GT.10000) STOP
c     IF(LASTIT.GE.1) CALL RACTF2(NA,RAD)
      IF (LASTIT.GE.1) GOTO 220
      DO I = 1, NRR
c     Do I = 1, NB*NB
       H(I) = FT(I)
      END DO
      IF (IT.GE.1) THEN
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
c     IF (IT.eq.1) THEN
       DO I = 1, NRR
         F(I) = BB(I)
         FA(I) = AA(I)
         Dbl_MB(iPRA+I-1) = RA(I)
         Dbl_MB(iPRB+I-1) = RB(I)
       END DO
c     ELSE
c       DO I = 1, NRR
c        FA(I) = BB(I)
c        F(I) = AA(I)
c        Dbl_MB(iPRA+I-1) = RB(I)
c        Dbl_MB(iPRB+I-1) = RA(I)
c       END DO
c     END IF
      ELSE
       Do I = 1, NRR
         F(I) = BB(I)
         BB(I)=SH(I)
       END DO
      END IF
      END IF
C
C   
C     ---------------
      IF (ISWANN .EQ. IZERO) THEN
C     CLOSED SHELL RESTRICTED
c      WRITE(6,*) ' SCFIT BEFORE DIAGR'
c      WRITE(6,*) " SCFIT BEFORE DIAGr, F",(F(I),I=1,100)
c      write(6,*)"FA before diagr",(FA(I),I=1,100)
c      write(6,*)"H before diagr",(H(I),I=1,100)
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
       CALL DIAGR(rtdb,BB,F,H,FA,RAD,ITMAX,IT,PREVS,SCFTOL,MIM,ISW3,
     1  NB,EIG,FOP,ISWANN,NEL,MULT,Dbl_MB(iW),JSLIT,NRR,Q,NIN,
     2  NCORE,NASS)
       call ga_dgop(Msg_diagr,F,NB*NB,'+')!Propagate status

c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "diagr cpu time", total_cpu, "total time",
c    1 total_wall
c     end if
       
c       write(6,*)"JSLIT in scfit",JSLIT
C      CALL SLITEN(0,'DIAGR')
c      write(6,*)"BB after diagr",(BB(I),I=1,nrr)
c      write(6,*)"H after diagr",(H(I),I=1,nrr)
c      write(6,*)"F after diagr",(F(I),I=1,nrr)
      ENDIF
c      WRITE(6,*) " SCFIT BEFORE DIAGU, F",(F(I),I=1,nrr)
c      write(6,*)"FA before diagu",(FA(I),I=1,nrr)
c      write(6,*)"H before diagu",(H(I),I=1,nrr)
      IF ((ISWANN .EQ. 1) .OR. (ISWANN .EQ. 2)) then
C     OPEN SHELL UNRESTRICTED
       CALL DIAGU(BB,F,H,FA,FZ,RAD,ITMAX,IT,PREVS,SCFTOL,MIM,ISW3,NB,
     1 EIG,FOP,ISWANN,NEL,MULT,Dbl_MB(iW),JSLIT,NRR,Q,NIN,NU
     1 ,NTYP,NDN,RA,RB,NCORE,NASS)
c      write(6,*)"FA after diagu",(FA(I),I=1,nrr)
c      call ga_dgop(Msg_diagu,F,NB*NB,'+')!Propagate status
c      IF (IT.EQ.0) then
c      call ga_dgop(Msg_diagu,FA,NB*NB,'+')
c      END IF
      END IF
c     IF (ISWANN.EQ.3)
C     OPEN SHELL RESTRICTED
c    .   CALL DIAGO(BB,F,H,FA,FZ,RAD,ITMAX,IT,PREVS,SCFTOL,IXS,DRYRUN)
C        WRITE(6,*) ' SCFIT AFTER DIAGO'
C     ---------------
  220 CONTINUE
      IT = IT + 1
C
C     TO CALCULATE THE GASPHASE ENERGY OF THE SOLVATED DENSITY
C     THE SOLVATION PART CONTRIBUTED TO CORE HAMILTONIAN AND FOCK
C     MATRIX SHOULD BE DELETED
C
C
C     CALL SLITEN(0,'SING1')
c     write(6,*)"JSLIT before sing",JSLIT
c      write(6,*)"FA after diagr",(FA(I),I=1,nrr) 
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec() 
       CALL SING(BB,AA,F,H,FA,FZ,SVEC,GAMMA,SC,RAD,E1,BH,
     .  IT,ITGL,ITMAX,ISW1,ISWANN,NDN,NB,NA,ISW3,NTYP,NU,Dbl_MB(iW),
     .  JSLIT,FOP,SUMELEL,NG,QQ,NRR,LABI,LABJ,LABK,LABL,LL,Q,
     .  REP,D12,NIN,Dbl_MB(iPNOW),RA,RB,Dbl_MB(iPRA), Dbl_MB(iPRB),
     . Dbl_MB(iXBACK1),Dbl_MB(iXBACK2),Dbl_MB(iPBACK1),slowscf)
C     CALL SLITEN(0,'SING2')
c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "sing cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

      CALL SCONS(D12,D13,IT,ITMAX,NA,SCFTOL,E1,JSLIT)

c    CALL SLITEN(0,'SING ')
C     
  230 CONTINUE
      DO 240 I = 1, NA
       D13(I) = ZERO
  240 CONTINUE
      J = 0
      DO 250 I=1,NB
       K      = NU(I)
       J      = J + I
       D13(K) = D13(K) + F(J)
  250 CONTINUE
      IF (IT.EQ.1.AND.ISWANN.EQ.1.AND.ga_nodeid().eq.0) 
     .WRITE(6,9200)
      IF (IT .NE. 1) GO TO 260
      IF (IQMMM .EQ. IZERO.AND.IPRINT.NE.-2 .and.ga_nodeid().eq.0) 
     .WRITE (6,9010) IT, E1
      GO TO 270

  260 DE1 = (PREVS-E1)
      IF (IQMMM .EQ. IZERO.AND.IPRINT.NE.-2 .and.ga_nodeid().eq.0) 
     .WRITE (6,9010) IT,E1,DE1
  270 PREVS = E1
      CALL SLITET(1,I11,JSLIT)
C     CALL SLITEN(0,'TEST1')

C     END SCF IF SENSE LIGHT ON ***************************
C     ELSE LOOP BACK TO LABEL 200 FOR ANOTHER SCF CYCLE
C
C AVM2010OCT 
CLI
C     I11 = 1 IS ON
c     write(6,*)"LASTIT",LASTIT
      IF (LASTIT.GE.1) GOTO 275
c     write(6,*)"I11",I11
      IF (I11 .NE. 1) GO TO 200
CJiabo
C     DO an addtional iteration for construct the reaction field for
C     possible VEM type calculation
 275  LASTIT=LASTIT+1
      IF(LASTIT.LE.1) GOTO 200
CLi
      IF (ISWANN .EQ. IZERO) THEN
c     if (ga_nodeid().eq.0) then
c     do i = 1, NB
c     write(6,*)"EIG = ", EIG(I)
c      end do
        DO  280 I=1,NRR
c       DO 280 I =1,NB*NB
         F(I) = BB(I)
c        write(6,*)"F(I)",F(I)
  280    CONTINUE
c       end if 
      ELSE
        DO 290 I=1,NRR
        FA(I)  = RB(I)+ RA(I)
  290   CONTINUE
      END IF
      if (.not.ma_pop_stack(lPBACK1))
     &   call errquit('scfit: cannot pop stack PNOW',3, MA_ERR)
      if (.not.ma_pop_stack(lXBACK2))
     &   call errquit('scfit: cannot pop stack XBACK2',3, MA_ERR)
      if (.not.ma_pop_stack(lXBACK1))
     &   call errquit('scfit: cannot pop stack XBACK1',3, MA_ERR)
      if (.not.ma_pop_stack(lPNOW))
     &   call errquit('scfit: cannot pop stack PNOW',3, MA_ERR)
c     if (.not.ma_pop_stack(lQ))
c    &   call errquit('scfit: cannot pop stack Q',3, MA_ERR)
      if (.not.ma_pop_stack(lPRB))
     &   call errquit('scfit: cannot pop stack lPRB',3, MA_ERR)
      if (.not.ma_pop_stack(lPRA))
     &   call errquit('scfit: cannot pop stack lPRA',3, MA_ERR)
      if (.not.ma_pop_stack(lW))
     &   call errquit('scfit: cannot pop stack lW',3, MA_ERR)
c     IF (USEUNIT0) WRITE (0,*) 'SCFIT out, SCF iterations done'
 5000 RETURN
 9010 FORMAT (7X,'CYCLE ',I4,6X,'E(SCF) =',F16.8,' au',4X,'DIFF.     =',
     .        1PD13.5)
 9020 FORMAT (8F10.6)
 9030 FORMAT (7X,'CYCLE=',I4,34X,'DENS. TOL.=',F10.6)
 9040 FORMAT (16I5)
 9050 FORMAT (8F10.6)
 9060 FORMAT (16F5.2)
 9070 FORMAT (F10.6,3 (I5,F10.6))
 9200 FORMAT (/,' UNRESTRICTED HARTREE-FOCK CALC.- ',/
     .       )
      END


