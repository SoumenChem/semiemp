c       
      logical function orbital_qmnumber(valence,ntype,nz, n_qmnum,
     1 l_qmnum,nbasis,basis_type)
c This function calculte n and l quantum numbers for different
c electrons.
      implicit none
c
c
c
      double precision valence
      character*16 AOS(32)
      character*8 orbital_sym(int(valence))
      integer ntype
      integer n_qmnum(*) !N quantum numbers
      integer l_qmnum(*) !L quantum numbers
      integer jj,i,k,j,KK
      integer l
      integer btype(32)
      integer nlqn(32)
      integer basis(32), nbasis !number of basis function
      integer basis_type(32) !s=0, s=1, px=2,py=3 ....
      integer nz  !atomic number
      data btype/2*0,1,2,3,1,2,3,4,5,6,7,8,4,5,6,7,8,9
     X ,10,11,12,13,14,15,9,10,11,12,13,14,15/
      data basis/2*1,2,3,4*4,5,6,7,8,6*9,10,11,12,13,14,15,
     X 8*16/
      data nlqn / 2*0,6*1,10*2,14*3/
       data   AOS/' S      ',' S      ',
     X    ' P(X)   ', ' P(Y)   ', ' P(Z)   ',
     X   ' P(X)   ', ' P(Y)   ', ' P(Z)   ',
     X    ' D(Z2)  ', 'D(X2-Y2)', ' D(XY)  ', ' D(XZ)  ',
     X   ' D(YZ)  ', ' D(Z2)  ', 'D(X2-Y2)', ' D(XY)  ', 
     X    ' D(XZ)  ',' D(YZ)  ', ' F(Z3)  ', ' F(XZ2) ', 
     X   ' F(YZ2) ','Z(X2-Y2)', ' F(XYZ) ', ' X3-3XY2', 
     X   ' 3YX2-Y3',' F(Z3)  ', ' F(XZ2) ', ' F(YZ2) ', 
     X   'Z(X2-Y2)',' F(XYZ) ', ' X3-3XY2', ' 3YX2-Y3'/

       orbital_qmnumber=.false.
       nbasis= basis(int(valence))
      do i=1, int(valence)  
c      write(6,*)"calculating qm" 
       basis_type(i)= btype(i)
       l_qmnum(i)= nlqn(i)
       orbital_sym(i)= AOS(i)
        do j= 1, i
          k = 1
         IF (nz.LE.2) GO TO 10
          k = 2
         IF (nz.LE.10) GO TO 10
          k = 3
         IF (nz.LE.18) GO TO 10
          k = 4
         IF (nz.LE.36) GO TO 10
          k = 5
         IF (nz.LE.54) GO TO 10
          k = 6
         IF (nz.LE.86) GO TO 10
          k = 7
C        CORRECT PRINC Q. NO.
C        FOR 4S 4P 3D (NTYP=4): nS, nP, (n-1)d
C        FOR 6S 6P 5D 4F (NTYP=5): ns, np, (n-1)d, (n-2)f
C        
  10     KK = K
c        IF (ntype.GE.4 .AND. basis_type(j).GT.3) KK = k - 1
c        IF (ntype.EQ.5 .AND. basis_type(j).GT.8) KK = k - 2
C        POLARIZATION FUNCTIONS FOR FIRST ROW ARE 3D
c        IF (nz.LE.10 .AND. basis_type(j).GT.3) KK = 3
C        POLARIZATION FUNCTIONS FOR H ARE 2S, 2P
c        IF (NZ.LE.2 .AND. basis_type(j).GT.0) KK = 2
         n_qmnum(j) = KK
        end do
       end do
       end
c---------------------------------------------------------------------------------
c This function calculates core energy of the system. 
c--------------------------------------------------------------------------------
          logical function  core_energy(nat,valence,co,corenergy)


          integer nat
          double precision co(nat*3), a(3), b(3)
          double precision valence(nat)
          logical status
          double precision r, corenergy
          logical  geom_calc_distance
          integer i,j,ij,k
c         parameter (BOHRS = 5.29177249D-1)

           core_energy= .false.
           corenergy= 0.0
           do i=1, nat-1
            do j= i+1, nat
             do k= 1, 3
              a(k)= co((i-1)*3+k)
              b(k)= co((j-1)*3+k)
             end do
             status= geom_calc_distance(a,b,r)
             corenergy= corenergy+(valence(i)*valence(j)/r)
            end do
           end do
           end
c---------------------------------------------------------------------------------
c This function calculate the rad matrix. It contains interatomic
c distance netween atoms of the molecule.
c--------------------------------------------------------------------------------
          logical function radmat(nat,co,radmatrix)
 
          double precision radmatrix(nat,nat)
          integer nat
          double precision co(nat*3), a(3), b(3)
          logical status
          double precision ab
          logical  geom_calc_distance
          integer i,j,ij,k
c         parameter (BOHRS = 5.29177249D-1)

           radmat= .false.
           do i=1, nat
            do j= 1, nat
             do k= 1, 3
              a(k)= co((i-1)*3+k)
              b(k)= co((j-1)*3+k)
             end do
             status= geom_calc_distance(a,b,ab)
             radmatrix(i,j)= ab
             ab=0.0
            end do
           end do
           end   
c-----------------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------------          
        SUBROUTINE GEOMAT(NA,F,RAD,CO)
C-----------------------------------------------------------------------
C     THIS SUBROUTINE IS CONCERNED WITH SETTING UP VARIOUS
C     EULERIAN TRANSFORMATION GEOMETRICAL MATRICES
C-----------------------------------------------------------------------
C     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER  (ZERO = 0.0D0, ONE = 1.0D0)
      PARAMETER    (THRESH=1.0D-7)
C     CO(1to3,I) ARE THE X, Y, Z COORDINATES OF ATOM I IN ANGSTROM
C     RAD(IJ) ARE THE INTERATOMIC DISTANCES IN ANGSTROM
      DOUBLE PRECISION F(*)
      DOUBLE PRECISION  RAD(NA,NA)
      DOUBLE PRECISION CO(NA*3)
      INTEGER NA,I,J,IJ
C     F is 3*(NA**2), RAD is NA-triangle
C-----------------------------------------------------------------------
c     write(6,*)"i am in geomat"
      NA1 = NA*NA
      NA2 = 2*NA1
      ITJ = 0
      DO 20 J = 1, NA
         DO 10 I = 1, NA
            ITJ = ITJ + 1
c           write(6,*) "F(ITJ)", F(ITJ)
            IF (I.NE.J .AND. RAD(I,J).LT.1.0D-4) THEN
              write(6,*)"*"
               F(ITJ) = ONE
               F(ITJ+NA1) = ONE
               IF (F(ITJ).LT.ZERO) F(ITJ+NA1) = -ONE
               F(ITJ+NA2) = ZERO
            ELSE
               IF (I.NE.J) THEN
C                 COS (THETA)
                  F(ITJ) = (CO((J-1)*3+3)-CO((I-1)*3+3))/RAD(I,J)
c                write(6,*)"F(ITJ)",F(ITJ),RAD(I,J)
C                 SIN(THETA)**2
                  G = ONE - F(ITJ)*F(ITJ)
                  IF ((THRESH-G).GE.ZERO) THEN
                     F(ITJ+NA1) = ONE
                     IF (F(ITJ).LT.ZERO) F(ITJ+NA1) = -ONE
                     F(ITJ+NA2) = ZERO
                  ELSE
                     G = ONE/SQRT(G)
C                    COS (PHI)
                   F(ITJ+NA1) = (CO((J-1)*3+1)-CO((I-1)*3+1))*G/RAD(I,J)
c                    write(6,*)"F(ITJ+NA1)", F(ITJ+NA1), RAD(I,J)
C                    SIN (PHI)
                   F(ITJ+NA2) = (CO((J-1)*3+2)-CO((I-1)*3+2))*G/RAD(I,J)
c                    write(6,*)"F(ITJ+NA2)",F(ITJ+NA2),RAD(I,J)
                  END IF
                ELSE
C                 THIS IS DIAGONAL I = J
                  F(ITJ) = ONE
                  F(ITJ+NA1) = ONE
                  F(ITJ+NA2) = ZERO
c                 write(6,*) F(ITJ), F(ITJ+NA1), F(ITJ+NA2)
c                 write(6,*)"**"
               END IF
            END IF
   10    CONTINUE
   20 CONTINUE
      RETURN

      END
c-----------------------------------------------------------------------------------
c-----------------------------------------------------------------------------------
c This function make an array with number of basis functions for each
c atom.
c----------------------------------------------------------------------------------     
     
      logical function basis_function(nat,nz,valence,NB,NP,ntype)
c This function calculte n and l quantum numbers for different
c electrons.
c     implicit none
      integer nat,k,AN
      integer nz(nat)
      integer NP(nat)
      integer ntype(nat)
      double precision valence(nat)
      integer bfn(5), NB(nat) !number of basis function
      data bfn/1,4,9,9,16/
      
      basis_function=.false.
      do i = 1, nat
         k= valence(i)
c       write(6,*)"Valence",Valence(i)
         NB(i)= bfn(ntype(i))
c       write(6,*)"NB",NB(i)
         AN=nz(i)
          N = 1
         IF (AN.LE.2) GO TO 10
          N = 2
         IF (AN.LE.10) GO TO 10
          N = 3
         IF (AN.LE.18) GO TO 10
          N = 4
         IF (AN.LE.36) GO TO 10
          N = 5
         IF (AN.LE.54) GO TO 10
          N = 6
         IF (AN.LE.86) GO TO 10
          N = 7
  10      KK=N
          NP(i)=KK 
      end do
      END
c----------------------------------------------------------------------------------------------
c
c----------------------------------------------------------------------------------------------
      SUBROUTINE ZEXIT(IE)
      INTEGER MAXERR,IE
      INTEGER ERR_NORMAL, ERR_DIM, ERR_INPUT, ERR_OUTPUT, ERR_MEMORY,
     .     ERR_ZIO
      PARAMETER (MAXERR=10)
      PARAMETER (ERR_NORMAL=1)
      PARAMETER (ERR_DIM=2)
C     3 has been abused in the past and is given no name.
      PARAMETER (ERR_INPUT=4)
      PARAMETER (ERR_OUTPUT=5)
      PARAMETER (ERR_MEMORY=6)
      PARAMETER (ERR_ZIO=7)
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/errno.cmn'

      CHARACTER*40  ERROR(MAXERR)
      DATA ERROR /
     1     'Normal exit                             ',
     2     'Dimensioning error                      ',
     3     'Old error code; no further information  ',
     4     'Input error                             ',
     5     'Output error                            ',
     6     'Memory allocation error                 ',
     7     'In-core file I/O subsystem error        ',
     8     'Option not available this version       ',
     9     'Illegal error code                      ',
     .     'Illegal error code                      '/


      IF (IE.LE.0 .OR. IE.GT.MAXERR) THEN
         WRITE (6,8990) IE
 8990    FORMAT(/' ZEXIT CALLED WITH BAD ERROR CODE:',I8)
      ELSE IF (IE.EQ.1) THEN
         WRITE (6,9000)
 9000    FORMAT (/' End of run. Normal exit.')
      ELSE
         WRITE (6,9020) IE,ERROR(IE)
 9020    FORMAT (/' End of run. Error code',i3,' which means: ',A)
      END IF

      STOP

      END
c---------------------------------------------------------------------
C      REMOVE THE ONE CENTER KINETIC ENERGY FOR THIS IS INCLUDED
C      THROUGH PARAMETERS FOR THE CORE INTEGRAL U
c--------------------------------------------------------------------
      Subroutine ONEKIE(ISW1,NB,FA,NU)
      integer NB,NU(NB),ISW1
      double precision FA(*)
      integer II,KK,MA,IK
      parameter (ZERO = 0.0d0)
      
c      write(6,*)"i am in onekie"
       IK = 0 
       DO II = 1, NB
        KA     = NU(II)
        DO  KK=1,II
         MA = NU(KK)
         IK = IK + 1
         IF(MA.EQ.KA) FA(IK) = ZERO
        end do
         end do
       IF(ISW1.EQ.8) THEN
C
       IK = 0
C      ZERO THE DIAGONAL KINETIC ENERGY FOR XINDO
       DO II = 1, NB
        KA = NU(II)
        DO  KK=1,II
         MA = NU(KK)
         IK = IK + 1
         IF(MA.EQ.KA) FA(IK) = ZERO
        end do
      end do
      end if 
       return
       end
c----------------------------------------------------------------------------------
C Create N*N symmetric matrix from a trigonal matrix
C------------------------------------------------------------------------------------
      Subroutine mat_sym(M,N,ITAG)
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
      integer N,it,lt,i,j,L,ITAG
      double precision M(N*N)

        if (.not.MA_Push_Get(MT_Dbl,n*n,'T',lT,iT))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)
 
        L =0     
      do i =1, n*n
        Dbl_MB(iT+L) = M(i)
        M(i) =0.0d0
c      write(6,*)"fock",Dbl_MB(iT+L)
        L = L+1
      end do
      if (ITAG.eq.1) then
      call sym(M,Dbl_MB(iT),n)
       else
      call antisym(M,Dbl_MB(iT),n)
      end if
c     L = 0      
c     do j = 0, n-1   
c      do i = j, n
c       M(L) = Dbl_MB(iT+j+i)
c       L = L + 1
c      end do
c     end do
c     do i =1, n*n
c      write(6,*)" after fock",M(i)
c     end do
      
      if (.not.ma_pop_stack(lT))
     &   call errquit('diag_para: cannot pop stack T',3, MA_ERR)
      return
      end
C-------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------
      Subroutine sym(M,P,N)
      integer N
      double precision M(N,N),P(*)

       L = 0
      do i = 1, n
       do j = 1, i
        L = L + 1
        M(i,j) = P(L)
       if (i.ne.j) M(j,i) = P(L)
       end do
      end do
      return
      end
C----------------------------------------------------------------------------------------
C
C----------------------------------------------------------------------------------------
      Subroutine antisym(M,P,N)
      integer N
      double precision M(N,N),P(*)

       L = 0
      do i = 1, n
       do j = 1, i
        L = L + 1
        M(i,j) = P(L)
       if (i.ne.j) M(j,i) = -P(L)
       end do
      end do

      return
      end
C--------------------------------------------------------------------------------------
C Symmetric matric lower triangular matrix
C--------------------------------------------------------------------------------------
      Subroutine mat_triangular(M,N)
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
      integer N,it,lt,i,j,L
      double precision M(N*N)

        if (.not.MA_Push_Get(MT_Dbl,n*n,'T',lT,iT))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)

        L =0
      do i =1, n*n
        Dbl_MB(iT+L) = M(i)
        M(i) =0.0d0
c      write(6,*)"fock",Dbl_MB(iT+L)
        L = L+1
      end do
      call triangular(M,Dbl_MB(iT),n)
c     L = 0
c     do j = 0, n-1
c      do i = j, n
c       M(L) = Dbl_MB(iT+j+i)
c       L = L + 1
c      end do
c     end do
c     do i =1, n*n
c      write(6,*)" after fock",M(i)
c     end do

      if (.not.ma_pop_stack(lT))
     &   call errquit('diag_para: cannot pop stack T',3, MA_ERR)
      return
      end
C-------------------------------------------------------------------------------------
C-------------------------------------------------------------------------------------
      Subroutine triangular(P,M,N)
      integer N
      double precision M(N,N),P(*)

       L = 0
      do i = 1, n
       do j = 1, i
        L = L + 1
        P(L) = M(i,j)
       end do
      end do

      return
      end
C---------------------------------------------------------------------------------------
C Rebuliding FOCK matrix after scaling
C--------------------------------------------------------------------------------------
      subroutine build_fock(F,D,E,NB)
      integer NB
      double precision F(NB,NB),D(NB),E(NB)
c     write(6,*)"i am in build_fock"    
      F=0.0d0
c     do i = 1, NB
c      do j = 1, NB   
c        F(i,j) = 0.0d0
c       end do
c     end do

      do i = 1, NB
        F(i,i) = D(i)
      end do
      do i = 2, NB
       F(i,i-1) = E(i)
       F(i-1,i) = E(i)
      end do
      return
      end
c-------------------------------------------------------------------------------------------
c parallel diagonalization
c------------------------------------------------------------------------------------------
      subroutine diag_para(F,NRR,NB,BB,EIG)

#include "cdft.fh"
#include "eaf.fh"
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
      integer NRR,NB,g_F,g_alpha
      integer iD,lD,iE,lE,iE2,lE2,iA,lA
      double precision F(*),BB(*),EIG(*)
      double precision total_wall, total_cpu     

c      write(6,*)"i am in diag_para",nb,nrr

             if (.not.MA_Push_Get(MT_Dbl,nb,'D',lD,iD))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)
c      write(6,*)"i am here"
             if (.not.MA_Push_Get(MT_Dbl,nb,'E',lE,iE))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)
c        write(6,*)"i am here 1"
             if (.not.MA_Push_Get(MT_Dbl,nb,'E2',lE2,iE2))
     &   call errquit('diag_para: cannot allocate E2',0,
     1    MA_ERR)
c          write(6,*)"i am here 2"
             if (.not.MA_Push_Get(MT_Dbl,nrr,'A',lA,iA))
     &   call errquit('diag_para: cannot allocate A',0,
     1    MA_ERR)


       do i = 1, NRR
c      write(6,*)"F(I)",F(I)
          Dbl_MB(iA+I-1)= F(I)
       end do
      call tred3_sg(NB,NRR,Dbl_MB(iA),Dbl_MB(iD),Dbl_MB(iE),Dbl_MB(iE2))
       call build_fock(F,Dbl_MB(iD),Dbl_MB(iE),NB)
c      write(6,*)"after build_fock"
c      write(6,*)"here"
c     call ga_dgop(Msg_buildfock,F,NB*NB,'+')!Propagate status 
c      write(6,*)"F",(F(I),I=1,nrr)
          if (.not. ga_create(mt_dbl, nb, nb, 'F', 0, 0, g_F))
     & call errquit('Could not create g_F', 0, GA_ERR)
c        write(6,*)"created g_F"
         if (.not.ga_create(mt_dbl, nb, nb, 'alpha', 0, 0, g_alpha))
     & call errquit('Could not create g_alpha', 0, GA_ERR)
c      write(6,*)"before ga_put"
       call ga_put(g_F,1,nb,1,nb,F,nb)

c     iotal_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c      write(6,*)"ga put done"
C--------------------------------------------------------------------------------------
C DIAGONALIZATION OF THE MATRIX
C---------------------------------------------------------------------------------------
c      call ga_print(g_F)
       call ga_diag_std(g_F,g_alpha,eig)
C--------------------------------------------------------------------------------------
c      write(6,*)"ga diag_std"
c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "diag_para cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

       call ga_get(g_alpha,1,nb,1,nb,BB,nb)
            if (.not. ga_destroy(g_F))
     & call errquit('Could not destroy g_F', 0, GA_ERR)
c           if (.not. ga_destroy(g_alpha))
c    & call errquit('Could not destroy g_alpha', 0, GA_ERR)
c      call ga_dgop(Msg_diag_para,F,NB*NB,'+')!Propagate status
c      call ga_dgop(Msg_diag_para,BB,NB*NB,'+')
       call TRBAK3(nb,nb,nrr,Dbl_MB(iA),nb,BB)
       call ga_put(g_alpha,1,nb,1,nb,BB,nb)

c      call ga_print(g_alpha)
            if (.not. ga_destroy(g_alpha))
     & call errquit('Could not destroy g_alpha', 0, GA_ERR)
c       do i =1, nb*nb
c        write(6,*)"BB",BB(I)
c       end do
      if (.not.ma_pop_stack(lA))
     &   call errquit('diag_para: cannot pop stack A',3, MA_ERR)
      if (.not.ma_pop_stack(lE2))
     &   call errquit('scfit: cannot pop stack E2',3, MA_ERR)
      if (.not.ma_pop_stack(lE))
     &   call errquit('scfit: cannot pop stack E',3, MA_ERR)
      if (.not.ma_pop_stack(lD))
     &   call errquit('scfit: cannot pop stack D',3, MA_ERR)
       RETURN
       END



C---------------------------------------------------------------------------------------
C Unitary trabsformation of parse matrix
C--------------------------------------------------------------------------------------
      subroutine sparse_mat(F,NRR,NB,FA)
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
      integer NRR,NB
      integer iD,lD,iE,lE,iE2,lE2,iA,lA
      double precision F(*),FA(*)


c      write(6,*)"i am in diag_para"

             if (.not.MA_Push_Get(MT_Dbl,nb,'D',lD,iD))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nb,'E',lE,iE))
     &   call errquit('diag_para: cannot allocate E',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nb,'E2',lE2,iE2))
     &   call errquit('diag_para: cannot allocate E2',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nrr,'A',lA,iA))
     &   call errquit('diag_para: cannot allocate A',0,
     1    MA_ERR)
       do i = 1, NRR
          Dbl_MB(iA+I-1)= F(I)
       end do
      call tred3_sg(NB,NRR,Dbl_MB(iA),Dbl_MB(iD),Dbl_MB(iE),Dbl_MB(iE2))
       call build_fock(F,Dbl_MB(iD),Dbl_MB(iE),NB)
       do i = 0, nrr-1
        write(6,*)"A",Dbl_MB(iA+i)
       end do
c      call TRBAK3(nb,nb,nrr,Dbl_MB(iA),nb,FA)

      if (.not.ma_pop_stack(lA))
     &   call errquit('diag_para: cannot pop stack A',3, MA_ERR)
      if (.not.ma_pop_stack(lE2))
     &   call errquit('scfit: cannot pop stack E2',3, MA_ERR)
      if (.not.ma_pop_stack(lE))
     &   call errquit('scfit: cannot pop stack E',3, MA_ERR)
      if (.not.ma_pop_stack(lD))
     &   call errquit('scfit: cannot pop stack D',3, MA_ERR)
       RETURN
       END
C-------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------
      SUBROUTINE XTHX(NIN,H,NH1,X,NM,IREAL)
C***********************************************************************
C
C     SUBROUTINE XTHX CALCULATES THE PRODUCT (X)(H)(X)
C     WHERE H IS THE REAL (IREAL = 0) OR THE IMAGINARY PART (IREAL = 1)
C     OF  A HERMITAN MATRIX AND X IS REAL TRIANGULAR.
C     ONLY THE LOWER TRIANGLE OF THE MATRICES IS STORED.
C
C     WRITTEN FEBRUARY 1990 BY M.KOTZIAN AND N.ROESCH, GARCHING
C
C     CALLS FUNCTION IDIJM
C
C----------------------------------------------------------------------
c     INCLUDE '../include/param.cmn'
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
c     PARAMETER ZERO = 0.0d0
      DIMENSION    NIN(*)
      DIMENSION    X(*), H(*), SC(NH1*NH1), SCR(NH1*NH1), SCRM(NH1*NH1)
C
      ZERO = 0.0D0
      IF (IREAL.EQ.0) THEN
         HSIGN = 1.0D0
      ELSE
         HSIGN = -1.0D0
      END IF
C
      DO 70 J = 1, NH1
         JMD = MOD(J,NM)
         DO 10 L = J, NH1
            IF (MOD(L,NM).EQ.JMD) THEN
               JL = IDIJM(NIN,L,J,NM)
               SC(L) = X(JL)
            ELSE
               SC(L) = ZERO
            END IF
   10    CONTINUE
         DO 40 I = 1, NH1
            SCR(I) = ZERO
            IF (J.LE.I) THEN
               DO 20 K = 1, I
                  IK = NIN(I) + K
                  SCR(I) = SCR(I) + H(IK)*SC(K)
   20          CONTINUE
               KST = I + 1
            ELSE
               KST = J
            END IF
C
            DO 30 K = KST, NH1
               IK = NIN(K) + I
               SCR(I) = SCR(I) + H(IK)*SC(K)*HSIGN
   30       CONTINUE
C
   40    CONTINUE
C

         DO 60 I = J, NH1
            IJ = NIN(I) + J
            SS = ZERO
            IMD = MOD(I,NM)
            DO 50 K = 1, I
               IF (MOD(K,NM).EQ.IMD) THEN
                  IK = IDIJM(NIN,I,K,NM)
                  SS = SS + X(IK)*SCR(K)
               END IF
   50       CONTINUE
            SCRM(IJ) = SS
   60    CONTINUE
   70 CONTINUE
C

      IJ = 0
      DO 90 I = 1, NH1
         DO 80 J = 1, I
            IJ = IJ + 1
            H(IJ) = SCRM(IJ)
   80    CONTINUE
   90 CONTINUE
C
      RETURN

      END
C-------------------------------------------------------
C***********************************************************************
      INTEGER FUNCTION IDIJM(NIN,I,J,M)
C***********************************************************************
C
      INTEGER      NIN(*)
      INTEGER      I, J, M
C
      IDI = (I-1)/M + 1
      IDJ = (J-1)/M + 1
      IDIJM = NIN(IDI) + IDJ
C
      RETURN

      END

