      subroutine semiemp_cs_fock_build(g_zdens,SUMELEL,LL,NIN,nbf,Q,
     1 NG,GAMMA,g_zfock,ITAG,GPARAM,emin,g_c)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
c#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      integer  LL, NIN(*), NG(*)
      double precision SUMELEL(*), Q(*), GAMMA(*), GPARAM, emin
      integer g_zdens, g_zfp, g_c
      integer g_densre          !real part of dens mat
      integer g_densim          !imag part of dens mat
      integer g_fockre          !real part of fock mat
      integer g_fockim          !imag part of fock mat
      integer g_zfock, g_tmp1, g_tmp2, g_tmp3
      logical status
      double precision TEST
      integer i,itag
      integer nbf,nrr
      integer idensre,ldensre, lfockre, ifockre, lfockim, ifockim
      integer ihdense, lhdense, idensim, ldensim, lfock, ifock
      integer leig, ieig, lcomut, icomut
      integer msg_diag_para
      double precision total_wall, total_cpu
      integer homo
      DATA TEST /1.0D-5/

         nrr = nbf*(nbf+1)/2
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c     write(6,*)"i am in fock build"
       
      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'densre',ldensre,idensre))
     &   call errquit('semiemp_cs_fock_build: cannot allocate densre',0,
     1    MA_ERR)
     
      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockre',lfockre,ifockre))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockre',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'hdense',lhdense,ihdense))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockre',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'densim',ldensim,idensim))
     &   call errquit('semiemp_cs_fock_build: cannot allocate densim',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockim',lfockim,ifockim))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockim',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fock',lfock,ifock))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fock',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf,'EIG',leig,ieig))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fock',0,
     1    MA_ERR)


      if (.not. ga_create(mt_dbl, nbf, nbf, "Re[Dens]", 0, 0, g_densre))
     $     call errquit ("couldnt create g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_densim, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockre, "Re[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockim, "Im[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_tmp1, "tmp1"))
     $     call errquit ("couldnt duplicate g_tmp1", 0, GA_ERR)
      if (.not. ga_duplicate(g_zdens, g_tmp2, "tmp2"))
     $     call errquit ("couldnt duplicate g_tmp2", 0, GA_ERR)
      call ga_zero(g_fockre)
      call ga_zero(g_fockim)
      call ga_zero(g_densre)
      call ga_zero(g_densim)
c     call ga_print(g_zdens)
      call convert_z2d (g_zdens, g_densre, g_densim)
c     GO TO 100
c     call ga_get(g_densre,1,nbf,1,nbf,Dbl_MB(idensre),nbf) 
c     call ga_get(g_densim,1,nbf,1,nbf,Dbl_MB(idensim),nbf)


cSG   call mat_triangular(Dbl_MB(idensre),nbf)


cSG   call mat_triangular(Dbl_MB(idensim),nbf)

      call ga_trimat(g_densre,dbl_mb(idensre),nbf)   
      call ga_trimat(g_densim,dbl_mb(idensim),nbf) 

      CALL FILLF(Dbl_MB(ifockre),Dbl_MB(ihdense),Dbl_MB(idensre),Q,GAMMA
     1 ,SUMELEL,NG,NRR,LL,NIN,1,nbf)
      

c     do i = 0, nbf*nbf-1
c      write(6,*)"fock",Dbl_MB(ifockre+i)
c     end do

      CALL FILLF(Dbl_MB(ifockim),Dbl_MB(ihdense),Dbl_MB(idensim),Q,GAMMA
     1 ,SUMELEL,NG,NRR,LL,NIN,2,nbf)
c     CALL MATPAC(nbf,nrr,nrr,Dbl_MB(idensre),Dbl_MB(ifockre),NIN,
c    1 Dbl_MB(icomut))

c     do i = 0, nbf*(nbf+1)/2-1
c      write(6,*)"fockim",Dbl_MB(ifockim+i)
c      Dbl_MB(ifockim+i) = 0.0d0
c     end do


cSG   call mat_sym(Dbl_MB(idensre),nbf,1)
cSG   call ga_put(g_densre,1,nbf,1,nbf,Dbl_MB(idensre),nbf)
c 100 CONTINUE

c     call ga_fillf(g_densre,g_fockre,nbf,LABI,LABJ,LABK,LABL,QQ,
c    1 g_core,LL,1,NIN)

c     call ga_fillf(g_densim,g_fockim,nbf,LABI,LABJ,LABK,LABL,QQ,
c    1 g_core,LL,2,NIN)

c     call ga_print(g_fockre)
      if (itag.eq.1) then

cSG   call ga_get(g_fockre,1,nbf,1,nbf,Dbl_MB(ifock),nbf)
cSG   call mat_triangular(Dbl_MB(ifock),nbf)
       do i =0 , nrr-1
        Dbl_MB(ifock+i) = Dbl_MB(ifockre+i)
c       if (ga_nodeid().eq.0) write(6,*) "fock",Dbl_MB(ifock+i)
       end do

      call diag_para(Dbl_MB(ifock),NRR,nbf,Dbl_MB(idensre),Dbl_MB(ieig))
      call ga_put(g_c,1,nbf,1,nbf,Dbl_MB(idensre),nbf)
c     call ga_print(g_core)
c     call ga_sqmat(g_core,dbl_mb(idensre),nbf,1)
c     write(6,*)"HOMO",dbl_MB(ieig+homo-1)
c     write(6,*)"HOMO-1",dbl_MB(ieig+homo-2)
c     write(6,*)"HOMO-7",dbl_MB(ieig+homo-8)

      GPARAM = (Dbl_MB(ieig+nbf-1)-Dbl_MB(ieig+0))
      emin = (Dbl_MB(ieig+0))
      if (ga_nodeid().eq.0) write(6,*)"Scaling factor",GPARAM      
      end if
c     GO TO 150         
cSG   call mat_sym(Dbl_MB(ifockre),nbf,1)

cSG   call mat_sym(Dbl_MB(ifockim),nbf,2)

cSG   call ga_put(g_fockre,1,nbf,1,nbf,Dbl_MB(ifockre),nbf)


cSG   call ga_put(g_fockim,1,nbf,1,nbf,Dbl_MB(ifockim),nbf)
      call ga_sqmat(g_fockre,dbl_mb(ifockre),nbf,1)
      call ga_sqmat(g_fockim,dbl_mb(ifockim),nbf,2)
c 150 CONTINUE
c     call ga_print(g_fockre)
      call convert_d2z (1d0, g_fockre, 1d0, g_fockim, g_zfock)

      call convert_d2z (1d0, g_densre, 1d0, g_densim, g_zdens)
c     if (itag.eq.1) then
c     call ga_zgemm ("N", "N", nbf, nbf,
c    $     nbf, z1, g_zfock, g_zdens, z0, g_tmp1)

c     call ga_print(g_zfock)

c     call ga_print(g_tmp1)
c     call ga_print(g_zdens)
c     call ga_dgemm ("N", "N", nbf, nbf,
c    $     nbf, z1, g_zdens, g_zfock, z0, g_tmp2)

c     call ga_add (z1,g_tmp1, -z1, g_tmp2, g_zfp)

c      if ( .not.mat_is_same (g_tmp1,g_tmp2,test))
c    &   call errquit('prop_chebyshev: not stationary density',0,
c    1    MA_ERR)
c      write(6,*)"density is stationary"
c      end if

c     call ga_print(g_zfock)
c     call ga_print(g_zdens)

c     total_wall = total_wall +  util_wallsec()
c     total_cpu  =  total_cpu + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "fock build cpu time", total_cpu, "total time",
c    1 total_wall
c     end if

      status=.true.
      status=status.and.ga_destroy(g_tmp2)
      status=status.and.ga_destroy(g_tmp1)
      status=status.and.ga_destroy(g_densre)
      status=status.and.ga_destroy(g_densim)
      status=status.and.ga_destroy(g_fockre)
      status=status.and.ga_destroy(g_fockim)


      if (.not.ma_pop_stack(leig))
     &   call errquit('semiemp_cs_fock_build: cannot pop stack fock',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfock))
     &   call errquit('semiemp_cs_fock_build: cannot pop stack fock',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfockim))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockim',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(ldensim))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack densim',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lhdense))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockre',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfockre))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockre',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(ldensre))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack densre',3,
     1 MA_ERR)


      end subroutine


c---------------------------------------------------------------------------------
C
C---------------------------------------------------------------------------------
      SUBROUTINE MATPAC(NN,NR1,NR2,AA,BB,NIN,CC)
C     IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     COMMUTATOR OF TWO TRIANGULAR MATRICES.
C
C      C  = A*B - B*A  = A*B - (A*B)T  = -(C)T
C
c     INCLUDE '../include/param.cmn'
c     INCLUDE '../include/units.cmn'
      INTEGER NN, NR1, NR2,NIN(*)
      DOUBLE PRECISION    AA(*), BB(*), CC(*)
C
c     INCLUDE '../include/z.cmn'
C

C     WRITE(6,*) ' MATPAC',NN,NR1,NR2
      DO 10 J = 1, NR2
         CC(J) = 0.0D0
   10 CONTINUE
C
C
      DO 40 J = 1, NN
         DO 30 I = 1, NN
            IC1 = NIN(MAX(I,J)) + MIN(I,J)
            FAC = 1.0D0
            IF (I.LT.J) FAC = -1.0D0*FAC
            SUM = 0.0D0
            DO 20 K = 1, NN
               IA = NIN(MAX(I,K)) + MIN(I,K)
               IB = NIN(MAX(J,K)) + MIN(J,K)
               SUM = SUM + BB(IA)*AA(IB)
   20       CONTINUE
            CC(IC1) = CC(IC1) + SUM*FAC
   30    CONTINUE
   40 CONTINUE
C
C
      DO 50 J = 1, NN
         IDI = (J* (J+1))/2
         CC(IDI) = 0.0D0
   50 CONTINUE

      write(6,*)"CC",(CC(I),I=1,nr1)
C
C     RETURN TO UDIIS
C
C     WRITE(6,*) ' LEAVING MATPAC'

      RETURN

      END

