      subroutine semiemp_os_fock_build(g_zdens,g_zdensa,SUMELEL,QQ,LABI
     1 ,LABJ,LABK,LABL,LL,NIN,nbf,Q,NG,GAMMA,g_zfocka,g_zfockb,ITAG,
     1 GPARAMA,GPARAMB)

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

      integer LABI(*), LABJ(*), LABK(*), LABL(*), LL, NIN(*), NG(*)
      double precision QQ(*), SUMELEL(*), Q(*), GAMMA(*), GPARAMA 
      double precision GPARAMB
      integer g_zdens, g_zfp, g_core
      integer g_zdensa
      integer g_densre          !real part of dens mat
      integer g_densim          !imag part of dens mat
      integer g_fockare          !real part of fock mat
      integer g_fockaim          !imag part of fock mat
      integer g_fockbre          !real part of fock mat
      integer g_fockbim          !imag part of fock mat
      integer g_densare
      integer g_densaim
      integer g_zfocka
      integer g_zfockb
      integer g_zfock, g_tmp1, g_tmp2, g_tmp3
      integer g_ca, g_cb, g_zca, g_zcb 
      logical status
      double precision TEST
      integer i,itag
      integer nbf,nrr
      integer idensre,ldensre, lfockare, ifockare, lfockaim, ifockaim
      integer lfockbre, ifockbre, lfockbim, ifockbim
      integer ihdense, lhdense, idensim, ldensim, lfock, ifock
      integer leig, ieig, lcomut, icomut
      integer msg_diag_para
      double precision emina, eminb
      double precision total_wall, total_cpu
      DATA TEST /1.0D-10/

         nrr = nbf*(nbf+1)/2
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()


      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'densre',ldensre,idensre))
     &   call errquit('semiemp_cs_fock_build: cannot allocate densre',0,
     1    MA_ERR)
     
      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockre',lfockare,ifockare))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockre',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'hdense',lhdense,ihdense))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockre',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'densim',ldensim,idensim))
     &   call errquit('semiemp_cs_fock_build: cannot allocate densim',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockaim',lfockaim,ifockaim))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockim',0,
     1    MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockbre',lfockbre,ifockbre))
     &   call errquit('semiemp_cs_fock_build: cannot allocate fockim',0,
     1    MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'fockbim',lfockbim,ifockbim))
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
      if (.not. ga_duplicate(g_densre, g_densare, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_densaim, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockare, "Re[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockaim, "Im[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockbre, "Re[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_fockbim, "Im[F]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      call ga_zero(g_densre)
      call ga_zero(g_densim)
      call ga_zero(g_densare)
      call ga_zero(g_densaim)
      call ga_zero(g_fockare)
      call ga_zero(g_fockaim)
      call ga_zero(g_fockbre)
      call ga_zero(g_fockbim)

      call convert_z2d (g_zdens, g_densre, g_densim)
      call convert_z2d (g_zdensa, g_densare, g_densaim) 
c     GO TO 100

c     call ga_get(g_densre,1,nbf,1,nbf,Dbl_MB(idensre),nbf) 
c     call ga_get(g_densim,1,nbf,1,nbf,Dbl_MB(idensim),nbf)


cSG   call mat_triangular(Dbl_MB(idensre),nbf)


cSG   call mat_triangular(Dbl_MB(idensim),nbf)

      call ga_trimat(g_densre,dbl_mb(idensre),nbf)   
      call ga_trimat(g_densare,dbl_mb(ihdense),nbf) 


      CALL FILLFU(Dbl_MB(ifockare),Dbl_MB(ifockbre),Dbl_MB(ihdense),
     1 Dbl_MB(idensre),Q,GAMMA,SUMELEL,NG,QQ,NRR,LABI,LABJ,LABK,LABL,LL,
     2 NIN,1)

      call ga_trimat(g_densim,dbl_mb(idensim),nbf)
      call ga_trimat(g_densaim,dbl_mb(ihdense),nbf)

c     do i = 0, nbf*nbf-1
c      write(6,*)"fock"!,Dbl_MB(ifockre+i)
c     end do


      CALL FILLFU(Dbl_MB(ifockaim),Dbl_MB(ifockbim),Dbl_MB(ihdense),
     1 Dbl_MB(idensim),Q,GAMMA,SUMELEL,NG,QQ,NRR,LABI,LABJ,LABK,LABL,LL,
     2 NIN,2)


      if (itag.eq.1) then

       do i =0 , nrr-1
        Dbl_MB(ifock+i) = Dbl_MB(ifockare+i)
c       if (ga_nodeid().eq.0) write(6,*) "fock"!,Dbl_MB(ifock+i)
       end do

      call diag_para(Dbl_MB(ifock),NRR,nbf,Dbl_MB(idensre),Dbl_MB(ieig))

c     call ga_put(g_ca,1,nbf,1,nbf,Dbl_MB(idensre),nbf)
c     GPARAMA = Dbl_MB(ieig+nbf-1)
      emina = Dbl_MB(ieig+0)
      GPARAMA = (Dbl_MB(ieig+nbf-1)-Dbl_MB(ieig+0))
      if (ga_nodeid().eq.0) write(6,*)"Scaling factor ALPHA",GPARAMA

       do i =0 , nrr-1
        Dbl_MB(ifock+i) = Dbl_MB(ifockbre+i)
c       if (ga_nodeid().eq.0) write(6,*) "fock",Dbl_MB(ifock+i)
       end do

      call diag_para(Dbl_MB(ifock),NRR,nbf,Dbl_MB(idensre),Dbl_MB(ieig))

c     call ga_put(g_cb,1,nbf,1,nbf,Dbl_MB(idensre),nbf)
c     GPARAMB = Dbl_MB(ieig+nbf-1)
      eminb = Dbl_MB(ieig+0)
      GPARAMB = (Dbl_MB(ieig+nbf-1)-Dbl_MB(ieig+0))     
      if (ga_nodeid().eq.0) write(6,*)"Scaling factor BETA",GPARAMB 
      end if
      



      call ga_sqmat(g_fockare,dbl_mb(ifockare),nbf,1)
      call ga_sqmat(g_fockaim,dbl_mb(ifockaim),nbf,2)
      call ga_sqmat(g_fockbre,dbl_mb(ifockbre),nbf,1)
      call ga_sqmat(g_fockbim,dbl_mb(ifockbim),nbf,2)

c 150 CONTINUE
      call convert_d2z (1d0, g_fockare, 1d0, g_fockaim, g_zfocka)
      call convert_d2z (1d0, g_fockbre, 1d0, g_fockbim, g_zfockb)

      call convert_d2z (1d0, g_densre, 1d0, g_densim, g_zdens)
      call convert_d2z (1d0, g_densare, 1d0, g_densaim, g_zdensa)


      status=.true.
      status=status.and.ga_destroy(g_densare)
      status=status.and.ga_destroy(g_densaim)
      status=status.and.ga_destroy(g_densre)
      status=status.and.ga_destroy(g_densim)
      status=status.and.ga_destroy(g_fockbre)
      status=status.and.ga_destroy(g_fockbim)
      status=status.and.ga_destroy(g_fockare)
      status=status.and.ga_destroy(g_fockaim)


      if (.not.ma_pop_stack(leig))
     &   call errquit('semiemp_cs_fock_build: cannot pop stack fock',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfock))
     &   call errquit('semiemp_cs_fock_build: cannot pop stack fock',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfockbim))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockbim',3
     1 ,MA_ERR)

      if (.not.ma_pop_stack(lfockbre))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockibre',
     1 3,MA_ERR)

      if (.not.ma_pop_stack(lfockaim))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockim',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(ldensim))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack densim',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lhdense))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockre',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(lfockare))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack fockre',3,
     1 MA_ERR)

      if (.not.ma_pop_stack(ldensre))
     &  call errquit('semiemp_cs_fock_build: cannot pop stack densre',3,
     1 MA_ERR)


      end subroutine


c---------------------------------------------------------------------------------
C
C---------------------------------------------------------------------------------

