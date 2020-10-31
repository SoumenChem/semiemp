      subroutine semiemp_cis(rtdb,nbf,eig,FT)
*
* $Id: movecs_pr_anal.F 27889 2015-12-04 03:24:50Z edo $
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"

c
      integer rtdb
      integer iBB, lBB, iW, lW, ieig, leig
      integer nbf, nmo
      integer lwiajb, iwiajb, lwijab, iwijab
      integer MIM(2)
      integer g_coeff,g_ft, g_temp
      integer ndocc, nvirt, ndim
      integer i, a, j, b, ia, jb
      integer lham, iham
      double precision eig(*), FT(*)

       if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'BB',lBB,iBB))
     &   call errquit('mo_print: cannot allocate BB',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:coefficient',mt_dbl,nbf*nbf,
     1 Dbl_MB(iBB)))
     1 call errquit('mo_print: failed to write coefficient',0, RTDB_ERR)

c      if (.not.ga_create(mt_dbl, nbf, nbf, 'Coeff', 0, 0, g_coeff))
c     & call errquit('Could not create g_coff', 0, GA_ERR)
c
c      if (.not.ga_create(mt_dbl, nbf, nbf, 'FT', 0, 0, g_ft))
c     & call errquit('Could not create g_ft', 0, GA_ERR)
c
c      if (.not.ga_create(mt_dbl, nbf, nbf, 'Temp', 0, 0, g_temp))
c     & call errquit('Could not create g_temp', 0, GA_ERR)
c
c       call ga_put(g_coeff,1,nbf,1,nbf,Dbl_MB(iBB),nbf)
c
c      call ga_sqmat(g_ft,FT,nbf,1)
c
c      call ga_dgemm('N','N',nbf,nbf,nbf,1.0d0,g_coeff,g_ft,0.0d0,
c     1 g_temp)
c
c      call ga_dgemm('N','N',nbf,nbf,nbf,1.0d0,g_temp,g_ft,0.0d0,
c     1 g_coeff)
c      
c      call ga_get(g_coeff,1,nbf,1,nbf,Dbl_MB(iBB),nbf) 
c
c      if (.not. ga_destroy(g_temp))
c     & call errquit('Could not destroy g_temp', 0, GA_ERR)       
c      if (.not. ga_destroy(g_ft))
c     & call errquit('Could not destroy g_ft', 0, GA_ERR)       
c      if (.not. ga_destroy(g_coeff))
c     & call errquit('Could not destroy g_coeff', 0, GA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:MIM',mt_int,2,MIM))
     1  call errquit('semiemp: failed to write MIM',0, RTDB_ERR)

        ndocc = MIM(1)
        nvirt = nbf-ndocc
        ndim = ndocc*nvirt
        write(6,*)"ndocc,nvirt,ndim",ndocc,nvirt,ndim

      if (.not.MA_Push_Get(MT_Dbl,ndim*ndim,'WIAJB',lwiajb,
     1 iwiajb))
     & call errquit('semiemp_energy: cannot allocate WIAJB',0,
     1    MA_ERR)

      if (.not.MA_Push_Get(MT_Dbl,ndim*ndim,'WIJAB',lwijab,
     1 iwijab))
     & call errquit('semiemp_energy: cannot allocate WIAJB',0,
     1    MA_ERR)

       call semiemp_ao2mo(rtdb,Dbl_MB(iBB),nbf,ndocc,Dbl_MB(iwiajb),
     1 Dbl_MB(iwijab))

        

c       if (.not.MA_Push_Get(MT_Dbl,ndim*ndim,'Ham',lham,iham))
c    &  call errquit('semiemp_energy: cannot allocate Ham',0,MA_ERR)

        call semiemp_singlet_cis(nbf,ndocc,nvirt,eig,
     1  Dbl_MB(iwijab),Dbl_MB(iwiajb))

c       if (.not.ma_pop_stack(lham))
c    & call errquit('semiemp_energy: cannot pop stack ham',3, MA_ERR)   
        if (.not.ma_pop_stack(lwijab))
     & call errquit('semiemp_energy: cannot pop stack wijab',3, MA_ERR)
        if (.not.ma_pop_stack(lwiajb))
     & call errquit('semiemp_energy: cannot pop stack wiajb',3, MA_ERR)     
       if (.not.ma_pop_stack(lBB))
     &   call errquit('mo_print: cannot pop stack BB',3, MA_ERR)

       return
       end
C-------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------
       subroutine semiemp_singlet_cis(nbf,ndocc,nvirt,eig,WIJAB,
     1 WIAJB)

       implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"
        
        integer rtdb
        integer ndocc, nvirt,nbf
c       double precision Ham(ndocc*nvirt,ndocc*nvirt)
        double precision eig(*)
        double precision Ham, E1
        double precision WIJAB(ndocc,ndocc,nvirt,nvirt)
        double precision WIAJB(ndocc,nvirt,ndocc,nvirt)
        integer i, j, a, b, ia, jb, ndim
        integer iocc, jocc, avirt, bvirt
        integer g_ham, g_alpha
        integer ihameig, lhameig, ihamalpha, lhamalpha
        integer nrr


        ndim = nvirt*ndocc

        if (.not.ga_create(mt_dbl, ndim, ndim, 'Ham', 0, 0, g_ham))
     & call errquit('Could not create g_ham', 0, GA_ERR)
       if (.not.ga_create(mt_dbl,ndim,ndim, 'Alpha', 0, 0, g_alpha))
     & call errquit('Could not create g_alpha', 0, GA_ERR)


        call ga_zero(g_ham)
        call ga_zero(g_alpha)
c       call ga_put(g_ham,1,1,1,1,0.0d0,1)
c       Ham = 0.0d0
c       write(6,*)"ndocc,nvirt in cis",ndocc,nvirt
        do i = 1, ndocc
         do a = 1, nvirt
          ia = (i-1)*nvirt + a
           do j = 1, ndocc
            do b = 1, nvirt
             jb = (j-1)*nvirt + b
                Ham = 0.0d0
c          Ham  = -WIJAB(i,j,a,b)
           Ham = -1.0d0*WIJAB(i,j,a,b) + 2.0d0*WIAJB(i,a,j,b)
               if (ia.eq.jb) Ham = Ham + eig(ndocc+a)-eig(i)
c             write(6,*)"a,b,ia,ib",a,b,ia,jb
              call ga_put(g_ham,ia,ia,jb,jb,Ham,1)
              end do 
            end do
          end do
         end do

      if (.not.MA_Push_Get(MT_Dbl,ndim,'hameig',lhameig,ihameig))
     & call errquit('semiemp_energy: cannot allocate HAMEIG',0,
     1    MA_ERR)


c      NRR= ndim*(ndim+1)/2
c      call ga_print(g_ham)
c          write(6,*)"done Hamiltonian"
       call ga_diag_std(g_ham,g_alpha,dbl_mb(ihameig))
c          write(6,*)"after diagonalization"

       write(6,*) "SINGLET EXCITED STATE"
       do i = 1, ndim 
         write(6,*)"EIG=", Dbl_MB(ihameig+i-1)
       end do 

        call ga_zero(g_ham)
        call ga_zero(g_alpha)
c       call ga_put(g_ham,1,1,1,1,0.0d0,1)

        do i = 1, ndocc
         do a = 1, nvirt
          ia = (i-1)*nvirt + a 
           do j = 1, ndocc
            do b = 1, nvirt
             jb = (j-1)*nvirt + b 
                Ham = 0.0d0
           Ham  = -WIJAB(i,j,a,b)
               if (ia.eq.jb) Ham = Ham + eig(ndocc+a)-eig(i)
c             write(6,*)"a,b,ia,ib",a,b,ia,jb
              call ga_put(g_ham,ia,ia,jb,jb,Ham,1)
              end do
            end do
          end do
         end do

         call ga_diag_std(g_ham,g_alpha,dbl_mb(ihameig))
c          write(6,*)"after diagonalization"
       write(6,*)"TRIPLET EXCITED STATE"
       do i = 1, ndim  
         write(6,*)"EIG=", Dbl_MB(ihameig+i-1)
       end do



       if (.not.ma_pop_stack(lhameig))
     & call errquit('semiemp_energy: cannot pop stack hameig',3, MA_ERR)     

      if (.not. ga_destroy(g_alpha))
     & call errquit('Could not destroy g_alpha', 0, GA_ERR)
       if (.not. ga_destroy(g_ham))
     & call errquit('Could not destroy g_ham', 0, GA_ERR)

         return
         end        
C-------------------------------------------------------------------------------
C AO TO MO TRANFORMATION
C-------------------------------------------------------------------------------

       subroutine semiemp_ao2mo(rtdb,BB,nbf,ndocc,WIAJB,WIJAB)
       
       implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"

       integer rtdb, nbf, nintgl,ndocc, nvirt
       integer ilabi,llabi,ilabj,llabj,ilabk,llabk,ilabl,llabl
       integer imoi, lmoi, imoj, lmoj, imok, lmok, imol, lmol
       integer iqq, lqq, iro, lro
       integer lwijkl, iwijkl
       integer iwiqrs, lwiqrs, iwijrs, lwijrs, iwijks, lwijks
       double precision BB(nbf,nbf), WIJAB(*),WIAJB(*)

        if (.not.rtdb_get(rtdb,'semiemp:numint',mt_int,1,nintgl))
     &  call errquit('semiemp failed read to inttyp',0,RTDB_ERR)

          if (.not.MA_Push_Get(MT_INT,NINTGL,'LABI',lLABI,
     1   iLABI))
     &   call errquit('semiemp_energy: cannot allocate LABI',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABJ',lLABJ,
     1   iLABJ))
     &   call errquit('semiemp_energy: cannot allocate LABJ',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABK',lLABK,
     1   iLABK))
     &   call errquit('semiemp_energy: cannot allocate LABK',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,NINTGL,'LABL',lLABL,
     1   iLABL))
     &   call errquit('semiemp_energy: cannot allocate LABL',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,NINTGL,'QQ',lQQ,iQQ))
     &   call errquit('semiemp_energy: cannot allocate QQ',0,
     1    MA_ERR)

           
        call semiemp_integral_load(NINTGL,INT_MB(iLABI),INT_MB(iLABJ),
     A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ))
        
          write(6,*) "NINTGL", NINTGL
c          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIQRS',lwiqrs,
c     1    iwiqrs))
c     &   call errquit('semiemp_energy: cannot allocate WIqrs',0,
c     1    MA_ERR)
c
c          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIJRS',lwijrs,
c     1    iwijrs))
c     &   call errquit('semiemp_energy: cannot allocate WIJRS',0,
c     1    MA_ERR)
c
c          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIJKS',lwijks,
c     1    iwijks))
c     &   call errquit('semiemp_energy: cannot allocate WIJKS',0,
c     1    MA_ERR)
          nvirt = nbf - ndocc
         call semiemp_mo_trans(NINTGL,INT_MB(iLABI),INT_MB(iLABJ),
     A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ),WIJAB,WIAJB,BB,nbf,
     A   ndocc,nvirt)

c        if (.not.ma_pop_stack(lwijks))
c     & call errquit('semiemp_energy: cannot pop stack wijks',3, MA_ERR) 
c        if (.not.ma_pop_stack(lwijrs))
c     & call errquit('semiemp_energy: cannot pop stack wijrs',3, MA_ERR) 
c        if (.not.ma_pop_stack(lwiqrs))
c     & call errquit('semiemp_energy: cannot pop stack wiqrs',3, MA_ERR)          
        if (.not.ma_pop_stack(lQQ))
     & call errquit('semiemp_energy: cannot pop stack QQ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABL))
     & call errquit('semiemp_energy: cannot pop stack LABL',3, MA_ERR)
         if (.not.ma_pop_stack(lLABK))
     & call errquit('semiemp_energy: cannot pop stack LABK',3, MA_ERR)
         if (.not.ma_pop_stack(lLABJ))
     & call errquit('semiemp_energy: cannot pop stack LABJ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABI))
     & call errquit('semiemp_energy: cannot pop stack LABI',3, MA_ERR)

       return 
       end
C--------------------------------------------------------------------------------
C
C----------------------------------------------------------------------------------
       subroutine semiemp_mo_trans(NINTGL,LABI,LABJ,LABK,LABL,QQ,WIJAB,
     1 WIAJB,BB,nbf,ndocc,nvirt)

       implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "bas.fh"
#include "global.fh"
#include "inp.fh"
#include "cscfps.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"

       integer NINTGL, LABI(*), LABJ(*), LABK(*), LABL(*)
       integer nbf, LL, p, q, r , s, i, j, k, l, ndocc, nvirt
       integer a, b, avirt, bvirt
       double precision BB(nbf,nbf), QQ(*), V
       double precision WIAJB(ndocc,nvirt,ndocc,nvirt)
       double precision WIJAB(ndocc,ndocc,nvirt,nvirt)
c      double precision WIQRS(nbf,nbf,nbf,nbf), WIJRS(nbf,nbf,nbf,nbf)
c      double precision WIJKS(nbf,nbf,nbf,nbf), V
       
c      WIQRS = 0.0d0
c      WIJRS = 0.0d0
c      WIJKS = 0.0d0
c      WIJKL = 0.0d0
c      WIAJB = 0.0d0
c      WIJAB = 0.0d0
       do i = 1, ndocc
        do j = 1, ndocc
         do  a = 1, nvirt
          do b = 1, nvirt
         WIJAB(i,j,a,b) = 0.0d0
         WIAJB(i,a,j,b) = 0.0d0
          end do
         end do
        end do
       end do

         write(6,*)"i am mo_trans"
       do i = 1, ndocc
        do avirt = 1 , nvirt
         do j = 1, ndocc
          do bvirt = 1, nvirt
           a = ndocc + avirt
           b = ndocc + bvirt
           do LL = 1, nintgl
            p = LABI(LL)
            q = LABJ(LL)
            r = LABK(LL)
            s = LABL(LL)
            V = QQ(LL)
        WIAJB(i,avirt,j,bvirt)=WIAJB(i,avirt,j,bvirt)
     1         +BB(p,i)*BB(q,a)*BB(r,j)*BB(s,b)*V  
           end do
          end do
         end do
        end do
       end do

       do i = 1, ndocc
        do j = 1, ndocc
         do avirt = 1 , nvirt       
          do bvirt = 1, nvirt
           a = ndocc + avirt
           b = ndocc + bvirt          
           do LL = 1, nintgl
            p = LABI(LL)
            q = LABJ(LL)
            r = LABK(LL)
            s = LABL(LL)
            V = QQ(LL)
        WIJAB(i,j,avirt,bvirt)=WIJAB(i,j,avirt,bvirt)
     1       +BB(p,i)*BB(q,j)*BB(r,a)*BB(s,b)*V  
           end do
          end do
         end do
        end do
       end do
       
c       do LL = 1, nintgl
c        p = LABI(LL)
c        q = LABJ(LL)
c        r = LABK(LL)
c        s = LABL(LL)
c        V = QQ(LL)
c        do i = 1, nbf
c         WIQRS(i,q,r,s) = WIQRS(i,q,r,s) + BB(p,i)*V
c         do j = 1, nbf
c          WIJRS(i,j,r,s) = WIJRS(i,j,r,s)+BB(q,j)*WIQRS(i,q,r,s)
c          do k = 1, nbf
c           WIJKS(i,j,k,s)=WIJKS(i,j,k,s)+BB(r,k)*WIJRS(i,j,r,s)
c           do l = 1, nbf
c            WIJKL(i,j,k,l)=WIJKL(i,j,k,l)+BB(s,l)*WIJKS(i,j,k,s)
cc           write(6,*)i,j,k,l,WIJKL(i,j,k,l)
c            end do
c          end do
c         end do
c        end do
c       end do
       
c      do i = 1, nbf
c       do LL = 1, nintgl
c         p = LABI(LL)
c         q = LABJ(LL)
c         r = LABK(LL)
c         s = LABL(LL)
c         V = QQ(LL)
c         WIQRS(i,q,r,s) = WIQRS(i,q,r,s)+ BB(i,p)*V
c       end do
c       do j = 1, nbf
c        do LL = 1, nintgl
c         p = LABI(LL)
c         q = LABJ(LL)
c         r = LABK(LL)
c         s = LABL(LL)
c         WIJRS(i,j,r,s) = WIJRS(i,j,r,s)+ BB(j,q)*WIQRS(i,q,r,s)
c        end do
c        do k = 1, nbf
c          do LL = 1, nintgl
c           p = LABI(LL)
c           q = LABJ(LL)
c           r = LABK(LL)
c           s = LABL(LL)
c           WIJKS(i,j,k,s) = WIJKS(i,j,k,s)+ BB(k,r)*WIJRS(i,j,r,s)
c          end do
c          do l = 1, nbf
c           do LL = 1, nintgl
c            p = LABI(LL)
c            q = LABJ(LL)
c            r = LABK(LL)
c            s = LABL(LL)
c            WIJKL(i,j,k,l) = WIJKL(i,j,k,l)+ BB(l,s)*WIJKS(i,j,k,s)             
c            end do
c            write(6,*)i,j,k,l,WIJKL(i,j,k,l)
c          end do 
c         end do
c        end do
c       end do 

       return
       end

               

       
