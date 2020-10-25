      subroutine semiemp_cis(rtdb,nbf,eig)
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
      integer lwijkl, iwijkl
      integer MIM(2)
      integer ndocc, nvirt
      integer i, a, j, b, ia, jb
      integer lham, iham
      double precision eig(*)

       if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'BB',lBB,iBB))
     &   call errquit('mo_print: cannot allocate BB',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:coefficient',mt_dbl,nbf*nbf,
     1 Dbl_MB(iBB)))
     1 call errquit('mo_print: failed to write coefficient',0, RTDB_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIJKL',lwijkl,
     1 iwijkl))
     & call errquit('semiemp_energy: cannot allocate WIJKL',0,
     1    MA_ERR)
      
       call semiemp_ao2mo(rtdb,Dbl_MB(iBB),nbf,Dbl_MB(iwijkl))

        if (.not.rtdb_get(rtdb,'semiemp:MIM',mt_int,2,MIM))
     1  call errquit('semiemp: failed to write MIM',0, RTDB_ERR)
        
        ndocc = MIM(1)
        nvirt = nbf-ndocc

        write(6,*)"ndocc,nvirt",ndocc,nvirt

c       if (.not.MA_Push_Get(MT_Dbl,ndocc*nvirt*ndocc*nvirt,'Ham',lham,
c    1  iham))
c    &  call errquit('semiemp_energy: cannot allocate Ham',0,MA_ERR)

        call semiemp_singlet_cis(nbf,ndocc,nvirt,eig,Dbl_MB(iwijkl))

c       if (.not.ma_pop_stack(lham))
c    & call errquit('semiemp_energy: cannot pop stack ham',3, MA_ERR)   
        if (.not.ma_pop_stack(lwijkl))
     & call errquit('semiemp_energy: cannot pop stack wijkl',3, MA_ERR)     
       if (.not.ma_pop_stack(lBB))
     &   call errquit('mo_print: cannot pop stack BB',3, MA_ERR)

       return
       end
C-------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------
       subroutine semiemp_singlet_cis(nbf,ndocc,nvirt,eig,WIJKL)

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
        
        integer ndocc, nvirt,nbf
c       double precision Ham(ndocc*nvirt,ndocc*nvirt)
        double precision eig(*)
        double precision Ham
        double precision WIJKL(nbf,nbf,nbf,nbf)
        integer i, j, a, b, ia, jb, ndim
        integer iocc, jocc, avirt, bvirt
        integer g_ham

        ndim = nvirt*ndocc

        if (.not.ga_create(mt_dbl, ndim, ndim, 'Ham', 0, 0, g_ham))
     & call errquit('Could not create g_ham', 0, GA_ERR)

        call ga_zero(g_ham)

        write(6,*)"ndocc,nvirt in cis",ndocc,nvirt
        do iocc = 1, ndocc
         do avirt = 1, nvirt
          a = ndocc+avirt
          ia = (iocc-1)*nvirt+avirt
          write(6,*)"ia",ia
           do jocc = 1, ndocc
            do bvirt = 1, nvirt
             b = bvirt + ndocc
             jb = (jocc-1)*nvirt + bvirt
             write(6,*)"jb", jb
              Ham = 0.0d0
              Ham = 2.0*(WIJKL(iocc,a,jocc,b)-WIJKL(iocc,jocc,a,b))
              write(6,*)"Ham",Ham
              if (ia.eq.jb) Ham =Ham + eig(a)-eig(iocc)
              write(6,*)"Ham = ", Ham
              call ga_put(g_ham,ia,ia,jb,jb,Ham,1)
             end do 
            end do 
          end do
         end do

           write(6,*)"done Hamiltonian"
           call ga_print(g_ham)

       if (.not. ga_destroy(g_ham))
     & call errquit('Could not destroy g_alpha', 0, GA_ERR)

         return
         end        
C-------------------------------------------------------------------------------
C AO TO MO TRANFORMATION
C-------------------------------------------------------------------------------

       subroutine semiemp_ao2mo(rtdb,BB,nbf,WIJKL)
       
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

       integer rtdb, nbf, nintgl
       integer ilabi,llabi,ilabj,llabj,ilabk,llabk,ilabl,llabl
       integer imoi, lmoi, imoj, lmoj, imok, lmok, imol, lmol
       integer iqq, lqq, iro, lro
       integer lwijkl, iwijkl
       integer iwiqrs, lwiqrs, iwijrs, lwijrs, iwijks, lwijks
       double precision BB(nbf,nbf), WIJKL(*)

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
        

          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIQRS',lwiqrs,
     1    iwiqrs))
     &   call errquit('semiemp_energy: cannot allocate WIqrs',0,
     1    MA_ERR)

          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIJRS',lwijrs,
     1    iwijrs))
     &   call errquit('semiemp_energy: cannot allocate WIJRS',0,
     1    MA_ERR)

          if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf*nbf,'WIJKS',lwijks,
     1    iwijks))
     &   call errquit('semiemp_energy: cannot allocate WIJKS',0,
     1    MA_ERR)
          
         call semiemp_mo_trans(NINTGL,INT_MB(iLABI),INT_MB(iLABJ),
     A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ),WIJKL,BB,nbf
     B   ,Dbl_MB(iwiqrs),Dbl_MB(iwijrs),Dbl_mb(iwijks))

        if (.not.ma_pop_stack(lwijks))
     & call errquit('semiemp_energy: cannot pop stack wijks',3, MA_ERR) 
        if (.not.ma_pop_stack(lwijrs))
     & call errquit('semiemp_energy: cannot pop stack wijrs',3, MA_ERR) 
        if (.not.ma_pop_stack(lwiqrs))
     & call errquit('semiemp_energy: cannot pop stack wiqrs',3, MA_ERR)          
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
       subroutine semiemp_mo_trans(NINTGL,LABI,LABJ,LABK,LABL,QQ,WIJKL,
     1 BB,nbf,WIQRS,WIJRS,WIJKS)

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
       integer nbf, LL, p, q, r , s, i, j, k, l
       double precision BB(nbf,nbf), WIJKL(nbf,nbf,nbf,nbf), QQ(*)
       double precision WIQRS(nbf,nbf,nbf,nbf), WIJRS(nbf,nbf,nbf,nbf)
       double precision WIJKS(nbf,nbf,nbf,nbf), V

       do i = 1, nbf
        do LL = 1, nintgl
          p = LABI(LL)
          q = LABJ(LL)
          r = LABK(LL)
          s = LABL(LL)
          V = QQ(LL)
          WIQRS(i,q,r,s) = WIQRS(i,q,r,s)+ BB(i,p)*V
        end do
        do j = 1, nbf
         do LL = 1, nintgl
          p = LABI(LL)
          q = LABJ(LL)
          r = LABK(LL)
          s = LABL(LL)
          WIJRS(i,j,r,s) = WIJRS(i,j,r,s)+ BB(j,q)*WIQRS(i,q,r,s)
         end do
         do k = 1, nbf
           do LL = 1, nintgl
            p = LABI(LL)
            q = LABJ(LL)
            r = LABK(LL)
            s = LABL(LL)
            WIJKS(i,j,k,s) = WIJKS(i,j,r,s)+ BB(k,r)*WIJRS(i,j,r,s)
           end do
           do l = 1, nbf
            do LL = 1, nintgl
             p = LABI(LL)
             q = LABJ(LL)
             r = LABK(LL)
             s = LABL(LL)
             WIJKL(i,j,k,l) = WIJKL(i,j,k,s)+ BB(l,s)*WIJKS(i,j,r,s)
            end do
           end do 
          end do
         end do
        end do 


       return
       end

               

       
