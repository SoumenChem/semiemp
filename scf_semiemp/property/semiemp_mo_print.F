      subroutine semiemp_mo_print(rtdb,atnum,nbf)
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
      double precision thresh   ! [input] Print coeffs with absval >= thresh
      integer g_vecs
      logical oevals            ! [input] If true print eigenvalues
      logical oocc              ! [input] If true print occupations
c
c     Print a summary of the MO vectors in the specified range.
c
      integer atnum(*)
      integer l_vecs, k_vecs, i, j, k_list, l_list
      integer iatom, latom
      integer l_aos,k_aos
      integer n, k, m, klo, khi, ibuf
      character buf*80
      double precision cur_thresh
c 
      integer iBB, lBB, iW, lW, ieig, leig
      integer type, nbf, nmo
      double precision maxv
c
c      
        thresh = 0.025d0
        call ga_sync()
  
        if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'BB',lBB,iBB))
     &   call errquit('mo_print: cannot allocate BB',0,MA_ERR)

        if (.not.MA_Push_Get(MT_Dbl,nbf,'W',lW,iW))
     &   call errquit('mo_print: cannot allocate W',0,MA_ERR)

        if (.not.MA_Push_Get(MT_Dbl,nbf,'EIG',leig,ieig))
     &   call errquit('mo_print: cannot allocate EIG',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:coefficient',mt_dbl,nbf*nbf,
     1 Dbl_MB(iBB)))
     1 call errquit('mo_print: failed to write coefficient',0, RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:occupation',mt_dbl,nbf,
     1 Dbl_MB(iW)))
     1 call errquit('mo_print: failed to write occupation',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:eig',mt_dbl,nbf,Dbl_MB(ieig)))
     1 call errquit('mo_print: failed to write EIG',0,RTDB_ERR)
 
      if (.not. ga_create(mt_dbl, nbf, nbf, 'ground state coefficient',
     &        0, 0, g_vecs))
     1 call errquit('mo_print: failed to write coefficient',0, RTDB_ERR)

            if (ga_nodeid().eq.0) write(LuOut,20)
 20         format(20x,('  Final Molecular Orbital analysis '))
       
       call ga_put(g_vecs,1,nbf,1,nbf,Dbl_MB(iBB),nbf)
c
c
         if (.not. ma_push_get(mt_dbl,nbf,'coefficient',l_vecs,k_vecs))
     $        call errquit('mo_print: ', nbf, MA_ERR)
         if (.not. ma_push_get(mt_int,nbf,'bf number',l_list,k_list))
     $        call errquit('mo_print: ma 2 failed', nbf, MA_ERR)
         if (.not. ma_push_get(mt_byte,nbf*8,'atom',latom,iatom))
     $        call errquit('mo_print: ma 3 failed', nbf, MA_ERR)
         if (.not. ma_push_get(mt_byte,nbf*8,'aos',l_aos,k_aos))
     $        call errquit('mo_print: ma 3 failed', nbf, MA_ERR)

c
      if (.not.rtdb_get(rtdb,'semiemp:atom',mt_byte,nbf*8,
     1 byte_mb(iatom)))
     1 call errquit('mo_print: failed to write ATOMn',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:aso',mt_byte,nbf*8,
     1 byte_mb(k_aos)))
     1 call errquit('mo_print: failed to write AOS',0,RTDB_ERR)

c 
c if there is no title sent, don't print it.
c
 110     continue
         do i = 1, nbf
c
            call ga_get(g_vecs, 1, nbf, i, i, dbl_mb(k_vecs), nbf)
c
c     Identify significant coefficients and sort by size
c
            n = 0
            cur_thresh = thresh
 111        do j = 0, nbf-1
               if (abs(dbl_mb(k_vecs+j)).ge.cur_thresh) then
                  int_mb(k_list + n) = j
                  n = n + 1
               endif
            enddo
            if (n.eq.0) then
c
c     found nothing. get max +10% of lowe
c
               maxv=1d-99
c     max
               do j = 0, nbf-1
                  if (abs(dbl_mb(k_vecs+j)).ge.maxv) then
                     int_mb(k_list) = j
                     maxv=abs(dbl_mb(k_vecs+j))
                  endif
               enddo
c     get within 10% range
               maxv=0.9d0*maxv
               do j = 0, nbf-1
                  if (abs(dbl_mb(k_vecs+j)).ge.maxv) then
                     int_mb(k_list + n) = j
                     n=n+1
                  endif
               enddo
            endif
            do j = 0, n-1
               do k = 0, j
                  if ((abs(dbl_mb(k_vecs+int_mb(k_list+k)))*0.995.lt.
     $                 abs(dbl_mb(k_vecs+int_mb(k_list+j)))).and.
     $                (abs(dbl_mb(k_vecs+int_mb(k_list+k)))*1.005.lt.
     $                 abs(dbl_mb(k_vecs+int_mb(k_list+j))))) then
                     m = int_mb(k_list+j)
                     int_mb(k_list+j) = int_mb(k_list+k)
                     int_mb(k_list+k) = m
                  elseif((abs(dbl_mb(k_vecs+int_mb(k_list+k)))*0.995.lt.
     $                    abs(dbl_mb(k_vecs+int_mb(k_list+j))))
     &                   .and..not.
     $                   (abs(dbl_mb(k_vecs+int_mb(k_list+k)))*1.005.lt.
     $                    abs(dbl_mb(k_vecs+int_mb(k_list+j))))) then
                     if (int_mb(k_list+j).lt.int_mb(k_list+k)) then
                       m = int_mb(k_list+j)
                       int_mb(k_list+j) = int_mb(k_list+k)
                       int_mb(k_list+k) = m
                     endif
                  else
c                    the ordering is OK so do nothing
                  endif
               enddo
            enddo
c
c     Construct optional output line
c
            ibuf = 1
            buf = ' '
            if (ga_nodeid().eq.0)
     1      write(buf(ibuf:),'(''Occ='',1p,d12.6)') Dbl_MB(iW+i-1)
               ibuf = ibuf + 18
            if (ga_nodeid().eq.0)
     1      write(buf(ibuf:),'(''E='',1p,d13.6)') Dbl_MB(ieig+i-1)
               ibuf = ibuf + 17
            if (ga_nodeid().eq.0)
     1      write(LuOut,1) i, buf(1:max(inp_strlen(buf),1))
 1          format(' Vector',i5,2x,a) 
c
c     Output the analysis
c
            
            if (ga_nodeid().eq.0) write(LuOut,22)
 22         format(1x,2('  atn.  Coefficient  Atom+Function  ',5x))
            if (ga_nodeid().eq.0) write(LuOut,23)
 23         format(1x,2(' ----- ------------  ---------------',5x))
            do klo = 0, min(n-1,9), 2
               khi = min(klo+1,n-1)
            if (ga_nodeid().eq.0)   write(LuOut,2) (
     $              atnum(int_mb(k_list+k)+1), 
     $              dbl_mb(k_vecs+int_mb(k_list+k)),
c    $              atnum(int_mb(k_list+k)),
     $              (byte_mb(iatom+int_mb(k_list+k)*8+m),m=0,7),
     $              (byte_mb(k_aos+int_mb(k_list+k)*8+m),m=0,7),
     $              k = klo,khi)
 2             format(1x,2(i5,2x,f12.6,2x,16a1,4x))
            enddo
            if (ga_nodeid().eq.0) write(LuOut,*)
         enddo
         call util_flush(LuOut)
c
      if (.not.ma_pop_stack(l_aos))
     &   call errquit('mo_print: cannot pop stack aos',3, MA_ERR)

      if (.not.ma_pop_stack(latom))
     &   call errquit('mo_print: cannot pop stack atom',3, MA_ERR)

      if (.not.ma_pop_stack(l_list))
     &   call errquit('mo_print: cannot pop stack list',3, MA_ERR)

      if (.not.ma_pop_stack(l_vecs))
     &   call errquit('mo_print: cannot pop stack vecs',3, MA_ERR)

      if (.not.ma_pop_stack(leig))
     &   call errquit('mo_print: cannot pop stack eig',3, MA_ERR)

      if (.not.ma_pop_stack(lW))
     &   call errquit('mo_print: cannot pop stack W',3, MA_ERR)

      if (.not.ma_pop_stack(lBB))
     &   call errquit('mo_print: cannot pop stack BB',3, MA_ERR)

c
      end



