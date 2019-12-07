      subroutine rt_semiemp_mullpop(rtdb,g_mullpop,tag)
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
      integer g_mullpop
      integer tag
c
c     Print a summary of the MO vectors in the specified range.
c
      integer num_bf
      integer lps, ips, i, j, k_list, l_list
      integer iaopop, laopop, lnumbf, inumbf
      integer icharge, lcharge
      integer iatom, latom
      integer l_aos,i_aos
      integer iformal, lformal
      integer iz, lz
      integer n, k, m, klo, khi, ibuf
      character buf*80
      double precision cur_thresh, term
c 
      integer iBB, lBB, iW, lW, ieig, leig
      integer  nbf, nat
      double precision maxv
c
c  
      if (.not.rtdb_get(rtdb,'semiemp:nbf',mt_int,1,nbf)) 
     1 call errquit('mullpop: failed to write nbf',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'semiemp:natom',mt_int,1,nat))
     1 call errquit('mullpop: failed to write nat',0,RTDB_ERR)
 
          if (.not.MA_Push_Get(MT_Dbl,nbf,'PS',lps,ips))
     &   call errquit('mullpop: cannot allocate PS',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nbf,'ao_pop',laopop,iaopop))
     &   call errquit('mullpop: cannot allocate aopop',0,MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'at_charge',lcharge,icharge))
     &   call errquit('mullpop: cannot allocate charge',0,MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'formal charge',lformal,
     1   iformal))
     &   call errquit('mullpop: cannot allocate charge',0,MA_ERR)
          if (.not.MA_Push_Get(mt_int,nat,'numbf',lnumbf,inumbf))
     &   call errquit('mullpop: cannot allocate numbf',0,
     1    MA_ERR)
        if (.not. ma_push_get(mt_byte,nbf*8,'atom',latom,iatom))
     $        call errquit('mullpop: ma 3 failed', nbf*8, MA_ERR)
         if (.not. ma_push_get(mt_byte,nbf*8,'aos',l_aos,i_aos))
     $        call errquit('mullpop: ma 3 failed', nbf*8, MA_ERR)
          if (.not.MA_Push_Get(mt_dbl,nat,'Z',lZ,iZ))
     &   call errquit('mullpop: cannot allocate Z',0,
     1    MA_ERR)
c
      if (.not.rtdb_get(rtdb,'semiemp:atom',mt_byte,nbf*8,
     1 byte_mb(iatom)))
     1 call errquit('mullpop: failed to write ATOM',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:aso',mt_byte,nbf*8,
     1 byte_mb(i_aos)))
     1 call errquit('mullpop: failed to write AOS',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:numbf',mt_int,nat,int_mb(inumbf)))
     1 call errquit('mullpop: failed to write numbf',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:DK',mt_dbl,nat,dbl_mb(iz)))
     1 call errquit('mullpop: failed to write Z',0,RTDB_ERR)
c 
c if there is no title sent, don't print it.
c
         do i = 1, nbf
c
            call ga_get(g_mullpop, 1, nbf, i, i, dbl_mb(ips), nbf)
           term = 0.0d0

           do j = 1, nbf
             term = term + dbl_mb(ips+j-1)
           end do
           dbl_mb(iaopop+i-1) = term
         end do

        if (ga_nodeid().eq.0.and.tag.eq.0) write(6,9030)
         num_bf = 0
         do i = 1, nat
          do j = 1, int_mb(inumbf+i-1)
           num_bf = num_bf + 1
          if (ga_nodeid().eq.0.and.tag.eq.0)
     1    WRITE (6,9020) num_bf,(byte_mb(iatom+(num_bf-1)*8+m),m=0,7),i,
     1  (byte_mb(i_aos+(num_bf-1)*8+m),m=0,7), dbl_mb(iaopop+num_bf-1)
          end do
         end do

        if (ga_nodeid().eq.0.and.tag.eq.0)   WRITE (6,9000)
         num_bf = 0
         do i = 1, nat
          term = 0.0d0
          do j = 1, int_mb(inumbf+i-1)
           num_bf = num_bf + 1
           term = term + dbl_mb(iaopop+num_bf-1)
           end do
           dbl_mb(icharge+i-1) = term
           dbl_mb(iformal+i-1) = dbl_mb(iZ+i-1) - term
        if (ga_nodeid().eq.0.and.tag.eq.0)
     1  WRITE (6,9010) I,(byte_mb(iatom+(num_bf-1)*8+m),m=0,7), 
     1  dbl_mb(iformal+i-1)  ,dbl_mb(icharge+i-1)
          end do

      if (.not.rtdb_put(rtdb,'mullpop:aopop',mt_dbl,nbf,dbl_mb(iaopop)))
     1 call errquit('mullpop: failed to write aopop',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'mullpop:charge',mt_dbl,nat,
     1 dbl_mb(icharge)))
     1 call errquit('mullpop: failed to write charge',0,RTDB_ERR)
      if (.not.rtdb_put(rtdb,'mullpop:formal',mt_dbl,nat,
     1 dbl_mb(iformal)))
     1 call errquit('mullpop: failed to write formal',0,RTDB_ERR)

c
c     Identify significant coefficients and sort by size
c
 9000 FORMAT (1X,/,9X,'ATOM',3X,'TYPE',8X,'FORMAL',5X,'TOTAL'/28X,
     x       'CHARGE',5X,'CHARGE',/)
 9010 FORMAT (9X,I4,3X,8A1,2F11.3)
 9020 FORMAT (1X,I4,8A1,I4,'->',8A1,F9.3)
 9030 FORMAT (10X,'ORBITAL POPULATIONS',//,8X,'MULLIKEN',/)

c        if (.not. ma_chop_stack(l_aos)) call errquit
c    $        ('bas_vec_info: ma pop?', 0, BASIS_ERR)

      if (.not.ma_pop_stack(lz))
     &   call errquit('mullpop: cannot pop stack Z',3, MA_ERR)
      if (.not.ma_pop_stack(l_aos))
     &   call errquit('mullpop: cannot pop stack AOS',3, MA_ERR)
      if (.not.ma_pop_stack(latom))
     &   call errquit('mullpop: cannot pop stack atom',3, MA_ERR)
      if (.not.ma_pop_stack(lnumbf))
     &   call errquit('mullpop: cannot pop stack numbf',3, MA_ERR)
      if (.not.ma_pop_stack(lformal))
     &   call errquit('mullpop: cannot pop stack formal',3, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('mullpop: cannot pop stack charge',3, MA_ERR)
      if (.not.ma_pop_stack(laopop))
     &   call errquit('mullpop: cannot pop stack aopop',3, MA_ERR)
      if (.not.ma_pop_stack(lps))
     &   call errquit('mullpop: cannot pop stack PS',3, MA_ERR)

c
      end



