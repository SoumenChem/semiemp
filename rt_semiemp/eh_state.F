       Subroutine eh_state(rtdb,FA,nb)
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"
#include "inp.fh"

       character(255) fname_restart
       integer, parameter      :: unitno = 7813259    !XXX TERRIBLE:HARDCODED UNIT #               
       integer rtdb, nb
       double precision FA(*)
       integer lBB, iBB, iW, lW, norb
       integer numorb(100)
       double precision occuorb(10)
       integer ehinfo(5)
       integer soccup(2), i, j, iorb
       integer ocor, uocor  
       integer numfo
       logical lmoocc
c       write(6,*)"i am at elh_state"
        NRR = NB*(NB+1)/2

c       write(6,*)"i am in eh_state"
        if (.not.MA_Push_Get(MT_Dbl,nb*nb,'BB',lBB,iBB))
     &   call errquit('eh_state: cannot allocate BB',0,MA_ERR)

        if (.not.MA_Push_Get(MT_Dbl,nb,'W',lW,iW))
     &   call errquit('eh_state: cannot allocate BB',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:coefficient',mt_dbl,nb*nb,
     1 Dbl_MB(iBB)))
     1 call errquit('eh_state: failed to write coefficient',0, RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:occupation',mt_dbl,nb,Dbl_MB(iW)))
     1 call errquit('eh_state: failed to write ocuupation',0,RTDB_ERR)    


       if (.not.rtdb_get(rtdb,'rt_semiemp:norb',mt_int,1,norb))
     1 call errquit('eh_state: failed to get norb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:numorb',mt_int,100,numorb))
     1 call errquit('eh_state: failed to get numorb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:occuorb',mt_dbl,10,occuorb))
     1 call errquit('eh_state: failed to get occuorb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rtsemiemp:ocor',mt_int,1,ocor))
     1 call errquit('eh_state: failed to get ocor',0,RTDB_ERR)

c      I = 1
c      J = 2
       uocor = 0
c      write(6,*)"HERE",norb

       numfo = 0

       do iorb = 1, norb
        I = numorb(iorb)
        if (Dbl_MB(iw+i-1).eq.2.0d0) numfo = numfo +1
        Dbl_MB(iw+i-1) = occuorb(iorb)
c      write(6,*) "occupation", occuorb(iorb)
       end do
      if (.not.rtdb_get (rtdb, "rt_semiemp:moocc",
     $     mt_log, 1, lmoocc))
     $     lmoocc = .false.

      if (lmoocc) then
         call util_file_name ("rt_moocc",
     $        .false., .false., fname_restart)
         open (unitno, status="old", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open: "
     $        //trim(fname_restart), 0, 0)
           do i = 1, nb
            read (unitno, fmt="(2e20.12)") Dbl_MB(iW+i-1)
           enddo
c        enddo
         close (unitno)

         if (ios .ne. 0) call errquit("failed to close: "
     $        //trim(fname_restart), 0, 0)

      end if

       CALL DENSIT(Dbl_MB(iBB),FA,Dbl_MB(iW),NB,NB,NRR)
c     write(6,*)"numfo,numorb",numfo,numorb(1), numorb(numfo) 
      ehinfo(3) = numfo
      ehinfo(1) = numorb(1)
      ehinfo(2) = numorb(numfo)
      ehinfo(4) = numorb(norb)
      ehinfo(5) = numorb(norb) - numorb(numfo) 
c     write(6,*)"ehinfo", ehinfo
       ocor = ocor + uocor


       if (.not.rtdb_put(rtdb,'rtsemiemp:ehinfo',mt_int,5, ehinfo))
     1 call errquit('eh_state: failed to put ehinfo',0,RTDB_ERR)


      if (.not.ma_pop_stack(lW))
     &   call errquit('el-h state: cannot pop stack W',3, MA_ERR)

      if (.not.ma_pop_stack(lBB))
     &   call errquit('el-h_state: cannot pop stack BB',3, MA_ERR)

       end


c---------------------------------------------------------------------------
c   eh state for UHF
c--------------------------------------------------------------------------  
       Subroutine eh_state_u(rtdb,RA,RB,nb)
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"
#include "inp.fh"

       integer rtdb, nb
       double precision RA(*), RB(*)
       integer lAA, iAA, iWA, lWA, naorb
       integer lBB, iBB, iWB, lWB, nborb
       integer numorba(100), numorbb(100)
       double precision occuorba(10), occuorbb(10)
       integer ehinfo(5)
       integer soccup(2), i, j, iorb
       integer ocor, uocor
       integer numfoa, numfob

        NRR= NB*(NB+1)/2
c       write(6,*)"i am in eh_state_u"
          if (.not.MA_Push_Get(MT_Dbl,nb*nb,'AA',lAA,iAA))
     &   call errquit('eh_state: cannot allocate AA',0,MA_ERR)

        if (.not.MA_Push_Get(MT_Dbl,nb,'WA',lWA,iWA))
     &   call errquit('eh_state: cannot allocate WA',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:acoefficient',mt_dbl,nb*nb,
     1 Dbl_MB(iAA)))
     1 call errquit('eh_state: failed to write alpha coefficient',0, 
     1 RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:aoccupation',mt_dbl,nb,Dbl_MB(iWA)
     1 ))
     1 call errquit('eh_state: failed to write alpha ocuupation',0,
     1 RTDB_ERR)

      if (.not.MA_Push_Get(MT_Dbl,nb*nb,'BB',lBB,iBB))
     &   call errquit('eh_state: cannot allocate BB',0,MA_ERR)

        if (.not.MA_Push_Get(MT_Dbl,nb,'WB',lWB,iWB))
     &   call errquit('eh_state: cannot allocate WA',0,MA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:bcoefficient',mt_dbl,nb*nb,
     1 Dbl_MB(iBB)))
     1 call errquit('eh_state: failed to write beta coefficient',0,
     1 RTDB_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:boccupation',mt_dbl,nb,
     1 Dbl_MB(iWB)))
     1 call errquit('eh_state: failed to write beta ocuupation',0,
     1 RTDB_ERR)

      if (.not.rtdb_get(rtdb,'rt_semiemp:naorb',mt_int,1,naorb))
     1 call errquit('eh_state: failed to get naorb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:numorba',mt_int,100,numorba))
     1 call errquit('eh_state: failed to get numorba',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:occuorba',mt_dbl,10,occuorba))
     1 call errquit('eh_state: failed to get occuorba',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'rt_semiemp:nborb',mt_int,1,nborb))
     1 call errquit('eh_state: failed to get nborb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:numorbb',mt_int,100,numorbb))
     1 call errquit('eh_state: failed to get numorbb',0,RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:occuorbb',mt_dbl,10,occuorbb))
     1 call errquit('eh_state: failed to get occuorba',0,RTDB_ERR)


       numfoa =0 
       do iorb = 1, naorb
        I = numorba(iorb)
        if (Dbl_MB(iwa+i-1).eq.1.0d0) numfoa = numfoa +1
        Dbl_MB(iwa+i-1) = occuorba(iorb)
c      write(6,*) "occupation", occuorb(iorb)
       end do
      CALL DENSIT(Dbl_MB(iAA),RA,Dbl_MB(iWA),NB,NB,NRR)

       numfob =0
       do iorb = 1, nborb
        I = numorbb(iorb)
        if (Dbl_MB(iwb+i-1).eq.1.0d0) numfob = numfoa +1
        Dbl_MB(iwb+i-1) = occuorbb(iorb)
c      write(6,*) "occupation", occuorb(iorb)
       end do
      CALL DENSIT(Dbl_MB(iBB),RA,Dbl_MB(iWB),NB,NB,NRR)
      
       if (.not.rtdb_put(rtdb,'rtsemiemp:ehinfo',mt_int,5, ehinfo))
     1 call errquit('eh_state_u: failed to put ehinfo',0,RTDB_ERR)

      if (.not.ma_pop_stack(lWB))
     &   call errquit('el-h state: cannot pop stack WB',3, MA_ERR)

      if (.not.ma_pop_stack(lBB))
     &   call errquit('el-h_state: cannot pop stack BB',3, MA_ERR)

      if (.not.ma_pop_stack(lWA))
     &   call errquit('el-h state: cannot pop stack WA',3, MA_ERR)

      if (.not.ma_pop_stack(lAA))
     &   call errquit('el-h_state: cannot pop stack AA',3, MA_ERR)
      
      end       
