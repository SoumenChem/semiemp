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

       integer rtdb, nb
       double precision FA(*)
       integer lBB, iBB, iW, lW, norb
       integer numorb(10)
       double precision occuorb(10)
       integer ehinfo(5)
       integer soccup(2), i, j, iorb
       integer ocor, uocor  
       integer numfo

c       write(6,*)"i am at elh_state"
        NRR = NB*(NB+1)/2

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

       if (.not.rtdb_get(rtdb,'rt_semiemp:numorb',mt_int,10,numorb))
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
     1 call errquit('eh_state: failed to get ehinfo',0,RTDB_ERR)


      if (.not.ma_pop_stack(lW))
     &   call errquit('el-h state: cannot pop stack W',3, MA_ERR)

      if (.not.ma_pop_stack(lBB))
     &   call errquit('el-h_state: cannot pop stack BB',3, MA_ERR)

       end 
