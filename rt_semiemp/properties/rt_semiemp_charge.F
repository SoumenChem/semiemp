      Subroutine rt_semiemp_charge(rtdb,g_zdens_ao,g_s,g_ft,nb,tt)

      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
C#include "cdft.fh"
#include "matutils.fh"
#include "zora.fh"
#include "rt_tddft.fh"

      integer rtdb
      integer g_zdens_ao, g_temp
      integer g_densre, g_densim, g_s, g_scr, g_ft
      integer NEL,NA,ISW1,NB,NRR
      integer i,me,it,nrestarts
      integer frag_charge(2)
      integer iformal, lformal
      integer nat
      double precision tt, pstrace, ga_trace_diag
      double precision term
      logical status

       me = ga_nodeid()
       
c       write(6,*)"i am in charge"!,(KAN(I),I=1,na)

      if (.not. ga_create(mt_dbl, nb, nb, "Re[Dens]", 0, 0, g_densre))
     $     call errquit ("couldnt create g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_densim, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_scr, "temp"))
     $     call errquit ("couldnt duplicate g_scr", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_temp, "temp1"))
     $     call errquit ("couldnt duplicate g_scr", 0, GA_ERR)

      if (.not.rtdb_get(rtdb,'semiemp:natom',mt_int,1,nat))
     1 call errquit('rt_semiemp_charge: failed to write nat',0,RTDB_ERR)

       if (.not.MA_Push_Get(MT_Dbl,nat,'formal charge',lformal,iformal))
     & call errquit('rt_semiemp_charge:cannot allocate formal',0,MA_ERR)

      call convert_z2d(g_zdens_ao,g_densre,g_densim)

c     call ga_print(g_densre)
      call ga_dgemm('N','N',nb,nb,nb,z1,g_densre,g_ft,z0,
     1 g_scr)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_ft,g_scr,z0,
     1 g_temp)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_s,g_temp,z0,
     1 g_scr)

       pstrace = ga_trace_diag(g_scr)
       if (me.eq.0) write(6,*)"pstrace",pstrace

      call ga_elem_multiply(g_temp,g_s,g_scr)
      call rt_semiemp_mullpop(rtdb,g_scr,1)

      if (.not.rtdb_get(rtdb,'rt_semiemp:frag_charge',mt_int,2,
     1  frag_charge))
     1 call errquit('rt_semiemp_charge: failed to read frag',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'mullpop:formal',mt_dbl,nat,
     1 dbl_mb(iformal)))
     1 call errquit('rt_semiemp_charge:failed to read formal',0,
     1 RTDB_ERR)

       term = 0.0d0
      do i = frag_charge(1), frag_charge(2)
       term = term + dbl_mb(iformal+i-1)
      end do

      call convert_d2z(z1,g_densre,z1,g_densim,g_zdens_ao)


      if (ga_nodeid().eq.0) 
     1 write (luout, "(a, 3x, 1f11.5, 1es22.12e3, 5x, a)")
     $        "fragment charge", tt, term, "system"

         call util_flush (luout)
      status=.true.
      status = status.and.ga_destroy(g_temp)
      status = status.and.ga_destroy(g_scr)
      status=status.and.ga_destroy(g_densim)
      status=status.and.ga_destroy(g_densre)

      if (.not.ma_pop_stack(lformal))
     & call errquit('rt_semiemp_charge:cannot pop stack formal',3,
     1 MA_ERR)


      end subroutine


