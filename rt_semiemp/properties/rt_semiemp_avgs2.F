      Subroutine rt_semiemp_avgs2(rtdb,g_zdensa,g_zdensb,g_densa,g_densb
     1 ,g_scr,g_zscr,n,nel,tt)

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
C     == Parameters ==
      character(*), parameter :: pname =
     $     "rt_semiemp_avgs2: "

      integer g_zdensa, g_zdensb, g_zscr
      integer g_densa, g_densb, g_scr
      integer  nel,n
      double precision tt
      double precision pptrace, avgs2, mult, ga_trace_diag
      
      call ga_zero(g_densa)
      call ga_zero(g_densb)
      call ga_zero(g_scr)  
c     write(6,*)"pptrace,avgs2",pptrace,avgs2 
      if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $       call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)

      call convert_z2d(g_zdensa,g_densa,g_scr)
      call convert_z2d(g_zdensb,g_densb,g_scr)
      
      call ga_zero(g_scr)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_densa,g_densb,0d0,
     1 g_scr)

       pptrace=ga_trace_diag(g_scr)
c      write(6,*)"mult,nel,pptrace",mult,nel,pptrace
       avgs2 = ((mult-1)/2)**2+dble(nel)/2-pptrace

      if (ga_nodeid().eq.0) write(6,*) "Time", tt, "<S2>", avgs2

      end subroutine

