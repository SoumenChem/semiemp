        subroutine semiemp_init(rtdb)!,scftype,nel,infa,mult,itmax,scftol,
c    1  III,IA,gfac,inttyp)
c
            implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
#include "inp.fh"
#include "geom.fh"

c     
                                      
      integer rtdb
      character*20 test
      character*20 scftype,apx
      logical rhf,uhf
      logical cndo,indo
      integer maxiter,inttyp,gfac
      integer nel,geom
      double precision mult, infa(6)
      integer charge
      double precision scftol
      integer III,IA

c     if (ga_nodeid().eq.0) 
c    &        write(luout,*) "i am in semiemp_init" 
c
c      if (.not.geom_create(geom,'geometry'))
c    1  call errquit('semiemp_energy: geom_create failed',0, GEOM_ERR)
c      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
c    1  call errquit('semiemp_energy: no geometry found',0, RTDB_ERR)

       if (.not.rtdb_cget(rtdb,'semiemp:scftype',1,scftype)) then
       scftype='rhf'
        if (.not.rtdb_cput(rtdb,'semiemp:scftype',1,scftype)) 
     1 call errquit('semiemp_init: failed to write scftype',0, RTDB_ERR)
      endif
      if (.not.rtdb_get(rtdb,'semiemp:nel',mt_int,1,nel)) then
       write(luout,*)"No valence electron data in input"
       endif
      if (.not.rtdb_get(rtdb,'semiemp:infa',mt_dbl,6,infa)) then
          infa(1) = 1.0
          infa(2) = 1.0
          infa(3) = 1.0
          infa(4) = 1.0
          infa(5) = 1.0
          infa(6) = 1.0
       if (.not.rtdb_put(rtdb,'semiemp:infa',mt_dbl,6,infa)) 
     1  call errquit('semiemp_init: failed to write infa',0, RTDB_ERR)
       endif
c       if (.not.rtdb_get(rtdb,'semiemp:charge',mt_int,1,charge)) then
c        charge=0
c        if (.not.rtdb_put(rtdb,'semiemp:charge',mt_int,1,charge))
c    1  call errquit('semiemp_init: failed to write charge',0, RTDB_ERR)
c      endif
        if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult)) then
         mult=1.0
         if (.not.rtdb_put(rtdb,'semiemp:mult',mt_dbl,1,mult))
     1  call errquit('semiemp_dbl: failed to write mult',0, RTDB_ERR)
       endif
        if (.not.rtdb_get(rtdb,'semiemp:maxiter',mt_int,1,maxiter)) then
         maxiter=50
          if (.not.rtdb_put(rtdb,'semiemp:maxiter',mt_int,1,maxiter))
     1 call errquit('semiemp_init: failed to write maxiter',0, RTDB_ERR)
       endif
        if (.not.rtdb_get(rtdb,'semiemp:scftol',mt_dbl,1,scftol)) then
        scftol=1.0D-8
         if (.not.rtdb_put(rtdb,'semiemp:scftol',mt_dbl,1,scftol))
     1  call errquit('semiemp: failed to write scftol',0, RTDB_ERR)
       endif
       if (.not.rtdb_get(rtdb,'semiemp:III',mt_int,1,III)) then
         III=0
          if (.not.rtdb_put(rtdb,'semiemp:III',mt_int,1,III))
     1  call errquit('semiemp: failed to write III',0, RTDB_ERR)
       endif 
       if (.not.rtdb_get(rtdb,'semiemp:IA',mt_int,1,IA)) then
        IA=0
         if (.not.rtdb_put(rtdb,'semiemp:IA',mt_int,1,IA))
     1  call errquit('semiemp: failed to write IA',0, RTDB_ERR)
       endif
       if (.not.rtdb_get(rtdb,'semiemp:gfac',mt_int,1,gfac)) then
         gfac=0
          if (.not.rtdb_put(rtdb,'semiemp:gfac',mt_int,1,gfac))
     1 call errquit('semiemp: failed to write gfac',0, RTDB_ERR)
         endif 
         if (.not.rtdb_cget(rtdb,'semiemp:apx',1,apx)) then
         apx='indo/1'
          if (.not.rtdb_cput(rtdb,'semiemp:apx',1,apx))
     1  call errquit('semiemp: failed to write apx',0, RTDB_ERR)
        endif
         if (.not.rtdb_get(rtdb,'semiemp:inttyp',mt_int,1,inttyp)) then
         inttyp=0
         if (.not.rtdb_put(rtdb,'semiemp:inttyp',mt_int,1,inttyp))
     1  call errquit('semiemp: failed to write inttyp',0, RTDB_ERR)
        endif 
        return
      end
