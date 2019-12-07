      logical function semiemp_gradient(rtdb)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "rtdb.fh"
#include "bas.fh"
#include "geom.fh"
#include "stdio.fh"
#include "msgids.fh"
#include "sym.fh"
#include "util.fh"
#include "dftps.fh"


       integer rtdb
       logical semiemp_energy
       external semiemp_energy
       logical status
       
      semiemp_gradient = .true.

      status = semiemp_energy(rtdb)
      return
      end

















