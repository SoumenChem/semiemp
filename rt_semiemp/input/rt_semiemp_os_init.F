C
C     rt_semiemp_cs_init.F
C
C     Initialize closed shell RT-TDDFT.
C
      subroutine rt_semiemp_os_init (params)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "matutils.fh"
c#include "cdft.fh"
#include "rt_tddft.fh"


C     == In/out ==
      type(rt_params_t), intent(inout) ::  params


C     == Parameters ==
      character(*), parameter :: pname = "rt_tddft_cs_init: "


C     == Variables ==
      integer me
      double precision nuclear_charge


      me = ga_nodeid ()


      params%theory = "OS"


      return
      end
