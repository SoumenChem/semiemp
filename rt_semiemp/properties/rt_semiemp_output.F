
      subroutine rt_semiemp_print_initialization (params)
      implicit none

#include "stdio.fh"
#include "global.fh"
#include "rt_tddft.fh"
#include "geom.fh"


C     == Inputs ==
      type(rt_params_t), intent(in) :: params


C     == Parameters ==
      character(len=*), parameter :: pname =
     $     "rt_semiemp_print_initialization: "


C     == Variables ==
      integer me
      integer active_geom_indx


      me = ga_nodeid ()

C
C     Print brief summary to stdout.
C


      if (me.eq.0) then
         write (luout, *)
         write (luout, "(a)")
c    $        "=============== Geometry Initialization ==============="

c        write (luout, "(2x,a)") 'Active geometry: "'//
c    $        trim (params%geom(active_geom_indx)%name)//'"'

         if (params%nbf_ao .eq. params%nbf_mo) then
            write (luout, "(2x,a)") "No linear dependencies"
         elseif (params%nbf_ao .gt. params%nbf_mo) then
            write (luout, "(2x,a,i0,a)") "Found ",
     $           params%nbf_ao - params%nbf_mo, " linear dependencies."
         else
            call errquit (pname//"found nbf_mo > nbf_ao ?!", 0, 0)
         endif

         write (luout, "(2x,a,i0)") "Number of atomic orbitals    = ",
     $        params%nbf_ao
         write (luout, "(2x,a,i0)") "Number of molecular orbitals = ",
     $        params%nbf_mo

         write (luout, "(a)")
     $        "======================================================="
         write (luout, *)
         call util_flush (luout)
      endif


      end subroutine









