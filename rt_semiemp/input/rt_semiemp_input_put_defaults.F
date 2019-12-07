      subroutine rt_semiemp_input_put_defaults (rtdb)
      implicit none

#include "rt_tddft.fh"
#include "errquit.fh"
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "stdio.fh"

      
C     == Inputs ==
      integer, intent(in) :: rtdb


C     == Parameters ==
      character(len=*), parameter ::
     $     pname ="rt_semiemp_input_put_defaults: "
      double precision, parameter :: tol_zero = 1d-8
      double precision, parameter :: tol_interpol = 1d-7  ! 10% default schwartz screening
      double precision, parameter :: tol_series = 1d-10
      integer, parameter          :: terms_series = 3
      integer, parameter          :: terms_interpol = 1

      integer, parameter          :: prop_method = 3 ! 2nd order magnus w/ interpolation
      integer, parameter          :: exp_method = 1  ! pseries

C      integer, parameter          :: nsnapshots = 0
      integer, parameter          :: nrestarts = 100
      integer, parameter          :: nchecks = 10

      integer, parameter          :: nprints =  999999999 ! dummy # when using "*" as # prints
      integer, parameter          :: checklvl = 1   ! default to base checks only

      logical, parameter          :: prof = .false.
      logical, parameter          :: noprop = .false.
      logical, parameter          :: static = .false.
      logical, parameter          :: dplot_do = .false. ! do postproc
      logical, parameter          :: subgs = .false.  !sub ground state dens w/ dplot (works?)
      logical, parameter          :: viz = .false.
      logical, parameter          :: nodisk = .false.   !XXX DISABLED NOW--ie, just use DFT settings
      logical, parameter          :: matrix_checks = .false. 
      logical, parameter          :: use_dmat = .false. 
      logical, parameter          :: mocap = .false. 


      character(len=*), parameter :: tag = "<rt_semiemp>: "

C
C     Default output tag
C
      if (.not. rtdb_cput (rtdb, "rt_semiemp:tag",
     $     1, tag)) call errquit (pname//
     $     "failed to put 'tag' target into rtdb",0,0)



C
C     Print options by default: dipole, field, energy, convergence, s2
C     
      if (.not.rtdb_put(rtdb,'rt_semiemp:print_dipole',
     $     mt_log,1,.true.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_quad',
     $     mt_log,1,.false.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_field',
     $     mt_log,1,.true.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_moocc',
     $     mt_log,1,.false.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_energy',
     $     mt_log,1,.true.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_cputime',
     $     mt_log,1,.false.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_charge',
     $     mt_log,1,.false.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_convergence',
     $     mt_log,1,.true.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:print_s2',
     $     mt_log,1,.true.))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)




C
C     Default tolerances
C
      if (.not.rtdb_put(rtdb,'rt_semiemp:tol_zero',
     $     mt_dbl,1,tol_zero))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:tol_series',
     $     mt_dbl,1,tol_series))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:tol_interpol',
     $     mt_dbl,1,tol_interpol))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:terms_series',
     $     mt_int,1,terms_series))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:terms_interpol',
     $     mt_int,1,terms_interpol))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Default propagator: magnus, interpol2, pseries 
C
c     if (.not.rtdb_put(rtdb,'rt_semiemp:prop_method',
c    $     mt_int,1,prop_method))
c    $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      
c     if (.not.rtdb_put(rtdb,'rt_semiemp:exp_method',
c    $     mt_int,1,exp_method))
c    $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Number of prints, snapshots, etc
C
c$$$      if (.not.rtdb_put(rtdb,'rt_tddft:nsnapshots',
c$$$     $     mt_int,1,nsnapshots))
c$$$     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:nrestarts',
     $     mt_int,1,nrestarts))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:nchecks',mt_int,1,nchecks))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:nprints',mt_int,1,nprints))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:checklvl',
     $     mt_int,1,checklvl))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Vizualization
C
      if (.not.rtdb_put(rtdb,'rt_semiemp:viz_subgs', mt_log,1,subgs))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:viz_dplot', mt_log,1,dplot_do))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:viz_active', mt_log,1,viz))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Misc flags
C      
      if (.not.rtdb_put(rtdb,'rt_semiemp:prof',mt_log,1,prof))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:noprop',mt_log,1,noprop))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      
      if (.not.rtdb_put(rtdb,'rt_semiemp:static',mt_log,1,static))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      
c     if (.not.rtdb_put(rtdb,'rt_semiemp:nodisk', mt_log,1,nodisk))
c    $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:matrix_checks',
     $     mt_log,1,matrix_checks))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:use_dmat',mt_log,1,use_dmat))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:mocap_active',mt_log,1,mocap))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


      end subroutine


c $Id: rt_tddft_input_put_defaults.F 28234 2016-05-21 05:33:52Z klopata $
