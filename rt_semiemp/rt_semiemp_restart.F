C
C     rt_semiemp_restart.F
C
C     Routines for dumping and reading in restart point info.
C
C
      subroutine rt_semiemp_restart_save (it,n,tt,g_zdens_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "inp.fh"
#include "rt_tddft.fh"

C     == Inputs ==
c     type(rt_params_t), intent(in)     :: params
      integer               :: it
      double precision      :: tt
      integer              :: nmats
      integer               :: g_zdens_ao

      
C     == Parameters ==
      character(*), parameter :: pname = "rt_semiemp_restart_save: "


C     == Variables ==
      integer n
      character(255) fname_restart
      integer, parameter      :: unitno = 6457251    !XXX TERRIBLE: HARDCODED UNIT #
      integer ios, me, lmat, kmat
      logical ok

      integer i, j, imat
      double complex zval
      double precision checksum
      double precision elapsed
      character(len=50) date


      me = ga_nodeid()
      
c     if (params%prof) call prof_start (elapsed)

c     if ((nmats .ne. 1).and.(nmats .ne. 2))
c    $     call errquit(pname//
c    $     "only valid fornmats = 1 or nmats = 2", 0, 0)

c     checksum = 0d0
c     do imat = 1, nmats
c        checksum = checksum + zmat_checksum (g_zdens_ao(imat))
c     enddo

c     n = params%ns_ao   !alias for clarity

C     note: date returns a newline at the end, so we strip it off
      call util_date (date)
      date = date(:inp_strlen(date)-1) 

      if (me.eq.0) then
         call util_file_name ("rt_restart",
     $        .false., .false., fname_restart)
         
         open (unitno, status="replace", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open/create: "
     $        //trim(fname_restart), 0, 0)
         
         write (unitno, "(a)")        "RT-Semiemp restart file"
         write (unitno, "(a,a)")      "created   ", trim(date)
c        write (unitno, "(a,i0)")     "nmats     ", nmats
         write (unitno, "(a,i0)")     "nbf_ao    ", n
         write (unitno, "(a,i0)")     "it        ", it
         write (unitno, "(a,e20.12)") "t         ", tt
c        write (unitno, "(a,e20.12)") "checksum  ", checksum

c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  call ga_get (g_zdens_ao, i, i, j, j, zval, 1)
                  write (unitno, advance="no", fmt="(2e20.10)") zval
               enddo
            enddo
c        enddo
         close (unitno)
         
         if (ios .ne. 0) call errquit("failed to close: "
     $        //trim(fname_restart), 0, 0)

      endif
c     call ga_sync()

c     call rt_tddft_print_msg (params, tt, "Saved restart file")

c     if (params%prof) call prof_end (elapsed, "Saving restart point")

      end subroutine




      subroutine rt_semiemp_restart_load (it,n,tt,g_zdens_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "inp.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      integer              :: it
      double precision     :: tt
c     integer, intent(in)               :: nmats
      integer              :: g_zdens_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_semiemp_restart_load: "


C     == Variables ==
      integer n
      character(255) fname_restart
      integer, parameter      :: unitno = 7813259    !XXX TERRIBLE:HARDCODED UNIT #
      integer ios, me, lmat, kmat
      logical ok

      integer i, j, imat
      double precision val_re, val_im
      double complex zval
      integer nmats_in, n_in, it_in
      double precision tt_in
      character(len=32) str_long_in
      character(len=10) str_in        !note fixed length 10
      character(len=50) date_in
      double precision checksum_in, checksum_calc
      double precision elapsed
      integer msglen
      character(len=50) date_now


      me = ga_nodeid()

c     n = params%ns_ao   !alias for clarity

c     if (params%prof) call prof_start (elapsed)

c     if ((nmats .ne. 1).and.(nmats .ne. 2))
c    $     call errquit(pname//
c    $     "only valid fornmats = 1 or nmats = 2", 0, 0)
      
      if (me.eq.0) then

         call util_file_name ("rt_restart",
     $        .false., .false., fname_restart)

         open (unitno, status="old", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open: "
     $        //trim(fname_restart), 0, 0)


C
C     Read in header info and check
C
         read (unitno, fmt="(a)") str_long_in
         if (trim(str_long_in) .ne. "RT-Semiemp restart file")
     $        call errquit (pname//"failed to find header line",0,0)

         read (unitno, fmt="(a,a)") str_in, str_long_in
         if (trim(str_in) .ne. "created")
     $        call errquit (pname//"failed to find creation date",0,0)
         date_in = str_long_in

c        read (unitno, fmt="(a10,i10)") str_in, nmats_in
c        if (str_in .ne. "nmats")
c    $        call errquit (pname//"failed to find 'nmats'",0,0)
c        if (nmats_in .ne. nmats)
c    $        call errquit (pname//"inconsistent nmats",0,0)

         read (unitno, fmt="(a10,i10)") str_in, n_in
         if (str_in .ne. "nbf_ao")
     $        call errquit (pname//"failed to find 'nbf_ao'",0,0)
         if (n_in .ne. n)
     $        call errquit (pname//"inconsistent nbf_ao",0,0)

         read (unitno, fmt="(a10,i10)") str_in, it_in
         if (str_in .ne. "it")
     $        call errquit (pname//"failed to find 'it'",0,0)

         read (unitno, fmt="(a10,e20.12)") str_in, tt_in
         if (str_in .ne. "t")
     $        call errquit (pname//"failed to find 't'",0,0)

c        read (unitno, fmt="(a10,e20.12)") str_in, checksum_in
c        if (str_in .ne. "checksum")
c    $        call errquit (pname//"failed to find 'checksum'",0,0)


C
C     Read in complex density matrix (or matrices)
C
c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  read (unitno, fmt="(e20.12)",advance="no") val_re
                  read (unitno, fmt="(e20.12)",advance="no") val_im
                  zval = val_re + zi*val_im

                  call ga_put (g_zdens_ao, i, i, j, j, zval, 1)
               enddo
            enddo
c        enddo
         close(unitno)

      endif
c     call ga_sync()

C
C     Broadcast it and tt to all procs
C
c     msglen = MA_sizeof(MT_INT,1,MT_BYTE)
c     call ga_brdcst(100001, it_in, msglen, 0)

c     msglen = MA_sizeof(MT_DBL,1,MT_BYTE)
c     call ga_brdcst(100002, tt_in, msglen, 0)

c     msglen = MA_sizeof(MT_DBL,1,MT_BYTE)
c     call ga_brdcst(100003, checksum_in, msglen, 0)


C     confirm checksum
c     checksum_calc = 0d0
c     do imat = 1, nmats
c        checksum_calc = checksum_calc + zmat_checksum(g_zdens_ao(imat))
c     enddo

c     if (abs(checksum_calc - checksum_in) > 1d-4)
c    $     call errquit(pname//"checksums inconsistent",0,0)

C     note: date returns a newline at the end, so we strip it off
      call util_date (date_now)
      date_now = date_now(:inp_strlen(date_now)-1)


      if (me.eq.0) then
         write (luout,*) ""
         write (luout,"(a)") "========= Restart Information =========="
         write (luout,"(x,a,a)")      "file name     ",
     $        trim(fname_restart)
         write (luout,"(x,a,a)")      "creation time ", date_in
         write (luout,"(x,a,a)")      "import time   ", date_now
c        write (luout,"(x,a,i0)")     "nmats         ", nmats_in
         write (luout,"(x,a,i0)")     "nbf_ao        ", n_in
         write (luout,"(x,a,i0)")     "it            ", it_in
         write (luout,"(x,a,e20.12)") "t             ", tt_in
         write (luout,"(a)") "========================================"
         write (luout,*) ""
      endif

C
C     Overwrite it and tt for output.  g_zdens_ao is already
C     overwritten.
C
      it = it_in
      tt = tt_in

c     call ga_sync ()

      write(6,*)"HERE in load"
c     call rt_tddft_print_msg (params, tt, "Loaded restart file")

      end subroutine
C------------------------------------------------------------------------
C     rt_semiemp_restart.F
C
C     Routines for dumping and reading in restart point info.
C
C
      subroutine rt_semiemp_os_restart_save(it,n,tt,g_zdensa_ao,
     1 g_zdensb_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "inp.fh"
#include "rt_tddft.fh"

C     == Inputs ==
c     type(rt_params_t), intent(in)     :: params
      integer               :: it
      double precision      :: tt
      integer              :: nmats
      integer               :: g_zdensa_ao, g_zdensb_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_semiemp_restart_save: "


C     == Variables ==
      integer n
      character(255) fname_restart
      integer, parameter      :: unitno = 6457251    !XXX TERRIBLE:HARDCODED UNIT #
      integer ios, me, lmat, kmat
      logical ok

      integer i, j, imat
      double complex zval
      double precision checksum
      double precision elapsed
      character(len=50) date


      me = ga_nodeid()
c     write(6,*)"save density"
c     if (params%prof) call prof_start (elapsed)

c     if ((nmats .ne. 1).and.(nmats .ne. 2))
c    $     call errquit(pname//
c    $     "only valid fornmats = 1 or nmats = 2", 0, 0)

c     checksum = 0d0
c     do imat = 1, nmats
c        checksum = checksum + zmat_checksum (g_zdens_ao(imat))
c     enddo

c     n = params%ns_ao   !alias for clarity

C     note: date returns a newline at the end, so we strip it off
      call util_date (date)
      date = date(:inp_strlen(date)-1)

      if (me.eq.0) then
         call util_file_name ("rt_alpha_restart",
     $        .false., .false., fname_restart)

         open (unitno, status="replace", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open/create: "
     $        //trim(fname_restart), 0, 0)

         write (unitno, "(a)")        "RT-Semiemp restart file"
         write (unitno, "(a,a)")      "created   ", trim(date)
c        write (unitno, "(a,i0)")     "nmats     ", nmats
         write (unitno, "(a,i0)")     "nbf_ao    ", n
         write (unitno, "(a,i0)")     "it        ", it
         write (unitno, "(a,e20.12)") "t         ", tt
c        write (unitno, "(a,e20.12)") "checksum  ", checksum

c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  call ga_get (g_zdensa_ao, i, i, j, j, zval, 1)
                  write (unitno, advance="no", fmt="(2e20.10)") zval
               enddo
            enddo
c        enddo
         close (unitno)

         if (ios .ne. 0) call errquit("failed to close: "
     $        //trim(fname_restart), 0, 0)

      endif
c     call ga_sync()
      if (me.eq.0) then
         call util_file_name ("rt_beta_restart",
     $        .false., .false., fname_restart)

         open (unitno, status="replace", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open/create: "
     $        //trim(fname_restart), 0, 0)

         write (unitno, "(a)")        "RT-Semiemp restart file"
         write (unitno, "(a,a)")      "created   ", trim(date)
c        write (unitno, "(a,i0)")     "nmats     ", nmats
         write (unitno, "(a,i0)")     "nbf_ao    ", n
         write (unitno, "(a,i0)")     "it        ", it
         write (unitno, "(a,e20.12)") "t         ", tt
c        write (unitno, "(a,e20.12)") "checksum  ", checksum

c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  call ga_get (g_zdensb_ao, i, i, j, j, zval, 1)
                  write (unitno, advance="no", fmt="(2e20.10)") zval
               enddo
            enddo
c        enddo
         close (unitno)

         if (ios .ne. 0) call errquit("failed to close: "
     $        //trim(fname_restart), 0, 0)

      endif
c     call rt_tddft_print_msg (params, tt, "Saved restart file")

c     if (params%prof) call prof_end (elapsed, "Saving restart point")

      end subroutine

C-----------------------------------------------------------------------------
      subroutine rt_semiemp_os_restart_load (it,n,tt,g_zdensa_ao,
     1 g_zdensb_ao)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "rtdb.fh"
#include "util.fh"
#include "cdft.fh"
#include "inp.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

C     == Inputs ==
      integer              :: it
      double precision     :: tt
c     integer, intent(in)               :: nmats
      integer              :: g_zdensa_ao, g_zdensb_ao


C     == Parameters ==
      character(*), parameter :: pname = "rt_semiemp_restart_load: "


C     == Variables ==
      integer n
      character(255) fname_restart
      integer, parameter      :: unitno = 7813259    !XXXTERRIBLE:HARDCODED UNIT #
      integer ios, me, lmat, kmat
      logical ok

      integer i, j, imat
      double precision val_re, val_im
      double complex zval
      integer nmats_in, n_in, it_in
      double precision tt_in
      character(len=32) str_long_in
      character(len=10) str_in        !note fixed length 10
      character(len=50) date_in
      double precision checksum_in, checksum_calc
      double precision elapsed
      integer msglen
      character(len=50) date_now


      me = ga_nodeid()

c     n = params%ns_ao   !alias for clarity

c     if (params%prof) call prof_start (elapsed)

c     if ((nmats .ne. 1).and.(nmats .ne. 2))
c    $     call errquit(pname//
c    $     "only valid fornmats = 1 or nmats = 2", 0, 0)

      if (me.eq.0) then

         call util_file_name ("rt_alpha_restart",
     $        .false., .false., fname_restart)

         open (unitno, status="old", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open: "
     $        //trim(fname_restart), 0, 0)


C
C     Read in header info and check
C
         read (unitno, fmt="(a)") str_long_in
         if (trim(str_long_in) .ne. "RT-Semiemp restart file")
     $        call errquit (pname//"failed to find header line",0,0)

         read (unitno, fmt="(a,a)") str_in, str_long_in
         if (trim(str_in) .ne. "created")
     $        call errquit (pname//"failed to find creation date",0,0)
         date_in = str_long_in

c        read (unitno, fmt="(a10,i10)") str_in, nmats_in
c        if (str_in .ne. "nmats")
c    $        call errquit (pname//"failed to find 'nmats'",0,0)
c        if (nmats_in .ne. nmats)
c    $        call errquit (pname//"inconsistent nmats",0,0)

         read (unitno, fmt="(a10,i10)") str_in, n_in
         if (str_in .ne. "nbf_ao")
     $        call errquit (pname//"failed to find 'nbf_ao'",0,0)
         if (n_in .ne. n)
     $        call errquit (pname//"inconsistent nbf_ao",0,0)

         read (unitno, fmt="(a10,i10)") str_in, it_in
         if (str_in .ne. "it")
     $        call errquit (pname//"failed to find 'it'",0,0)

         read (unitno, fmt="(a10,e20.12)") str_in, tt_in
         if (str_in .ne. "t")
     $        call errquit (pname//"failed to find 't'",0,0)

c        read (unitno, fmt="(a10,e20.12)") str_in, checksum_in
c        if (str_in .ne. "checksum")
c    $        call errquit (pname//"failed to find 'checksum'",0,0)


C
C     Read in complex density matrix (or matrices)
C
c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  read (unitno, fmt="(e20.12)",advance="no") val_re
                  read (unitno, fmt="(e20.12)",advance="no") val_im
                  zval = val_re + zi*val_im

                  call ga_put (g_zdensa_ao, i, i, j, j, zval, 1)
               enddo
            enddo
c        enddo
         close(unitno)
      endif 

      if (me.eq.0) then

         call util_file_name ("rt_beta_restart",
     $        .false., .false., fname_restart)

         open (unitno, status="old", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open: "
     $        //trim(fname_restart), 0, 0)


C
C     Read in header info and check
C
         read (unitno, fmt="(a)") str_long_in
         if (trim(str_long_in) .ne. "RT-Semiemp restart file")
     $        call errquit (pname//"failed to find header line",0,0)

         read (unitno, fmt="(a,a)") str_in, str_long_in
         if (trim(str_in) .ne. "created")
     $        call errquit (pname//"failed to find creation date",0,0)
         date_in = str_long_in

c        read (unitno, fmt="(a10,i10)") str_in, nmats_in
c        if (str_in .ne. "nmats")
c    $        call errquit (pname//"failed to find 'nmats'",0,0)
c        if (nmats_in .ne. nmats)
c    $        call errquit (pname//"inconsistent nmats",0,0)

         read (unitno, fmt="(a10,i10)") str_in, n_in
         if (str_in .ne. "nbf_ao")
     $        call errquit (pname//"failed to find 'nbf_ao'",0,0)
         if (n_in .ne. n)
     $        call errquit (pname//"inconsistent nbf_ao",0,0)

         read (unitno, fmt="(a10,i10)") str_in, it_in
         if (str_in .ne. "it")
     $        call errquit (pname//"failed to find 'it'",0,0)

         read (unitno, fmt="(a10,e20.12)") str_in, tt_in
         if (str_in .ne. "t")
     $        call errquit (pname//"failed to find 't'",0,0)

c        read (unitno, fmt="(a10,e20.12)") str_in, checksum_in
c        if (str_in .ne. "checksum")
c    $        call errquit (pname//"failed to find 'checksum'",0,0)


C
C     Read in complex density matrix (or matrices)
C
c        do imat = 1, nmats
            do i = 1, n
               do j = 1, n
                  read (unitno, fmt="(e20.12)",advance="no") val_re
                  read (unitno, fmt="(e20.12)",advance="no") val_im
                  zval = val_re + zi*val_im

                  call ga_put (g_zdensb_ao, i, i, j, j, zval, 1)
               enddo
            enddo
c        enddo
         close(unitno)

      endif
c     call ga_sync()

C
C     Broadcast it and tt to all procs
C
c     msglen = MA_sizeof(MT_INT,1,MT_BYTE)
c     call ga_brdcst(100001, it_in, msglen, 0)

c     msglen = MA_sizeof(MT_DBL,1,MT_BYTE)
c     call ga_brdcst(100002, tt_in, msglen, 0)

c     msglen = MA_sizeof(MT_DBL,1,MT_BYTE)
c     call ga_brdcst(100003, checksum_in, msglen, 0)


C     confirm checksum
c     checksum_calc = 0d0
c     do imat = 1, nmats
c        checksum_calc = checksum_calc + zmat_checksum(g_zdens_ao(imat))
c     enddo

c     if (abs(checksum_calc - checksum_in) > 1d-4)
c    $     call errquit(pname//"checksums inconsistent",0,0)

C     note: date returns a newline at the end, so we strip it off
      call util_date (date_now)
      date_now = date_now(:inp_strlen(date_now)-1)


      if (me.eq.0) then
         write (luout,*) ""
         write (luout,"(a)") "========= Restart Information =========="
         write (luout,"(x,a,a)")      "file name     ",
     $        trim(fname_restart)
         write (luout,"(x,a,a)")      "creation time ", date_in
         write (luout,"(x,a,a)")      "import time   ", date_now
c        write (luout,"(x,a,i0)")     "nmats         ", nmats_in
         write (luout,"(x,a,i0)")     "nbf_ao        ", n_in
         write (luout,"(x,a,i0)")     "it            ", it_in
         write (luout,"(x,a,e20.12)") "t             ", tt_in
         write (luout,"(a)") "========================================"
         write (luout,*) ""
      endif

C
C     Overwrite it and tt for output.  g_zdens_ao is already
C     overwritten.
C
      it = it_in
      tt = tt_in

c     call ga_sync ()

      write(6,*)"HERE in load"
c     call rt_tddft_print_msg (params, tt, "Loaded restart file")

      end subroutine

