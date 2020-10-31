      subroutine semiemp_integral_save (LL,LABI,LABJ,LABK,LABL,QQ,n)
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

      
C     == Parameters ==
      character(*), parameter :: pname = "semiemp_integral_save: "


C     == Variables ==
      integer n
      character(255) fname_restart
      integer, parameter      :: unitno = 6457251    !XXX TERRIBLE: HARDCODED UNIT #
      integer ios, me, lmat, kmat
      logical ok

      integer i, j, imat
      integer LL
      integer LABI(*),LABJ(*), LABK(*),LABL(*)
      double precision checksum
      double precision QQ(*)
      double precision elapsed
      character(len=50) date


      me = ga_nodeid()
      
C     note: date returns a newline at the end, so we strip it off
      call util_date (date)
      date = date(:inp_strlen(date)-1) 

      if (me.eq.0) then
         call util_file_name ("aoint",
     $        .false., .false., fname_restart)
         
         open (unitno, status="replace", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open/create: "
     $        //trim(fname_restart), 0, 0)
         
         write (unitno, "(a)")        "AO integral file"
         write (unitno, "(a,a)")      "created   ", trim(date)
         write (unitno, "(a,i0)")     "nbf_ao    ", n
         write (unitno, "(a,i0)")     "numint    ", LL

         do i = 1, LL
           write (unitno, fmt="(i5,i5,i5,i5,2e20.10)") 
     1  LABI(i),LABJ(i),LABK(i),LABL(i),QQ(i)
         enddo
         close (unitno)
         
         if (ios .ne. 0) call errquit("failed to close: "
     $        //trim(fname_restart), 0, 0)

      endif


      end subroutine


      subroutine semiemp_integral_load (LL,LABI,LABJ,LABK,LABL,QQ)
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
      integer              :: LL
      double precision     :: QQ(*)
c     integer, intent(in)               :: nmats
      integer LABI(*), LABJ(*), LABK(*), LABL(*)


C     == Parameters ==
      character(*), parameter :: pname = "semiemp_integral_load: "


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
      integer ll_in
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

         call util_file_name ("aoint",
     $        .false., .false., fname_restart)

         open (unitno, status="old", file=fname_restart, iostat=ios)
         if (ios .ne. 0) call errquit("failed to open: "
     $        //trim(fname_restart), 0, 0)


C
C     Read in header info and check
C
         read (unitno, fmt="(a)") str_long_in
         if (trim(str_long_in) .ne. "AO integral file")
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

         read (unitno, fmt="(a10,i10)") str_in, ll_in
         if (str_in .ne. "numint")
     $  call errquit (pname//"failed to find 'number of integrals'",0,0)
         if (ll_in .ne. LL)
     $        call errquit (pname//"inconsistent LL",0,0)


c        read (unitno, fmt="(a10,e20.12)") str_in, checksum_in
c        if (str_in .ne. "checksum")
c    $        call errquit (pname//"failed to find 'checksum'",0,0)

C
C     Read in complex density matrix (or matrices)
C
         do i = 1, ll_in
            read (unitno, fmt="(i5,i5,i5,i5,2e20.10)")
     1       LABI(i),LABJ(i),LABK(i),LABL(i),QQ(i)
         enddo

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


c     if (me.eq.0) then
c        write (luout,*) ""
c        write (luout,"(a)") "========= Integral Read =========="
c        write (luout,"(x,a,a)")      "file name     ",
c    $        trim(fname_restart)
c        write (luout,"(x,a,a)")      "creation time ", date_in
c        write (luout,"(x,a,a)")      "import time   ", date_now
c        write (luout,"(x,a,i0)")     "nmats         ", nmats_in
c        write (luout,"(x,a,i0)")     "nbf_ao        ", n_in
c        write (luout,"(x,a,i0)")     "number of integral ", LL
c        write (luout,"(a)") "========================================"
c        write (luout,*) ""
c     endif

C
C     Overwrite it and tt for output.  g_zdens_ao is already
C     overwritten.
C

c     call ga_sync ()

c     call rt_tddft_print_msg (params, tt, "Loaded restart file")

      end subroutine      

