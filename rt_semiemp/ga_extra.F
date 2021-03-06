C---------------------------------------------------------------
C    Lower traingular matrix to square matrix
C----------------------------------------------------------------
      Subroutine ga_sqmat(g_a,p,n,itag)
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
#include "inp.fh"
      integer n, l, ilo, ihi, jlo, jhi, numi, numj
      integer g_a, iT, lT
      integer me, istart, iend
      integer m0, m
      double precision p(*)


        if (.not.MA_Push_Get(MT_Dbl,n*n,'T',lT,iT))
     &   call errquit('ga_sqmat: cannot allocate E',0,
     1    MA_ERR)

      if (ITAG.eq.1) then
      call sym(Dbl_MB(iT),P,n)
       else
      call antisym(Dbl_MB(iT),P,n)
      end if


      call ga_put(g_a,1,n,1,n,Dbl_MB(iT),n)

      if (.not.ma_pop_stack(lT))
     &   call errquit('ga_sqmat: cannot pop stack T',3, MA_ERR)

      return
      end
C---------------------------------------------------------------
C    Squre matrix to lower traingular matrix 
C----------------------------------------------------------------
      Subroutine ga_trimat(g_a,p,n)
#include "errquit.fh"
#include "mafdecls.fh"
#include "tcgmsg.fh"
#include "global.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "sym.fh"
#include "util.fh"
#include "msgids.fh"
#include "stdio.fh"
#include "dra.fh"
#include "inp.fh"
      integer n, l
      integer g_a, it, lt
      double precision p(*)

        if (.not.MA_Push_Get(MT_Dbl,n*n,'T',lT,iT))
     &   call errquit('ga_trimat: cannot allocate E',0,
     1    MA_ERR)


      
      call ga_get(g_a, 1, n, 1, n, dbl_mb(it),n)

      call triangular(p,dbl_mb(it),n)

      if (.not.ma_pop_stack(lT))
     &   call errquit('ga_trimat: cannot pop stack T',3, MA_ERR)

      return
      end
C--------------------------------------------------------------------
C     Projection from complex AO to MO basis
C--------------------------------------------------------------------
      Subroutine mo_zproj(g_za,g_zb,n,g_zc)
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

      integer g_za, g_zb, g_zc, g_ztemp1
      integer n
      logical ok

c     write(6,*)"i am in mo_zproj"
      call ga_zero(g_zc)
      ok = .true.
      ok=ok.and.ga_create(mt_dcpl,n,n,"complex temp", 0,0, g_ztemp1)

      call ga_zgemm('T','N',n,n,n,z1,g_zb,g_za,z0,g_ztemp1)

      call ga_zgemm('N','N',n,n,n,z1,g_ztemp1,g_zb,z0,g_zc)
      
      ok = .true.
      ok = ok .and. ga_destroy(g_ztemp1)
      return
      end
C-----------------------------------------------------------------------------
C    Projection from  AO to MO basis
C---------------------------------------------------------------------------
      Subroutine mo_proj(g_a,g_b,n,g_c)
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

      integer g_a, g_b, g_c, g_temp1
      integer n
      logical ok

      ok = .true.
      ok=ok.and.ga_create(mt_dbl,n,n,"complex temp", 0,0, g_temp1)
      call ga_zero(g_temp1)
      call ga_zero(g_c)

      call ga_dgemm('T','N',n,n,n,z1,g_b,g_a,z0,g_temp1)

      call ga_dgemm('N','N',n,n,n,z1,g_temp1,g_b,z0,g_c)

      ok = .true.
      ok = ok .and. ga_destroy(g_temp1)
      return
      end
C-------------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------------
      Subroutine rt_activespace_os(rtdb,n)
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

        integer leig, ieig
        integer iocc, iuocc, ocor
        double precision ecut
        integer  orbspace(2), MIM(2)
        integer aorb(3)
        integer i, j, n

             if (.not.MA_Push_Get(MT_Dbl,n,'EIG',leig,ieig))
     &   call errquit('semiemp_energy: cannot allocate EIG',0,
     1    MA_ERR)


       if (.not.rtdb_get(rtdb,'semiemp:eig',mt_dbl,n,Dbl_MB(ieig)))
     1  call errquit('rt_activespace: failed to write eig',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1 continue

       if (.not.rtdb_get(rtdb,'semiemp:MIM',mt_int,2,MIM))
     1 call errquit('semiemp: failed to read MIM',0, RTDB_ERR)


       if (ga_nodeid().eq.0)
     1 write(6,*)"orbitals",orbspace(1), orbspace(2)

        aorb(1) = MIM(1) - orbspace(1)+1
        aorb(2) = orbspace(2) - MIM(1)

        if (.not.rtdb_put(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1  call errquit('rt_activespace: failed to write orbspace',0,
     2   RTDB_ERR)

       if (.not.rtdb_put(rtdb,'rt_semiemp:betaactsize',mt_int,3,aorb))
     1  call errquit('rt_activespace: failed to write betaactsize',0,
     2   RTDB_ERR)

        aorb(1) = MIM(2) - orbspace(1)+1
        aorb(2) = orbspace(2) - MIM(2)

        if (.not.rtdb_put(rtdb,'rt_semiemp:alphaactsize',mt_int,3,aorb))
     1  call errquit('rt_activespace: failed to write alphaactsize',0,
     2   RTDB_ERR)
       
           if (.not.ma_pop_stack(leig))
     &  call errquit('semiemp_energy: cannot pop stack EIG',3, MA_ERR)
      return
      end

C-------------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------------
      Subroutine rt_activespace(rtdb,n)
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

        integer leig, ieig
        integer iocc, iuocc, ocor, nel
        double precision ecut
        integer  orbspace(2)
        integer aorb(3)
        integer i, j, n

             if (.not.MA_Push_Get(MT_Dbl,n,'EIG',leig,ieig))
     &   call errquit('semiemp_energy: cannot allocate EIG',0,
     1    MA_ERR)


       if (.not.rtdb_get(rtdb,'semiemp:eig',mt_dbl,n,Dbl_MB(ieig)))
     1  call errquit('rt_activespace: failed to write eig',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'semiemp:nel',mt_int,1,nel))
     1  call errquit('rt_activespace: failed to write nel',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1 continue

       if (.not.rtdb_get(rtdb,'rtsemiemp:ocor',mt_int,1,ocor))
     1  call errquit('rt_activespace: failed to write ocor',0, RTDB_ERR)

       if (.not.rtdb_get(rtdb,'rt_semiemp:ecut',mt_dbl,1,ecut))
     1  GO TO 20
c    1  call errquit('rt_activespace: failed to write ecut',0, RTDB_ERR)


       if (ga_nodeid().eq.0)  write(6,*)"ecut",ecut
c      write(6,*)"orbitals",orbspace(1), orbspace(2)
       iocc = 0
       iuocc = 0
       ecut = ecut/27.2113961
       aorb(1) = 0
       aorb(2) = 0
       aorb(3) = 0
c      write(6,*)"ecut",ecut

       do i = 1, ocor
        do j = ocor+1, n
        ediff = dbl_mb(ieig+j-1) - dbl_mb(ieig+i-1)
c         write(6,*)"ediff",ediff
          if (ediff.le.ecut) then
c         write(6,*)"ediffi,i,j",ediff,i,j
          if (iocc.eq.0.or.iocc.gt.i)  iocc = i
          if (iuocc.lt.j) iuocc = j
          else
           go to 10
          end if

        end do
  10   continue
       end do

        aorb(1) = ocor - iocc + 1
        aorb(2) = iuocc - ocor
c       write(6,*)"iocc,iuocc",iocc,iuocc
       orbspace(1) = iocc
       orbspace(2) = iuocc
  20   continue
       if (ga_nodeid().eq.0) 
     1 write(6,*)"orbitals",orbspace(1), orbspace(2)

        aorb(1) = ocor - orbspace(1)+1
        aorb(2) = orbspace(2) - ocor
       if (.not.rtdb_put(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1  call errquit('rt_activespace: failed to write orbspace',0,
     2   RTDB_ERR)

       if (.not.rtdb_put(rtdb,'rt_semiemp:actsize',mt_int,3,aorb))
     1  call errquit('rt_activespace: failed to write orbspace',0,
     2   RTDB_ERR)

           if (.not.ma_pop_stack(leig))
     &  call errquit('semiemp_energy: cannot pop stack EIG',3, MA_ERR)
      return
      end
C--------------------------------------------------------------------------------------
C
C--------------------------------------------------------------------------------------
      Subroutine fock_extrapolate(g_b,g_a,g_c)
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

       integer g_a, g_b, g_c

       call ga_add(2*z1,g_a,-z1,g_b,g_c)


      return
      end
C-------------------------------------------------------------------------------------------------------
C
C-----------------------------------------------------------------------------------------------------
      Subroutine fock_intrapolate(g_b,g_a,g_c)
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

       integer g_a, g_b, g_c

       call ga_add(0.5*z1,g_a,0.5*z1,g_b,g_c)


      return
      end

C====================================================================
C
C     Compares two matricies to see if they differ.
C     Works for mixed real and/or complex matricies.
C
      logical function mat_is_same_semiemp (g_a, g_b, tol)
      implicit none

#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"

C     == Inputs ==
      integer, intent(in)          :: g_a, g_b
      double precision, intent(in) :: tol


C     == Variables ==
      integer dtype_a, dtype_b
      integer n1_a, n2_a, n1_b, n2_b
      double precision norm
      logical same
      integer g_tmp
      integer me

C     == Parameters ==
      character(*), parameter      :: pname = "mat_is_same: "
      double complex, parameter    :: zone = (1d0, 0d0)
      double complex, parameter    :: nzone = (-1d0, 0d0)
      
C     == External ==
      double precision, external :: mat_norm

      me = ga_nodeid ()
      
C     Check GAs and get their data types
      call ga_check_handle (g_a, "1st arg of "//pname//"not valid GA")
      call ga_check_handle (g_b, "2nd arg of "//pname//"not valid GA")

      call ga_inquire (g_a, dtype_a, n1_a, n2_a)
      call ga_inquire (g_b, dtype_b, n1_b, n2_b)

      if ((n1_a .ne. n1_b).or.(n2_a .ne. n2_b))
     $     call errquit (pname//"mats must be same size", 0, 0)

      same = .false.
      
      if (dtype_a .eq. dtype_b)  then !they are the same data type

         if (dtype_a .eq. mt_dbl) then ! real data
            if (.not. ga_create (mt_dbl,n1_a,n1_a, "tmp", 0, 0, g_tmp))
     $           call errquit (pname//"cannot create g_tmp",0,0)

            call ga_zero (g_tmp)
            call ga_dadd (1d0, g_a, -1d0, g_b, g_tmp)
            norm = mat_norm (g_tmp)

            if (norm < tol) then
               same = .true.
            else
               same = .false.
               if (me.eq.0) write (luout, *)
     $              pname//"|| A - B || = ", norm
               call util_flush (luout)
            endif

            if (.not. ga_destroy (g_tmp))
     $           call errquit (pname//"cannot destroy g_tmp",0,0)

         elseif (dtype_a .eq. mt_dcpl) then ! complex data
            if (.not. ga_create (mt_dcpl,n1_a,n1_a, "tmp", 0, 0, g_tmp))
     $           call errquit (pname//"cannot create g_tmp",0,0)

            call ga_zero (g_tmp)
C            call ga_zadd (zone, g_a, nzone, g_b, g_tmp)
            call ga_add (zone, g_a, nzone, g_b, g_tmp)
            norm = mat_norm (g_tmp)

            if (norm < tol) then
               same = .true.
            else
               same = .false.
c              if (me.eq.0) write (luout, *)
c    $              pname//"|| A - B || = ", norm
               call util_flush (luout)
            endif

            if (.not. ga_destroy (g_tmp))
     $           call errquit (pname//"cannot destroy g_tmp",0,0)

      else
         call errquit (pname//"unsupported data type", 0, 0)
      endif                     ! ends real/complex data check
            
      else                      ! they are different data types
         call errquit (pname//"both must be same data type", 0, 0)
         
      endif                     ! ends same data check
      
      mat_is_same_semiemp = same

      end function

C====================================================================
C
C     Computes Gaussian function for windowing pulse.
C
      double precision function rt_semiemp_gaussian_window (frequency, 
     $ center,width,tt)
      implicit none

#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
#include "rt_tddft.fh"

C     == Inputs ==

      double precision tt, center, frequency, width, phase 

C     == Parameters ==
      character(*) ,parameter  :: pname = "rt_semiemp_gaussian_window: "

C     == Variables ==
      double precision arg, gout, val
      phase = 0.0d0
      rt_semiemp_gaussian_window = 0.0d0 
      val = 1.0d0
      gout = 0.0d0
      arg = 0.0d0

      if (frequency.ne.0d0) val =  sin(frequency*tt + phase) 
      arg = -1d0 * (tt - center)*(tt - center) /
     $     (2d0*width*width)

c$$$  dd          c = field%max / sqrt (2d0 * dpi) / s  !normalization

      gout = exp (arg)
c     write(6,*)" val, gout =", val, gout
      rt_semiemp_gaussian_window = gout*val

      end function

C------------------------------------------------------------------------------------
C
C------------------------------------------------------------------------------------
      Subroutine rt_semiemp_anhi(rtdb,g_densa,g_densb,g_scr,
     $ n,nel)

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

C     == Parameters ==
      character(*), parameter     :: pname       = "rt_semiemp_anhi:"
      integer rtdb 
      integer g_zdensa, g_zdensb, g_zscr
      integer g_densa, g_densb, g_scr
      integer g_ab, g_ba, g_anticom, g_correc
      integer  nel,n
      double precision tt, xterm
      double precision unel, norma, normb
      double precision pptrace, avgs2, mult, ga_trace_diag
      logical ok
      
c     call ga_zero(g_densa)
c     call ga_zero(g_densb)
      call ga_zero(g_scr)
c     write(6,*)"pptrace,avgs2",pptrace,avgs2 
      if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     $       call errquit(pname//'failed to read mult rtdb',0,RTDB_ERR)

c     call convert_z2d(g_zdensa,g_densa,g_scr)
c     call convert_z2d(g_zdensb,g_densb,g_scr)

      ok = .true.
      ok=ok.and.ga_create(mt_dbl,n,n,"AB" , 0, 0, g_ab)
      ok=ok.and.ga_create(mt_dbl,n,n,"BA" , 0, 0, g_ba)
      ok=ok.and.ga_create(mt_dbl,n,n,"BA" , 0, 0, g_anticom)
      ok=ok.and.ga_create(mt_dbl,n,n,"BA" , 0, 0, g_correc)
      
      call ga_zero(g_ab)
      call ga_zero(g_ba)
      call ga_zero(g_anticom)
      call ga_zero(g_correc)

      if (.not.ok) call errquit(pname//"GA allocation failed", 0,GA_ERR)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_densa,g_densb,0.0d0,
     1 g_scr)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_densa,g_densb,0d0,
     1 g_ab)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_densb,g_densa,0d0,
     1 g_ba)

       pptrace=ga_trace_diag(g_scr)
       unel = 2*((mult-1)/2)
       xterm = unel+2-(DBLE(nel)-unel)/2+pptrace  
       write(6,*)"unel,pptrace,multi,xterm",unel,pptrace,mult,xterm

       call ga_add(1.0d0,g_ab,1.0d0,g_ba,g_anticom)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_ab,g_densa,0d0,
     1 g_correc)
      
      call ga_add(1.0,g_correc,-0.5,g_anticom,g_scr) 

      call ga_add(1.0,g_densa,-1/xterm,g_scr,g_correc)

      call ga_copy(g_correc,g_densa)
      

      call ga_zero(g_correc)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_ba,g_densb,0d0,
     1 g_correc)

      call ga_add(1.0,g_correc,-0.5,g_anticom,g_scr)

      call ga_add(1.0,g_densb,-1/xterm,g_scr,g_correc)

      call ga_copy(g_correc,g_densb)
      

c     call ga_norm_infinity(g_densa,norma)

c     call ga_norm_infinity(g_densb,normb)

c     call ga_scale(g_densa,1/norma)
c     call ga_scale(g_densb,1/normb)

      call ga_dgemm('N','N',n,n,n,1.0d0,g_densa,g_densb,0d0,
     1 g_scr)

       pptrace=ga_trace_diag(g_scr)
c      write(6,*)"mult,nel,pptrace",mult,nel,pptrace
       avgs2 = ((mult-1)/2)**2+dble(nel)/2-pptrace

      if (ga_nodeid().eq.0) write(6,*) "<S2>", avgs2

      
       ok = .true.
       ok = ok .and. ga_destroy(g_correc)
       ok = ok .and. ga_destroy(g_anticom)       
       ok = ok .and. ga_destroy(g_ba)
       ok = ok .and. ga_destroy(g_ab)


      end subroutine 
