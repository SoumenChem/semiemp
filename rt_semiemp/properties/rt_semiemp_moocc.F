C====================================================================
C
C     Calculates the time-dependent occupations of the molecular orbitals.
C
C     n_k(t) = C'_k^+ P'(t) C'_k
C
C     where C'_k is the k^th eigenvector of the ground state Fock
C     matrix, and P' is the density matrix in the MO basis.  Note that P
C     is complex, but C' is real since it is from the SCF.  This means
C     that we can just use the real part of the dens mat, and take
C     transpose of C' instead of conjg transpose.
C
C     Note, you can send this either the full dens mat, or just the
C     alpha or beta spin part.
C
C     Note, can also compute using matrix multiplications (XXX double check):
C     n_k(t) = [C'^+ P'(t) C']_kk
C
      subroutine rt_semiemp_moocc(rtdb,g_densre_mo,nbf,g_movecs_gs,tt,
     1 tag)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "geom.fh"
#include "util.fh"
#include "rt_tddft.fh"
#include "matutils.fh"      
#include "rtdb.fh"
#include "inp.fh"

C     == Inputs ==
      integer rtdb
      integer g_densre_mo  !re part of density matrix in MO basis
      integer g_movecs_gs  !ground state movecs
      integer tag

C     == Outputs ==
      double precision moocc(10) !MO occupations

      
C     == Parameters ==
      character(*),parameter :: pname = "rt_tddft_moocc: "


C     == Variables ==
      integer me
      integer i
      integer nbf
      integer lveck, iveck      !handle and index for kth eigenvector
      double precision occk     !occupation of orbital k
      integer g_veck, g_tmp, g_veckt
      integer k
      integer norb, numorb(100)
      double precision dval, tt
      character*20 scftype



      me = ga_nodeid()

       if (.not.rtdb_cget(rtdb,'semiemp:scftype',1,scftype)) then
      call errquit('semiemp_energy: failed to read scftype',0, RTDB_ERR)
       endif

      if (inp_compare(.false.,scftype,'rhf'))  then

      if (.not.rtdb_get(rtdb,'rt_semiemp:norb',mt_int,1,norb))
     1 call errquit('rt_semiemp_moocc: failed to get norb',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'rt_semiemp:numorb',mt_int,100,numorb))
     1 call errquit('rt_semiemp_moocc: failed to get numorb',0,RTDB_ERR)

       elseif (inp_compare(.false.,scftype,'uhf'))  then
         if (tag.eq.1) then
          if (.not.rtdb_get(rtdb,'rt_semiemp:naorb',mt_int,1,norb))
     1 call errquit('rt_semiemp_moocc: failed to get naorb',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'rt_semiemp:numorba',mt_int,100,numorb))
     1 call errquit('rt_semiemp_moocc: failed to get numorba',0,
     1 RTDB_ERR)
          
        elseif (tag.gt.1) then
       if (.not.rtdb_get(rtdb,'rt_semiemp:nborb',mt_int,1,norb))
     1 call errquit('rt_semiemp_moocc: failed to get nborb',0,RTDB_ERR)

      if (.not.rtdb_get(rtdb,'rt_semiemp:numorbb',mt_int,100,numorb))
     1 call errquit('rt_semiemp_moocc: failed to get numorbb',0,
     1 RTDB_ERR)
                   
          end if 
         end if 
      
      if (.not. ga_create(mt_dbl, nbf, 1,
     $     "k^th evec", 0, 0, g_veck))
     $     call errquit ("failed to create veck", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_tmp, "moocc tmp"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

      if (.not. ga_duplicate(g_veck, g_veckt, "col of transpose C"))
     $     call errquit ("failed to create g_tmp", 0, GA_ERR)

C
C     Load k^th evec in g_veck and k^th column of C in g_vectk.
C
      do i = 1, norb
CXXX  [KAL]: its redundant to have two g_veck
        k = numorb(i)
         call ga_zero (g_veck)
         call ga_zero (g_veckt)
         call ga_zero (g_tmp)

         call ga_copy_patch ("N",
     $        g_movecs_gs, 1, nbf, k, k,
     $        g_veck, 1, nbf, 1, 1)
         call ga_copy_patch ("T",
     $        g_movecs_gs, 1, nbf, k, k,
     $        g_veckt, 1, nbf, 1, 1)


C     Compute P'(t) C'_k.
         call ga_dgemm ("N", "N", nbf, 1, nbf,
     $        1d0, g_densre_mo, g_veck, 0d0, g_tmp)

         
C     Compute n_k = C'_k^T P'(t) C'_k.
         
         moocc(i) = ga_ddot (g_veckt, g_tmp)

      if (me.eq.0) then
         write (luout, "(a, i4, 3x, 1f11.5, 1es22.12e3, 5x, a)")
     $        "MO Occupation", k ,tt, moocc(i), "system"
       end if

      enddo


C     == Clean up ==
      if (.not. ga_destroy(g_veck))
     $     call errquit ("failed to destroy g_veck", 0, GA_ERR)
      if (.not. ga_destroy(g_tmp))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)
      if (.not. ga_destroy(g_veckt))
     $     call errquit ("failed to destroy g_tmp", 0, GA_ERR)

   10 continue
      end subroutine


C====================================================================
