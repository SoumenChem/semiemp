C
C     prop_euler.F
C
C     Propagate MO density matrix/matricies using Euler method.
C
      subroutine prop_semiemp_euler(g_zdens_ao,g_zfock_ao,n,dt)
      implicit none

#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "msgids.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "util.fh"
#include "cdft.fh"
#include "matutils.fh"
#include "rt_tddft.fh"

      
C     == Inputs ==


C     == Parameters ==
      character(*), parameter     :: pname      = "prop_semiemp_euler: "


C     == Variables ==
      integer g_zdens_ao
      integer g_zfock_ao
      integer g_ztmp1, g_ztmp2
      logical ok
      integer imat, n
      double complex zidt
      double complex zdamp1, zdamp2
      double complex zjunk
      double precision dt
      integer me
      integer i

      double precision kcen, kwid, k1, k2, kval, kw, kk, kmax

      


c     me = ga_nodeid()
      
C
C     P'(t+dt) = P'(t) + dP'/dt * dt  (all mats in MO basis)
C
C     where, dP'/dt = -i [F', P']   (in atomic units)
C
C     so, P'(t+dt) = P'(t) - i*dt*[F', P']
C                  = P'(t) - i*dt*(F'P' - P'F')
C

      
C
C     Allocation
C     

         if (.not.ga_create(mt_dcpl, n, n,
     $        "temp1", 0, 0, g_ztmp1))
     $        call errquit ("failed to create ztemp1", 0, GA_ERR)

         if (.not.ga_create(mt_dcpl, n, n,
     $        "ztemp2", 0, 0, g_ztmp2))
     $        call errquit ("failed to create ztemp2", 0, GA_ERR)


C

         
C
C     Compute (dP'/dt)*dt = -i*dt*[F', P'] = -i*dt(F'P' - P'F'); ends up stored in g_ztmp1_mo
C

         call ga_zgemm ("N", "N",
     $        n, n, n, z1, g_zdens_ao, g_zfock_ao,
     $        z0, g_ztmp1)
         
         call ga_zgemm ("N", "N",
     $        n, n, n, z1, g_zfock_ao, g_zdens_ao,
     $        z0, g_ztmp2)
         
         call ga_zadd (z1, g_ztmp2, zn1,
     $        g_ztmp1, g_ztmp1)


C     
C     Step density matrix forward and replace input dens mat.
C     
C     zidt = zni*dcmplx (params%dt, 0d0)  ! -i*dt
         zidt = dcmplx (0d0, -1d0*dt)
         
         call ga_zadd (z1, g_zdens_ao, zidt,
     $        g_ztmp1, g_zdens_ao)



C
C     Clean up
C
         if (.not. ga_destroy (g_ztmp2))
     $        call errquit (pname//"failed to destroy zfock_ao", 0, 0)

         if (.not. ga_destroy (g_ztmp1))
     $        call errquit (pname//"failed to destroy zdens_ao", 0, 0)



      end subroutine
c $Id: prop_euler.F 24344 2013-06-22 06:15:21Z d3y133 $
