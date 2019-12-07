      subroutine semiemp_canorg( nmo, svals, g_svecs, g_xij,nbf)
c
c     Canonical Orthogonalization
c     Szabo & Ostlund page 144 (Dover Ed.)
c
c     $Id: dft_canorg.F 27479 2015-09-08 22:47:26Z edo $
c
      implicit none
#include "errquit.fh"
c
      double precision svals(*) !  S evals [input]
      integer g_svecs ! GA handle for S evecs [input]
      integer g_xij  ! GA [ouput]
      integer nmo, nbf! [out]
c
#include "mafdecls.fh"
#include "global.fh"
c#include "cdft.fh"
c
      integer me, nproc, i
      integer ltmpm, itmpm,ndep
      double precision aaz
      integer ilo,ihi,jlo,jhi,numi,numj,icount
      integer jfirst,jlast,ptrfirst
c
      me = ga_nodeid()
      nproc = ga_nnodes()
      ndep=nbf-nmo
c
c
       write(6,*)"here 1",nbf,nmo
c        Build Xij
c
      call ga_zero(g_xij)
      call ga_distribution(g_svecs,
     .           ga_nodeid(), ilo, ihi, jlo, jhi)
      if (ilo.gt.0 .and. ilo.le.ihi) then
         numi =  ihi-ilo+1
         numj =  jhi-jlo+1
         if (numi.gt.0 .and. numj.gt.0) then
            if (.not. MA_Push_Get(MT_Dbl, numi*numj, 'tmpm',
     T           ltmpm, itmpm)) call errquit(
     T           'dftcano: failed to alloc tmpm',0, MA_ERR)
            call ga_get(g_svecs,ilo,ihi,jlo,jhi,
     .           dbl_mb(itmpm),numi)
            icount=itmpm
            jfirst=99999
            jlast=-1
            ptrfirst=0
            do i=jlo,jhi
               if(i.gt.ndep) then
                  if(ptrfirst.eq.0) ptrfirst=icount
                  jfirst=min(i-ndep,jfirst)
                  jlast=max(i-ndep,jlast)
                  aaz = svals(i)
                  aaz = 1.d0/sqrt(aaz)
                  call dscal(numi, aaz, DBL_MB(icount), 1)
               endif
               icount=icount+numi
            enddo
            write(6,*)"here 2",jlast
            if (jlast.ne.-1)
     C           call ga_put(g_xij,ilo,ihi,jfirst,jlast,
     .           dbl_mb(ptrfirst),numi)
           write(6,*)"here 4"
            if(.not.MA_Pop_Stack(ltmpm)) call errquit(
     .           'failed popstack in dftcanorg ',0, 0)
         endif
      endif
c
      return
      end










