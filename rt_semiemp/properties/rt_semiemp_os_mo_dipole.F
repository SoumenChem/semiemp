      Subroutine rt_semiemp_os_mo_dipole(rtdb,g_zdens,g_zdensa,g_zdensb
     1 ,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,NB,NRR,ISW1,Z12,D12,NEL,E,
     1 adipx,adipy,adipz,bdipx,bdipy,bdipz,tag,tt,g_s,g_ft,it,nrestarts
     1 ,got_orbspace)

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

      integer rtdb
      integer g_dipx, g_dipy, g_dipz, g_zdens, g_temp
      integer g_zdensa, g_zdensb
      integer g_densre, g_densim, g_s, g_scr, g_ft
      double precision RAD(*),Z12(*),D12(*),H(*),BB(*)
      double precision ALPHA(*),CO(*),E(*)
      double precision ADIPX(*), ADIPY(*), ADIPZ(*)
      double precision BDIPX(*), BDIPY(*), BDIPZ(*)
      integer NU(*),NW(*),NP(*),KAN(*)
      integer LUMO(1000), HOMO(1000), MIM(2)
      integer NEL,NA,ISW1,NB,NRR
      integer NIN(*),i,me,it,nrestarts
      integer ldens, idens, ldipx, idipx, ldipy, idipy, ldipz, idipz
      integer orbspace(2), aorb(3)
      integer iPME, lPME, ihomo, lhomo, ilumo, llumo 
      integer iocc, locc, iunocc, lunocc
      double precision DA(3),DB(3), tt, pstrace, ga_trace_diag
      integer Msg_mat_triangular, Msg_dipole
      integer  iorb, norb
      integer ehinfo(5)
      double precision PME(3,505)
      logical status, got_orbspace, got_elh
      character*20 tag


       got_elh = .false.
       me = ga_nodeid()
       
c       write(6,*)"i am in semiemp_dipole"!,(KAN(I),I=1,na)
             if (.not.MA_Push_Get(MT_Dbl,nb*nb,'Densre',ldens,idens))
     &   call errquit('rtsemiemp: cannot allocate svecs',0,
     1    MA_ERR)

      if (.not. ga_create(mt_dbl, nb, nb, "Re[Dens]", 0, 0, g_densre))
     $     call errquit ("couldnt create g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_densim, "Im[Dens]"))
     $     call errquit ("couldnt duplicate g_densre", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_scr, "temp"))
     $     call errquit ("couldnt duplicate g_scr", 0, GA_ERR)
      if (.not. ga_duplicate(g_densre, g_temp, "temp1"))
     $     call errquit ("couldnt duplicate g_scr", 0, GA_ERR)

       call ga_zero(g_densre)
       call ga_zero(g_densim)
       call ga_zero(g_scr)
       call ga_zero(g_temp)

c     call ga_get (g_dipx, 1, nb, 1, nb, Dbl_MB(idipx),nb)
c     call ga_get (g_dipy, 1, nb, 1, nb, Dbl_MB(idipy),nb)
c     call ga_get (g_dipz, 1, nb, 1, nb, Dbl_MB(idipz),nb)


       if (.not.rtdb_get(rtdb,'semiemp:MIM',mt_int,2,MIM))
     1 call errquit('semiemp: failed to read MIM',0, RTDB_ERR)

       if (got_orbspace) then
c       write(6,*)"here"
        if (.not.rtdb_get(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1  call errquit('rt_semiemp_dipole: failed to read orbspace',0,
     2   RTDB_ERR)
        if (.not.rtdb_get(rtdb,'rt_semiemp:alphaactsize',mt_int,3,aorb))
     1  call errquit('rt_semiemp_dipole: failed to read alphaaorb',0,
     2   RTDB_ERR)

        IORB= aorb(1)*aorb(2)

        if (.not.MA_Push_Get(MT_Dbl,3*IORB,'PME',lPME,iPME))
     &   call errquit('rt_Semiemp_dipole: cannot allocate PME',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Int,IORB,'IHOMO',locc,iocc))
     &   call errquit('rt_Semiemp_dipole: cannot allocate IHOMO',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Int,IORB,'ILUMO',lunocc,
     1  iunocc))
     &   call errquit('rt_semiemp_dipole: cannot allocate ILUMO',0,
     1    MA_ERR)


c       write(6,*)"in dipole",aorb(1),aorb(2)!orbspace(1),orbspace(2)
      end if      
       call convert_z2d(g_zdensa,g_densre,g_densim)

c     call ga_print(g_densre)
      if (mod(it,nrestarts).eq.0) then
      call ga_dgemm('N','N',nb,nb,nb,z1,g_densre,g_ft,z0,
     1 g_scr)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_ft,g_scr,z0,
     1 g_temp)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_s,g_temp,z0,
     1 g_scr)

       pstrace = ga_trace_diag(g_scr)
       if (me.eq.0) write(6,*)"alpha electrons=",pstrace
       end if

c     call ga_elem_multiply(g_temp,g_s,g_scr)
c     call rt_semiemp_mullpop(rtdb,g_scr,0)

      call ga_trimat(g_densre,Dbl_MB(idens),nb)
c     do i = 0, 5
c      write(6,*)"dens",Dbl_MB(idens+i)
c     end do
c     if (got_elh) then
c      IORB =  aorb(1)*aorb(2) + ehinfo(1)*(aorb(1)-ehinfo(1))
c     else
c      IORB = aorb(1)*aorb(2)
c     end if

c     write(6,*)"IORB",IORB, aorb(1), aorb(2)

      call DIPOLE(rtdb,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,NB,NRR,
     1 ISW1,Z12,D12,Dbl_MB(idens),NEL,DA,adipx,adipy,adipz,2,
     1 Dbl_mb(iPME),int_mb(iocc),int_mb(iunocc),IORB,MIM(2))


c     do i = 0, nrr-1
c      write(6,*)"dens after"!,int_mb(ilumo+iorb-1),Dbl_MB(idens+i)
c     end do
c     write(6,*)"me",me
      if (me.eq.0) then
c     write(6,*)"IORB",IORB
       write (luout, "(a, a, 2x, 1f11.5, 3es22.12e3, 5x, a)")
     $ trim(tag),"totalalpha", tt, DA(1), DA(2), DA(3), 
     $ "# "//"Dipole moment [system]"
      if (got_orbspace) then
       write(6,*)"IROB" , IORB
       DO i = 1, iorb
c      write(luout,*)"Occupied-unoccpied",HOMO(i),"-",LUMO(i)
       write(luout,*)"Occupied-unoccpied",int_mb(iocc+i-1),"-",
     1  int_mb(iunocc+i-1)
       write(luout,"(a, 2x, i3, 1f11.5, 3es22.12e3, 6x, a)")
     1 trim(tag), i, tt, Dbl_MB(iPME+3*(i-1)), Dbl_MB(iPME+3*(i-1)+1), 
     1 Dbl_MB((iPME+3*(i-1))+2),"# "//"Alpha MO pair
     1 [system]"
       END DO  
      endif
         call util_flush (luout)

         call util_flush (luout)
      endif


      if (got_orbspace) then
        if (.not.ma_pop_stack(lunocc))
     &   call errquit('rt_semiem_dipole: cannot pop stack unocc',3,
     1   MA_ERR)
        if (.not.ma_pop_stack(locc))
     &   call errquit('rt_semiemp_dipole: cannot pop stack occ',3,
     1   MA_ERR)
        if (.not.ma_pop_stack(lpme))
     &   call errquit('rt_semiep_dipole: cannot pop stack lPME',3,
     1   MA_ERR)
       end if

       if (got_orbspace) then
c       write(6,*)"here"
        if (.not.rtdb_get(rtdb,'rt_semiemp:orbspace',mt_int,2,orbspace))
     1  call errquit('rt_semiemp_dipole: failed to read orbspace',0,
     2   RTDB_ERR)
        if (.not.rtdb_get(rtdb,'rt_semiemp:alphaactsize',mt_int,3,aorb))
     1  call errquit('rt_semiemp_dipole: failed to read alphaaorb',0,
     2   RTDB_ERR)

        IORB= aorb(1)*aorb(2)

        if (.not.MA_Push_Get(MT_Dbl,3*IORB,'PME',lPME,iPME))
     &   call errquit('rt_Semiemp_dipole: cannot allocate PME',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Int,IORB,'IHOMO',locc,iocc))
     &   call errquit('rt_Semiemp_dipole: cannot allocate IHOMO',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Int,IORB,'ILUMO',lunocc,
     1  iunocc))
     &   call errquit('rt_semiemp_dipole: cannot allocate ILUMO',0,
     1    MA_ERR)


c       write(6,*)"in dipole",aorb(1),aorb(2)!orbspace(1),orbspace(2)
      end if

       call ga_zero(g_densre)
       call ga_zero(g_densim)
       call ga_zero(g_scr)
       call ga_zero(g_temp)

      call convert_z2d(g_zdensb,g_densre,g_densim)
       
      if (mod(it,nrestarts).eq.0) then
      call ga_dgemm('N','N',nb,nb,nb,z1,g_densre,g_ft,z0,
     1 g_scr)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_ft,g_scr,z0,
     1 g_temp)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_s,g_temp,z0,
     1 g_scr)

       pstrace = ga_trace_diag(g_scr)
       if (me.eq.0) write(6,*)"beta electrons =",pstrace
       end if

c     call ga_elem_multiply(g_temp,g_s,g_scr)
c     call rt_semiemp_mullpop(rtdb,g_scr,0)

      call ga_trimat(g_densre,Dbl_MB(idens),nb)
c     do i = 0, 5
c      write(6,*)"dens",Dbl_MB(idens+i)
c     end do
c     if (got_elh) then
c      IORB =  aorb(1)*aorb(2) + ehinfo(1)*(aorb(1)-ehinfo(1))
c     else
c      IORB = aorb(1)*aorb(2)
c     end if

c     write(6,*)"IORB",IORB, aorb(1), aorb(2)

      call DIPOLE(rtdb,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,NB,NRR,
     1 ISW1,Z12,D12,Dbl_MB(idens),NEL,DB,bdipx,bdipy,bdipz,2,
     1 Dbl_mb(iPME),int_mb(iocc),int_mb(iunocc),IORB,MIM(1))


      if (me.eq.0) then
c     write(6,*)"IORB",IORB
       write (luout, "(a, a, 2x, 1f11.5, 3es22.12e3, 5x, a)")
     $ trim(tag),"totalbeta", tt, DB(1), DB(2), DB(3),
     $ "# "//"Dipole moment [system]"
      if (got_orbspace) then
       write(6,*)"IROB" , IORB
       DO i = 1, iorb
c      write(luout,*)"Occupied-unoccpied",HOMO(i),"-",LUMO(i)
       write(luout,*)"Occupied-unoccpied",int_mb(iocc+i-1),"-",
     1  int_mb(iunocc+i-1)
       write(luout,"(a, 2x, i3, 1f11.5, 3es22.12e3, 6x, a)")
     1 trim(tag), i, tt, Dbl_MB(iPME+3*(i-1)), Dbl_MB(iPME+3*(i-1)+1),
     1 Dbl_MB((iPME+3*(i-1))+2),"# "//"Beta MO pair
     1 [system]"
       END DO
      endif
         call util_flush (luout)

         call util_flush (luout)
      endif
      

       if (got_orbspace) then
        if (.not.ma_pop_stack(lunocc))
     &   call errquit('rt_semiem_dipole: cannot pop stack unocc',3,
     1   MA_ERR)
        if (.not.ma_pop_stack(locc))
     &   call errquit('rt_semiemp_dipole: cannot pop stack occ',3,
     1   MA_ERR)
        if (.not.ma_pop_stack(lpme))
     &   call errquit('rt_semiep_dipole: cannot pop stack lPME',3,
     1   MA_ERR)
       end if


       call ga_zero(g_densre)
       call ga_zero(g_densim)
       call ga_zero(g_scr)
       call ga_zero(g_temp)

      call convert_z2d(g_zdens,g_densre,g_densim)

      if (mod(it,nrestarts).eq.0) then
      call ga_dgemm('N','N',nb,nb,nb,z1,g_densre,g_ft,z0,
     1 g_scr)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_ft,g_scr,z0,
     1 g_temp)

      call ga_dgemm('N','N',nb,nb,nb,z1,g_s,g_temp,z0,
     1 g_scr)

       pstrace = ga_trace_diag(g_scr)
       if (me.eq.0) write(6,*)"total electrons =",pstrace
       end if

c     call ga_elem_multiply(g_temp,g_s,g_scr)
c     call rt_semiemp_mullpop(rtdb,g_scr,0)

c     call ga_trimat(g_densre,Dbl_MB(idens),nb)

c     call DIPOLE(rtdb,RAD,KAN,NA,CO,NU,NW,NP,NIN,ALPHA,H,BB,NB,NRR,
c    1 ISW1,Z12,D12,Dbl_MB(idens),NEL,D,dipx,dipy,dipz,2,Dbl_mb(iPME),
c    1 int_mb(iocc),int_mb(iunocc),IORB)


      if (me.eq.0) then
c     write(6,*)"IORB",IORB
       write (luout, "(a, a, 2x, 1f11.5, 3es22.12e3, 5x, a)")
     $ trim(tag),"total", tt, DA(1)+DB(1), DA(2)+DB(2), DA(3)+DB(3),
     $ "# "//"Dipole moment [system]"
         call util_flush (luout)

         call util_flush (luout)
      endif
      
      status=.true.
      status = status.and.ga_destroy(g_temp)
      status = status.and.ga_destroy(g_scr)
      status=status.and.ga_destroy(g_densim)
      status=status.and.ga_destroy(g_densre)

      if (.not.ma_pop_stack(ldens))
     & call errquit('rt_semiemp_mo_dipole:cannot pop stack dens',3, 
     & MA_ERR)


      end subroutine


