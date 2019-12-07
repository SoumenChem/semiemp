      subroutine semiemp_input(rtdb)
c
      implicit none
c
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
c
      integer rtdb
      character*20 test
      character*8 scftype     !type of hf calculation
      character*8 apx          !approximate hamiltonina type
      logical rhf,uhf
      logical cndo,indo
      integer maxiter           !maximum iteration
      integer inttyp            !type of integral
      integer gfac              !component of ci part
      logical semiemp_energy    
      double precision scftol   !convergence criteria
      integer nel                !number of valence electron
      integer charge
      double precision mult              !multiplicity
      double precision infa(6)            !interaction factor
      integer III                !type of calculation (e.g energy)
      integer IA            !types of calculation (e.g energy,TS)
      integer geom
      logical status 
      logical  geom_rtdb_ncent  !for getting number of atoms
      integer nat               !number of atoms
      
      if (ga_nodeid().eq.0)
     &     write(luout,*) "building semiemp"
c
c      Create geometry object and load in geometry from rtdb
c      if (.not.geom_create(geom,'geometry'))
c    1  call errquit('semiemp_init: geom_create failed',0, GEOM_ERR)
c      if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
c    1      call errquit('semiemp_init: no geometry found',0, RTDB_ERR)

      scftype= 'rhf'
      if (.not.rtdb_cput(rtdb,'semiemp:scftype',1,scftype))
     &  call errquit('semiemp_input: failed writing to rtdb',0,RTDB_ERR)
       charge=0
      if (.not.rtdb_put(rtdb,'semiemp:charge',mt_int,1,charge))
     &  call errquit('semiemp_input: failed writing to rtdb',0,RTDB_ERR)
       mult=1.0D0
      if (.not.rtdb_put(rtdb,'semiemp:mult',mt_dbl,1,mult))
     &  call errquit('semiemp_input: failed writing to rtdb',0,RTDB_ERR)
       scftol=1.0D-8
      if (.not.rtdb_put(rtdb,'semiemp:scftol',mt_dbl,1,scftol))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
       apx='INDO/1'
      if (.not.rtdb_cput(rtdb,'semiemp:apx',1,apx))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
      inttyp=0 
      if (.not.rtdb_put(rtdb,'semiemp:inttyp',mt_int,1,inttyp))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
       nel=0
       if (.not.rtdb_put(rtdb,'semiemp:nel',mt_int,1,nel))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
       III=000
        if (.not.rtdb_put(rtdb,'semiemp:III',mt_int,1,III))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
        IA=0 
        if (.not.rtdb_put(rtdb,'semiemp:IA',mt_int,1,IA))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
        infa(1)=1.0D0
        infa(2)=1.0D0
        infa(3)=1.0D0
        infa(4)=1.0D0
        infa(5)=1.0D0
        infa(6)=1.0D0
         if (.not.rtdb_put(rtdb,'semiemp:infa',mt_dbl,6,infa))
     &  call errquit('smeiemp_input: failed writing to rtdb',0,RTDB_ERR)
 
c----------
c Read input
c ----------
c
 10     continue
        if (.not. inp_read()) 
     &  call errquit('semiemp_input: failed reading input',0, INPUT_ERR)
           if (.not. inp_a(test)) 
     &  call errquit('semiemp_input: failed reading input',0, INPUT_ERR)
         if (inp_compare(.false.,test,'scftype')) then
         if (.not.inp_a(scftype)) then
          write(LuOut,*) 'semiemp_input: SCFTYPE value not found; ',
     1      'default value of rhf will be used'
          scftype='rhf'
          endif
          if (.not.rtdb_cput(rtdb,'semiemp:scftype',1,scftype))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         else if (inp_compare(.false.,test,'apx')) then
         if (.not.inp_a(apx)) then
          write(LuOut,*) 'semiemp_input: apx value not found; ',
     1      'default value of INDO/1 will be used'
          apx='INDO/1'
          endif
          if (.not.rtdb_cput(rtdb,'semiemp:apx',1,apx))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
          else if (inp_compare(.false.,test,'mult')) then
          if (.not.inp_f(mult)) then
          write(LuOut,*) 'semiemp_input:mult value not found; ',
     1      'default value of 1 will be used'
          mult=1.0D0
          endif
          if (.not.rtdb_put(rtdb,'semiemp:mult',mt_dbl,1,mult))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
           else if (inp_compare(.false.,test,'scftol')) then
         if (.not.inp_f(scftol)) then
          write(LuOut,*) 'semiemp_input:scftol value not found; ',
     1      'default value of .0001 will be used'
          scftol=1.0D-7
          endif
          if (.not.rtdb_put(rtdb,'semiemp:scftol',mt_dbl,1,scftol))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        else if (inp_compare(.false.,test,'maxiter')) then
        if (.not.inp_i(maxiter)) then
          write(LuOut,*) 'semiemp_input:maxiter value not found; ',
     1      'default value of 1000 will be used'
          maxiter=50
         endif
          if (.not.rtdb_put(rtdb,'semiemp:maxiter',mt_int,1,maxiter))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        else if (inp_compare(.false.,test,'nel')) then
        if (.not.inp_i(nel)) then
          write(luout,*) 'NEL value not found' 
         endif
          if (.not.rtdb_put(rtdb,'semiemp:nel',mt_int,1,nel))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        else if (inp_compare(.false.,test,'inttyp')) then
        if (.not.inp_i(inttyp)) then
          write(LuOut,*) 'semiemp_input: INTTYP value not found; ',
     1      'default value of 0 will be used'
          inttyp=0
        endif
          if (.not.rtdb_put(rtdb,'semiemp:inttyp',mt_int,1,inttyp))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         else if (inp_compare(.false.,test,'gfac')) then
        if (.not.inp_i(gfac)) then
          write(LuOut,*) 'semiemp_input: GFAC value not found; ',
     1      'default value of 0 will be used'
          gfac=0
        endif
          if (.not.rtdb_put(rtdb,'semiemp:gfac',mt_int,1,gfac))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
        else if (inp_compare(.false.,test,'III')) then
        if (.not.inp_i(III)) then
          write(LuOut,*) 'semiemp_input: III value not found; ',
     1      'default value of 0000 will be used'
          III=0000
        endif
          if (.not.rtdb_put(rtdb,'semiemp:III',mt_int,1,III))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)  
         else if (inp_compare(.false.,test,'IA')) then
        if (.not.inp_i(IA)) then
          write(LuOut,*) 'semiemp_input: IA value not found; ',
     1      'default value of 0 will be used'
          IA=0
        endif
          if (.not.rtdb_put(rtdb,'semiemp:IA',mt_int,1,IA))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
           else if (inp_compare(.false.,test,'nat')) then
        if (.not.inp_i(nat)) then
          write(LuOut,*) 'semiemp_input: nat value not found; ',
     1      'it will be calculated'
          status= geom_rtdb_ncent(rtdb,"geometry",nat) 
        endif
          if (.not.rtdb_put(rtdb,'semiemp:nat',mt_dbl,1,nat))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       RTDB_ERR)
         else if (inp_compare(.false.,test,'infa')) then
        if (.not.inp_i(infa(1))) infa(1)=1.0
         if (.not.inp_i(infa(2))) infa(2)=1.0
          if (.not.inp_i(infa(3))) infa(3)=1.0
         if (.not.inp_i(infa(4))) infa(4)=1.0
         if (.not.inp_i(infa(5))) infa(5)=1.0
         if (.not.inp_i(infa(6))) infa(6)=1.0
          if (.not.rtdb_put(rtdb,'semiemp:infa',mt_dbl,6,infa))
     1    call errquit('semiemp_input: failed writing to rtdb',0,
     &       INPUT_ERR)
            else if (inp_compare(.false.,test,'end')) then
          goto 20
           else
           call errquit('semiemp_input: unknown directive',0, INPUT_ERR)
          endif
          goto 10
c
c ------
c Return
c ------
c
  20       return
           end
