C
C     rt_tddft_input.F
C     
C     Parses input deck for rt-tddft parameters.
C     
C
      subroutine rt_semiemp_input(rtdb)
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
      character(*), parameter :: pname = "rt_semiemp_input: "
      integer, parameter      :: maxnum = 999999999 ! dummy # when using "*" as # prints

C     == Varibles ==
      logical done
      character*255 test, curr_popt, loadopts, loadtarget
      logical got_opt, got_elh

      double precision tmin, tmax, dt
      
      logical prof, noprop, static, nodisk, matrix_checks
C      logical dplot_do, subgs
C      character*20 dplot_opts
      integer nchecks, nprints, nrestarts
C      integer nsnapshots
      character*20 num_str
      integer checklvl
      double precision field_max, center, width, frequency
      double precision ecut
      integer orbitals(2), orbspace(2)
      integer frag_charge(2)
      integer norb, numorb(100)
      integer naorb, nborb, numorba(100), numorbb(100)
      double precision occuorb(10), occuorba(10), occuorbb(10)
      character*255 field_type
      character*18  polarization, spin
      character*255 tag_in, method
      integer prop_method, exp_method
      logical use_dmat

      
C     (parameters that must be supplied--no defaults)
      logical got_tmin, got_tmax, got_dt, got_fmax, got_ecut
      logical got_orbitals, got_orbspace
      integer nfields, nexcites
      integer i 
      

C     == External ==
      integer, external :: atoi


C     (values which have no defaults and must be supplied)
c     write(6,*)"i am in rt_semiemp_input"
      got_tmin = .false.
      got_tmax = .false.
      got_dt = .false.
      got_fmax = .false.
      got_orbspace = .false.
      got_elh = .false.
     
      nfields = 0
      nexcites = 0

      if (.not.rtdb_put(rtdb,'rt_semiemp:nfields',
     $     mt_int,1,nfields))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      if (.not.rtdb_put(rtdb,'rt_semiemp:nexcites',
     $     mt_int,1,nexcites))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

         if (.not.rtdb_put(rtdb,'rt_semiemp:got_orbspace',
     $           mt_log,1,got_orbspace))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
      
         if (.not.rtdb_put(rtdb,'rt_semiemp:got_elh',
     $           mt_log,1,got_elh))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Dump defaults into rtdb (will be replaced later if values found in
C     input deck).
C
c     call rt_semiemp_input_put_defaults (rtdb)


      
C     
C     Parse the input; we will put in rtdb later after checking.
C     
      done = .false.
      do while (.not. done)

         
         if (.not. inp_read()) 
     $        call errquit(pname//'Read failed input',0, INPUT_ERR)
         if (.not. inp_a(test)) 
     $        call errquit(pname//'Read failed keyword',0, INPUT_ERR)



         if (inp_compare(.false.,test,"polarization")) then
            if (.not.inp_a(polarization)) then
               call errquit (pname//"failed to parse polarization",0,0)
          endif

            if (.not. rtdb_cput (rtdb, "rt_semiemp:polarization",
     $           1, polarization))
     $           call errquit (pname//
     $           "failed to put 'polarization' target into rtdb",0,0)

       elseif (inp_compare(.false.,test,'norb')) then
            if (.not.inp_i(norb))
     $     call errquit (pname//"norb takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:norb',mt_int,1,
     $     norb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      elseif (inp_compare(.false.,test,'naorb')) then
            if (.not.inp_i(naorb))
     $     call errquit (pname//"naorb takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:naorb',mt_int,1,
     $     naorb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      elseif (inp_compare(.false.,test,'nborb')) then
            if (.not.inp_i(nborb))
     $     call errquit (pname//"nborb takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:nborb',mt_int,1,
     $     nborb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
           
       elseif (inp_compare(.false.,test,'numorb')) then
        do i = 1, norb
            if (.not.inp_i(numorb(i)))
     $     call errquit (pname//"numorb takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:numorb',mt_int,100,
     $     numorb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

c         got_elh=.true.
c        if (.not.rtdb_put(rtdb,'rt_semiemp:got_elh',mt_log,1,got_elh))
c    $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
           
       elseif (inp_compare(.false.,test,'numorba')) then
        do i = 1, naorb
            if (.not.inp_i(numorba(i)))
     $     call errquit (pname//"numorba takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:numorba',mt_int,100,
     $     numorba))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

       elseif (inp_compare(.false.,test,'numorbb')) then
        do i = 1, nborb
            if (.not.inp_i(numorbb(i)))
     $     call errquit (pname//"numorbb takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:numorbb',mt_int,100,
     $     numorbb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
           
       elseif (inp_compare(.false.,test,'occuorb')) then
        do i = 1, norb
            if (.not.inp_f(occuorb(i)))
     $     call errquit (pname//"occuorb takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:occuorb',mt_dbl,10,
     $     occuorb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

         got_elh=.true.
         if (.not.rtdb_put(rtdb,'rt_semiemp:got_elh',mt_log,1,got_elh))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

      elseif (inp_compare(.false.,test,'occuorba')) then
        do i = 1, naorb
            if (.not.inp_f(occuorba(i)))
     $     call errquit (pname//"occuorba takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:occuorba',mt_dbl,10,
     $     occuorba))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

       elseif (inp_compare(.false.,test,'occuorbb')) then
        do i = 1, nborb
            if (.not.inp_f(occuorbb(i)))
     $     call errquit (pname//"occuorbb takes a float", 0, 0)
        end do
           if (.not.rtdb_put(rtdb,'rt_semiemp:occuorbb',mt_dbl,10,
     $     occuorbb))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
           
          got_elh=.true.
         if (.not.rtdb_put(rtdb,'rt_semiemp:got_elh',mt_log,1,got_elh))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Tag/title for the run.
C
        elseif (inp_compare(.false.,test,"tag")) then
            if (.not.inp_a(tag_in)) then
               call errquit (pname//"failed to parse tag",0,0)
            endif
            
            if (len_trim(tag_in) .gt. 24)
     $           call errquit (pname//"tag too long, max length is 24",
     $           0,0)
            
            if (.not. rtdb_cput (rtdb, "rt_semiemp:tag",
     $           1, trim(tag_in)))
     $           call errquit (pname//
     $           "failed to put 'tag' target into rtdb",0,0)

        elseif (inp_compare(.false.,test,"field_type")) then
            if (.not.inp_a(field_type)) then
               call errquit (pname//"failed to parse field_type",0,0)
            endif
            if (.not. rtdb_cput (rtdb, "rt_semiemp:field_type",
     $           1, field_type))
     $           call errquit (pname//
     $           "failed to put 'field_type' target into rtdb",0,0)


C
C     center the field (only for gaussian and Hann)
C
         elseif (inp_compare(.false.,test,'center')) then
            if (.not.inp_f(center)) call errquit (pname//
     $           "center takes a float >= 0", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_semiemp:center',mt_dbl,1,
     $    center))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

C
C     width the field (only for gaussian and Hann)
C
         elseif (inp_compare(.false.,test,'width')) then
            if (.not.inp_f(width)) call errquit (pname//
     $           "width takes a float >= 0", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_semiemp:width',mt_dbl,1,
     $    width))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

C
C     frequency the field (only for gaussian and cw)
C
         elseif (inp_compare(.false.,test,'frequency')) then
            if (.not.inp_f(frequency)) call errquit (pname//
     $           "frequency takes a float >= 0", 0, 0)
           write(6,*)"frequency",frequency
            if (.not.rtdb_put(rtdb,'rt_semiemp:frequency',mt_dbl,1,
     $    frequency))
     $         call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)   
C
C     field max
C

       elseif (inp_compare(.false.,test,'field_max')) then
            if (.not.inp_f(field_max)) call errquit (pname//
     $           "field_max takes a float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_semiemp:field_max',mt_dbl,1,
     $    field_max))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_fmax = .true.
            
C
C
C
C     Look for starting state to load.
C
         elseif (inp_compare(.false.,test,"load")) then
            if (.not.inp_a(loadopts)) then
               call errquit (pname//"failed to parse load option",0,0)
            endif
            
            if (loadopts .eq. "vectors") then
               if (.not.inp_a(loadtarget)) then
                  call errquit (pname//
     $                 "failed to parse 'load vectors' target",0,0)
               endif
               if (.not. rtdb_cput (rtdb, "rt_semiemp:init_movecs",
     $              1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_movecs' target into rtdb",0,0)

            elseif (loadopts .eq. "density") then
               use_dmat=.true.
               if (.not.rtdb_put(rtdb,'rt_semiemp:use_dmat',mt_log,
     $              1,use_dmat))
     $            call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c              Read total (if CS) or alpha (if OS) dmat
               if (.not.inp_a(loadtarget)) then
                  call errquit (pname//
     $                 "failed to parse 'load density' target",0,0)
               endif
               if (.not. rtdb_cput (rtdb, "rt_semiemp:init_dmat1",
     $              1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_density' target into rtdb",0,0)
c              Read beta dmat if given (i.e. doing OS)
               if (inp_a(loadtarget)) then
                 if (.not. rtdb_cput (rtdb, "rt_semiemp:init_dmat2",
     $               1, loadtarget)) call errquit (pname//
     $              "failed to put 'init_density' target into rtdb",0,0)
               endif

            elseif (loadopts .eq. "scf") then
C     (no need to do anything--will use SCF vectors as starting point if nothing in rtdb)
               
            elseif (loadopts .eq. "restart") then ! will look for *.rt_restart file
               if (.not.rtdb_put (rtdb, "rt_semiemp:restart",
     $              mt_log, 1, .true.))
     $              call errquit(pname//
     $              'Write failed to rtdb',0,RTDB_ERR)
               
            else
               call errquit (pname//"invalid 'load' option: "
     $              //loadopts,0,0)
            endif
                  
            
C
C     applied field sub-block
C
c        elseif (inp_compare(.false.,test,'field')) then
c           
c           if (.not. inp_a (field_name))
c    $           call errquit (pname//"failed to read field name",0,0)

c           nfields = nfields + 1
c           
c         if (.not.rtdb_put(rtdb,'rt_semiemp:nfields',mt_int,1,nfields))
c    $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

c           call rt_semiemp_input_field (rtdb, field_name, nfields)


            

C
C     visualization sub-block
C
c        elseif (inp_compare(.false.,test,'visualization')) then
c           call rt_tddft_input_visualization (rtdb)


C
C     MO CAP sub-block
C
c        elseif (inp_compare(.false.,test,'mocap')) then
c           call rt_tddft_input_mocap (rtdb)
            
            
C     
C     tmin
C     
C     XXX HARDCODED TO -5*dt
C
c$$$         elseif (inp_compare(.false.,test,'tmin')) then
c$$$            if (.not.inp_f(tmin)) call errquit (pname//
c$$$     $           "tmin takes a float", 0, 0)
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:tmin',mt_dbl,1,tmin))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$            got_tmin = .true.


C
C     active space
C
            
        elseif (inp_compare(.false.,test,'ecut')) then
            if (.not.inp_f(ecut)) call errquit (pname//
     $           "ecut takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:ecut',mt_dbl,1,
     $    ecut))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_orbspace = .true.
        
         if (.not.rtdb_put(rtdb,'rt_semiemp:got_orbspace',
     $           mt_log,1,got_orbspace))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

                elseif (inp_compare(.false.,test,'orbitals')) then
            if (.not.inp_i(orbitals(1))) 
     $     call errquit (pname//"orbitals takes a float", 0, 0)
            if (.not.inp_i(orbitals(2)))
     $     call errquit (pname//"orbitals takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:orbitals',mt_int,2,
     $     orbitals))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            got_orbspace = .true.
         if (.not.rtdb_put(rtdb,'rt_semiemp:got_orbspace',
     $           mt_log,1,got_orbspace))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

                elseif (inp_compare(.false.,test,'orbspace')) then
          if (.not.inp_i(orbspace(1)))
     $    call errquit (pname//
     $           "orbspace takes a float", 0, 0)
          if (.not.inp_i(orbspace(2)))
     $    call errquit (pname//
     $           "orbspace takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:orbspace',mt_int,2,
     $    orbspace))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_orbspace = .true.
           write(6,*)"orbspace",got_orbspace

         if (.not.rtdb_put(rtdb,'rt_semiemp:got_orbspace',
     $           mt_log,1,got_orbspace))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

C
C  fragment charge
C
         elseif (inp_compare(.false.,test,'frag_charge')) then
          if (.not.inp_i(frag_charge(1)))
     $    call errquit (pname//
     $           "frag_charge takes a float", 0, 0)
          if (.not.inp_i(frag_charge(2)))
     $    call errquit (pname//
     $           "frag_charge takes a float", 0, 0)
           if (.not.rtdb_put(rtdb,'rt_semiemp:frag_charge',mt_int,2,
     $    frag_charge))
     $      call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
C     
C     tmax
C     
         elseif (inp_compare(.false.,test,'tmax')) then
            if (.not.inp_f(tmax)) call errquit (pname//
     $           "tmax takes a float", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_semiemp:tmax',mt_dbl,1,tmax))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_tmax = .true.
            
            
C     
C     dt
C     
         elseif (inp_compare(.false.,test,'dt')) then
            if (.not.inp_f(dt)) call errquit (pname//
     $           "dt takes a float >= 0", 0, 0)
            if (.not.rtdb_put(rtdb,'rt_semiemp:dt',mt_dbl,1,dt))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            got_dt = .true.


C     
C     checklvl
C     
         elseif (inp_compare(.false.,test,'checklvl')) then
            if (.not.inp_i(checklvl)) call errquit (pname//
     $           "checklvl takes a value of 1, 2, or 3", 0, 0)
            
            if (.not.rtdb_put(rtdb,'rt_semiemp:checklvl',
     $           mt_int,1,checklvl))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)





c        elseif (inp_compare(.false.,test,'tolerances')) then
c           call rt_tddft_input_tolerances (rtdb)
            
C     
C     nchecks
C     
         elseif (inp_compare(.false.,test,'nchecks')) then
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nchecks takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nchecks = maxnum
            else
               nchecks = atoi (num_str)
               if (nchecks < 0)
     $              call errquit (pname//
     $              "nchecks takes an int >= 0 (or *)", 0, 0)
            endif

          if (.not.rtdb_put(rtdb,'rt_semiemp:nchecks',mt_int,1,nchecks))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)




C     
C     nprints (note * means every step, uses large dummy number)
C     
         elseif (inp_compare(.false.,test,'nprints')) then
            
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nprints takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nprints = maxnum
            else
               nprints = atoi (num_str)
               if (nprints < 0)
     $              call errquit (pname//
     $              "nprints takes an int >= 0 (or *)", 0, 0)
            endif

          if (.not.rtdb_put(rtdb,'rt_semiemp:nprints',mt_int,1,nprints))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)



C     
C     nsnapshots
C     
c$$$         elseif (inp_compare(.false.,test,'nsnapshots')) then
c$$$            if (.not.inp_a(num_str))
c$$$     $           call errquit (pname//
c$$$     $           "nsnapshots takes an int >= 0 (or *)", 0, 0)
c$$$
c$$$            if ( trim(num_str) .eq. "*") then
c$$$               nsnapshots = maxnum
c$$$            else
c$$$               nsnapshots = atoi (num_str)
c$$$               if (nsnapshots < 0)
c$$$     $              call errquit (pname//
c$$$     $              "nsnapshots takes an int >= 0 (or *)", 0, 0)
c$$$            endif
c$$$            
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:nsnapshots',
c$$$     $           mt_int,1,nsnapshots))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     
C     nrestarts
C     
         elseif (inp_compare(.false.,test,'nrestarts')) then
            if (.not.inp_a(num_str))
     $           call errquit (pname//
     $           "nrestarts takes an int >= 0 (or *)", 0, 0)

            if ( trim(num_str) .eq. "*") then
               nrestarts = maxnum
            else
               nrestarts = atoi (num_str)
               if (nrestarts < 0)
     $              call errquit (pname//
     $              "nrestarts takes an int >= 0 (or *)", 0, 0)
            endif

            if (.not.rtdb_put(rtdb,'rt_semiemp:nrestarts',
     $           mt_int,1,nrestarts))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            


C
C     perturbing alpha or beta density
C
         elseif (inp_compare(.false.,test,'spin')) then

            if (.not.inp_a(spin)) then
               call errquit (pname//
     $              "failed to parse spin",0,0)
            endif
            if (.not. rtdb_cput (rtdb, "rt_semiemp:spin",
     $           1, spin))
     $           call errquit (pname//
     $           "failed to put 'spin' target into rtdb",0,0)


C
C     Propagator
C
         prop_method = 0
         elseif (inp_compare(.false.,test,'propagator')) then

            if (.not.inp_a(method)) then
               call errquit (pname//
     $              "failed to parse propagator method",0,0)
            endif

            if (method.eq."euler") then
               prop_method = 1  ! euler
            elseif (method.eq."chebyshev") then
               prop_method = 2   
            else
               call errquit (pname//"invalid propagator: "//method)
            endif

            if (.not.rtdb_put(rtdb,'rt_semiemp:prop_method',
     $           mt_int,1,prop_method))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)



C
C     Exponentiation method
C
         elseif (inp_compare(.false.,test,'exp')) then
            
            if (.not.inp_a(method)) then
               call errquit (pname//
     $              "failed to parse exponentiation method",0,0)
            endif

            if (method.eq."diag") then
               exp_method = 2  ! diagonalization

            elseif (method.eq."pseries") then
               exp_method = 1  ! power series
               
            elseif (method.eq."magnus") then
               exp_method = 3   ! baker-campbell-hausdorff (disabled) !!since assumes e^X A e^-X
               call errquit (pname//"BCH exp disabled"//method)
               
            else
               call errquit (pname//
     $              "invalid exponentiation method: "//method)
            endif
            
            if (.not.rtdb_put(rtdb,'rt_semiemp:exp_method',
     $           mt_int,1,exp_method))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            


C
C     profiling
C
         elseif (inp_compare(.false.,test,'prof')) then
            prof=.true.

            if (.not.rtdb_put(rtdb,'rt_semiemp:prof',mt_log,1,prof))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     profiling
C
         elseif (inp_compare(.false.,test,'matrix_checks')) then
            matrix_checks=.true.

            if (.not.rtdb_put(rtdb,'rt_semiemp:matrix_checks',
     $           mt_log,1,matrix_checks))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C     
C     Override DFT settings so that there is *absolutely* no I/O.  If
C     "usedisk" is present, we will use DFT rtdb options, otherwise we
C     disable everything (default).
C     
         elseif (inp_compare(.false.,test,'usedisk')) then
            call errquit (pname//"usedisk disabled",0,0)

c$$$            nodisk=.false.
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:nodisk',mt_log,1,nodisk))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
            

C
C     post process using dplot? (i.e. make density plots)
C
c$$$         elseif (inp_compare(.false.,test,'dplot')) then
c$$$            dplot_do=.true.
c$$$            
c$$$            if (inp_a(dplot_opts)) then
c$$$               if (trim(dplot_opts).eq."subgs") then
c$$$                  subgs = .true.
c$$$                  if (.not.rtdb_put(rtdb,'rt_tddft:subgs',
c$$$     $                 mt_log,1,subgs))
c$$$     $                 call errquit(pname//'Write failed to rtdb',
c$$$     $                 0,RTDB_ERR)
c$$$               else
c$$$                  call errquit (pname//"invalid dplot option: "
c$$$     $                 //trim(dplot_opts), 0, 0)
c$$$               endif
c$$$            endif
c$$$
c$$$            if (.not.rtdb_put(rtdb,'rt_tddft:dplot', mt_log,1,dplot_do))
c$$$     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)
c$$$
            
            
C
C     noprop directive (dont propagate)
C
         elseif (inp_compare(.false.,test,'noprop')) then
            noprop=.true.
            if (.not.rtdb_put(rtdb,'rt_semiemp:noprop',mt_log,1,noprop))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     static directive (never rebuild Fock matrix)
C
         elseif (inp_compare(.false.,test,'static')) then
            static = .true.
            if (.not.rtdb_put(rtdb,'rt_semiemp:static',mt_log,1,static))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

            
            
C
C     Rules for applying fields to geometries
C
         elseif (inp_compare(.false.,test,"excite")) then
            nexcites = nexcites + 1
            
            if (.not.rtdb_put(rtdb,'rt_semiemp:nexcites',
     $           mt_int,1,nexcites))
     $           call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)

c           call rt_semiemp_input_excite (rtdb, nexcites, nfields)


C
C     List of time-dependent system properties to print.
C
         elseif (inp_compare(.false.,test,'print')) then
            call rt_semiemp_input_print (rtdb)

C     
C     end
C     
         else if (inp_compare(.false.,test,'end')) then
            done = .true.
         else
            call errquit(pname//'Unknown directive: '//trim(test),
     $           0, INPUT_ERR)
         endif
      enddo  !end main parsing loop

      

C      if (.not.got_tmin) call errquit(pname//"must supply tmin",0,0)
      if (.not.got_tmax) call errquit(pname//"must supply tmax",0,0)
      if (.not.got_dt) call errquit(pname//"must supply dt",0,0)


C     HARDCODED TMIN
C      tmin = -5.0*dt
      tmin = 0d0


      if (.not.rtdb_put(rtdb,'rt_semiemp:tmin',mt_dbl,1,tmin))
     $     call errquit(pname//'Write failed to rtdb',0,RTDB_ERR)


C
C     Check that parameters are valid and compatible with each other.
C     
C      if (tmin.lt.0d0) call errquit(pname//"tmin must be > 0", 0, 0)
      if (tmax.lt.0d0) call errquit(pname//"tmax must be > 0", 0, 0)
      if (tmax.lt.tmin) call errquit(pname//"tmax must be > tmin", 0, 0)
      if (dt.lt.0d0) call errquit(pname//"dt must be > 0", 0, 0)
      
      end subroutine
      
c $Id: rt_tddft_input.F 28240 2016-05-22 00:18:37Z klopata $
