      logical function semiemp_energy(rtdb)
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
#include "bas.fh"
#include "schwarz.fh"
#include "dftps.fh"
#include "bgj.fh"
c

      integer rtdb             ! Runtime database handle
      character*20 test
      character*20 scftype
      character*20 apx
      character*30 operation
      logical rhf,uhf
      logical cndo,indo
      integer itmax
      integer inttyp
      integer gfac
      integer g_F
      double precision  charge
      double precision mult
      double precision scftol,E1
      integer infa(5)
      double precision ME
      double precision total_wall, total_cpu
      logical numerical 
      integer III,L
      integer IA
      integer MIM(2) ! NUM. OF ORBITALS IN EACH SHELL I.
      double precision FOP(31) ! THE NUM. OF ELECTRONS IN EACH SHELL I. CLOSED
                      ! SHELL FIRST  THEN GREATEST OPENSHELL FIRST, ETC.
      character*8 atom_type
      integer ihomo(505),ilumo(505)
      integer atom_num
      integer iorb, aorb1, aorb2
      integer lnin,inin !NIN(6000)
      integer iddmat,lddmat,ifdmat,lfdmat,ifdmat2,lfdmat2,lDK,iDK
      integer id12,ld12
      integer JSLIT(4)
      integer lQQ,iQQ,iFZ,lFZ
      double precision sexpo(105,6) !s orbital exponent
      double precision pexpo(105,6) !p orbital exponenet
      double precision dexpo(105,6) !d orbital exponent
      double precision fexpo(105,6) !f orbital exponent
      double precision UCORE(4,86) ! for AM1,PM3 and MNDO
      double precision  EIS(105),EIP(105),EID(105),EIF(105) !for beta 
c     double precision PO(9,105), AALP(105), DAD(9,105)
c     double precision NB(nat)      ! 
      integer lPO,iPO,lAALP,iAALP,iDAD,lDAD
      integer nat,kan       !number of atoms and atomic numbers 
      double precision q  !charges
      logical geom_rtdb_ncent  !subroutine to get numbet of atoms
      logical status
      integer i,ii,j,ielem,iexpo,k,kk,X
      double precision c  !coordinates
      integer geom    !geometry
      integer l_xyz, i_xyz
      integer iH,lH  !Overlap matrix
      integer iSH, lSH !Weighted overlap matrix
      integer iW,lW
      integer iFA,lFA
      integer iFT,lFT !for fock transformation matrix
      integer iF,lF
      integer iBB,lBB,iAA,lAA,ibkb,ibka,lbka,lbkb
      integer iETS,lETS,iETP,lETP,iETD,lETD,iETF,lETF
      integer lgrad, igrad,irho,lrho,irhoa,lrhoa,irhob,lrhob
      integer iRA,lRA,iRB,lRB
      integer ieig,leig !Eigen value matrix
      integer iwcor,lwcor  !for wcor matrix
      integer lz, iz
      integer lcharge, icharge
      integer ltags, itags
c     logical  geom_cart_get2 !for getting coordinate and atomic number
      integer geometry, iat
      integer atnum       !atomic numbers
      integer NEL  !number of total electrons
      integer NCORE
      double precision  valence  !valence electrons
      integer lval,ival       !valence electrons
      double precision cx,cy,cz !coordinates
      double precision energy
      double precision GSS(86),GPP(105),GSP(105),GP2(86),HSP(86),HPP(86)
      integer  ntype
      integer lntype,intype   !type of orbitals in an atom
      integer lnqm,inqm    !n quantum number
      integer llqm,ilqm    !l quantum number
      integer lnb,inb      !number of basis for each atom
      integer lnu,inu      !type of atom array
      integer lalpha,ialpha  !Exponent for basis function
      logical orbital_qmnumber  !to get the qm number fo a orbital
      integer lcor,icor   !for core electrons
      integer inp,lnp     !principle qunatum number
      integer lpqn,ipqn   !Principle quantum number array
      integer lban,iban   !atomic number array for basis function
      integer lng, ing    !Pointer for gama
      integer lnw,inw     !type of basis function array
      integer lgama,igama !two electron integral
      integer lbeta,ibeta,lbetb,ibetb !beta values
      integer iLABI,iLABJ,iLABK,iLABL,lLABI,lLABJ,lLABK,lLABL
      integer idipx,idipy,idipz,ldipx,ldipy,ldipz, lnfirst,infirst
      integer iatnum, latnum
      integer iaso,laso
      integer iatom, latom
      double precision CORE
      logical core_energy
      logical basis_function
      double precision REP
      double precision nrep
      integer nbasis !number of basis functions for an atom
      integer tot_basis,nbf ! Total number of basis functions
      integer isumelel,lsumelel
      integer lbtype,ibtype !basis type
      integer lrad,irad !interatomic distance
      logical radmat
      integer ieul,leul !Eulerian transformation
      integer K0,K1,K4,K9,K16 !basis function type
      integer NGDD       !number of basis function type
      integer NRR
      integer NF       !size of F matrix
      integer LL !number of atomic integral
      integer ISW1,ISW2,ISW3,ISWANN
      integer NDN
      double precision COFSQ2(105,5), E(3), PME(3,505)
      integer iSCG1SP,iSCF2PP,
     1 iSCG2SD,iSCG1PD,iSCF2PD,iSCG3PD,iSCF2DD,iSCF4DD,iSCG3SF,iSCG2PF,
     2 iSCG4PF,iSCG1DF,iSCG3DF,iSCG5DF,iSCF2PF,iSCF2DF,iSCF4DF,iSCF2FF,
     3 iSCF4FF,iSCF6FF,iR1SPPD,iR2SDDD,iR2SDPP,iR1PDDF,iR1SPDF,iR2PFDD,
     4 iR2PFFF,iR2PPPF,iR2SDFF,iR2SDPF,iR3PDDF,iR3SFDF,iR3SFPD,iR4PFDD,
     5 iR4PFFF
      integer lSCG1SP,lSCF2PP,
     1 lSCG2SD,lSCG1PD,lSCF2PD,lSCG3PD,lSCF2DD,lSCF4DD,lSCG3SF,lSCG2PF,
     2 lSCG4PF,lSCG1DF,lSCG3DF,lSCG5DF,lSCF2PF,lSCF2DF,lSCF4DF,lSCF2FF,
     3 lSCF4FF,lSCF6FF,lR1SPPD,lR2SDDD,lR2SDPP,lR1PDDF,lR1SPDF,lR2PFDD,
     4 lR2PFFF,lR2PPPF,lR2SDFF,lR2SDPF,lR3PDDF,lR3SFDF,lR3SFPD,lR4PFDD,
     5 lR4PFFF
       integer g_S, g_alpha, g_A, g_c
       logical rt_semiemp
       external rt_semiemp
       logical docis
c    
       total_cpu = 0.0d0
       total_wall = 0.0d0 
c      if (ga_nodeid().eq.0) 
c    .      write(luout,*) " i am in semiemp_energy"
c
       semiemp_energy=.true. 
c
       if (.not.rtdb_put(rtdb,'semiemp:energy',mt_dbl,1,0.0d0))
     1  call errquit('semiemp: failed to write energy',0, RTDB_ERR)

c      Create geometry object and load in geometry from rtdb
       if (.not.geom_create(geom,'geometry'))
     1  call errquit('semiemp_energy: geom_create failed',0, GEOM_ERR)
       if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     1  call errquit('semiemp_energy: no geometry found',0, RTDB_ERR)
        if(.not.geom_print(geom)) call errquit(
     .       'geom print failed',0, GEOM_ERR) 
c      Get the atomic information from the geometry 
      if (.not. geom_ncent(geom, nat))
     &   call errquit('semiemp_energy: geom_ncent failed',geom,GEOM_ERR)
      if (.not.MA_push_get(MT_Dbl,nat*3,'xyz',l_xyz,i_xyz))
     &   call errquit('semiemp_energy: cannot allocate xyz',0, MA_ERR)
      if (.not.MA_Push_Get(MT_int,nat,'atnums',lz,iz))
     &   call errquit('semiemp_energy: cannot allocate atns',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nat,'charge',lcharge,icharge))
     &   call errquit('semiemp_enrgy: cannot allocate charge',0, MA_ERR)
       if (.not.MA_Push_Get(MT_Dbl,nat,'valence',lval,ival))
     &   call errquit('semiemp_energy: cannot allocate valence',0, 
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_int,nat,'NTYPE',lntype,intype))
     &   call errquit('semiemp_energy: cannot allocate NTYPE',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'CORE',lcor,icor))
     &   call errquit('semiemp_energy: cannot allocate CORE',0,
     1    MA_ERR)
               if (.not.MA_Push_Get(MT_int,nat,'NB',lnb,inb))
     &   call errquit('semiemp_energy: cannot allocate NB',0,
     1    MA_ERR)
               if (.not.MA_Push_Get(MT_int,nat,'NP',lnp,inp))
     &   call errquit('semiemp_energy: cannot allocate NP',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nat*nat,'rad',lrad,irad))
     &   call errquit('semiemp_energy: cannot allocate RAD',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'DK',lDK,iDK))
     &   call errquit('semiemp_energy: cannot allocate T',0,
     1    MA_ERR)
      if (.not.MA_Push_Get(MT_Byte,nat*16,'tags',ltags,itags))
     &   call errquit('semiemp_energy: cannot allocate tags',0, MA_ERR)

c      Initialize semiempirical parameters
       call semiemp_init(rtdb)

c

C   ISW1 specifies model Hamiltonian:   (sw_hamiltonian)
C
C   ISW1  =   0      1       2      3      4      5      6     7
C            EHT   CNDO/1 CNDO/2  INDO/1 INDO/2 NDDO/1 NDDO/2 PPP
C         =  51     52      53
C           MNDO    AM1    PM3
       if (.not.rtdb_cget(rtdb,'semiemp:apx',1,apx)) then
       call errquit('oneint_semiemp: failed to read apx',0, RTDB_ERR)
       endif
       if (inp_compare(.false.,apx,'EHT'))  then
        ISW1 = 0
       elseif (inp_compare(.false.,apx,'cndo/1'))  then
        ISW1 = 1
       elseif (inp_compare(.false.,apx,'cndo/2'))  then
        ISW1 = 2
       elseif (inp_compare(.false.,apx,'indo/1'))  then
        ISW1 = 3
       elseif (inp_compare(.false.,apx,'indo/2'))  then
        ISW1 = 4
       elseif (inp_compare(.false.,apx,'nndo/1'))  then
        ISW1 = 5
       elseif (inp_compare(.false.,apx,'nndo/2'))  then
        ISW1 = 6
       elseif (inp_compare(.false.,apx,'ppp'))  then
        ISW1 = 7
       elseif (inp_compare(.false.,apx,'mndo'))  then
        ISW1 = 51
       elseif (inp_compare(.false.,apx,'am1'))  then
        ISW1 = 52
       elseif (inp_compare(.false.,apx,'pm3'))  then
        ISW1 = 53
      end if 
        ISW2 = 4
        ISW3 = 0 
      if (.not.rtdb_get(rtdb,'semiemp:scftol',mt_dbl,1,scftol))
     &  call errquit('semiemp failed read to scftol',0,RTDB_ERR)
c     end if
         if (.not.rtdb_get(rtdb,'semiemp:inttyp',mt_int,1,inttyp))
     &  call errquit('semiemp failed read to inttyp',0,RTDB_ERR)
      if (.not.rtdb_get(rtdb,'semiemp:maxiter',mt_int,1,itmax))
     &  call errquit('scfit failed read to itmax',0,RTDB_ERR)
c      end if
c     write(6,*)"itmax",itmax

      if (.not.rtdb_get(rtdb,'semiemp:mult',mt_dbl,1,mult))
     &  call errquit('scfit failed read to mult',0,RTDB_ERR)
c      end if
       if (.not.rtdb_put(rtdb,'semiemp:natom',mt_int,1,nat))
     1  call errquit('semiemp: failed to write natom',0, RTDB_ERR)


C   ISWANN = 0  for RHF
C          = 1  for UHF
C          = 2  for UHF + projection
C          = 3  for ROHF
       if (.not.rtdb_cget(rtdb,'semiemp:scftype',1,scftype)) then
      call errquit('semiemp_energy: failed to read scftype',0, RTDB_ERR)
       endif
c      write(6,*) "SCFTYPE", scftype
       if (inp_compare(.false.,scftype,'rhf'))  then
        ISWANN = 0
       elseif (inp_compare(.false.,scftype,'uhf'))  then
        ISWANN = 1
       elseif (inp_compare(.false.,scftype,'rohf'))  then
        ISWANN = 3
       end if
      NDN = (ISWANN + 1) / 2 + 1
c
      CALL SLITE(0,JSLIT)
c     write(6,*)"JSLIT",JSLIT
      if (.not. geom_cart_get2(geom, nat, Byte_MB(itags),
     &        Dbl_MB(i_xyz), Dbl_MB(icharge), int_mb(iz)))
     &  call errquit('semiemp_energy: geom_cart_get failed',74, 
     1  GEOM_ERR)

c    calculating number of  valence electrons
         call core_valence_charge(nat,int_mb(iz),Dbl_MB(ival),
     1    Dbl_MB(icor),int_mb(intype),dbl_mb(iDK))
      if (.not. rtdb_get(rtdb, 'charge', MT_DBL, 1, charge))then
        charge = 0.0d0
      endif

c    make basis function array
       status= basis_function(nat,int_mb(iz),dbl_mb(ival),int_mb(inb),
     1  int_mb(inp),int_mb(intype))
       nbf= 0
       cx = 0.0
       cy = 0.0
       cz = 0.0
       K0=0
       K1=0
       K4=0
       K9=0
       K16=0
       NEL=0
       NCORE=0
       CORE=0
       valence = 0.0
       do iat = 0, nat-1
         valence =valence +dbl_mb(ival+iat)
        end do
c       write(6,*)"Number of valence electrons =", valence
         if (.not.MA_Push_Get(MT_int,INT(valence),'nqm',lnqm,inqm))
     &   call errquit('semiemp_energy: cannot allocate nqm',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_int,INT(valence),'lqm',llqm,ilqm))
     &   call errquit('semiemp_energy: cannot allocate lqm',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_int,INT(valence),'basis',lbtype,
     1   ibtype))
     &   call errquit('semiemp_energy: cannot allocate BTYPE',0,
     1    MA_ERR)
      do iat = 0,nat-1
        atnum =int_mb(iz+iat)
        cx =  dbl_mb(i_xyz+3*iat )
        cy =  dbl_mb(i_xyz+3*iat + 1 )
        cz =  dbl_mb(i_xyz+3*iat + 2 )
        valence= dbl_mb(ival+iat)
        CORE= dbl_mb(icor+iat)
        NCORE = NCORE + CORE/2
        ntype= int_mb(intype+iat)
c  populating n and l quantum numbers for different electrons for the
c  atoms
      status= orbital_qmnumber(valence,ntype,atnum, int_mb(inqm),
     1 int_mb(ilqm), nbasis,int_mb(ibtype))
      if (ga_nodeid().eq.0) then
c         write(luout,*) "cx: ", cx
c         write(luout,*) "cy: ", cy
c         write(luout,*) "cz: ", cz
c         write(luout,*) "atnum: ",atnum
          NEL=NEL+valence
c         write(luout,*) "valence electron",valence
c         write(6,*) "core electrons", CORE
c         write(luout,*) "atom type", ntype
          IF(ntype .LE. 0) THEN
          K0 = K0 + 1
          ELSEIF(ntype .LE. 1) THEN
          K1 = K1 + 1
          ELSEIF(ntype .LE. 2) THEN
          K4 = K4 + 1
          ELSEIF(ntype .LE. 4) THEN
          K9 = K9 + 1
          ELSE
          K16 = K16 + 1
          ENDIF
c         write(luout,*) "basis fn", int_mb(inb+iat)
c         do j=0, valence-1
c         write(luout,*) "n qm",int_mb(inp+iat)
c         write(luout,*) "l qm",int_mb(ilqm+j)
c         write(luout,*) "orbital_sym"    
c         end do  
c         write(luout,*)
        end if
        nbf = nbf+ int_mb(inb+iat)
      enddo
c      if (ga_nodeid().eq.0)
c      write(luout,*) "total number of basis fn", nbf
c      write(6,*)"charge",charge

       IF(ISWANN .NE. 3) THEN
       NEL=NEL- INT(charge)
       MIM(2) = (NEL + MULT- 1 ) / 2
       MIM(1) = MIM(2) - MULT + 1
       FOP(2) = REAL(MIM(2))
       FOP(1) = REAL(MIM(1))
c      write(6,*)"MIM",MIM(1),MIM(2)
c      write(6,*)"NEL=",NEL
       IF(ISWANN .EQ. 0) FOP(1) = REAL(NEL)

        if (.not.rtdb_put(rtdb,'semiemp:MIM',mt_int,2,MIM))
     1  call errquit('semiemp: failed to write MIM',0, RTDB_ERR)

c     ELSE
c      CALL RHFST(MULT)
      ENDIF

      if (ga_nodeid().eq.0)  then
       write(luout,*) "APPROXIMATE HAMILTONIAN", "=",  apx
       write(luout,*)"SCFTOL = ",scftol
       write(luout,*)"MULT = ",mult
       write(luout,*) "Total Number of Basis Functions =", nbf
      endif
c      write(6,*)"MIM",MIM        
c        write(6,*)"Ks"
c        write(6,*)K0,K1,K4,K9,K16
           NGDD = K1 + 2*K4 + 3*K9 + 4*K16
           NRR = (nbf*(nbf+1)/2)
           NF = MAX(3 * nat*nat, (nbf*(nbf + 1)) / 2)

       if (.not.ma_pop_stack(lbtype))
     &  call errquit('semiemp_energy: cannot pop stack btype',3,MA_ERR)
       if (.not.ma_pop_stack(llqm))
     &   call errquit('semiemp_energy:cannot pop stack lqm',3, MA_ERR)
       if (.not.ma_pop_stack(lnqm))
     &   call errquit('semiemp_energy: cannot pop stack nqm',3, MA_ERR)

       if (.not.rtdb_put(rtdb,'semiemp:nbf',mt_int,1,nbf))
     1  call errquit('semiemp: failed to write nbf',0, RTDB_ERR)

       if (.not.MA_Push_Get(MT_int,nbf,'lqm',llqm,ilqm))
     &   call errquit('semiemp_energy: cannot allocate nqm',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,3*nat*nat,'eul',leul,ieul))
     &   call errquit('semiemp_energy: cannot allocate EUL',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Int,nbf,'NW',lnw,inw))
     &   call errquit('semiemp_energy: cannot allocate NW',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Int,nbf,'NU',lnu,inu))
     &   call errquit('semiemp_energy: cannot allocate NU',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Int,nbf,'N',lpqn,ipqn))
     &   call errquit('semiemp_energy: cannot allocate N',0,
     1    MA_ERR)
               if (.not.MA_Push_Get(MT_Int,nbf,'BAN',lban,iban))
     &   call errquit('semiemp_energy: cannot allocate BAN',0,
     1    MA_ERR)
               if (.not.MA_Push_Get(MT_Int,nbf,'NG',lng,ing))
     &   call errquit('semiemp_energy: cannot allocate NG',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'GAMA',
     1    lgama,igama))
     &   call errquit('semiemp_energy: cannot allocate GAMA',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*6,'Alpha',lalpha,ialpha))
     &   call errquit('semiemp_energy: cannot allocate ALPHA',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,105*10,'BETA',lbeta,ibeta))
     &   call errquit('semiemp_energy: cannot allocate BETA',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,105*4,'BETB',lbetb,ibetb))
     &   call errquit('semiemp_energy: cannot allocate BETB',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'H',lH,iH))
     &   call errquit('semiemp_energy: cannot allocate H',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'SH',lSH,iSH))
     &   call errquit('semiemp_energy: cannot allocate SH',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'FA',lFA,iFA))
     &   call errquit('semiemp_energy: cannot allocate FA',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'FT',lFT,iFT))
     &   call errquit('semiemp_energy: cannot allocate FT',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf,'EIG',leig,ieig))
     &   call errquit('semiemp_energy: cannot allocate EIG',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'WCOR',lwcor,iwcor))
     &   call errquit('semiemp_energy: cannot allocate WCOR',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_int,nbf*nbf,'NIN',lnin,inin))
     &   call errquit('semiemp_energy: cannot allocate NIN',0, MA_ERR)
c------------------------------------------------------------------------------
c     integrals
c------------------------------------------------------------------------------
             if (.not.MA_Push_Get(MT_Dbl,nat,'SCG1SP',lscg1sp,iscg1sp))
     &   call errquit('semiemp_energy: cannot allocate SCG1SP',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2PP',lscf2pp,iscf2pp))
     &   call errquit('semiemp_energy: cannot allocate SCF2PP',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCG2SD',lscg2sd,iscg2sd))
     &   call errquit('semiemp_energy: cannot allocate SCG2SD',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCG1PD',lscg1pd,iscg1pd))
     &   call errquit('semiemp_energy: cannot allocate SCF1PD',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2PD',lscf2pd,iscf2pd))
     &   call errquit('semiemp_energy: cannot allocate SCF2PD',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCG3PD',lscg3pd,iscg3pd))
     &   call errquit('semiemp_energy: cannot allocate SCG3PD',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2DD',lscf2dd,iscf2dd))
     &   call errquit('semiemp_energy: cannot allocate SCF2DD',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF4DD',lscf4dd,iscf4dd))
     &   call errquit('semiemp_energy: cannot allocate SCF4DD',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCG3SF',lscg3sf,iscg3sf))
     &   call errquit('semiemp_energy: cannot allocate SCG3SF',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCG2PF',lscg2pf,iscg2pf))
     &   call errquit('semiemp_energy: cannot allocate SCG2PF',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'SCG4PF',lscg4pf,iscg4pf))
     &   call errquit('semiemp_energy: cannot allocate SCG4PF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCG1DF',lscg1df,iscg1df))
     &   call errquit('semiemp_energy: cannot allocate SCG1DF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCG3DF',lscg3df,iscg3df))
     &   call errquit('semiemp_energy: cannot allocate SCG3DF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCG5DF',lscg5df,iscg5df))
     &   call errquit('semiemp_energy: cannot allocate SCG5DF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2PF',lscf2pf,iscf2pf))
     &   call errquit('semiemp_energy: cannot allocate SCF2PF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2DF',lscf2df,iscf2df))
     &   call errquit('semiemp_energy: cannot allocate SCF2DF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF4DF',lscf4df,iscf4df))
     &   call errquit('semiemp_energy: cannot allocate SCF4DF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF2FF',lscf2ff,iscf2ff))
     &   call errquit('semiemp_energy: cannot allocate SCF2FF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF4FF',lscf4ff,iscf4ff))
     &   call errquit('semiemp_energy: cannot allocate SCF4FF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'SCF6FF',lscf6ff,iscf6ff))
     &   call errquit('semiemp_energy: cannot allocate SCF6FF',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'R1SPPD',lr1sppd,ir1sppd))
     &   call errquit('semiemp_energy: cannot allocate R1SPPD',0,
     1    MA_ERR)
          if (.not.MA_Push_Get(MT_Dbl,nat,'R2SDDD',lr2sddd,ir2sddd))
     &   call errquit('semiemp_energy: cannot allocate R2SDDD',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R2SDPP',lr2sdpp,ir2sdpp))
     &   call errquit('semiemp_energy: cannot allocate R2SDPP',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R1PDDF',lr1pddf,ir1pddf))
     &   call errquit('semiemp_energy: cannot allocate R1PDDF',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R1SPDF',lr1spdf,ir1spdf))
     &   call errquit('semiemp_energy: cannot allocate R1SPDF',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R2PFDD',lr2pfdd,ir2pfdd))
     &   call errquit('semiemp_energy: cannot allocate R2PFDD',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R2PFFF',lr2pfff,ir2pfff))
     &   call errquit('semiemp_energy: cannot allocate R2PFFF',0,
     1    MA_ERR)
         if (.not.MA_Push_Get(MT_Dbl,nat,'R2PPPF',lr2pppf,ir2pppf))
     &   call errquit('semiemp_energy: cannot allocate R2PPPF',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R2SDFF',lr2sdff,ir2sdff))
     &   call errquit('semiemp_energy: cannot allocate R2SDFF',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R2SDPF',lr2sdpf,ir2sdpf))
     &   call errquit('semiemp_energy: cannot allocate R2SDPF',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R3PDDF',lr3pddf,ir3pddf))
     &   call errquit('semiemp_energy: cannot allocate R3PDDF',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R3SFDF',lr3sfdf,ir3sfdf))
     &   call errquit('semiemp_energy: cannot allocate R3SFDF',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R3SFPD',lr3sfpd,ir3sfpd))
     &   call errquit('semiemp_energy: cannot allocate R3SFPD',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R4PFDD',lr4pfdd,ir4pfdd))
     &   call errquit('semiemp_energy: cannot allocate R4PFDD',0,
     1    MA_ERR)
        if (.not.MA_Push_Get(MT_Dbl,nat,'R4PFFF',lr4pfff,ir4pfff))
     &   call errquit('semiemp_energy: cannot allocate R4PFFF',0,
     1    MA_ERR)
C--------------------------------------------------------------------------------------
C
C-------------------------------------------------------------------------------------
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'F',lF,iF))
     &   call errquit('semiemp_energy: cannot allocate F',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'BB',lBB,
     1    iBB))
     &   call errquit('semiemp_energy: cannot allocate BB',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,3*nat,'FZ',lFZ,iFZ))
     &   call errquit('semiemp_energy: cannot allocate FZ',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf,'SUMELEL',lsumelel,
     1         isumelel))
     &   call errquit('semiemp_energy: cannot allocate SUMELEL',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf,'DDMAT',lddmat,iddmat))
     &   call errquit('semiemp_energy: cannot allocate D',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf,'fdmat',lfdmat,
     1    ifdmat))
     &   call errquit('semiemp_energy: cannot allocate E',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf,'fdmat2',lfdmat2,
     1    ifdmat2))
     &   call errquit('semiemp_energy: cannot allocate E',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nat,'D12',ld12,id12))
     &   call errquit('semiemp_energy: cannot allocate T',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'AA',lAA,
     1    iAA))
     &   call errquit('semiemp_energy: cannot allocate AA',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf,'BKA',lbka,
     1    ibka))
     &   call errquit('semiemp_energy: cannot allocate BKA',0,
     1    MA_ERR)
              if (.not.MA_Push_Get(MT_Dbl,nbf,'BKB',lbkb,
     1    ibkb))
     &   call errquit('semiemp_energy: cannot allocate BKB',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'RA',lRA,iRA))
     &   call errquit('semiemp_energy: cannot allocate RA',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'RB',lRB,iRB))
     &   call errquit('semiemp_energy: cannot allocate RB',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nat,'ETS',lETS,iETS))
     &   call errquit('semiemp_energy: cannot allocate RB',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nat,'ETP',lETP,iETP))
     &   call errquit('semiemp_energy: cannot allocate RB',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nat,'ETD',lETD,iETD))
     &   call errquit('semiemp_energy: cannot allocate RB',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,nat,'ETF',lETF,iETF))
     &   call errquit('semiemp_energy: cannot allocate ETF',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,3*nat,'GRAD',lgrad,igrad))
     &   call errquit('semiemp_energy: cannot allocate GRAD',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'Rho',lrho,irho))
     &   call errquit('semiemp_energy: cannot allocate Rho',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'FX',ldipx,idipx))
     &   call errquit('semiemp_energy: cannot allocate FX',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'FY',ldipy,idipy))
     &   call errquit('semiemp_energy: cannot allocate PRA',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,nbf*nbf,'FZ',ldipz,idipz))
     &   call errquit('semiemp_energy: cannot allocate PRB',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_INT,nat,'NFIRST',lnfirst,infirst))
     &   call errquit('semiemp_energy: cannot allocate nfirst',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,9*105,'PO',lPO,iPO))
     &   call errquit('semiemp_energy: cannot allocate PO',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_byte,nbf*8,'ASO',laso,iaso))
     &   call errquit('semiemp_energy: cannot allocate ASO',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_byte,nbf*8,'ATOM',latom,iatom))
     &   call errquit('semiemp_energy: cannot allocate ATOM',0, MA_ERR)
             if (.not.MA_Push_Get(MT_Int,nbf,'ATOM_NUM',latnum,iatnum))
     &   call errquit('semiemp_energy: cannot allocate ATOM_NUM',0,
     1    MA_ERR)
c-----------------------------------------------------------------------------------
c
       if (.not.rtdb_put(rtdb,'semiemp:DK',mt_dbl,nat,dbl_mb(iDK)))
     1  call errquit('semiemp: failed to write DK',0, RTDB_ERR)

C----------------------------------------------------------------------------------
C
C----------------------------------------------------------------------------------       
           IF (ISW1.GT.50) THEN
           if (.not.MA_Push_Get(MT_INT,nbf*nbf*nbf,'LABI',lLABI,
     1   iLABI))
     &   call errquit('semiemp_energy: cannot allocate LABI',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,nbf*nbf*nbf,'LABJ',lLABJ,
     1   iLABJ))
     &   call errquit('semiemp_energy: cannot allocate LABJ',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,nbf*nbf*nbf,'LABK',lLABK,
     1   iLABK))
     &   call errquit('semiemp_energy: cannot allocate LABK',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_INT,nbf*nbf*nbf,'LABL',lLABL,
     1   iLABL))
     &   call errquit('semiemp_energy: cannot allocate LABL',0,
     1    MA_ERR)
           if (.not.MA_Push_Get(MT_Dbl,nbf*nbf*nbf,'QQ',lQQ,iQQ))
     &   call errquit('semiemp_energy: cannot allocate QQ',0,
     1    MA_ERR)
            if (.not.MA_Push_Get(MT_Dbl,9*105,'DAD',lDAD,iDAD))
     &   call errquit('semiemp_energy: cannot allocate DAD',0,
     1    MA_ERR)
             if (.not.MA_Push_Get(MT_Dbl,105,'AALP',lAALP,iAALP))
     &   call errquit('semiemp_energy: cannot allocate AALP',0,
     1    MA_ERR)
          ENDIF

c      write(luout,*)"memory allocation done"

c
c
c     Initialize exponents
      do ielem = 1,105
           EIS(ielem)= 0.0d0
           EIP(ielem)= 0.0d0
           EID(ielem)= 0.0d0
           EIF(ielem)= 0.0d0
        do iexpo = 1,6   
           sexpo(ielem,iexpo) = 0.d0
           pexpo(ielem,iexpo) = 0.d0
           dexpo(ielem,iexpo) = 0.d0
           fexpo(ielem,iexpo) = 0.d0
        end do 
       end do
C-----------------------------------------------------------------------
      DO  I=0,nbf*nbf-1
       Dbl_MB(iFA+I) = 0.0d0
       Dbl_MB(iBB+I) = 0.0d0
       Dbl_MB(iF+I) = 0.0d0
       Dbl_MB(igama+I) = 0.0d0
      END DO

c       write(luout,*)"initialization"

c
c      Get atomic data: basis functions
       call atom_data(int_mb(iz),nat,sexpo,pexpo,dexpo,fexpo) 
c---------------------------------------------------------------       
c
        L=0
        do I=0,nbf*nbf-1
           INT_MB(iNIN+I)= L
c        write(6,*)"NIN"!,INT_MB(iNIN+I)
         L=L+I+1
        end do  

       call expjam(int_mb(iz),nat,nbf,int_mb(inb),sexpo,pexpo,
     1   dexpo,fexpo,Dbl_MB(ialpha),Int_MB(inw),Int_MB(inu),Int_MB(inp),
     2   Int_MB(ipqn),Int_MB(iban),int_mb(ilqm),INT_MB(infirst),ISW1,
     3   INT_MB(intype),byte_mb(iaso),byte_mb(iatom),int_mb(iatnum))

       if (.not.rtdb_put(rtdb,'semiemp:aso',mt_byte,nbf*8,
     1  byte_mb(iaso)))
     1  call errquit('semiemp: failed to write ASO',0, RTDB_ERR)

       if (.not.rtdb_put(rtdb,'semiemp:atom',mt_byte,nbf*8,
     1  byte_mb(iatom)))
     1  call errquit('semiemp: failed to write ATOM',0, RTDB_ERR)
       if (.not.rtdb_put(rtdb,'semiemp:numbf',mt_int,nat,int_mb(inb)))
     1  call errquit('semiemp: failed to write numbf',0, RTDB_ERR)


       IF (ISW1.GT.50.AND.ISW1.LT.54)
     1  CALL AMDATA(UCORE,ISW1,NAT,NBF,Dbl_MB(ialpha),INT_MB(ilqm),
     1  INT_MB(inu),INT_MB(iban),Dbl_MB(iPO),Dbl_MB(iAALP),Dbl_MB(iDAD)
     1 ,GSS,GPP,GSP,GP2,HSP,HPP)

c      Put the energy into the rtdb
c      calculating core energy
       nrep = 0.0
       rep = 0.0
       status= geom_nuc_rep_energy( geom, nrep)
       status= core_energy(nat,Dbl_MB(ival),Dbl_MB(i_xyz),REP)

c      calculating radius matrix
       status= radmat(nat,Dbl_MB(i_xyz),Dbl_MB(irad))
c----------------------------------------------------------------------------------------------
c     constructing various Eulerian transformation matrix
c----------------------------------------------------------------------------------------------
      call geomat(nat,Dbl_MB(ieul),Dbl_MB(irad),Dbl_MB(i_xyz))     
c---------------------------------------------------------------------------------------------
c      calculating two electron integrals
c------------------------------------------------------------------------------------------
       call twoint_semiemp(nat,nbf,int_mb(inb),Dbl_MB(ialpha),
     1 Dbl_MB(irad),Dbl_MB(igama),Int_MB(inw),Int_MB(ipqn),Int_MB(iban),
     2 Int_MB(inu),Int_MB(ing),ISW1,inttyp,Dbl_MB(iPO),INT_MB(iNIN))

       IF ((ISW1.GT.50.AND.ISW1.LT.54).OR.INTTYP.ne.0)
     1 CALL NUCREP(nat,Dbl_MB(irad),INT_MB(iz),Dbl_MB(ival),REP,
     1 Dbl_MB(iPO),Dbl_MB(iAALP),ISW1,INT_MB(inb),INT_MB(ing),
     1 Dbl_MB(igama),inttyp)
       if (ga_nodeid().eq.0) then
       write(luout,*) "NUCLEAR REPULSION ENERGY = ", nrep
       write(luout,*) "CORE REPULSION ENERGY = ", REP
       endif

c-------------------------------------------------------------------------------------------
c    getting beta values
c-------------------------------------------------------------------------------------------
       call para(nat,int_mb(iz),int_mb(intype),EIS,EIP,EID,EIF,ISW1 
     1  ,ISW2,COFSQ2)
        IBETA = 0
       IF(ISW1.GT.50.AND.ISW1.LT.54) GO TO 10
c       IF (IBETA.EQ.6) THEN
c      BETA
c this has to be filled in later      
       call paramb(nat,int_mb(iz),int_mb(intype),Dbl_MB(ibeta),inttyp)
       call paramd(Dbl_MB(ibetb))
   10  Continue
c----------------------------------------------------------------------------------------------
c      calculating one electron integrals  
c---------------------------------------------------------------------------------------------
C-----------------------------------------------------------------------
c     total_wall = -util_wallsec()
c     total_cpu  = -util_cpusec()
c     constructing the one electron matrix element of the Hamiltonina   
       call oneint_semiemp(rtdb,nat,nbf,Dbl_MB(ialpha),Int_Mb(ipqn),
     1 Int_MB(inw),int_mb(iban),Int_MB(inu),Dbl_MB(irad),Dbl_MB(ieul),
     2 Dbl_MB(iFA),Dbl_MB(iH),Dbl_MB(iSH),INT_MB(inin),ISW1,INTTYP,
     3 Dbl_MB(iFT))
c      write(6,*)"after oneint_semiemp"
c     total_wall = total_wall + util_wallsec()
c     total_cpu  = total_cpu  + util_cpusec()
c     if (ga_nodeid() .eq. 0) then
c      write(LuOut,*) "oneint cpu time", total_cpu, "total time", 
c    1 total_wall
c     end if
C-----------------------------------------------------------------------
C     SINCE SH IS A SUB-BLOCK OF BB MOVE SH INTO F (BB WAS USED IN
C     TRNMAT - TT)
C-----------------------------------------------------------------------
      DO  I=0,NRR-1
c      write(6,*)"SH" ,Dbl_MB(iSH+I)
       Dbl_MB(iFT+I) = Dbl_MB(iSH+I)
c      Dbl_MB(iSYO+I)= Dbl_MB(iSH+I)
      END DO
C--------------------------------------------------------------------------------------
C DIAGONALIZATION OF THE OVERLAP MATRIX
C---------------------------------------------------------------------------------------
c     IF (nbf.LE.2) THEN
c      call symorth(nbf,Dbl_MB(iBB),Dbl_MB(iSYO),Dbl_MB(iQQ),
c    1 INT_MB(iNIN))
c        II = 0
c        DO  I = 1, nbf
c           II = II + I
c           Dbl_MB(ieig+I-1) = Dbl_MB(iFT+II-1)
c        END DO
c       GO TO 15
c      end if

c      write(6,*)"before diag_para"
       call diag_para(Dbl_MB(iFT),NRR,nbf,Dbl_MB(iBB),Dbl_MB(iEIG))
   15    CONTINUE
C------------------------------------------------------------------------------------------  
C GET THE FOCK TRANSFORMATION MATRIX IN TRNMAT
C-------------------------------------------------------------------------------------------
c      write(6,*)"before trnmat"
       call TRNMAT(Dbl_MB(iFT),Dbl_MB(ieig),nbf,NRR,Dbl_MB(iBB))

c     DO  I=0,nrr-1
c      write(6,*)"FT"!,Dbl_MB(iFT+I) 
c     END DO

c      CALL XTHX(INT_MB(iNIN),Dbl_MB(iFT),nbf,Dbl_MB(iFT),nbf,0)
c      do i = 0, nrr-1
c       write(6,*)"SH",Dbl_MB(iFT+i)
c       end do
       if (ga_nodeid().eq.0) WRITE (6,9050)
 9050 FORMAT (/,10X,'FOCK TRANSFORMATION MATRIX FORMED (S**(-1/2))')
c      write(6,*)"OVERLAP"
C     HERE FT  SAVE (OVERLAP)**(-1/2)

c     DO I = 0, NRR-1
c     write(6,*)"FA after tranmat",Dbl_MB(iFA+I)
c     Dbl_MB(iFT+I) = Dbl_MB(iF+I)
c     END DO
C-------------------------------------------------------------------------------------------
C      REMOVE THE ONE CENTER KINETIC ENERGY FOR THIS IS INCLUDED
C      THROUGH PARAMETERS FOR THE CORE INTEGRAL U
C--------------------------------------------------------------------------------------------
c      write(6,*)"before onekie"
       call ONEKIE(ISW1,nbf,Dbl_MB(iFA),Int_MB(inu))
c--------------------------------------------------------------------------------------------
C     THIS PORTION OF THE PROGRAM SETS UP THE HAMILTONIAN, H(IJ),
C     AND INITIALIZES THE ONE-ELECTRON MATRIX FOR THE SCF:
C
C     PARAM - INSERT THE DIAGONAL PARAMETERS, EIS, EIP, EID, AND
C     THE OFF-DIAGONAL PARAMETERS, BETA.
C     EXHUCK- FILL IN THE DIAGONAL ELEMENTS OF THE HAMILTONIAN
C     AND SCALE THE OFF-DIAGONAL ELEMENTS.
C     CORK  - ? NOT IMPLEMENTED AT THIS TIME.
C
C     FETCH - INSERT SLATER-CONDON FACTORS AND S, P, D, F R-TYPE
C     INTEGRALS?
C     FAHUCK- INITIALIZE FA TO THE EXTENDED HUCKEL HAMILTONIAN, H.
C
C     IMPROV- IMPROVE THE INITIAL GUESS FOR THE SCF.
C
C     HCORE - INCLUDE EFFECTIVE CORE POTENTIAL W(B)
C
C-----------------------------------------------------------------------
c      write(6,*)"before exhuck"
       call exhuck(Dbl_MB(iH),Dbl_MB(iF),EIS,EIP,EID,EIF,Dbl_MB(irad),
     1 nat,Dbl_MB(iwcor),NRR,INT_MB(ilqm),nbf,INT_MB(iNIN),
     2 INT_MB(iban),INT_MB(inw),INT_MB(inu),Dbl_MB(ibeta),Dbl_MB(ibetb),
     3 ISW1,INT_MB(ipqn),Dbl_MB(ialpha),Dbl_MB(ibka),Dbl_MB(ibkb)) 
c
c      write(6,*)"before fetch"
       call
     1 fetch(nbf,INT_MB(inw),INT_MB(inu),INT_MB(ipqn),INT_MB(iz),
     2 Dbl_MB(ialpha),int_mb(intype),Dbl_MB(iSCG1SP),Dbl_MB(iSCF2PP),
     3 Dbl_MB(iSCG2SD),Dbl_MB(iSCG1PD),Dbl_MB(iSCF2PD),Dbl_MB(iSCG3PD),
     4 Dbl_MB(iSCF2DD),Dbl_MB(iSCF4DD),Dbl_MB(iSCG3SF),Dbl_MB(iSCG2PF),
     5 Dbl_MB(iSCG4PF),Dbl_MB(iSCG1DF),Dbl_MB(iSCG3DF),Dbl_MB(iSCG5DF),
     6 Dbl_MB(iSCF2PF),Dbl_MB(iSCF2DF),Dbl_MB(iSCF4DF),Dbl_MB(iSCF2FF),
     7 Dbl_MB(iSCF4FF),Dbl_MB(iSCF6FF),Dbl_MB(iR1SPPD),Dbl_MB(iR2SDDD),
     8 Dbl_MB(iR2SDPP),Dbl_MB(iR1PDDF),Dbl_MB(iR1SPDF),Dbl_MB(iR2PFDD),
     9 Dbl_MB(iR2PFFF),Dbl_MB(iR2PPPF),Dbl_MB(iR2SDFF),Dbl_MB(iR2SDPF),
     A Dbl_MB(iR3PDDF),Dbl_MB(iR3SFDF),Dbl_MB(iR3SFPD),Dbl_MB(iR4PFDD),
     A Dbl_MB(iR4PFFF)) 
C
      call fahuck(nbf,NRR,INT_MB(inw),INT_MB(inu),Dbl_MB(iFA),
     1 Dbl_MB(iH),Dbl_MB(iF),ISW1,UCORE)    
      call improv(nbf,Dbl_MB(inu),DBL_MB(ival),Dbl_MB(iH),ISW1) 
c--------------------------------------------------------------------------------------------
c         FORM THE ONE-ELECTRON MATRIX
c--------------------------------------------------------------------------------------------
c      write(6,*)"uform"
c      do i=0, nrr
c      write(6,*)"F before scfit", Dbl_MB(iF+i)
c      write(6,*)"FA before uform", Dbl_MB(iFA+i)
c      write(6,*)"H before scfit", Dbl_MB(iH+i)
c      write(6,*)"SH before scfit", Dbl_MB(iSH+i)
c     end do
      IF(ISW1.LT.50) THEN
      call uform(nat,nbf,INT_MB(inw),Dbl_MB(iF),Dbl_MB(iFA),
     1 Dbl_MB(iwcor),Dbl_MB(igama),INT_MB(iban),INT_MB(ilqm),
     2 INT_MB(inu),Dbl_MB(irad),ISW1,Int_MB(ing),Dbl_MB(ialpha),COFSQ2,
     3 INT_MB(inin),Int_MB(intype),Dbl_MB(ival),Dbl_MB(ival),ISW2,
     3 Int_MB(ipqn),Dbl_MB(iSCG1SP),Dbl_MB(iSCF2PP),
     3 Dbl_MB(iSCG2SD),Dbl_MB(iSCG1PD),Dbl_MB(iSCF2PD),Dbl_MB(iSCG3PD),
     4 Dbl_MB(iSCF2DD),Dbl_MB(iSCF4DD),Dbl_MB(iSCG3SF),Dbl_MB(iSCG2PF),
     5 Dbl_MB(iSCG4PF),Dbl_MB(iSCG1DF),Dbl_MB(iSCG3DF),Dbl_MB(iSCG5DF),
     6 Dbl_MB(iSCF2PF),Dbl_MB(iSCF2DF),Dbl_MB(iSCF4DF),Dbl_MB(iSCF2FF),
     7 Dbl_MB(iSCF4FF),Dbl_MB(iSCF6FF),Dbl_MB(iR1SPPD),Dbl_MB(iR2SDDD),
     8 Dbl_MB(iR2SDPP),Dbl_MB(iR1PDDF),Dbl_MB(iR1SPDF),Dbl_MB(iR2PFDD),
     9 Dbl_MB(iR2PFFF),Dbl_MB(iR2PPPF),Dbl_MB(iR2SDFF),Dbl_MB(iR2SDPF),
     A Dbl_MB(iR3PDDF),Dbl_MB(iR3SFDF),Dbl_MB(iR3SFPD),Dbl_MB(iR4PFDD),
     A Dbl_MB(iR4PFFF),Dbl_MB(isumelel),Dbl_MB(iETS),Dbl_MB(iETP),
     A Dbl_MB(iETD),Dbl_MB(iETF),inttyp)
        ELSEIF(ISW1.LT.54) THEN
         CALL UFORMA(nbf,nat,UCORE,DBL_MB(ieul),Dbl_MB(iFA),Dbl_MB(iRAD)
     1   ,INT_MB(ilqm),INT_MB(iNIN),INT_MB(iNU),INT_MB(iNW),INT_MB(iz) 
     2   ,Dbl_MB(iPO),Dbl_MB(iDAD),Dbl_MB(ival))
        ENDIF
c      ENDIF       
c     DO I = 0, nrr
c      write(6,*)"FA",Dbl_MB(iFA+I)
c     END DO
      LL = 0
c----------------------------------------------------------------------------------------------
c  FORM TWO-ELECTRON INTEGRALS FROM COULOUMB MATRIX, GAMMA FOR ALL ISW1
c----------------------------------------------------------------------------------------------
      IF (ISW1.GT.50.AND.ISW1.LT.54) THEN
       CALL GFORMAM(nbf,nat,Dbl_MB(igama),INT_MB(inu),
     1  INT_MB(iban),LL,INT_MB(iLABI),INT_MB(iLABJ),INT_MB(iLABK),
     2  INT_MB(iLABL),Dbl_MB(iQQ),GSS,GPP,GSP,GP2,HSP,HPP,Dbl_MB(iPO)
     3 ,Dbl_MB(iDAD),INT_MB(infirst),Dbl_MB(i_xyz))
       ELSE
        call gform(nbf,Int_MB(ing),Dbl_MB(igama),Dbl_MB(iFA),
     1 INT_MB(inu),INT_MB(inw),INT_MB(iban),INT_MB(intype),ISW1,
     2 Dbl_MB(iSCG1SP),Dbl_MB(iSCF2PP),
     3 Dbl_MB(iSCG2SD),Dbl_MB(iSCG1PD),Dbl_MB(iSCF2PD),Dbl_MB(iSCG3PD),
     4 Dbl_MB(iSCF2DD),Dbl_MB(iSCF4DD),Dbl_MB(iSCG3SF),Dbl_MB(iSCG2PF),
     5 Dbl_MB(iSCG4PF),Dbl_MB(iSCG1DF),Dbl_MB(iSCG3DF),Dbl_MB(iSCG5DF),
     6 Dbl_MB(iSCF2PF),Dbl_MB(iSCF2DF),Dbl_MB(iSCF4DF),Dbl_MB(iSCF2FF),
     7 Dbl_MB(iSCF4FF),Dbl_MB(iSCF6FF),Dbl_MB(iR1SPPD),Dbl_MB(iR2SDDD),
     8 Dbl_MB(iR2SDPP),Dbl_MB(iR1PDDF),Dbl_MB(iR1SPDF),Dbl_MB(iR2PFDD),
     9 Dbl_MB(iR2PFFF),Dbl_MB(iR2PPPF),Dbl_MB(iR2SDFF),Dbl_MB(iR2SDPF),
     A Dbl_MB(iR3PDDF),Dbl_MB(iR3SFDF),Dbl_MB(iR3SFPD),Dbl_MB(iR4PFDD),
     A Dbl_MB(iR4PFFF),INT_MB(inin),LL)
       ENDIF

       if (.not.rtdb_put(rtdb,'semiemp:numint',mt_int,1,LL))
     1  call errquit('semiemp: failed to write numint',0, RTDB_ERR)

c      call semiemp_integral_save (LL,INT_MB(iLABI),INT_MB(iLABJ),
c    A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ),nbf) 

c      call semiemp_integral_load(LL,INT_MB(iLABI),INT_MB(iLABJ),
c    A   INT_MB(iLABK),INT_MB(iLABL),Dbl_MB(iQQ))

c----------------------------------------------------------------------------------------------
c     ITS TIME TO DO SELF-CONSISTENT FIELD ITERATIONS
C     CALL SCFIT TO DO THE SELF-CONSISTENT FIELD ITERATIONS.
C---------------------------------------------------------------------------------------------
c      write(6,*)"here i am"
c      do i=0, nrr-1
c      write(6,*)"F before scfit", Dbl_MB(iF+i)
c      write(6,*)"FA before scfit", Dbl_MB(iFA+i)
c      write(6,*)"H before scfit", Dbl_MB(iH+i)
c      write(6,*)"SH before scfit", Dbl_MB(iSH+i)
c      end do
c      call mat_sym(Dbl_MB(iFT),nbf)
c      call mat_sym(Dbl_MB(iF),nbf)
c      call mat_sym(Dbl_MB(iFA),nbf)
c      call mat_sym(Dbl_MB(iH),nbf)
c      do i=1, nbf
c      write(6,*)"ASO,ATOM",(byte_mb(iaso+(i-1)*8+j),j=0,7),
c    1 (byte_mb(iatom+(i-1)*8+j),j=0,7)
c      end do
c      write(6,*)"before scfit"

       call scfit(rtdb,Dbl_MB(iH),Dbl_MB(iF),Dbl_MB(iFA),Dbl_MB(iSH),
     1 Dbl_MB(iFT),Dbl_MB(irad),Dbl_MB(igama),ISW1,nbf,nat,ISWANN,scftol
     2 ,itmax,MIM,ISW3,Dbl_MB(ieig),FOP,NEL,MULT,NRR,Dbl_MB(iBB),NDN
     3 ,INT_MB(intype),INT_MB(inu),INT_MB(ing),JSLIT,Dbl_MB(isumelel),
     4  REP,Dbl_MB(iFZ),INT_MB(iNIN),E1,Dbl_MB(id12),
     6 Dbl_MB(iAA),Dbl_MB(iRA),Dbl_MB(iRB),NCORE,Dbl_MB(irho),LL)
c       write(6,*)"energy",E1
c       E1 = E1/100
      total_wall = total_wall + util_wallsec()
      total_cpu  = total_cpu  + util_cpusec()

       if (.not.rtdb_put(rtdb,'semiemp:energy',mt_dbl,1,E1))
     1  call errquit('semiemp: failed to write energy',0, RTDB_ERR)
       if (.not.rtdb_put(rtdb,'semiemp:nel',mt_int,1,nel))
     1  call errquit('semiemp: failed to write nel',0, RTDB_ERR)

       CALL SLITE(0,JSLIT)
       CALL PURVIS(Dbl_MB(iF),nbf,Dbl_MB(irad),NAT,INT_MB(iz),
     1 INT_MB(intype),Dbl_MB(i_xyz),NCORE,ISW1)
       
       CALL BONDOR(Dbl_MB(irho),NRR,INT_MB(iNIN),INT_MB(inu),ISW1,NAT,
     1 nbf,ISWANN)

       IORB = 1
       CALL DIPOLE(rtdb,Dbl_MB(iRAD),INT_MB(iz),NAT,Dbl_MB(i_xyz),
     1 INT_MB(iNU),INT_MB(iNW),INT_MB(ipqn),INT_MB(iNIN),Dbl_MB(iALPHA),
     2 Dbl_MB(iH),Dbl_MB(iBB),nbf,nrr,ISW1,Dbl_MB(ival),Dbl_MB(id12),
     3 Dbl_MB(irho),NEL,E,Dbl_MB(idipx),Dbl_MB(idipy),Dbl_MB(idipz),1
     4 ,PME,IHOMO,ILUMO,IORB,MIM(1))

       if (.not.rtdb_put(rtdb,'semiemp:dipole',mt_dbl,3,E))
     1  call errquit('semiemp: failed to write dipole',0, RTDB_ERR)
       if (.not.rtdb_put(rtdb,'semiemp:eig',mt_dbl,nbf,Dbl_MB(ieig)))
     1  call errquit('semiemp: failed to write eig',0, RTDB_ERR)


       if (ISWANN.eq.0) then
        if (.not.rtdb_put(rtdb,'semiemp:eig',mt_dbl,nbf,Dbl_MB(ieig)))
     1  call errquit('semiemp: failed to write eig',0, RTDB_ERR)
       
       CALL semiemp_mo_print(rtdb,int_mb(iatnum),nbf)
      end if 
c      do i = 0 , nrr-1
c       write(6,*)"FZ",Dbl_MB(idipz+i)
c      end do
        if (.not.rtdb_get(rtdb,'semiemp:docis',mt_log,1,docis))
     1  call errquit('semiemp: failed to read docis',0, RTDB_ERR)
       
       if (docis) then
        Write(6,*)"CIS CALCULATION"
        call semiemp_cis(rtdb,nbf,Dbl_MB(ieig))       
        end if
       if (ISWANN.eq.1) 
     1 CALL semiemp_mo_uhf_print(rtdb,int_mb(iatnum),nbf)
C-----------------------------------------------------------------------------------
C THIS PART PERFORMS REAL TIME CALCULATION
C----------------------------------------------------------------------------------
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '
      if (inp_compare(.false., 'rt_semiemp', operation)) then  !CSG
         call rtsemiemp(rtdb,dbl_mb(iSH),nbf,ISWANN,Dbl_MB(irho),
     1 Dbl_MB(isumelel),LL,Dbl_MB(iFA),INT_MB(iNIN),
     3 INT_MB(ing),Dbl_MB(igama),Dbl_MB(idipx),Dbl_MB(idipy),
     4 Dbl_MB(idipz),Dbl_MB(iRAD),INT_MB(iz),nat,Dbl_MB(i_xyz),
     5 INT_MB(inu),INT_MB(inw),INT_MB(ipqn),Dbl_MB(iALPHA),
     2 Dbl_MB(iH),Dbl_MB(iBB),NEL,ISW1,Dbl_MB(ival),Dbl_MB(id12),
     3 Dbl_MB(iFT),Dbl_MB(iRA),Dbl_MB(iRB))                  
       end if
C-------------------------------------------------------------------------------------------
C THIS PART CALCULATES ANALYTICAL GRADIENT FOR THE SEMIEMPIRICAL METHOD
C AND STORE IT INTO RTDB
C-------------------------------------------------------------------------------------------
      if (.not. rtdb_get(rtdb, 'task:numerical', mt_log, 1,
     $     numerical)) numerical=.false.

      if (numerical) GO TO 20
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '
      if  ((inp_compare(.false.,'gradient',operation)).or.(inp_compare
     $ (.false.,'optimize',operation)).or.(inp_compare(.false.,'qmd',
     $ operation)))  then  
         call grad_semiemp(nat,nbf,Dbl_MB(ialpha),Int_Mb(ipqn),
     1 Int_MB(inw),int_mb(iban),Int_MB(inu),Dbl_MB(irad),Dbl_MB(ieul)
     2 ,INT_MB(iNIN),Dbl_MB(irho),Dbl_MB(iH),Dbl_MB(ibeta),Dbl_MB(ibetb)
     3 ,ISW1,Dbl_MB(ibka),Dbl_MB(ibkb),NRR,INT_MB(ing),Dbl_MB(iRA),
     4 Dbl_MB(iRB),Dbl_MB(ival),INT_MB(intype),ISWANN,Dbl_MB(i_xyz),
     A Dbl_MB(iETS),Dbl_MB(iETP),Dbl_MB(iETD),Dbl_MB(iETF),Dbl_MB(iFZ),
     A Dbl_MB(igrad),inttyp,NDN)
 
c      call grad_force(rtdb,nbf,geom)
c      CALL sym_grad_symmetrize(geom,dbl_mb(iFZ))
c      do i = 0 , nrr-1
c       write(6,*)"FZ",Dbl_MB(idipz+i)
c      end do
       if (.not.rtdb_put(rtdb,'semiemp:gradient',mt_dbl,3*nat,
     1 Dbl_MB(igrad)))
     1  call errquit('semiemp: failed to write gradient',0, RTDB_ERR)
       end if
  20   CONTINUE
C---------------------------------------------------------------------------------------------
       if (ga_nodeid().eq.0)
     & write(6,*) "back to semiemp" 
       if (.not. geom_destroy(geom)) call errquit
     $     ('geom_print_rtdb_ecce: geom destroy failed',0, GEOM_ERR)
c---------------------------------------------------------------------------------------------  
c     Clean memory
         IF (ISW1.LT.54.and.ISW1.GT.50) THEN
         if (.not.ma_pop_stack(lAALP))
     &   call errquit('semiemp_energy: cannot pop stack AALP',3, MA_ERR)
         if (.not.ma_pop_stack(lDAD))
     &   call errquit('semiemp_energy: cannot pop stack DAD',3, MA_ERR)
         if (.not.ma_pop_stack(lQQ))
     & call errquit('semiemp_energy: cannot pop stack QQ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABL))
     & call errquit('semiemp_energy: cannot pop stack LABL',3, MA_ERR)
         if (.not.ma_pop_stack(lLABK))
     & call errquit('semiemp_energy: cannot pop stack LABK',3, MA_ERR)
         if (.not.ma_pop_stack(lLABJ))
     & call errquit('semiemp_energy: cannot pop stack LABJ',3, MA_ERR)
         if (.not.ma_pop_stack(lLABI))
     & call errquit('semiemp_energy: cannot pop stack LABI',3, MA_ERR)
         END IF
      if (.not.ma_pop_stack(latnum))
     &   call errquit('semiemp_energy: cannot pop stack ATOM_NUM',3, 
     & MA_ERR)
      if (.not.ma_pop_stack(latom))
     &   call errquit('semiemp_energy: cannot pop stack ATOM',3, MA_ERR)
      if (.not.ma_pop_stack(laso))
     &   call errquit('semiemp_energy: cannot pop stack SYO',3, MA_ERR)
      if (.not.ma_pop_stack(lPO))
     &   call errquit('semiemp_energy: cannot pop stack PO',3, MA_ERR)
      if (.not.ma_pop_stack(lnfirst))
     & call errquit('semiemp_energy: cannot pop stack NFIRST',3, MA_ERR)
      if (.not.ma_pop_stack(ldipz))
     &   call errquit('semiemp_energy: cannot pop stack FZ',3, MA_ERR)
      if (.not.ma_pop_stack(ldipy))
     &   call errquit('semiemp_energy: cannot pop stack FY',3, MA_ERR)
      if (.not.ma_pop_stack(ldipx))
     &   call errquit('semiemp_energy: cannot pop stack FX',3, MA_ERR)
         if (.not.ma_pop_stack(lrho))
     & call errquit('semiemp_energy: cannot pop stack Rhob',3, MA_ERR)
         if (.not.ma_pop_stack(lgrad))
     & call errquit('semiemp_energy: cannot pop stack GRAD',3, MA_ERR)
         if (.not.ma_pop_stack(lETF))
     & call errquit('semiemp_energy: cannot pop stack ETF',3, MA_ERR)
         if (.not.ma_pop_stack(lETD))
     & call errquit('semiemp_energy: cannot pop stack ETD',3, MA_ERR)
         if (.not.ma_pop_stack(lETP))
     & call errquit('semiemp_energy: cannot pop stack ETP',3, MA_ERR)
         if (.not.ma_pop_stack(lETS))
     & call errquit('semiemp_energy: cannot pop stack ETS',3, MA_ERR)
         if (.not.ma_pop_stack(lRB))
     & call errquit('semiemp_energy: cannot pop stack RB',3, MA_ERR)
         if (.not.ma_pop_stack(lRA))
     & call errquit('semiemp_energy: cannot pop stack RA',3, MA_ERR)
         if (.not.ma_pop_stack(lbkb))
     & call errquit('semiemp_energy: cannot pop stack BKB',3, MA_ERR)
         if (.not.ma_pop_stack(lbka))
     & call errquit('semiemp_energy: cannot pop stack BKA',3, MA_ERR)
         if (.not.ma_pop_stack(lAA))
     & call errquit('semiemp_energy: cannot pop stack AA',3, MA_ERR)
         if (.not.ma_pop_stack(ld12))
     & call errquit('semiemp_energy: cannot pop stack D12',3, MA_ERR)
         if (.not.ma_pop_stack(lfdmat2))
     & call errquit('semiemp_energy: cannot pop stack E',3, MA_ERR)
         if (.not.ma_pop_stack(lfdmat))
     & call errquit('semiemp_energy: cannot pop stack E',3, MA_ERR)
         if (.not.ma_pop_stack(lddmat))
     & call errquit('semiemp_energy: cannot pop stack D',3, MA_ERR)
         if (.not.ma_pop_stack(lsumelel))
     & call errquit('semiemp_energy: cannot pop stack SUMELE',3, MA_ERR)
         if (.not.ma_pop_stack(lFZ))
     & call errquit('semiemp_energy: cannot pop stack FZ',3, MA_ERR)
         if (.not.ma_pop_stack(lBB))
     & call errquit('semiemp_energy: cannot pop stack BB',3, MA_ERR)
         if (.not.ma_pop_stack(lF))
     & call errquit('semiemp_energy: cannot pop stack F',3, MA_ERR)
         if (.not.ma_pop_stack(lr4pfff))
     & call errquit('semiemp_energy: cannot pop stack R4PFFF',3, MA_ERR)
         if (.not.ma_pop_stack(lr4pfdd))
     & call errquit('semiemp_energy: cannot pop stack R4PFDD',3, MA_ERR)
         if (.not.ma_pop_stack(lr3sfpd))
     & call errquit('semiemp_energy: cannot pop stack R3SFPD',3, MA_ERR)
         if (.not.ma_pop_stack(lr3sfdf))
     & call errquit('semiemp_energy: cannot pop stack R3SFDF',3, MA_ERR)
         if (.not.ma_pop_stack(lr3pddf))
     & call errquit('semiemp_energy: cannot pop stack R3PDDF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2sdpf))
     & call errquit('semiemp_energy: cannot pop stack R2SDPF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2sdff))
     & call errquit('semiemp_energy: cannot pop stack R2SDFF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2pppf))
     & call errquit('semiemp_energy: cannot pop stack R2PPPF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2pfff))
     & call errquit('semiemp_energy: cannot pop stack R2PFFF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2pfdd))
     & call errquit('semiemp_energy: cannot pop stack R2PFDD',3, MA_ERR)
         if (.not.ma_pop_stack(lr1spdf))
     & call errquit('semiemp_energy: cannot pop stack R1SPDF',3, MA_ERR)
         if (.not.ma_pop_stack(lr1pddf))
     & call errquit('semiemp_energy: cannot pop stack R1PDDF',3, MA_ERR)
         if (.not.ma_pop_stack(lr2sdpp))
     & call errquit('semiemp_energy: cannot pop stack R2SDPP',3, MA_ERR)
          if (.not.ma_pop_stack(lr2sddd))
     & call errquit('semiemp_energy: cannot pop stack R2SDDD',3, MA_ERR)
          if (.not.ma_pop_stack(lr1sppd))
     & call errquit('semiemp_energy: cannot pop stack R1SPPD',3, MA_ERR)
          if (.not.ma_pop_stack(lscf6ff))
     & call errquit('semiemp_energy: cannot pop stack SCF6FF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf4ff))
     & call errquit('semiemp_energy: cannot pop stack SCF4FF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2ff))
     & call errquit('semiemp_energy: cannot pop stack SCF2FF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf4df))
     & call errquit('semiemp_energy: cannot pop stack SCF4DF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2df))
     & call errquit('semiemp_energy: cannot pop stack SCF2DF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2pf))
     & call errquit('semiemp_energy: cannot pop stack SCF2PF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg5df))
     & call errquit('semiemp_energy: cannot pop stack SCG5DF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg3df))
     & call errquit('semiemp_energy: cannot pop stack SCG3DF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg1df))
     & call errquit('semiemp_energy: cannot pop stack SCG1DF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg4pf))
     & call errquit('semiemp_energy: cannot pop stack SCG4PF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg2pf))
     & call errquit('semiemp_energy: cannot pop stack SCG2PF',3, MA_ERR)
          if (.not.ma_pop_stack(lscg3sf))
     & call errquit('semiemp_energy: cannot pop stack SCG3SF',3, MA_ERR)
          if (.not.ma_pop_stack(lscf4dd))
     & call errquit('semiemp_energy: cannot pop stack SCF4DD',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2dd))
     & call errquit('semiemp_energy: cannot pop stack SCF2DD',3, MA_ERR)
          if (.not.ma_pop_stack(lscg3pd))
     & call errquit('semiemp_energy: cannot pop stack SCG3PD',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2pd))
     & call errquit('semiemp_energy: cannot pop stack SCF2PD',3, MA_ERR)
          if (.not.ma_pop_stack(lscg1pd))
     & call errquit('semiemp_energy: cannot pop stack SCG1PD',3, MA_ERR)
          if (.not.ma_pop_stack(lscg2sd))
     & call errquit('semiemp_energy: cannot pop stack SCG2SD',3, MA_ERR)
          if (.not.ma_pop_stack(lscf2pp))
     & call errquit('semiemp_energy: cannot pop stack SCF2PP',3, MA_ERR)
          if (.not.ma_pop_stack(lscg1sp))
     & call errquit('semiemp_energy: cannot pop stack SCG1SP',3, MA_ERR)
          if (.not.ma_pop_stack(lnin))
     &   call errquit('semiemp_energy: cannot pop stack lNIN',2, MA_ERR)
          if (.not.ma_pop_stack(lwcor))
     &  call errquit('semiemp_energy: cannot pop stack wcor',3, MA_ERR)
           if (.not.ma_pop_stack(leig))
     &  call errquit('semiemp_energy: cannot pop stack EIG',3, MA_ERR)
           if (.not.ma_pop_stack(lFT))
     &  call errquit('semiemp_energy: cannot pop stack FT',3, MA_ERR)
           if (.not.ma_pop_stack(lFA))
     &  call errquit('semiemp_energy: cannot pop stack FA',3, MA_ERR)
           if (.not.ma_pop_stack(lSH))
     &  call errquit('semiemp_energy: cannot pop stack SH',3, MA_ERR)
           if (.not.ma_pop_stack(lH))
     &  call errquit('semiemp_energy: cannot pop stack H',3, MA_ERR)  
           if (.not.ma_pop_stack(lbetb))
     &  call errquit('semiemp_energy: cannot pop stack betb',3, MA_ERR)
           if (.not.ma_pop_stack(lbeta))
     &  call errquit('semiemp_energy: cannot pop stack beta',3, MA_ERR)
            if (.not.ma_pop_stack(lalpha))
     &  call errquit('semiemp_energy: cannot pop stack alpha',3, MA_ERR)
            if (.not.ma_pop_stack(lgama))
     &  call errquit('semiemp_energy: cannot pop stack gama',3, MA_ERR)
            if (.not.ma_pop_stack(lng))
     &  call errquit('semiemp_energy: cannot pop stack NG',3, MA_ERR)
            if (.not.ma_pop_stack(lban))
     &  call errquit('semiemp_energy: cannot pop stack BAN',3, MA_ERR)
            if (.not.ma_pop_stack(lpqn))
     &  call errquit('semiemp_energy: cannot pop stack nu',3, MA_ERR)
           if (.not.ma_pop_stack(lnu))
     &  call errquit('semiemp_energy: cannot pop stack nu',3, MA_ERR)
           if (.not.ma_pop_stack(lnw))
     &  call errquit('semiemp_energy: cannot pop stack nw',3, MA_ERR)
           if (.not.ma_pop_stack(leul))
     &  call errquit('semiemp_energy: cannot pop stack eul',3, MA_ERR)
       if (.not.ma_pop_stack(llqm))
     &   call errquit('semiemp_energy:cannot pop stack lqm',3, MA_ERR)       
         if (.not.ma_pop_stack(ltags))
     &   call errquit('semiemp_energy: cannot pop stack tag',3, MA_ERR)         
          if (.not.ma_pop_stack(lDK))
     &  call errquit('semiemp_energy: cannot pop stack DK',3, MA_ERR)
          if (.not.ma_pop_stack(lrad))
     &  call errquit('semiemp_energy: cannot pop stack rad',3, MA_ERR)
           if (.not.ma_pop_stack(lnp))
     &  call errquit('semiemp_energy: cannot pop stack NP',3, MA_ERR)
           if (.not.ma_pop_stack(lnb))
     &  call errquit('semiemp_energy: cannot pop stack NB',3, MA_ERR)
          if (.not.ma_pop_stack(lcor))
     &   call errquit('semiemp_energy: cannot pop stack lcor',3, MA_ERR)
         if (.not.ma_pop_stack(lntype))
     &   call errquit('semiemp_energy: cannot pop stack ntyp',3, MA_ERR)
        if (.not.ma_pop_stack(lval))
     &   call errquit('semiemp_energy: cannot pop stack lval',3, MA_ERR)
      if (.not.ma_pop_stack(lcharge))
     &   call errquit('semiemp_energy: cannot pop stack lchr',3, MA_ERR)
      if (.not.ma_pop_stack(lz))
     &   call errquit('semiemp_energy: cannot pop stack lz',2, MA_ERR)
      if (.not.ma_pop_stack(l_xyz))
     &   call errquit('semiemp_energy: cannot pop stack lxyz',3, MA_ERR)

       semiemp_energy = .true.
c
       return
       end




