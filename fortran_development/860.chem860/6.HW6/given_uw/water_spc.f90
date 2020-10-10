  program water_spc
  use water_spc_mod
  implicit none
!  --------------------------------------------------------------------------
  integer            :: ncycle         !number of steps used used in Monte Carlo 
  integer            :: i,j,icount     !counters
  integer            :: istep          !counter of successful moves in mcmove
  integer            :: kstep,jstep    !used to calc fraction of successful steps in the volume change
  real(kind=8)       :: avg_density    !average density
  real(kind=8)       :: avg_energy     !average energy
  real(kind=8)       :: dens_cgs       !density in cgs units
  real(kind=8)       :: rnd            !a random number
  real(kind=8)       :: vpot           !the total potential energy
  real(kind=8)       :: kappa_check    !used to check if kappa is set ok
! ----------------------------------------------------------------
  open(16,file='water_spc.dat')
! ----------------------------------------------------------------
! read in the input parameters
  call input(ncycle)
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed()
! ----------------------------------------------------------------
! set up the initial positions 
  call init 
! ----------------------------------------------------------------
  first = 1          !must be recalled if the volume changes
  call ewald_setup   !ifirst=1 on the first call to setup
! ----------------------------------------------------------------
  call self          !calc intramolecular contrib. to coulomb energy
! ----------------------------------------------------------------
! check that the value of kappa in ewald sum is OK.
  kappa_check= 2.d0*erfc(0.5d0*kappa*box)/box
  if(kappa_check>1.d-5)then
   write(6,'(a,d12.6)')'2.d0*erfc(0.5d0*kappa*box)/box =',kappa_check
   write(6,*)' type cr if are you sure this value is OK?'
   read(5,*) 
  endif
! ----------------------------------------------------------------
! calculate the initial contributions to the pot energy
  call kwald      !long range ewald -  pot_lr
  call rwald      !short range ewald - pot_sr
  call potential  !6 12 potential    - pot_lj
  write(ntape,*)' the following energies are in kcal/mol'
  vpot = pot_lr+pot_sr+pot_lj-vself  !total pot
  write(ntape,*)' the initial lr potential is ',pot_lr/e_conv
  write(ntape,*)' the initial sr potential is ',pot_sr/e_conv
  write(ntape,*)' the initial lj potential is ',pot_lj/e_conv
  write(ntape,*)' the intramolecular contr is ',vself/e_conv
  write(ntape,*)' the initial tot coulomb  is ',(pot_lr+pot_sr-vself)/e_conv
  write(ntape,*)' the total potential is      ',vpot/e_conv
! ----------------------------------------------------------------
! set initial quantities to zero
  avg_density = 0.d0     
  avg_energy = 0.d0

  icount = 0  !counts total sampled values
  kstep = 0   !counts calls to mcvol
  jstep = 0   !counts successful volume moves
  istep = 0   !counts successful config moves
! ----------------------------------------------------------------
  write(ntape,*)'                  potential contributions (kcal/mol) '
  write(ntape,*)' density(g/cm^3)  pot_lr   pot_sr       pot_lj      vpot'
! ----------------------------------------------------------------
! begin the monte carlo cycle
  do i = 1,ncycle
     call random_number(rnd)                 !pick random number
     rnd = rnd*(nmol+1)
     if(rnd<=nmol+const_vol)then
      call mcmove(vpot,istep)                !attempt a config move
     else
      call mcvol(vpot,jstep)                 !attempt a volume move
      kstep = kstep + 1
     endif
      
     if(mod(i-1,nsamp)==0)then               !calculate running averages
       icount = icount + 1
       call center                           !maintain molecules in box
       dens_cgs = (nmol/volume)/av*(1.d24*cmass)  !conversion to g/cm^3 density
       avg_density = avg_density + dens_cgs 
       avg_energy = avg_energy + vpot
!      --------------------------------------------------------------
       write(ntape,'(6f12.5)')dens_cgs,pot_lr/e_conv,pot_sr/e_conv, &
             pot_lj/e_conv,vpot/e_conv
     endif
  enddo

  call center
  avg_density = avg_density/icount
  avg_energy = avg_energy/e_conv/icount
  write(ntape,*)' istep = ',istep
  write(ntape,*)' ncycle = ',ncycle
  write(ntape,*)' % accepted coord is        = ',100*dfloat(istep)/dfloat(ncycle-kstep)
  write(ntape,*)' % accepted volume moves is = ',100*dfloat(jstep)/dfloat(kstep)
  write(ntape,*)' the final volume is = ',volume
  write(ntape,'(a,f14.8)')' the average density is ',avg_density
  write(ntape,'(a,f14.8)')' the average energy is ',avg_energy
  write(ntape,'(a,f14.8)')' moles/liter = ',1000.d0*avg_density/cmass 

! --------------------------------------------------------------
! save final configuration, one file for xyz output
  open(18,file='old_config.ini')
  open(19,file='config_h2o.xyz')
    write(19,'(i5,a)')npart
    write(19,'(f12.6,a)')
    write(18,'(f13.7,a)')box,'  ! box length'
    write(18,*)x
    do i = 1,npart
     write(19,'(a4,3f10.3)')atom_name(ni(mod(i-1,natom)+1)),x(i,:)
    enddo
  close(18)
  close(19)
! --------------------------------------------------------------
    

  end program water_spc




