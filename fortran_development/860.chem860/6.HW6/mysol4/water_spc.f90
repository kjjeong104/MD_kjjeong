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
  real(kind=8)       :: rnd,tt,told            !a random number
  real(kind=8)       :: vpot           !the total potential energy
  real(kind=8)       :: kappa_check    !used to check if kappa is set ok
! ----------------------------------------------------------------
  open(16,file='spc_water.dat')
! ----------------------------------------------------------------
! read in the input parameters
  call input(ncycle)
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
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
   write(16,'(a,d13.6)')'2.d0*erfc(0.5d0*kappa*box)/box =',kappa_check
   write(16,*)' type cr if are you sure this value is OK?'
   read(5,*) 
  endif
! ----------------------------------------------------------------
! calculate the initial contributions to the pot energy
  call kwald      !long range ewald -  pot_lr
  call rwald      !short range ewald - pot_sr
  call potential  !6 12 potential    - pot_lj
  write(16,*)' the following energies are in kcal/mol'
  vpot = pot_lr+pot_sr+pot_lj-vself+e_lrc  !total pot
  write(16,*)' the initial lr potential is ',pot_lr/e_conv
  write(16,*)' the initial sr potential is ',pot_sr/e_conv
  write(16,*)' the initial lj potential is ',pot_lj/e_conv
  write(16,*)' the intramolecular contr is ',vself/e_conv
  write(16,*)' the initial tot coulomb  is ',(pot_lr+pot_sr-vself)/e_conv
  write(16,*)' the total potential is      ',vpot/e_conv
! ----------------------------------------------------------------
! set initial quantities to zero
  avg_density = 0.d0     
  avg_energy = 0.d0
  call cpu_time(told)
  icount = 0  !counts total sampled values
  kstep = 0   !counts calls to mcvol
  jstep = 0   !counts successful volume moves
  istep = 0   !counts successful config moves
  call gr(0,nmol,density,box,x) !initialize gr calc
! ----------------------------------------------------------------
  write(6,*)'                  potential contributions (kcal/mol) '
  write(6,*)' density(g/cm^3)  pot_lr   pot_sr       pot_lj      vpot    volume     e_lrc    ecut'
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
	do j = 1,nmol
		xoxy(j,1)=x(j*3-2,1)
		xoxy(j,2)=x(j*3-2,2)
		xoxy(j,3)=x(j*3-2,3)
	enddo
	call gr(1,nmol,density,box,xoxy)
       dens_cgs = (nmol/volume)/av*(1.d24*cmass)  !conversion to g/cm^3 density
       avg_density = avg_density + dens_cgs 
       avg_energy = avg_energy + vpot
!      --------------------------------------------------------------
       write(6,'(8f12.5)')dens_cgs,pot_lr/e_conv,pot_sr/e_conv, &
             pot_lj/e_conv,vpot/e_conv,volume,e_lrc/e_conv,ecut/e_conv
     endif
  enddo
  call cpu_time(tt)
  call center
  avg_density = avg_density/icount
  avg_energy = avg_energy/e_conv/icount
  write(16,*)' istep = ',istep
  write(16,*)' ncycle = ',ncycle
  write(16,*)' % accepted coord is        = ',100*dfloat(istep)/dfloat(ncycle-kstep)
  write(16,*)' % accepted volume moves is = ',100*dfloat(jstep)/dfloat(kstep)
  write(16,*)' the final volume is = ',volume
  write(16,'(a,f14.8)')' the average density is ',avg_density
  write(16,'(a,f14.8)')' the average energy is ',avg_energy
  write(16,'(a,f14.8)')' moles/liter = ',1000.d0*avg_density/cmass 
  write(16,*)'cpu time=',tt-told
! --------------------------------------------------------------
! save final configuration, one file for povray output
  call gr(2,nmol,density,box,xoxy)
  open(18,file='old_config.ini')
  open(19,file='config_pov.dat')
    write(19,'(i5,a)')nmol,','
    write(19,'(f12.6,a)')box,','
    write(18,'(f12.6,a)')box,'  ! box length'
    write(18,*)x
    write(19,99)((x(i,j),j=1,3),i=1,npart)
 99 format(3(f10.3,','))
  close(18)
  close(19)
! --------------------------------------------------------------
    

  end program water_spc



