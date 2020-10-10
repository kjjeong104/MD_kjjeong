  program lj_met_chain
  use lj_met_chain_mod
  implicit none
!  --------------------------------------------------------------------------
  integer                    :: i,j                      !counter variables
  integer                    :: irun,nequil,iflag    !number of steps used used in Monte Carlo 
  integer                    :: istep=0                  !counter of successful moves
  integer                    :: jstep=0,kstep=0          ! used to calc fraction of successful steps in the volume change
  real(kind=8)               :: dens_avg,vavg,pavg,PLRCavg            !average calculated density and potential,pressure
  real(kind=8)               :: rnd                      !random number
  real(kind=8)                                  :: vpot,PLJ, pres1,PLRC
  real(kind=8)                                  :: dens_std,pstd !standard deviation of the blocks
  real(kind=8), allocatable, dimension(:)       :: dens,pres_block,PLRC_block     !average density of a block 
! ----------------------------------------------------------------
  open(16,file='lj_met_chain.dat')
! read in the input parameters
  call input(irun,nequil)
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------
! set up the initial positions
  call init !set up an initial condition
  call potential   !calculate the potentials for the LJ potential
  vpot = pot_lj + pot_str + e_lrc
  write(6,*)' The initial potential (vpot,pot_lj,pot_str)is',vpot,pot_lj,pot_str
  write(*,*)' Start running equilibration.  Vpot and density will be plotted.'
  pause
! ----------------------------------------------------------------
        write(6,*)'vpot    density    pres1    box'
  do i = 1,nequil
     call random_number(rnd)        !pick angle to be moved
     rnd = rnd*(npart+1)
     if(rnd<=npart+const_vol)then
      call mcmove(vpot,istep)      !attempts to displace a particle
     else
      call mcvol(vpot,jstep,kstep) !attempts a volume move
     endif
     if(mod(i,nsamp)==1)then
       call center
       call potential    !calculate the potential energy
	vpot = pot_lj + pot_str + e_lrc
       call pressure(pres1,PLRC)
       write(6,*)vpot,density,pres1,box !PLJ(2,density)
     endif
  enddo
  write(*,*)' Equilibration is finished. '
  write(*,*)' istep = ',istep
  write(*,*)' nequil = ',nequil
  write(*,*)' % accepted config moves is = ',100*dfloat(istep)/dfloat(nequil-kstep)
  write(*,*)' % accepted volume moves is = ',100*dfloat(jstep)/dfloat(kstep)
  write(*,*)' These % can be changed by varying maximum moves allowed in input file.'
  write(*,*) 
  allocate (dens(nblocks))
  allocate (pres_block(nblocks))
  allocate (PLRC_block(nblocks))
  write(*,*)' Begin calculations of the block averages. '
  write(*,*)' The number of blocks is ',nblocks 
  write(*,*)' The number of points per block is ',irun
  write(*,*)' vpot and density are output to screen'
  pause
  istep = 0
  jstep = 0
  kstep = 0
  do j = 1,nblocks
  write(6,*)' density     vpot    volume    pressure     e_lrc    ecut'
    iflag = 0
    dens_avg = 0
    vavg = 0
    pavg = 0
    PLRCavg = 0
    do i = 1,irun
     call random_number(rnd)        !pick particle to be moved
     rnd = rnd*(npart+1)
     if(rnd<=npart+const_vol)then
      call mcmove(vpot,istep) !attempts to displace a particle
     else
!      call pressure(pres)	!
      call mcvol(vpot,jstep,kstep) !attempts to displace a particle
     endif
     if(mod(i,nsamp)==1)then
       call center
       call potential    !calculate the potential energy
	vpot = pot_lj + pot_str + e_lrc
       call pressure(pres1,PLRC)
       write(6,'(7f14.8)')density,vpot,volume,pres1,e_lrc,ecut
     endif
     vavg = vavg + vpot
     dens_avg = dens_avg + density
     pavg = pavg + pres1
     PLRCavg = PLRCavg + PLRC
    enddo
    dens_avg = dens_avg/irun
    vavg = vavg/irun
    pavg = pavg/irun
    PLRCavg = PLRCavg/irun
    iflag = 2
    write(*,'(a,f8.4)')' the temperature is ',tstar
    write(*,'(a,f8.4)')' the average density for this block is ',dens_avg
    write(*,'(a,f8.4)')' the average pressure for this block is ',pavg
    write(*,'(a,f8.4)')' the average LRC contribution on pressure for this block is ',PLRCavg
    dens(j) = dens_avg
    pres_block(j) = pavg
    PLRC_block(j) = PLRCavg
  enddo
  dens_avg = sum(dens)/dfloat(nblocks)
  dens_std= dsqrt(sum((dens-dens_avg)**2)/dfloat(nblocks-1))
  pavg = sum(pres_block)/dfloat(nblocks)
  pstd= dsqrt(sum((pres_block-pavg)**2)/dfloat(nblocks-1))
  PLRCavg = sum(PLRC_block)/dfloat(nblocks)
  write(*,*)
  write(*,'(a,f7.4)')' The average density of the blocks is ',dens_avg
  write(*,'(a,f7.4)')' The sdt of the the block densities is ',dens_std
  write(*,'(a,f7.4)')' The uncertainty of the the block densities is about ',dens_std/dsqrt(dfloat(nblocks))
  write(*,'(a,f7.4)')' The average pressure of the blocks is ',pavg
  write(*,'(a,f7.4)')' The sdt of the the block pressures is ',pstd
  write(*,'(a,f7.4)')' The uncertainty of the the block pressures is about ',pstd/dsqrt(dfloat(nblocks))
  write(*,'(a,f7.4)')' The average LRC contribution of pressure of the blocks is ',PLRCavg
!----------------------------------------------------------------
  write(16,*)' istep = ',istep
  write(16,*)' total run  = ',irun*nblocks
  write(16,*)' % accepted coord is        = ',100*dfloat(istep)/dfloat(irun*nblocks-kstep)
  write(16,*)' % accepted volume moves is = ',100*dfloat(jstep)/dfloat(kstep)
  write(16,*)' the final volume is = ',volume
  write(16,'(a,f14.8)')' the average density is ',dens_avg
  write(16,'(a,f14.8)')' the average pressure is ',pavg
  write(16,'(a,f14.8)')' the sdt of density is ',dens_std
  write(16,'(a,f14.8)')' The sdt of the the block pressures is ',pstd

! save final configuration, one file for povray output
  open(18,file='old_config.ini')
  open(19,file='config_pov.dat') 
     write(19,'i5')npart
     write(19,*)' '
     write(18,'(f12.6,a)')box,'  ! box length'
    write(18,*)x 
    write(19,99)((x(i,j),j=1,3),i=1,npart) 
 99 format('X    ',3(f10.3))
 close(18)
 close(19)
! ----------------------------------------------------------------    
  end program lj_met_chain

