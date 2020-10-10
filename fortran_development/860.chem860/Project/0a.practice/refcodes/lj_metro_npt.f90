  program lj_metro_npt
  use lj_metro_npt_mod
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
! read in the input parameters
  call input(irun,nequil)
  write(*,*)' irun = ',irun
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------
! set up the initial positions assuming fcc
  call init !set up an initial condition
  call potential_full(vpot)   !calculate the potentials for the LJ potential
  write(*,*)' The initial potential is',vpot
  write(*,*)' Start running equilibration.  Vpot and density will be plotted.'
  pause
! ----------------------------------------------------------------
  do i = 1,nequil
     call random_number(rnd)        !pick angle to be moved
     rnd = rnd*(npart+1)
     if(rnd<=npart)then
      call mcmove(vpot,istep)      !attempts to displace a particle
     else
!      call pressure(pres)        !
      call mcvol(vpot,jstep,kstep) !attempts a volume move
     endif
     if(mod(i,nsamp)==1)then
       call center
       call potential_full(vpot)    !calculate the potential energy
       call pressure(pres1,PLRC)
       write(6,'(5f14.8)')vpot,density,pres1,box !PLJ(2,density)
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
  do j = 1,nblocks
    iflag = 0
    dens_avg = 0
    vavg = 0
    pavg = 0
    PLRCavg = 0
    write(*,*)'      vpot        density   '
    do i = 1,irun
  
     call random_number(rnd)        !pick angle to be moved
     rnd = rnd*(npart+1)
     if(rnd<=npart)then
      call mcmove(vpot,istep) !attempts to displace a particle
     else
!      call pressure(pres)	!
      call mcvol(vpot,jstep,kstep) !attempts to displace a particle
     endif
     if(mod(i,nsamp)==1)then
       call center
       call potential_full(vpot)    !calculate the potential energy
       call pressure(pres1,PLRC)
       write(6,'(5f14.8)')vpot,density,pres1,box,PLRC
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

! ----------------------------------------------------------------

    
  end program lj_metro_npt




! ----------------------------------------------------------------
  subroutine init !set up an initial condition using fcc
  use lj_metro_npt_mod
  implicit none
  integer :: i,j,k,l,ipart,nt
  real*8  :: rnd,delta
! pick points on a fcc cubic lattice
   
  ipart  = (0.250001*dfloat(npart))**(1.0/3.0)
  delta = box/ipart
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nt = 0
 !put in the corners
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delta*i
     x(nt,2) = delta*j
     x(nt,3) = delta*k
    enddo
   enddo
  enddo
! put in the z faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delta*i+0.5d0*delta
     x(nt,2) = delta*j+0.5d0*delta
     x(nt,3) = delta*k
    enddo
   enddo
  enddo
! put in the x faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delta*i
     x(nt,2) = delta*j+0.5d0*delta
     x(nt,3) = delta*k+0.5d0*delta
    enddo
   enddo
  enddo
! put in the y faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delta*i+0.5d0*delta
     x(nt,2) = delta*j
     x(nt,3) = delta*k+0.5d0*delta
    enddo
   enddo
  enddo
  if(nt/=npart)then
     write(6,*)' pausing since nt is not equal to npart '
     read(*,'()')
  endif
! ------------------------------------------------------
  ecut=8.0D0*pi*density*((1/rcut)**9/3.0D0-(1/rcut)**3)/3.0D0  !cutoff value for the LJ
  end subroutine init


!-----------------------------------------------------------------------------------
  subroutine potential_full(eno)   !calculate the potentials for the LJ potential
  use lj_metro_npt_mod
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: eno,r2i,r6i
  real(kind=8)               :: xr,r2
  integer                    :: i,j    !counters

  eno = 0
  do i = 1,npart-1
   do j = i+1,npart
       dr = x(i,:)-x(j,:)
       dr = dr-box*nint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
       if(r2<=rc2)then
         r2i = 1/r2
         r6i=r2i**3
         eno =  eno+4*r6i*(r6i-1)
         endif
   enddo
  enddo
  eno = eno-ecut*dble(npart)
  end subroutine potential_full

!-----------------------------------------------------------------------------------
  subroutine potential(i,eno)   !calculate the potentials for the LJ potential
  use lj_metro_npt_mod
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: eno,r2i,r6i
  real(kind=8)               :: xr,r2
  integer                    :: i,j    !counters

  eno = 0
  do j = 1,npart
     if(j/=i)then
       dr = x(i,:)-x(j,:)
       dr = dr-box*nint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
       if(r2<=rc2)then
         r2i = 1/r2
         r6i=r2i**3
         eno =  eno+4*r6i*(r6i-1)
       endif
     endif
  enddo
  eno = eno-ecut*dble(npart)
  end subroutine potential

!-----------------------------------------------------------------------------------
  subroutine center   !keep particles in box
  use lj_metro_npt_mod
  implicit none
  integer                    :: i,j
  do i = 1,npart
   do j = 1,3
     if(x(i,j)>=box)then
       x(i,j) = x(i,j)-box
     endif
     if(x(i,j)<=0)then
       x(i,j) = x(i,j)+box
     endif
   enddo
  enddo
  end subroutine center


  subroutine pressure(p,PLRC)  !p is the pressure
  use lj_metro_npt_mod
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r2,r6i,p,rc3it,rc9it,PLRC
  integer                    :: i,j    !counters
  p = 0
  do i = 1,npart-1
   do j = i+1,npart
     dr = x(i,:)-x(j,:)
     dr = dr-box*nint(dr/box)  !periodic boundary conditions
     r2 = dot_product(dr,dr)
     if(r2<=rc2)then
       r6i=(1.d0/r2)**3
       p = p + r6i*(r6i-0.5)  !LJ potential
     endif
   enddo
  enddo
  rc3it = (1.d0/rcut)**3/2.d0
  rc9it = (1.d0/rcut)**9/3.d0
  PLRC = (32*pi*density*dble(npart)*(rc9it-rc3it)/3.0D0)/volume
  p = (npart*tstar+16*p)/volume+PLRC	!P LRC
  end subroutine pressure

  subroutine input(irun,nequil)
    use lj_metro_npt_mod
    implicit none
    integer      :: irun,nequil
    open(15,file='lj_metro_npt.inp')
!   ----------------------------------------------------------------
!   read in the data file
    read(15,*)npart     ! number of particles
    read(15,*)irun      ! number of points for the simulation
    read(15,*)nequil    ! number of points in equilibration
    read(15,*)nsamp     ! number of points between calls to boltmann or gr 
    read(15,*)delx      ! mc config step size
    read(15,*)vmax      ! mc volume step size
    read(15,*)tstar     ! the desired temperature
    read(15,*)pres      ! the value of the pressure 
    read(15,*)density   ! the initial density 
    read(15,*)seed      ! used for initializing random number generator
    read(15,*)nblocks   ! number of blocks used for averaging the density
    close (15)
    volume  = npart/density
    box = volume**(1.d0/3.d0)
    beta = 1.d0/tstar
    rcut = box/2             !the cut-off on the LJ forces is rcut
    ncalls = irun/nsamp         !number of calls to gr
    allocate(x(npart,3),xm(npart,3))
    rc2 = rcut*rcut
    write(6,'(a,f8.3)')' the desired temp is    ',tstar
    write(6,'(a,f8.3)')' the box length is      ',box
    write(6,'(a,i8)')  ' number of particles is ',npart
    write(6,'(a,f8.3)')' the initial density is ',density
    write(6,'(a,f8.3)')' the pressure is        ',pres
    write(6,'(a,f8.3)')' r cutoff is           ',rcut
    write(6,*)' hit enter to continue'
    read(*,'()')
! ----------------------------------------------------------------
  end subroutine input

!-----------------------------------------------------------------------------------
  subroutine gr(iflag)
  use lj_metro_npt_mod
  implicit none
  integer, parameter                 :: num_bins=200000 
  integer                            :: iflag,jbin,i,j
  real*8, save,  dimension(num_bins) :: gn
  real*8,  dimension(3)              :: dr
  real*8, save                       :: rs,rmin,rmax,rdel,ridel,rb,const,rinner,router



! --------------------------------------------------------------------
  select case(iflag)
! ------------------------------------------------
  case(0)
   gn = 0   !initialization
   rmin = 0.8d0
   rmax = 0.5d0*box
   ridel = num_bins/(rmax-rmin)
   rdel = 1.d0/ridel
   const = 4.d0/3.d0*pi*density
   iflag = 1
! ------------------------------------------------
  case(1)
   do i = 1,npart-1
    do j = i+1,npart
      dr = x(i,:)-x(j,:)
      dr = dr-box*nint(dr/box)  !periodic boundary conditions
      rs = dsqrt(dot_product(dr,dr))
      jbin =  1+(rs-rmin)*ridel
      if(jbin<=num_bins) gn(jbin) = gn(jbin) + 2
    enddo
   enddo
! ------------------------------------------------
  case(2)     !  carry out on the last call
    open(16,file='gr.out')
    gn = gn/(ncalls*npart)
    rb = rmin+0.5*rdel
    do i = 1,num_bins
     rb = rb + rdel
     rinner = rb-0.5d0*rdel 
     router = rinner+rdel
     gn(i) = gn(i)/(const*(router**3-rinner**3))
     write(16,'(f12.7,f13.8)')rb,gn(i)
    enddo
    close(16)
  end select
!  --------------------------------------------------------------------
  end subroutine gr

! ----------------------------------------------------------------
  subroutine mcmove(vpot,istep) !attempts to displace a particle
    use lj_metro_npt_mod
    implicit none
    integer                          :: j,o,istep
    real(kind=8), dimension(3)       :: xn
    real(kind=8)                     :: vpot,eno,enn,rnd

    call random_number(rnd)        !pick particle to be moved
    o=int(rnd*npart)+1
    call potential(o,eno)  !calculate the potential energy
   
    xn = x(o,:)                         !current position

    do j = 1,3
      call random_number(rnd)
      x(o,j) = x(o,j) + (rnd-0.5)*delx  !new position
    enddo

    call potential(o,enn)  !calculate the potential energy
    call random_number(rnd)

    if(rnd<exp(-beta*(enn-eno)))then  !decide if you want to make a move
           istep = istep + 1  !count on accepted moves.
    else
           x(o,:)=xn
    endif
  end subroutine mcmove


  subroutine mcvol(vpot,jstep,kstep) !attempts to displace a particle
    use lj_metro_npt_mod
    implicit none
    integer                          :: j,i,jstep,kstep
    real(kind=8)                     :: rnd       !random number
    real(kind=8)                     :: vpot,eno,enn,vo,lnvn,vn,boxo,arg
 !  -------------------------------------------------------------
 !  save the old values
    xm = x 
    boxo = box
    vo = box**3
 !  -------------------------------------------------------------
    call potential_full(eno)     !calculate the potential energy
    call random_number(rnd)        !pick molecule to be moved
    lnvn=dlog(vo)+(rnd-0.5d0)*vmax
    vn = dexp(lnvn)
    box = vn**(1.d0/3.d0)         !new box length
!   --------------------------------------------------------------
    x = box/boxo*xm               !scale the coordinates
!   --------------------------------------------------------------
    call potential_full(enn)    !calculate the potential energy
!   --------------------------------------------------------------
    arg = -beta*((enn-eno)+pres*(vn-vo)-(npart+1)*dlog(vn/vo)/beta)
    call random_number(rnd)       
    kstep = kstep + 1
    if(rnd>dexp(arg))then
     box = boxo
     volume = vo
     x = xm                        !move is rejected
    else
     density = npart/vn
     jstep = jstep + 1
     ecut=8.0D0*pi*density*((1/rcut)**9/3.0D0-(1/rcut)**3)/3.0D0  !cutoff value for the LJ

    endif
  end subroutine mcvol



