  program lj_metro
  use lj_metro_mod
  implicit none
!  --------------------------------------------------------------------------
  integer            :: irun,nequil,iflag  !number of steps used used in Monte Carlo 
  integer            :: i,j                      !counters
  integer            :: istep=0                  !counter of successful moves
  real(kind=8)       :: ptot        !ptot is instantaneous pressure
  real(kind=8)       :: pavg,vavg   !average calculated temp and pressure 

  real*8, dimension(npart,3) :: x
  real*8                     :: vpot,temp,pdel
  real*8, allocatable, dimension(:)       :: ppp   !average pressure
! ----------------------------------------------------------------
! read in the input parameters
  call input(irun,nequil)
  print*,' irun = ',irun
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------
! set up the initial positions assuming fcc
  call init(x) !set up an initial condition
! ----------------------------------------------------------------
  do i = 1,nequil
     call mcmove(x,vpot,istep) !attempts to displace a particle
     if(mod(i,nsamp)==1)then
       call center(x)
     endif
  enddo
  write(*,*)' istep = ',istep
  write(*,*)' nequil = ',nequil
  write(*,*)' the percentage of accepted moves is = ',dfloat(istep)/dfloat(nequil)
  write(*,*) ' type enter to continue'
  read(5,*) 
  allocate (ppp(5))
  do j = 1,5
    iflag = 0
!   call gr(x,iflag)  !initialize
    pavg = 0
    vavg = 0
    write(*,*)'      vpot           temp        pressure   '
    do i = 1,irun
  
     call mcmove(x,vpot,istep) !attempts to displace a particle
     if(mod(i,nsamp)==1)then
       call center(x)
       call pressure(ptot,x)
       write(6,'(2f14.8)')vpot,ptot
       pavg = pavg + ptot
     endif
!    if(mod(i,nsamp)==1)call gr(x,iflag)

     vavg = vavg + vpot
    enddo
    pavg = pavg/ncalls
    vavg = vavg/irun
    iflag = 2
!   call gr(x,iflag)  !finalize
    write(*,'(a,f8.4)')' the temperature is ',tstar
    write(*,'(a,f8.4)')' the average potential is ',vavg/2.d0
    write(*,'(a,f8.4)')' the average pressure is    ',pavg
    write(*,'(a,f10.4)')' PV = ', pavg*volume
    write(*,'(a,f8.4)')' nT = ', npart*tstar
    ppp(j) = pavg
  enddo
  write(*,'(5f10.5)')ppp
  pavg = sum(ppp)/5
  pdel = dsqrt(sum((ppp-pavg)**2/4))
  print*,pavg,pdel,irun

! ----------------------------------------------------------------
 ! iflag = 2
 ! call kebin(v,iflag)

    
  end program lj_metro




! ----------------------------------------------------------------
  subroutine init(x) !set up an initial condition using fcc
  use lj_metro_mod
  implicit none
  real*8, dimension(npart,3) :: x
  integer :: i,j,k,l,ipart,nt
  real*8  :: rnd,delta
! pick points on a fcc cubic lattice
   
  ipart  = (npart/4)**(1.0/3.0)
  delta = box(1)/ipart
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
!    write(6,*)' pausing since nt is not equal to npart '
!    read(*,'()')
  endif
! ------------------------------------------------------
  ecut = 4.0*(1/rc2**6-1/rc2**3)  !cutoff value for the LJ
  end subroutine init



!-----------------------------------------------------------------------------------
  subroutine potential(i,x,eno)   !calculate the potentials for the LJ potential
  use lj_metro_mod
  implicit none
  real*8, dimension(npart,3) :: x
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
         eno =  eno+4*r6i*(r6i-1)-ecut
       endif
     endif
  enddo
  end subroutine potential

!-----------------------------------------------------------------------------------
  subroutine center(x)   !keep particles in box
  use lj_metro_mod
  implicit none
  real*8, dimension(npart,3) :: x
  integer                    :: i,j
  do i = 1,npart
   do j = 1,3
     if(x(i,j)>=box(j))then
       x(i,j) = x(i,j)-box(j)
     endif
     if(x(i,j)<=0)then
       x(i,j) = x(i,j)+box(j)
     endif
   enddo
  enddo
  end subroutine center


  subroutine pressure(p,x)  !p is the pressure
  use lj_metro_mod
  implicit none
  real*8, dimension(npart,3) :: x
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r2,r6i,p
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
  p = (npart*tstar+16*p)/volume
  end subroutine pressure

  subroutine input(irun,nequil)
    use lj_metro_mod
    implicit none
    integer      :: irun,nequil
    open(15,file='lj_metro.inp')
!   ----------------------------------------------------------------
!   read in the data file
    read(15,*)irun      ! number of points for the simulation
    read(15,*)nequil    ! number of points in equilibration
    read(15,*)nsamp     ! number of points between calls to boltmann or gr 
    read(15,*)delx      ! integration step size
    read(15,*)tstar     ! the desired temperature
    read(15,*)density   ! the value of the density 
    read(15,*)seed      ! used for initializing random number generator
    close (15)
    volume  = npart/density
    box(1) = volume**(1.d0/3.d0)
    box(2) = box(1) 
    box(3) = box(1)
    beta = 1.d0/tstar
    rcut = box(1)/2             !the cut-off on the LJ forces is rcut
    ncalls = irun/nsamp         !number of calls to gr
    rc2 = rcut*rcut
    write(6,'(a,f8.3)')' the desired temp is    ',tstar
    write(6,'(a,f8.3)')' the box length is      ',box(1)
    write(6,'(a,i8)')  ' number of particles is ',npart
    write(6,'(a,f8.3)')' the density is         ',density
    write(6,'(a,f8.3)')' r cutoff is            ',rcut
    write(6,*)' hit enter to continue'
    read(*,'()')
! ----------------------------------------------------------------
  end subroutine input

!-----------------------------------------------------------------------------------
  subroutine gr(x,iflag)
  use lj_metro_mod
  implicit none
  integer, parameter                 :: num_bins=200000 
  integer                            :: iflag,jbin,i,j
  real*8, save,  dimension(num_bins) :: gn
  real*8,  dimension(npart,3)        :: x
  real*8,  dimension(3)              :: dr
  real*8, save                       :: rs,rmin,rmax,rdel,ridel,rb,const,rinner,router



! --------------------------------------------------------------------
  select case(iflag)
! ------------------------------------------------
  case(0)
   gn = 0   !initialization
   rmin = 0.8d0
   rmax = 0.5d0*box(1)
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
  subroutine mcmove(x,vpot,istep) !attempts to displace a particle
    use lj_metro_mod
    implicit none
    integer                          :: j,o,istep
    real(kind=8), dimension(3)       :: xn
    real(kind=8)                     :: vpot,eno,enn,rnd
    real(kind=8), dimension(npart,3) :: x

    call random_number(rnd)        !pick particle to be moved
    o=int(rnd*npart)+1
    call potential(o,x,eno)  !calculate the potential energy
   
    xn = x(o,:)                         !current position

    do j = 1,3
      call random_number(rnd)
      x(o,j) = x(o,j) + (rnd-0.5)*delx  !new position
    enddo

    call potential(o,x,enn)  !calculate the potential energy
    call random_number(rnd)
    if(rnd<exp(-beta*(enn-eno)))then  !decide if you want to make a move
           vpot = enn
           istep = istep + 1  !count on accepted moves.
    else
           x(o,:)=xn
    endif
  end subroutine mcmove
