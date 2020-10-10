  program lj
  use lj_ewald_mod
  implicit none
!  --------------------------------------------------------------------------
  integer            :: irun,itest,nequil,iflag  !time steps used used in MD 
  integer            :: i,j                      !counters
  real(kind=8)       :: ptot        !ptot is instantaneous pressure
  real(kind=8)       :: tavg,pavg   !average calculated temp and pressure 
  integer            :: ktemp       !for temperature bath
  real(kind=8)       :: lambda,tau  !used in temperature bath
  real(kind=8), dimension(npart,3) :: x,v,xm,f,fewald,frwald
  real*8, dimension(npart)   :: z 
  real*8                     :: ecut,vpot,etot,temp,vr,vk,vkp,vkm
  real*8                     :: delta !used for finite differences
! ----------------------------------------------------------------
! read in the input parameters
  call input(irun,itest,nequil)
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------
! set up the initial positions assuming fcc
  call init(x,v,xm,z,ecut) !set up an initial condition
  call ewald_setup
! ----------------------------------------------------------------
! Problem 3 - determine an appropriate value for kappa
  write(*,'(a,d12.6)')' erfc(2*rcut)/rcut     = ',erfc(2*rcut)/rcut
  write(*,'(a,d12.6)')' erfc(3*rcut)/rcut     = ',erfc(3*rcut)/rcut
  write(*,'(a,d12.6)')' erfc(kappa*rcut)/rcut = ',erfc(kappa*rcut)/rcut
! ----------------------------------------------------------------
! Problem 5 - check the results for the derivatives.
  delta = .002d0  ! variable used to test finite difference.  It is currently too big.
  x(12,2) = x(12,2)+.01d0  !move some arbitray atoms so we are not at a minimum.
  x(16,1) = x(16,1)+.01d0
  x(12,3) = x(16,3)+.02d0
  j = 3  !1-3 for x,y and z
  call kwald (x,z,vk,fewald)   ! long range
  call rwald (x,z,vr,frwald)   ! short range
  write(*,'(a,f10.3)')' the value of kappa is   ',kappa
  write(*,'(a,f12.6)')' the analytical force is ',frwald(12,j)+fewald(12,j)
  write(*,*)' These should be the same '
	write(*,'(a,3f12.6)')' LR,SR ewald pot, and LR+SR at initial : ', vk, vr, vk+vr
	write(*,'(a,f12.6)') 'Now, moving the particle no.12 which has charge ', z(12)
  x(12,j) = x(12,j)+delta   !move atoms for finite difference
  call kwald (x,z,vk,fewald)
  call rwald (x,z,vr,frwald)
  vkp = vr+vk
	write(*,'(a,f12.6)')' For positive displacement, ewald LR+SR=vkp is? ', vkp
	write(*,'(a,2f12.6)')'After positive displacement, (LR pot, SR pot) = ', vk, vr
  x(12,j) = x(12,j)-2.d0*delta
  call kwald (x,z,vk,fewald)   !problem 5 fewald is the force
  call rwald (x,z,vr,frwald)   !problem 5 frwald is the force
  vkm = vr+vk
        write(*,'(a,f12.6)')' For negative displacement, ewald LR+SR=vkm is? ', vkm
        write(*,'(a,2f12.6)')'After negative displacement, (LR pot, SR pot) = ', vk, vr
  write(*,'(a,f12.6)')' numerical force is      ',(vkp-vkm)/(2.d0*delta)
 

  end program lj




! ----------------------------------------------------------------
  subroutine init(x,v,xm,z,ecut) !set up an initial condition using fcc
  use lj_ewald_mod
  implicit none
  real*8, dimension(npart,3) :: x,v,xm
  real*8, dimension(npart)   :: z 
  integer :: i,j,k,l,ipart,nt
  real*8  :: delx,dely,delz,vx,norm,rnd,fs,ecut
! pick points on a bcc cubic lattice
   
  print*,' npart = ',npart

  ipart  = (float(npart/2)+.00001)**(1.0/3.0)
  print*,' the number of cells in each direction is ',ipart
! ------------------------------------------------------
! check that number is correct
  if((ipart**3)*2/=npart)then
    print*,2*ipart**3
    print*,' npart = ',npart
    write(*,*)'check the value of npart'
    stop
  endif
  delx = box/ipart
  dely = box/ipart
  delz = box/ipart
! ------------------------------------------------------
 nt = 0
 !put in the corners
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i
     x(nt,2) = dely*j
     x(nt,3) = delz*k
     z(nt) = 1.0
    enddo
   enddo
  enddo
! put in the box centers
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i+0.5d0*delx
     x(nt,2) = dely*j+0.5d0*dely
     x(nt,3) = delz*k+0.5d0*delz
     z(nt) = -1.d0
    enddo
   enddo
  enddo
  if(nt/=npart)then
     write(6,*)' pausing since nt is not equal to npart '
     read(*,'()')
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! here are the velocities
    do i = 1,npart
   do l = 1,3
     call random_number(rnd)
     v(i,l) = rnd-0.5d0   !give random velocities
   enddo
  enddo
! make sure average velocities are zero
  do l = 1,3
    vx = sum(v(:,l))/npart
    v(:,l) = v(:,l)-vx
  enddo 
! calculate the scale factor for adjusting ke
  fs = dsqrt(3*tstar*npart/(sum(v*v)))
  v = fs*v                           !rescale velocities
  xm = x - v*dt
! ------------------------------------------------------
  ecut = 4.0*(1/rc2**6-1/rc2**3)  !cutoff value for the LJ
  end subroutine init



!-----------------------------------------------------------------------------------
  subroutine force(x,f,vpot)   !calculate the forces for the LJ potential
  use lj_ewald_mod
  implicit none
  real*8, dimension(npart,3) :: x,f
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: ecut,vpot,r2i,r6i,ff
  real(kind=8)               :: xr,r2
  integer                    :: i,j    !counters
  vpot = 0
  f = 0
  do i = 1,npart-1
   do j = i+1,npart
     dr = x(i,:)-x(j,:)
     dr = dr-box*nint(dr/box)  !periodic boundary conditions
     r2 = dot_product(dr,dr)
     if(r2<=rc2)then
       r2i = 1/r2
       r6i=r2i**3
       ff=48*r2i*r6i*(r6i-0.5)  !LJ potential
       f(i,:) = f(i,:) + ff*dr
       f(j,:) = f(j,:) - ff*dr
       vpot =  vpot+4*r6i*(r6i-1)-ecut
     endif
   enddo
  enddo
  end subroutine force

!-----------------------------------------------------------------------------------
  subroutine integrate(f,x,xm,vi,etot,vpot,temp,lambda)   !integrate equations of motion
  use lj_ewald_mod
  implicit none
  real*8, dimension(npart,3) :: xx,x,f,xm,vi
  real*8                     :: sumv2, dt2, dit,temp, etot, vpot,lambda
  dt2 = dt*dt
  dit = 1/(2*dt)
  xx = 2*x-xm+dt2*f
  vi = (xx-xm)*dit*lambda
  xx = vi/dit+xm
  sumv2=sum(vi*vi)
  xm = x; x = xx;  temp = sumv2/(3*npart)
  etot=(vpot+0.5*sumv2)/npart   !total energy per particle
  end subroutine integrate




!-----------------------------------------------------------------------------------
  subroutine center(x,xm)   !keep particles in box
  use lj_ewald_mod
  implicit none
  real*8, dimension(npart,3) :: xm,x
  integer                    :: i,j
  do i = 1,npart
   do j = 1,3
     if(x(i,j)>=box)then
       x(i,j) = x(i,j)-box
       xm(i,j) = xm(i,j)-box
     endif
     if(x(i,j)<=0)then
       x(i,j) = x(i,j)+box
       xm(i,j) = xm(i,j)+box
     endif
   enddo
  enddo
  end subroutine center

!-----------------------------------------------------------------------------------
  subroutine kebin(v,tavg,iflag)   !bin the speeds
  use lj_ewald_mod
  implicit none
  integer                         :: iflag,jbin,i
  real*8, save,  dimension(nbins) :: rbin ,bn
  real*8,  dimension(npart,3)     :: v
  real*8,  dimension(npart)       :: rs
  real*8, save                    :: rmin,rmax,rdel,ridel
  real*8                          :: anorm,an,rr,tavg,sig
! --------------------------------------------------------------------
  select case(iflag)
   case(0)                                        !initializtion step
     bn = 0   
     rmin = 0.d0
     rmax = dsqrt(15.0d0*tstar)
     ridel = nbins/(rmax-rmin)
     rdel = 1.d0/ridel
     iflag = 1
   case(1)
     rs = dsqrt((v(:,1)**2+v(:,2)**2+v(:,3)**2))  !calculate the speed and bin
     do i = 1,npart
      jbin =  1+(rs(i)-rmin)*ridel
      if(jbin<=nbins) bn(jbin) = bn(jbin) + 1
     enddo
   case(2)                                        ! carry out on the last call
!    ------------------------------------------------
!    normalize and print the binned distribution
     bn = bn/(mcalls*npart*rdel)
     rbin(1) = rmin+0.5*rdel
     do i = 2,nbins
      rbin(i) = rbin(i-1) + rdel
     enddo
!    ------------------------------------------------
!    here is the Boltman distribution
     open(17, file='bolt.out')
     anorm = dsqrt(2.d0/(pi*tavg**3))
     write(*,*)'   Comparisons '
     write(*,*)'   Speed  Boltzmann    MD    Difference '
     sig = 0.d0
     do i = 1,nbins
       rr = rbin(i)
       an = anorm*rr*rr*dexp(-0.5d0*rr*rr/tstar)  
       write( 6,'(4f10.5)')rr,an,bn(i),an-bn(i)
       write(17,'(4f10.5)')rr,an,bn(i),an-bn(i)
       sig = sig + (an-bn(i))**2
     enddo
     sig = dsqrt(sig/(nbins-1)) !difference between MD and Boltzman
     write(*,*)' difference between Boltzman and MD is sig = ',sig
     write(*,*)' type enter to continue' 
     read(*,'()')
     close(17)
!    ------------------------------------------------
   end select
! --------------------------------------------------------------------
  end subroutine kebin

  subroutine pressure(p,temp,x)  !p is the pressure
  use lj_ewald_mod
  implicit none
  real*8, dimension(npart,3) :: x
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r2,r6i,temp,p
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
  p = (npart*temp+16*p)/volume
  end subroutine pressure

  subroutine input(irun,itest,nequil)
    use lj_ewald_mod
    implicit none
    integer      :: irun,itest,nequil
    real(kind=8) :: nkappa
    open(15,file='lj_ewald.inp')
!   ----------------------------------------------------------------
!   read in the data file
    read(15,*)trun      ! runtime of the simulation
    read(15,*)tequil    ! run_time of initial equilibration with T control.
    read(15,*)time_test ! run_time of the test for bolztmann statistics.
    read(15,*)tsample   ! time between calls to boltmann or gr 
    read(15,*)dt        ! integration step size
    read(15,*)tstar     ! the desired temperature
    read(15,*)density   ! the value of the density 
    read(15,*)nkappa    ! used to setup kappa below
	read(15,*)kmax
	read(15,*)KSQMAX
    read(15,*)seed      ! used for initializing random number generator
    close (15)
    volume  = npart/density
    box = volume**(1.d0/3.d0)
    kappa =  nkappa/box
    rcut = box/2             !the cut-off on the LJ forces is rcut
    irun = trun/dt              !time steps in the simulation
    itest = time_test/dt        !time steps for testing Boltzmann
    nequil = tequil/dt          !time steps for equilibraiton
    nsamp = tsample/dt          !number of steps before a call to gr
    ncalls = trun/tsample       !number of calls to gr
    mcalls = time_test/tsample  !number of calls to speedbin
    rc2 = rcut*rcut
    write(*,*)' HERE ARE THE INPUT VALUES'
    write(*,'(a,f8.3)')' the desired temp is    ',tstar
    write(*,'(a,f8.3)')' the box length is      ',box
    write(*,'(a,i8)')  ' number of particles is ',npart
    write(*,'(a,f8.3)')' the density is         ',density
    write(*,'(a,f8.3)')' r cutoff is            ',rcut
    write(*,'(a,f8.3)')' kappa is               ',kappa
    write(*,*)' For problem 5 you need to make sure you get an optimal kappa value.'
    write(*,*)' hit enter to continue'
    read(*,'()')
! ----------------------------------------------------------------
  end subroutine input

!-----------------------------------------------------------------------------------
  subroutine gr(x,iflag)
  use lj_ewald_mod
  implicit none
  integer, parameter                 :: num_bins=60 
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
     write(16,'(f10.5,f10.5)')rb,gn(i)
    enddo
    close(16)
  end select
!  --------------------------------------------------------------------
  end subroutine gr



