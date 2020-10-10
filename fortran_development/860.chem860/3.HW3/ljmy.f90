! program for lennard-jones particles in a box, developed for chem 860
! ----------------------------------------------------------------
  program lj
  use lj_mod
  use pressure_mod
  implicit none
  integer                                       :: i,j,k           !counter
  integer                                       :: iflag         !used in calls to kebin
  real(kind=8)                                  :: vpot,etot,sumv,temp,msd
  real(kind=8)                                  :: tavg !the average temperature
  real(kind=8)                                  :: vavg !the average potential 
! ----------------------------------------------------------------
! these varibles are used for pressure calculations
  real(kind=8), allocatable, dimension(:)       :: ppp   !average pressure
  real(kind=8)                                  :: pavg,volume,ptot,pdel
! ----------------------------------------------------------------
  call input   !read in the parameters that specify the simulation 
  call init    !set up an initial condition
! ----------------------------------------------------------------
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------
! lambda is used to get the bath to the desired temperature tstar
! tau is used to adjust the amount of the coupling
  tau = 2 
  temp = tstar   !first guess for temperature, used in lambda def.
! ----------------------------------------------------------------
! allow some equilibration
  do i = 1,nequil
   if(i>=nequil/2)tau=10   !reduce the size of the coupling half way through
   lambda = sqrt( 1.0 + (dt/tau)*(tstar/temp-1.0) )  
   call force(vpot)
   call integrate(etot,vpot,sumv,temp)   !integrate equations of motion
   if(mod(i,nsamp)==1.OR.nsamp==1)call center        !put molecules back into box if moved out of it
  enddo
  write(*,*)' finished initial equilibration'

! ---------------------------------------------------------------------------------------------------
! test for equilibration by checking distribution of velocities.
  write(*,*)' test for boltzmann statistics'
  iflag = 0
  tavg = 0
  tau = 100   !there is a weak thermostat
  lambda = 1.d0  !there is no thermostat
  call kebin(iflag)
  itest = ttest/dt
  do i = 1,itest
     call force(vpot)
     call integrate(etot,vpot,sumv,temp)   !integrate equations of motion
     tavg = tavg + temp
     if(mod(i,nsamp)==1.OR.nsamp==1)then
       call center(x,xm)
       call kebin(iflag)
     endif
  enddo
  tavg = tavg/itest
  iflag = 2
  beta = 1/tavg
  call kebin(iflag)  !bins kinetic energy 
  write(*,'(a,f10.5)')' The temperature is ',tavg
  write(*,*)' the value of beta is being reset to ',beta
! -----------------------------------------------------------------------------------------------------
! here I carry out nblock simulations in order to calculate block averages for the pressure
! the first do-loop is over the blocks.
  allocate (ppp(nblock))
  volume = box(1)*box(2)*box(3)

  do j = 1,nblock

   if(j < 10)then 
	open(17,file='vel_'//char(48+j)//'.out')
   else if(j>=10) then
	open(17,file='vel_'//char(48+j/10)//char(48+mod(j,10))//'.out')
   endif
   if(j < 10)then
	open(18,file='msd_'//char(48+j)//'.out')
   else if(j>=10) then
	open(18,file='msd_'//char(48+j/10)//char(48+mod(j,10))//'.out')
   endif
   write(17,'(2i,f14.8)')npart,INT(nrun/nsamp),tsample
   iflag = 0
   call gr(iflag,j,npart,density,box,x)  !initialize
	x0 = x !!!To calculate MSD, store the initial position of that block
   tavg = 0
   pavg = 0
   vavg = 0
   write(*,*)'      etot           temp        pressure   '
   ncalls = 0
   do i = 1,nrun
    call force(vpot)
    call integrate(etot,vpot,sumv,temp)   !integrate equations of motion
    if(mod(i,nsamp)==1.OR.nsamp==1)then
      ncalls = ncalls + 1
      call center
      call gr(iflag,j,npart,density,box,x)  
      call pressure(ptot,temp,x,box,volume,rc2)
	call msdcalc(msd)
      write(*,'(4f14.8)')etot,temp,ptot,etot
	do k=1,npart
	 write(17,'(3f14.8)')v(k,1),v(k,2),v(k,3)
	enddo
	write(18,'(2f14.8)')i*dt,msd
      pavg = pavg + ptot
    endif
    tavg = tavg + temp
    vavg = vavg + vpot
   enddo
   pavg = pavg/ncalls
   tavg = tavg/nrun
   iflag = 2
   call gr(iflag,j,npart,density,box,x)  !final call
   write(*,'(a,f8.4)')' the average temperature is ',tavg
   write(*,'(a,f8.4)')' the average potential is ',vpot/npart
   write(*,'(a,f8.4)')' the average pressure is    ',pavg
   write(*,'(a,f8.4)')' PV = ', pavg*volume
   write(*,'(a,f8.4)')' nT = ', npart*tavg
   ppp(j) = pavg
	close(17)
	close(18)
  enddo
 
  write(*,*)' here are the pressures for the blocks. '
  write(*,'(5f10.5)')ppp

  pavg = sum(ppp)/nblock
  write(*,*)' here is the uncertainty in the pressure '
  pdel = dsqrt(sum((ppp-pavg)**2/dfloat(nblock-1)))
  print*,pavg,pdel,nrun
    

! ----------------------------------------------------------------
    
  end program lj

  
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Set up an initial conditions. Pick points on a fcc cubic lattice
  subroutine init 
  use lj_mod
  implicit none
  integer                  :: i,j,k,l,ipart,nt
  real(kind=8)             :: delx,dely,delz,vx,norm,rnd,fs

  ipart  = (npart/4)**(1.0/3.0)   !four atoms per unit cell
  delx = box(1)/ipart
  dely = box(2)/ipart
  delz = box(3)/ipart
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 nt = 0
 !put in the corners
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i
     x(nt,2) = dely*j
     x(nt,3) = delz*k
    enddo
   enddo
  enddo
! put in the z faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i+0.5d0*delx
     x(nt,2) = dely*j+0.5d0*dely
     x(nt,3) = delz*k
    enddo
   enddo
  enddo
! put in the x faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i
     x(nt,2) = dely*j+0.5d0*dely
     x(nt,3) = delz*k+0.5d0*delz
    enddo
   enddo
  enddo
! put in the y faces
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     nt = nt + 1
     x(nt,1) = delx*i+0.5d0*dely
     x(nt,2) = dely*j
     x(nt,3) = delz*k+0.5d0*delz
    enddo
   enddo
  enddo
  if(nt/=npart)then
    write(*,*)' nt is not equal to npart'
    pause
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
    vx = sumx(v,l)
    v(:,l) = v(:,l)-vx
  enddo 
! calculate the scale factor for adjusting ke
  fs = dsqrt(3*npart/(sum(v*v)*beta))
  v = fs*v                           !rescale velocities
  xm = x - v*dt
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ecut = 4.0*(1/rc2**6-1/rc2**3)  !cutoff value for the LJ
  end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine force(vpot)
  use lj_mod
  implicit none
  integer                          :: i,j    !counters
  real(kind=8), dimension(3)       :: dr     !vector between atoms i and j
  real(kind=8)                     :: vpot,r2,r2i,r6i,ff
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
  subroutine integrate(etot,vpot,sumv,temp)   !integrate equations of motion
  use lj_mod
  implicit none
  real(kind=8)               :: sumv, sumv2, dt2, dit,temp, etot, vpot
  dt2 = dt*dt
  dit = 1/(2*dt)
  xx = 2*x-xm+dt2*f
  v = (xx-xm)*dit*lambda
  xx = v/dit+xm
  sumv = sum(v)
  sumv2=sum(v*v)
  xm = x
  x = xx
  temp = sumv2/(3*npart)
  etot=(vpot+0.5*sumv2)/npart   !total energy per particle
  end subroutine integrate

!-----------------------------------------------------------------------------------
  subroutine center   !keep particles in box
  use lj_mod
  implicit none
  integer                          :: i,j
  do i = 1,npart
   do j = 1,3
     if(x(i,j)>=box(j))then
       x(i,j) = x(i,j)-box(j)
       xm(i,j) = xm(i,j)-box(j)
     endif
     if(x(i,j)<=0)then
       x(i,j) = x(i,j)+box(j)
       xm(i,j) = xm(i,j)+box(j)
     endif
   enddo
  enddo
  end subroutine center

!-----------------------------------------------------------------------------------
  subroutine kebin(iflag)      !bins kinetic energy 
  use lj_mod
  implicit none
  integer                           ::  iflag,jbin,i,ntot
  real(kind=8)                      ::  rmin,rmax,rdel,ridel,anorm,an,rr,sig
  save rdel,rmin,ridel



  select case(iflag)

  case(0)
      bn = 0   !initialization
      rmin = 0.d0
      rmax = dsqrt(15.0d0/beta)
      ridel = nbins/(rmax-rmin)
      rdel = 1.d0/ridel
      iflag = 1
  case(1)
     rs = dsqrt(v(:,1)**2+v(:,2)**2+v(:,3)**2)  !calculate the speed
     do i = 1,npart
        jbin =  1+(rs(i)-rmin)*ridel
        if(jbin<=nbins) bn(jbin) = bn(jbin) + 1
     enddo
  case(2)
     bn = bn/(sum(bn)*rdel)
     rbin(1) = rmin+0.5*rdel
     do i = 2,nbins
        rbin(i) = rbin(i-1) + rdel
     enddo
!    --------------------------------------------------------------------
!    here is the Boltman distribution
     open(17, file='bolt.out')
     anorm = dsqrt(2.d0/(pi/beta**3))
     write(*,*)'   Comparisons '
     write(*,*)'   Speed  Boltzmann    MD    Difference '
     sig = 0.d0
     do i = 1,nbins
       rr = rbin(i)
       an = anorm*rr*rr*dexp(-0.5d0*rr*rr*beta)
       write( *,'(4f10.5)')rr,an,bn(i),an-bn(i)
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


  end subroutine kebin
!	--------
 subroutine msdcalc(msd)
 use lj_mod
 implicit none
 integer	:: i
 real(kind=8)		:: msd
 real(kind=8),dimension(3) :: drt
 msd = 0.d0
 drt = 0.d0
 do i=1,npart
	drt = x(i,:)-x0(i,:)
	msd = msd + dot_product(drt,drt)
 enddo
 msd = msd / npart
 end subroutine msdcalc
