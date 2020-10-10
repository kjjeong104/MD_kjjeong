  program fpl_lj
! ----------------------------------------------------------------SECTION A
! declare all the variables that will be used.
  implicit none
  integer, parameter                 :: npart=28,nint=27,nbins=100  !number of atoms and number internal dof
  integer                            :: ncycle,nsamp,nequil,iflag=0,istart=0
  integer                            :: i,j,ibinner             !counters
  real*8                             :: fii,req,delt,dt,etot,time_tot,vpot,temp,tavg,xmin,xmax,tequil
  real*8                             :: beta
  real*8,  dimension(npart)          :: x,xm,xx,v,f
  real*8,  dimension(nint)           :: r,fr
  real*8,  dimension(npart,nint)     :: drx
  real*8,  dimension(npart,npart)    :: umatrix          !used for normal modes
  real*8,  dimension(npart)          :: root,bk,q,pq,enm    !used for normal modes
  integer, dimension(2)              :: seed  !seed for random number generator
  real(kind=8),  dimension(nbins)    :: rbin ,bn         !used in binning the kinetic energy
! ----------------------------------------------------------------SECTION B
  open(26,file='coord.out')
! read in the date file
  open(15,file='fpl_lj.inp')
  read(15,*)time_tot  !total integration time
  read(15,*)dt        !time step
  read(15,*)delt      !sampling time
  read(15,*)fii       !force constant
  read(15,*)req       !equilibrium bond length
  read(15,*)beta      !inverse temperature for intital conditions
  read(15,*)xmin      !boundaries for the box
  read(15,*)seed      !used in a random number generator
  read(15,*)tequil    !steps used in equilibration
  read(15,*)ibinner	!option related to binning speed distribution
  close (15)
! ----------------------------------------------------------------SECTTION B1
! set up the random number generator
  call random_seed(put=seed)
! ----------------------------------------------------------------SECTION C
  req = 2.0**(1.0/6.0)   ! for L-J 
  xmax = npart*req-xmin
  fii = 4*(12*13/req**14-6*7/req**8)
  print*,' the wall boundaries are',xmin,xmax
! ----------------------------------------------------------------SECTION C
! set up derivatives of r wrt x
  drx = 0
  do i = 1,npart-1
   drx(i,i) = -1
   drx(i+1,i) = 1
  enddo
! ----------------------------------------------------------------SECTION D
  ncycle = time_tot/dt   !number of integration steps
  nsamp =  delt/dt       !sampling rate
  nequil = tequil/dt
! ----------------------------------------------------------------SECTION E
! calculate the normal modes of the linear chain
! root - are the eigenvalues
! umatrix - hessian on input, returned as transformation between normal modes and cartesian on return
  call normal_modes(npart,fii,umatrix,root,bk)
! ---------------------------------------------------------------SECTION F
  call init(npart,x,xm,v,dt,req,beta) !set up an initial condition
! ---------------------------------------------------------------SECTION F1
  do i = 1,nequil
   call force(npart,nint,x,f,r,fr,drx,fii,vpot,req,xmin,xmax) 
   call integrate(npart,dt,f,x,xx,xm,v ,etot,vpot,temp)   !integrate equations of motion
  enddo
! ---------------------------------------------------------------SECTION G  (binning has been added)
  tavg = 0  !used to calculate an average temperature
  bn = 0
  do i = 1,ncycle
   call force(npart,nint,x,f,r,fr,drx,fii,vpot,req,xmin,xmax) 
   call integrate(npart,dt,f,x,xx,xm,v ,etot,vpot,temp)   !integrate equations of motion
   tavg = tavg+temp 
   if(mod(i,nsamp)==1)call speedbin(npart,nbins,x,v,rbin,bn,temp,umatrix,root,iflag,ibinner)
 !  if(mod(i,nsamp)==1.and.i<=5000)write(26,'(9f10.5)')dt*i,(x(j),j=1,8)
    if(mod(i,nsamp*200)==1)write(6,'(9f10.5)')dt*i,etot
    if(x(1)<=-3.or.x(npart)>=xmax-2)then
       print*,' bonds have been broken'
       pause
    endif
    if(x(npart)>=xmax-2)print*,x(npart)
 !  if(mod(i,nsamp)==1)call nm_time(npart,req,x,v,q,pq,i,dt,root,bk,enm,umatrix,istart,etot)
  enddo
  temp = tavg/ncycle                !here is the average temperature

  iflag = 1
  call speedbin(npart,nbins,x,v,rbin,bn,temp,umatrix,root,iflag,ibinner)
  print*,'debug message: enm(0) is:'
  write(6,*) enm(0)
  end program fpl_lj


! ----------------------------------------------------------------
  subroutine init(npart,x,xm,v,dt,req,beta) !set up an initial condition
  implicit none
  real(kind=8), dimension(npart) :: x,xm,v
  integer                        :: i,npart
  real(kind=8)                   :: req,dt
  real(kind=8)                   :: rnd,beta,fs,vx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   This subroutine have been updated.
! here are the velocities
  do i = 1,npart
     call random_number(rnd)
     v(i) = rnd-0.5d0   !give random velocities
  enddo
! make sure average velocities are zero
  vx = sum(v)/npart
  v = v - vx
! calculate the scale factor for adjusting ke
  fs = dsqrt(npart/(sum(v*v)*beta))
  v = fs*v                           !rescale velocities
! put atoms at their equilibrium configuration
  x(1) = 0.d0
  do i = 2,npart
   x(i) = x(i-1)+req
  enddo
  xm = x - v*dt
! ------------------------------------------------------
  end subroutine init


!  ---------------------------------------------------------------------------------
!  Subroutine to calculate the forces.
!  Variables are: x cartesion coordinates; f forces wrt x; 
!                 r internals coordinates; f forces wrt r; 
!                 req equilbrium bond length; vpot potential energy
!  ------------------------------------------------------------------------------
  subroutine force(npart,nint,x,f,r,fr,drx,fii,vpot,req,xmin,xmax)
  implicit none
  real(kind=8), dimension(npart)      :: x,f
  real(kind=8), dimension(nint)       :: r,fr,ri6
  real(kind=8), dimension(npart,nint) :: drx 
  integer                             :: i,npart,nint   
  real(kind=8)         :: req,fii,vpot,x1m,xnm,xmin,xmax
   
! --------------------------------------------------------------------------------
! this is for the lennard jones potential
  do i = 1,npart-1  !calculate internals
    r(i) = x(i+1)-x(i)
  enddo
  ri6 = 1.d0/r**6
  fr = 48*(ri6*(ri6-0.5d0))/r
  vpot = 4*sum(ri6*(ri6-1.d0))
! --------------------------------------------------------------------------------
  f = matmul(drx,fr)                 !calculate cartesian forces
! --------------------------------------------------------------------------------
! these additional terms force the atoms to stay in a box with boundaries at xmin and xmax
! I have turned these terms off.  They can be included to study higher energies.
! x1m = 1.d0/(x(1)-xmin)**12
! xnm = 1.d0/(xmax-x(npart))**12
! f(1) = f(1) + 12*x1m/(x(1)-xmin)
! f(npart) = f(npart) - 12*xnm/(xmax-x(npart))
! vpot = vpot + x1m+xnm
! --------------------------------------------------------------------------------
  
  end subroutine force

!-----------------------------------------------------------------------------------
  subroutine integrate(npart,dt,f,x,xx,xm,v,etot,vpot,temp)   !integrate equations of motion
  implicit none
  real(kind=8), dimension(npart)   :: xx,x,f,xm,v
  real(kind=8)                     ::  dt,dt2, dit, etot, vpot,temp
  integer                          :: npart,mdim
    dt2 = dt*dt
    dit = 1/(2*dt)
    xx = 2*x-xm+dt2*f
    v = (xx-xm)*dit
    xm = x
    x = xx
    temp = sum(v*v)/(npart)
    etot=vpot+0.5*sum(v*v)   !total energy per particle
  end subroutine integrate

!-----------------------------------------------------------------------------------
  subroutine normal_modes(n,fii,u,root,bk)
  real(kind=8), dimension(n,n)     :: u
  real(kind=8), dimension(n)       :: bk,root
  real(kind=8)                     :: fii 
  integer                          :: n,i
! ---------------------------------------------------------
! set up the hessian
  u = 0
  do i = 1,n-1
   u(i,i) = 2.d0*fii
   u(i+1,i) = -fii
   u(i,i+1) = -fii
  enddo
  u(1,1) = fii
  u(n,n) = fii
! ---------------------------------------------------------
  call house(u,n,root,bk) !on return u is the transformation matrix to normal modes.
  write(6,*)' the normal mode force constants are shown below'
  write(6,*)root
  end subroutine normal_modes

!-----------------------------------------------------------------------------------
! subroutine converts from cartesian to normal coordinates
  subroutine nm_time(npart,req,x,v,q,pq,ic,dt,root,bk,enm,u,istart,etot)
  implicit none
  real(kind=8), dimension(npart,npart)    :: u
  real(kind=8), dimension(npart)          :: bk,root,x,v,q,pq,enm
  real(kind=8)                            :: dt,req,etot
  integer                                 :: i,ic,istart,npart
  if(istart==0)then
   open(55,file='nm_time.dat')
   print*,' root(3) = ',root(3)
   istart=1
  endif
  do i = 1,npart
    bk(i) = x(i) - req*(i-1)
  enddo
  q = matmul(transpose(u),bk)
  pq = matmul(transpose(u),v)
  enm = 0.5d0*(pq*pq+root*q*q)
  write(55,'(7f10.5)')q(3),pq(3),enm(3),etot  !write results for mode 3
  end subroutine nm_time


!---------------------------------------------------------------------------
  subroutine speedbin(npart,nbins,x,v,rbin,bn,temp,umatrix,root,iflag,ibinner)
  implicit none
  integer                          :: iflag,jbin,i,j,nbins,npart,ibinner
  real(kind=8),  dimension(nbins)   :: rbin ,bn
  real(kind=8),  dimension(npart)  :: x,v,pq ,q,enm,root                        !velocity
  real*8,  dimension(npart,npart)  :: umatrix          !used for normal modes
  real(kind=8)                     :: rmin,rmax,btot,rdel,ridel,anorm,anorm1,an,rr,temp
  
! -------------------------------------------------------------------------------
! One can bin various quantites.
! if ibinner = 1 !bin the energy in Cartestian speeds.
! if ibinner = 2 !bin the energy in NM speeds. 
! if ibinner > 2 bin in the energy in normal mode ibinner-2 
!  ibinner = 5
  pq = matmul(transpose(umatrix),v)
  q = matmul(transpose(umatrix),x)
  enm = 0.5d0*(pq*pq+root*q*q)
! -------------------------------------------------------------------------------
  rmin = 0.d0  !these determine the range of the binned variables
  rmax = (5*temp) 
  if(ibinner<=2)rmax = sqrt(rmax)
  rdel = (rmax-rmin)/nbins
  ridel = 1.d0/rdel
! -------------------------------------------------------------------------------
  if(iflag==0)then
     if(ibinner==1)then
       do i = 1,npart   
        jbin =  1+(dabs(v(i))-rmin)*ridel
        if(jbin<=nbins) bn(jbin) = bn(jbin) + 1
       enddo
     elseif(ibinner==2)then
       do i = 2,npart
        jbin =  1+(dabs(pq(i))-rmin)*ridel
        if(jbin<=nbins) bn(jbin) = bn(jbin) + 1
       enddo
     else
      jbin =  1+(dabs(enm(ibinner-2))-rmin)*ridel
      if(jbin<=nbins) bn(jbin) = bn(jbin) + 1
    endif
  else !only for last call
    btot = sum(bn) !normalize the array
    bn = bn/(btot*rdel)
    rbin(1) = rmin+0.5*rdel
    do i = 2,nbins
      rbin(i) = rbin(i-1) + rdel
    enddo
    open(16, file='bins.out')
    write(16,'(f10.5,f14.5)')(rbin(i),bn(i),i=1,nbins)
    close(16)
!   --------------------------------------------------------------------
!   here is the Boltman distribution
    open(17, file='bolt.out')
    if(ibinner<=2)then
      anorm = dsqrt(2.d0/(temp*dacos(-1.d0))) 
      do i = 1,nbins
        an = anorm*dexp(-.5d0*rbin(i)**2/temp)  
        write(17,'(2f16.5)')rbin(i),an
      enddo
    else
      anorm = 1.d0/temp 
      do i = 1,nbins
        an = anorm*dexp(-rbin(i)/temp)  
        write(17,'(2f16.5)')rbin(i),an
      enddo
   endif
   close(17)
!  --------------------------------------------------------------------
  endif

  end subroutine speedbin

