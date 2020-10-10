  program fpl
! ----------------------------------------------------------------SECTION A
! declare all the variables that will be used.
  implicit none
  integer, parameter                 :: npart=6,nint=5  !number of atoms and number internal dof
  integer                            :: ncycle,nsamp,istart=0
  integer                            :: i,j             !counter
  real*8                             :: fii,req,delt,dt,etot,time_tot,vpot
  real*8,  dimension(npart)          :: x,xm,xx,v,f
  real*8,  dimension(nint)           :: r,fr
  real*8,  dimension(npart,nint)     :: drx
  real*8,  dimension(npart,npart)    :: umatrix          !used for normal modes
  real*8,  dimension(npart)          :: root,bk,q,pq,enm    !used for normal modes

! ----------------------------------------------------------------SECTION A'
! this file is new.  You can you the output file init.xyz to make a jmol movie.
  open(25,file='init.xyz')   !jmol movie file
! ----------------------------------------------------------------SECTION B

! read in the date file
  open(15,file='fpl_simp.inp')
  read(15,*)time_tot  !total integration time
  read(15,*)dt        !time step
  read(15,*)delt      !sampling time
  read(15,*)fii       !force constant
  read(15,*)req       !equilibrium bond length
  close (15)
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
! ----------------------------------------------------------------SECTION E
! calculate the normal modes of the linear chain
! root - are the eigenvalues
! umatrix - hessian on input, returned as transformation between normal modes and cartesian on return
  call normal_modes(npart,fii,umatrix,root,bk)
! ---------------------------------------------------------------SECTION F
  call init(npart,x,xm,v,dt,req) !set up an initial condition
! ---------------------------------------------------------------SECTION G
  do i = 1,ncycle
   call force(npart,nint,x,f,r,fr,drx,fii,vpot,req) 
   call integrate(npart,dt,f,x,xx,xm,v ,etot,vpot)   !integrate equations of motion
   if(mod(i,10)==1)call nm_time(npart,req,x,v,q,pq,i,dt,root,bk,enm,umatrix,istart)
   if(mod(i,ncycle/40)==1)call jmol(npart,x)    !make a movie
  enddo
    
  end program fpl


  subroutine jmol(npart,x)
  implicit none
  integer                          :: npart,i
  character(len=4)                 :: atom=' C  '
  real(kind=8), dimension(npart)   :: x
  write(25,*)npart
  write(25,*)
  do i = 1,npart
    write(25,'(a4,3f10.3)')atom,x(i),0.d0,0.d0
  enddo
  end subroutine jmol

! ----------------------------------------------------------------
  subroutine init(npart,x,xm,v,dt,req) !set up an initial condition
  implicit none
  real(kind=8), dimension(npart) :: x,xm,v
  integer                        :: i,npart
  real(kind=8)                   :: req,dt
  v = 0
  x(1) = 0.d0
  do i = 2,npart
   x(i) = x(i-1)+req
  enddo
  x(1)=x(1)-1.d0                   !stretch one bond
  xm = x - v*dt
! ------------------------------------------------------
  end subroutine init


!  ---------------------------------------------------------------------------------
!  Subroutine to calculate the forces.
!  Variables are: x cartesion coordinates; f forces wrt x; 
!                 r internals coordinates; f forces wrt r; 
!                 req equilbrium bond length; vpot potential energy
!  ------------------------------------------------------------------------------
  subroutine force(npart,nint,x,f,r,fr,drx,fii,vpot,req)
  implicit none
  real(kind=8), dimension(npart)      :: x,f
  real(kind=8), dimension(nint)       :: r,fr
  real(kind=8), dimension(npart,nint) :: drx 
  integer                             :: i,npart,nint   
  real(kind=8)         :: req,fii,vpot
   
  do i = 1,npart-1  !calculate internals
    r(i) = x(i+1)-x(i)-req
  enddo

  fr = -fii*r                        !calculate internal forces
  vpot = 0.5d0*fii*dot_product(r,r)  !calculate potential
  f = matmul(drx,fr)                 !calculate cartesian forces
  
  end subroutine force

!-----------------------------------------------------------------------------------
  subroutine integrate(npart,dt,f,x,xx,xm,v,etot,vpot)   !integrate equations of motion
  implicit none
  real(kind=8), dimension(npart)   :: xx,x,f,xm,v
  real(kind=8)                     ::  dt,dt2, dit, etot, vpot
  integer                          :: npart,mdim
    dt2 = dt*dt
    dit = 1/(2*dt)
    xx = 2*x-xm+dt2*f
    v = (xx-xm)*dit
    xm = x
    x = xx
    etot=vpot+0.5*sum(v*v)   !total energy
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
! print out the hessian
  write(*,*)' if the number of atoms is less that 15 the Hessian is now printed.'
  if(n<=15)then
    do i = 1,n
      write(*,'(15f7.3)')u(i,:)
    enddo
  endif
! ---------------------------------------------------------
! find the normal modes by diagonalizing the Hessian
  call house(u,n,root,bk) !on return u is the transformation matrix to normal modes.

! print out the roots
  write(*,*)' here are the eigenvalues of the hessian.'
  write(*,*)root

  write(*,*)' If the number of atoms is less that 15 the normal modes are now printed.'
  write(*,*)' Each column corresponds to a normal mode vector.'
  if(n<=15)then
    do i = 1,n
      write(*,'(15f8.4)')u(i,:)
    enddo
  endif

  end subroutine normal_modes

!-----------------------------------------------------------------------------------
! subroutine converts from cartesian to normal coordinates
  subroutine nm_time(npart,req,x,v,q,pq,ic,dt,root,bk,enm,u,istart)
  implicit none
  real(kind=8), dimension(npart,npart)    :: u
  real(kind=8), dimension(npart)          :: bk,root,x,v,q,pq,enm
  real(kind=8)                            :: dt,req
  integer                                 :: i,ic,istart,npart
  if(istart==0)then
   open(55,file='nm_time.dat')
   write(55,*)' time       q(3)      pq(3)    energy(3) '
   print*,' root(3) = ',root(3)
   istart=1
  endif
  do i = 1,npart
    bk(i) = x(i) - req*(i-1)
  enddo
  q = matmul(transpose(u),bk)
  pq = matmul(transpose(u),v)
  enm = 0.5d0*(pq*pq+root*q*q)
  write(55,'(7f10.5)')ic*dt,q(3),pq(3),enm(3)  !write results for mode 3
  end subroutine nm_time


