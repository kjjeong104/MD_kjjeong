  program ljo
  implicit none
  integer :: i                !counter
  integer :: ncycle           !number time steps
  real*8  :: x,xm,xx,v,f      !positions, velocities v, and forces f
  real*8  :: time_tot = 62.82  !total run time : 10*T = 20*pi
  real*8  :: etot,vpot        !total and potential energies
  real*8  :: dt        !step size
  real*8  :: k = 1.0          !force constant (mass m=1)
! ----------------------------------------------------------------
  print *, 'Timestep = ?'
  read *, dt
  ncycle = time_tot/dt  !number of time steps
! ----------------------------------------------------------------
! set initial conditions
  x = 1.01		! let x slightly greater than 1, to avoid escape
  v = 0.0
  xm = 1.01-v*dt
  !xa = 1.0
  !va = 0.0
        open(15,File='STDOUT_LJ')

  do i = 1,ncycle
   call force(x,f,k,vpot)                      !calculate forces and potential
   call integrate(dt,f,x,xx,xm,v ,etot,vpot)   !integrate equations of motion
   !if(mod(i,10)==1)write(6,'(4f12.7)')i*dt,xm,v,etot
	write(15,'(4f12.7)')i*dt,xm,v,etot
  enddo
  end program ljo

!  ---------------------------------------------------------------------------------
!  Subroutine to calculate the forces.
!  Variables are: x cartesion coordinates; f forces wrt x; 
!                 k is the force constant; vpot is the potential energy
!  ------------------------------------------------------------------------------
  subroutine force(x,f,k,vpot)
  implicit none
  real(kind=8) :: x,f,k,vpot
  real*8	::	x6i
    x6i = 1/(x**6)
    f = 48/x * (x6i * (x6i - 0.5))
    vpot = 4 * (x6i * (x6i - 1.0))
	!write(6,'(4f12.7)') x,x6i,f,vpot
  end subroutine force

!-----------------------------------------------------------------------------------
  subroutine integrate(dt,f,x,xx,xm,v,etot,vpot)   !integrate equations of motion
  implicit none
  real(kind=8)      :: xx,x,f,xm,v
  real(kind=8)      :: dt,etot, vpot
    xx = 2*x-xm+dt*dt*f  !position and new time
    v = (xx-xm)/(2*dt)   !velocity at current time
    xm = x
    x = xx
    etot=vpot+0.5*v*v
  end subroutine integrate

