  program sho
  implicit none
  integer :: i                !counter
  integer :: ncycle           !number time steps
  real*8  :: x,xm,xx,v,f      !positions, velocities v, and forces f
  real*8  :: time_tot = 10.0  !total run time
  real*8  :: etot,vpot        !total and potential energies
  real*8  :: dt       =  0.02 !step size
  real*8  :: k = 1.0          !force constant (mass m=1)

! ----------------------------------------------------------------
  ncycle = time_tot/dt  !number of time steps
! ----------------------------------------------------------------
! set initial conditions
  x = 1.0
  v = 0.0
  xm = 1.0-v*dt

  do i = 1,ncycle
   call force(x,f,k,vpot)                      !calculate forces and potential
   call integrate(dt,f,x,xx,xm,v ,etot,vpot)   !integrate equations of motion
   if(mod(i,10)==1)write(6,'(4f12.7)')i*dt,xm,v,etot
  enddo
    
  end program sho

  

 
  
!  ---------------------------------------------------------------------------------
!  Subroutine to calculate the forces.
!  Variables are: x cartesion coordinates; f forces wrt x; 
!                 k is the force constant; vpot is the potential energy
!  ------------------------------------------------------------------------------
  subroutine force(x,f,k,vpot)
  implicit none
  real(kind=8) :: x,f,k,vpot  
    f = -k*x
    vpot = 0.5*k*x*x
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

