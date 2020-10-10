  program sho_simp
  implicit double precision (a-h,o-z)
! ----------------------------------------------------------------
  time_tot = 10    !total time of integration
  dt       =  0.02 !step size
  k = 1.0          !force constant (mass m=1)
  ncycle = time_tot/dt  !number of time steps
! ----------------------------------------------------------------
! set initial conditions
  x = 1.0
  v = 0.0
  xm = 1.0-v*dt
! ----------------------------------------------------------------
  do i = 1,ncycle
!   -------------------------------------------------------------
!   calculate the force
    f = -k*x    !determine the force
    vpot = 0.5*k*x*x !calculate the potential
!   -------------------------------------------------------------
!   integrate 
    xx = 2*x-xm+dt*dt*f  !position and new time
    v = (xx-xm)/(2*dt)   !velocity at current time
    xm = x
    x = xx
    etot=vpot+0.5*v*v    !calculate the total energy
!   -------------------------------------------------------------
   if(mod(i,10)==1)write(6,'(4f12.7)')i*dt,xm,v,etot
  enddo
  end program sho_simp

  

 
