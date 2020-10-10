! subroutine calculates corrlation functions directly and with FT.
! the diffusion constant is also calculated.
  program correlate
  use correlate_mod
  use fourier_mod
  implicit none
  integer :: i,j ,k ,iatom                         !counter
  integer :: isign
!  --------------------------------------------------------------------------
  real(kind=8), allocatable, dimension(:,:)   :: a
  real(kind=8), allocatable, dimension(:)     :: acf,norm
  real(kind=8), allocatable, dimension(:)     :: fft,ft,split,hold
  real(kind=8), allocatable, dimension(:,:,:) :: v
  real(kind=8)                                :: dt,tcor,t1,t2,integral,jntegral
  integer                                     :: npoints,nfour,nfour2,npart,icor
!  --------------------------------------------------------------------------
  open(15,file='velocity.dat')         !here you read in velocities from lj program.
  write(*,*)' enter the correlation time.  1.5 is a good guess.'
  read(*,*)tcor
  read(15,*)npart,npoints,dt     !This dt can be some multiple of the dt used in traj program
  icor = tcor/dt
! -----------------------------------------------------Section A
  allocate(v(npoints,npart,3))
  allocate(a(npoints,3),acf(icor+1),norm(icor+1))
! --------------------------------------------------------------
  do i = 1,npoints     !read in the velocities
   read(15,*)v(i,:,:)
  enddo
! -----------------------------------------------------Section B
! average over the atoms  !calcluate acf
   acf = 0
   norm = 0
   call cpu_time(t1)
   do i = 1,npart
     a = v(:,i,:)
    call corfun(a,acf,norm)
   enddo
   acf = acf/norm
   open(26,file='acf.out')
   integral = -acf(1)*0.5d0
   do i = 1,icor+1
     integral = acf(i)+integral
     write(26,*)(i-1)*dt,acf(i)
   enddo
   integral = integral*dt/3.d0
   write(*,*)' the diffusion constant is',integral
   close(26)
    write(*,*)'the diffusion constant is',(0.5d0*acf(1)+sum(acf(2:icor+1)))*dt/3.d0
 ! deallocate(a,acf)
   call cpu_time(t2)
   write(*,*)' cpu time for acf is ',t2-t1
! -----------------------------------------------------Section C
! determine the number of fourier points
  nfour = 2
  do i = 1,24
   if(nfour<npoints)nfour = nfour*2
  enddo
  nfour = nfour*2
  nfour2 = 2*nfour
  write(*,'(a,i7)')' number of fourier points is ',nfour
  allocate (fft(nfour),ft(nfour),split(nfour2),hold(nfour2))
! -----------------------------------------------------Section D
  call cpu_time(t1)
  fft =  0
  do iatom = 1,npart
   do j = 1,3                !cartesian coordinates 
    do i = 1,npoints
      ft(i) = v(i,iatom,j)
    enddo
    ft(npoints+1:nfour) = 0.0   !padding an array with zeros
    isign = 1
    call  fftet(ft,split,hold,isign)
! -----------------------------------------------------Section E
!   put real and imaginary parts of function back together
    k=1
    do i=1,nfour
         ft(i) = (hold(k)**2+hold(k+1)**2)
	 k=k+2
    enddo
    isign = -1
    call  fftet(ft,split,hold,isign)
    hold = dt*hold/nfour              
    k = 1
    isign = 1
    do i=1,nfour
        ft(i) = isign*hold(k)
        isign = -isign
        k=k+2
    enddo
    fft = ft + fft
   enddo
  enddo
  fft = fft/(npoints*dt*npart)
  call cpu_time(t2)
  write(*,*)' cpu time for acf calculated with FT is ',t2-t1

  open(26,file='corr.out')

! integral = -0.5d0*(fft(nfour/2+1)-acf(1))*dt
! jntegral = -0.5d0*fft(nfour/2+1)*dt

  do i = nfour/2,nfour-1
    integral = integral + (fft(i+1)-acf(i-nfour/2+1))*dt
!   jntegral = jntegral + fft(i+1)*dt
    write(26,'(4f12.4)')dt*(i-nfour/2),fft(i+1),acf(i-nfour/2+1),integral
  enddo
  write(*,*)'the diffusion constant is',sum(fft(nfour/2-icor:nfour/2+icor))*dt/6.d0
  write(*,*)'the diffusion constant is',sum(fft)*dt/6.d0
  close(26)
! ----------------------------------------------------------------

    
  end program correlate




