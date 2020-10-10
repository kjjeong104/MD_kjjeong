  subroutine kebin(v,ke,rbin,bn,beta,nn,iflag)
  implicit none
  integer                      :: iflag,jbin,i,j,ntot,nbins,npart
  real(kind=8),  dimension(:)  :: rbin ,bn
  real(kind=8),  dimension(:)  :: v                           !velocity
  real(kind=8),  dimension(:)  :: ke  !kinetic energy
  real(kind=8)                 :: rmin,rmax,rdel,ridel,anorm,an,rr,beta
  integer, dimension(:)        :: nn  !distribution of kinetic energies
  
  nbins = ubound(rbin,1)
  npart = ubound(ke,1)

  if(iflag==0)nn = 0   !initialization of numbers in bins

  ke = 0.5d0*v*v       ! calculate the kinetic energies

! --------------------------------------------------------------------
! set the bounds on the distribution
  rmin = 0.d0
  rmax = 2.5d0/beta
! --------------------------------------------------------------------
! set the size of the bins
  ridel = nbins/(rmax-rmin)
  rdel = 1.d0/ridel
! --------------------------------------------------------------------

  do i = 1,npart
   jbin =  1+(ke(i)-rmin)*ridel
   if(jbin<=nbins) nn(jbin) = nn(jbin) + 1
  enddo
  
! --------------------------------------------------------------------
! carry out on the last call
  if(iflag==2)then
   ntot = sum(nn) !normalize the array
   write(*,*)' ntot = ',ntot
   bn = nn/(ntot*rdel)
   rbin(1) = rmin+0.5*rdel
   do i = 2,nbins
    rbin(i) = rbin(i-1) + rdel
   enddo
!  --------------------------------------------------------------------
!  here is the Boltman distribution (perhaps).  This might need to be 
!  fixed
   open(17, file='bolt.out')
   anorm = 4.d0*beta
   do i = 1,nbins
    rr = rbin(i)
    an = anorm*dexp(-rr*4.d0*beta)  
    write(17,'(2f10.5)')rr,an
    write(*,'(2f10.5)')rr,an
   enddo

   pause
   close(17)
!  --------------------------------------------------------------------
   print*,' here are the final bins'
   open(16, file='bins.out')
   write(16,'(f10.5,f10.5)')(rbin(i),bn(i),i=1,nbins)
   write( *,'(f10.5,f10.5)')(rbin(i),bn(i),i=1,nbins)
   close(16)
  endif

  end subroutine kebin

