  program quadmc
  implicit none
  real(kind=8)  :: x, xn, eno, enn, delx, vpot, beta, xmin, xmax, rnd
  integer	:: j, nstep, istep, iflag
  integer, dimension(2) :: seed
! ---Read info
 open(15,file='quadmc.inp')
 read(15,*)beta !scaled temperature kT
 read(15,*)nstep !# of steps
 read(15,*)delx !spatial size of one maximum step
 read(15,*)xmin !minimum value of x for bin
 read(15,*)xmax !maximum value of x for bin
 read(15,*)seed !random seed
 call random_seed(put=seed)
! ---Initialization
  iflag=0
  istep=0
  x = 0
  call potential(x,eno)
! ----
 do j=1,nstep
  call random_number(rnd) !displacement variable
  xn = x + (rnd-0.5)*delx
  call potential(xn,enn)
  call random_number(rnd) !acceptance variable
  if(rnd<exp(-beta*(enn-eno)))then
	eno = enn
	istep = istep + 1
	x = xn
  endif
  call xbin(x,iflag,xmin,xmax,nstep,beta)
 enddo
  iflag = 2
  call xbin(x,iflag,xmin,xmax,nstep,beta)
	write(*,'(a,f12.7)')' The percentage of accepted moves is = ',dfloat(istep)/dfloat(nstep)

  end program quadmc
! --------------------
 subroutine potential(x,vpot)
 implicit none
 real*8	:: x, vpot
 vpot = -(x**2/2.0d0) + x**4/4.0d0	!HO potential. One more formula in xbin subroutine
 end subroutine potential
! --------------------
 subroutine xbin(x,iflag,xmin,xmax,nstep,beta)
 implicit none
 real*8		:: x, xmin, xmax, vpot, beta, sig
 integer	:: iflag, jbin, j, nstep
 integer, parameter	:: num_bins=1000
 real*8, save, dimension(num_bins) :: ppos, boltp
 real*8, save		:: xs,xdel,xidel,xb

 select case(iflag)
! --
  case(0)
  ppos = 0 	!initialization
  boltp = 0
  xidel = num_bins/(xmax-xmin)
  xdel = 1.d0/xidel
  iflag = 1
! --
  case(1)
   jbin = 1+(x-xmin)*xidel
   if(jbin<=num_bins) ppos(jbin) = ppos(jbin) + 1
! --
  case(2)	! carry out on the last call
   open(16,file='histo.out') ! will bin probability density
   ppos = ppos/(nstep*xdel)
   xb = xmin+0.5*xdel
   do j=1,num_bins
    vpot = -(xb**2/2.0d0) + xb**4/4.0d0
    boltp(j) = exp(-beta*vpot)
    write(16,'(2f14.8)')xb,ppos(j)
    xb = xb + xdel
   enddo
   
   sig=0.d0
   boltp = boltp /(sum(boltp)*xdel)
   open(17,file='bolt.out') !Boltzmann distribution
   xb = xmin+0.5*xdel
   do j=1,num_bins
    write(17,'(2f14.8)')xb,boltp(j)
    sig = sig + (ppos(j)-boltp(j))**2
    xb = xb + xdel
   enddo
   sig = dsqrt(sig/(num_bins-1))
   write(*,*)'difference between Boltzmann and MC is sig =',sig
   close(16)
   close(17)
 end select
! --

 end subroutine xbin
