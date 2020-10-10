  program rndarray
  implicit none
  integer	:: i !counter
  integer, dimension(2)	::	seed
  real(kind=8)	::	rnd	!a random number
 open(16,file='rndarray.dat')
 seed = (/ 53705, 1643 /)
 call random_seed(put=seed)
do i=1,50000
	call random_number(rnd)
	write(16,*)rnd
enddo
 close(16)
  end program rndarray
