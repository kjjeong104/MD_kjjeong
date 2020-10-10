  program integer_test
  implicit none
  integer                  :: i,j

  j = 65
  i = j**(1.0/3.0)+1

  write(*,*) 'i =',i

  end program integer_test
