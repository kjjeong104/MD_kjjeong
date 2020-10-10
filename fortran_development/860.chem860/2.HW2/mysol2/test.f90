  program test
  implicit none
  integer                   :: i,j
  integer, dimension(2,6)   :: iflag
  integer, dimension(2)     :: isum
  integer, dimension(6)     :: jsum
  do i = 1,2
   write(*,*)i
   do j = 1,6
     iflag(i,j) =10*i+j
   enddo
  enddo
  write(*,'(6I4)')transpose(iflag)
  isum = sum(iflag,2)
  jsum = sum(iflag,1)
  write(*,*)isum
  write(*,*)jsum
  write(*,*)sum(iflag(:,1))
  write(*,*)sum(iflag(:,6))
  end program test
