 Module correlate_mod
 implicit none

  contains

!   ------------------------------------------------------------------------------------------
    subroutine corfun(a,acf,norm)
     implicit none
     real(kind=8), dimension(:,:)   :: a
     real(kind=8), dimension(:)     :: acf,norm
     real(kind=8), dimension(3)     :: a0
     integer                        :: i,j,k,ttomax,irun,icor

       irun = ubound(a,1)
       icor = ubound(acf,1)-1
       do i = 1,irun
        a0 = a(i,:)
        ttomax = min(irun,i+icor)
        do j = i,ttomax
          k = j-i+1
          acf(k) = acf(k) + dot_product(a0,a(j,:))
          norm(k) = norm(k) + 1
        enddo
      enddo
   end subroutine corfun
!   ------------------------------------------------------------------------------------------
         
  END MODULE correlate_mod 
