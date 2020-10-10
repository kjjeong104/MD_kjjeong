    module pressure_mod
      implicit none

      contains

      subroutine pressure(p,temp,x,box,volume,rc2)
      implicit none
      real(kind=8), dimension(:,:) :: x
      real(kind=8), dimension(3)   :: dr,box
      real(kind=8)                 :: r2,r6i,rc2,temp,volume,p
      integer                      :: i,j    !counters
      integer                      :: npart  !number of particles

      p = 0
      npart = ubound(x,1)
      do i = 1,npart-1
       do j = i+1,npart
         dr = x(i,:)-x(j,:)
         dr = dr-box*nint(dr/box)  !periodic boundary conditions
         r2 = dot_product(dr,dr)
         if(r2<=rc2)then
           r6i=(1.d0/r2)**3
           p = p + r6i*(r6i-0.5)  !LJ potential
         endif
       enddo
      enddo
      p = (npart*temp+16*p)/volume
     end subroutine pressure

 
    end  module pressure_mod

