  MODULE lj_metro_npt_mod
      implicit none
      real*8, parameter ::clight=2.99792458D10,av=6.0221367D23,hb=1.05457266D-27 
      real*8, parameter ::mass_h=1.007825d0
      real*8, parameter :: pi=3.141592654
      integer, parameter:: nbins=50
      real*8            :: box
      real*8            :: density,tstar,pres 
      real(kind=8)          :: beta,delx,vmax
      real(kind=8)          :: rcut,rc2,ecut,volume
      integer, dimension(2) :: seed
      integer               :: nsamp,ncalls
      integer               :: nblocks                    ! number of blocks
      integer               :: npart                      ! number of particles
      real(kind=8), allocatable, dimension(:,:) :: x,xm   !coordinate
!   ------------------------------------------------------------------------------------------
          
    contains


    END MODULE lj_metro_npt_mod




