  MODULE lj_mod
      implicit none
      integer                                   :: npart  !number of particles 
      real(kind=8), dimension(3)                :: box       !size of the box
      real(kind=8)                              :: beta,dt   !beta=1/kT and dt is integration step size
      integer                                   :: nequil,itest,nrun,nsamp
      real(kind=8)                              :: tequil,ttest,trun,tsample
      real(kind=8)                              :: density  !the density of the gas
      real(kind=8), allocatable, dimension(:,:) :: x,xx,v,xm,f
      integer, dimension(2)                     :: seed       !seed used in random number generator
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     these variables are used in the potential cutoffs
      real(kind=8)                              :: rcut,rc2,ecut  !cut off squared
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     used in kebin
      integer                                   :: nbins  !number of particles 
      real(kind=8), allocatable, dimension(:)   :: rbin ,bn
      real(kind=8), allocatable, dimension(:)   :: rs  !this is the ojbect to be binned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=8), parameter ::clight=2.99792458D10,av=6.0221367D23,hb=1.05457266D-27 
      real(kind=8), parameter ::mass_h=1.007825d0
      real(kind=8), parameter :: pi=3.141592654
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine input ! read in the date file
      implicit none
      open(15,file='lj.inp')
      read(15,*)npart
      read(15,*)tequil
      read(15,*)ttest
      read(15,*)trun
      read(15,*)tsample
      read(15,*)dt
      read(15,*)beta
      read(15,*)box
      read(15,*)nbins
      read(15,*)seed

      nrun = trun/dt
      itest = ttest/dt
      nequil = tequil/dt
      nsamp = tsample/dt
      if(nsamp==0)nsamp = 1

      density = npart/(box(1)*box(2)*box(3))

      rcut = minval(box)/2.d0
      write(*,*)' the number of particles is ',npart
      write(*,*)' the density is             ',npart/(box(1)*box(2)*box(3))
      write(*,*)' the value of rcut is       ',rcut
      write(*,*)' the initial guess for T is ',1/beta
      write(*,*)' the step size is           ',dt
      allocate(x(npart,3),v(npart,3),xm(npart,3),xx(npart,3),f(npart,3))
      allocate(rbin(nbins),bn(nbins),rs(npart))
      rc2 = rcut**2
      close (15)
      end subroutine input

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=8) function sumx(v,icart)
      implicit none
      integer                      :: icart
      real(kind=8), dimension(:,:) :: v
       sumx = sum(v(:,icart))/ubound(v,1)
      end function sumx

    END MODULE lj_mod	
