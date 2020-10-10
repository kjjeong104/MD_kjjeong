!-----------------------------------------------------------------------------------
  subroutine gr(iflag,nblock,npart,density,box,x)
  implicit none
  integer                                  :: npart
  integer                                  :: nblock
  integer, parameter                       :: num_bins=200 
  integer                                  :: iflag,jbin,i,j
  integer, save                            :: ncalls
  real(kind=8), dimension(npart,3)         :: x
  real(kind=8), dimension(3)               :: dr,box
  real(kind=8)                             :: density,const,rinner,router
  real(kind=8), save,  dimension(num_bins) :: gn
  real(kind=8), save                       :: rs,rmin,rmax,rdel,ridel,rb



  
! --------------------------------------------------------------------
  select case(iflag)
! ------------------------------------------------
  case(0)
   gn = 0   !initialization
   rmin = 0.8d0
   rmax = 0.5d0*box(1)  
   ridel = num_bins/(rmax-rmin)
   rdel = 1.d0/ridel
   ncalls = 0
   iflag = 1
! ------------------------------------------------
  case(1)
   ncalls = ncalls + 1
   do i = 1,npart-1
    do j = i+1,npart
      dr = x(i,:)-x(j,:)
      dr = dr-box*nint(dr/box)  !periodic boundary conditions
      rs = dsqrt(dot_product(dr,dr))
      jbin =  1+(rs-rmin)*ridel
      if(jbin<=num_bins) gn(jbin) = gn(jbin) + 2
    enddo
   enddo
! ------------------------------------------------
  case(2)     !  carry out on the last call
    open(16,file='gr_'//char(48+nblock)//'.out')
    gn = gn/(ncalls*npart)
    rb = rmin+0.5*rdel
    const = 4.d0/3.d0*dacos(-1.d0)*density
    do i = 1,num_bins
     rb = rb + rdel
     rinner = rb-0.5d0*rdel 
     router = rinner+rdel
     gn(i) = gn(i)/(const*(router**3-rinner**3))
     write(16,'(f12.7,f13.8)')rb,gn(i)
    enddo
    close(16)
  end select
!  --------------------------------------------------------------------
  end subroutine gr
