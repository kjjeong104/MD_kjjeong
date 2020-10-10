  MODULE lj_met_chain_mod
      implicit none
      real*8, parameter ::clight=2.99792458D10,av=6.0221367D23,hb=1.05457266D-27 
      real*8, parameter ::mass_h=1.007825d0
      real*8, parameter :: pi=3.141592654, kbolt=8.31447d0
      integer, parameter:: nbins=50
      real*8            :: box, bl, kvib
      real*8            :: density,tstar,pres 
      real(kind=8)          :: beta,delx,vmax
      real(kind=8)          :: rcut,rc2,ecut,volume
      integer		    :: iread, const_vol
      integer, dimension(2) :: seed
      integer               :: nsamp,ncalls
      integer               :: nblocks                    ! number of blocks
      integer               :: npart,nmol,natom   ! number of particles,chains,atoms in a chain
      real(kind=8)		:: pot_lj,pot_str,e_lrc
      real(kind=8), allocatable, dimension(:,:) :: x,xm,old_mol   !coordinate
!   ------------------------------------------------------------------------------------------
          
    contains


!   -----------------------
  subroutine input(irun,nequil)
    implicit none
    integer      :: irun,nequil
    open(15,file='lj_met_chain.inp')
!   ----------------------------------------------------------------
!   read in the data file
    read(15,*)iread	!if iread=1 read in old_configuration
    read(15,*)const_vol !set to zero for NPT otherwise NVT
    read(15,*)nmol     ! number of chain molecules
    read(15,*)natom	! number of 'beads' in one chain
    read(15,*)bl	!bond length
    read(15,*)kvib	!bond harmonic force const
    npart = natom*nmol
    read(15,*)irun      ! number of points for the simulation
    read(15,*)nequil    ! number of points in equilibration
    read(15,*)nsamp     ! number of points between calls to boltmann or gr 
    read(15,*)delx      ! mc config step size
    read(15,*)vmax      ! mc volume step size
    read(15,*)tstar     ! the desired temperature
    read(15,*)pres      ! the value of the pressure 
    read(15,*)density   ! the initial density 
    read(15,*)seed      ! used for initializing random number generator
    read(15,*)nblocks   ! number of blocks used for averaging the density
    close (15)
    volume  = npart/density
    box = volume**(1.d0/3.d0)
    beta = 1.d0/tstar
    rcut = 2.5d0             !the cut-off on the LJ forces is rcut
    ncalls = irun/nsamp         !number of calls to gr
    allocate(x(npart,3),xm(npart,3),old_mol(natom,3))
    rc2 = rcut*rcut
    write(6,'(a,f8.5)')' the box length is     ',box
    write(6,'(a,f8.5)')' r cutoff is           ',rcut
	e_lrc = 8.0D0*pi*npart*density*(((1/rcut)**9)/3.0D0-(1/rcut)**3)/3.0D0 
    write(6,*)' e_lrc is              ',e_lrc
! ----------------------------------------------------------------
  end subroutine input


! ----------------------------------------------------------------
  subroutine init !set up an initial condition using cubic arrangement
  implicit none
  integer :: i,j,k,l,ipart,nt
  real*8  :: rnd,delta !bl: bond length btwn beads
! pick points on a cubic grid
   
  ipart  = (1.00001*dfloat(nmol))**(1.0/3.0)
  delta = box/dfloat(ipart)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!reading old config if needed
 if(iread==1)then
   open(18,file='old_config.ini')
   read(18,*) box
   write(16,'(a,f10.4)')' the current box length is ',box
   volume = box**3
   read(18,*)x
   close(18)
 else
!!!!!!!!!!!!!!!!!!
 nt = 0
  if(natom/=8)then
    write(16,*)' For this version, please try only 8 beads per chain. stopping'
    stop
  endif
 !put on a grid. in each grid, grow chain by small cube
  do i = 0,ipart-1
   do j = 0,ipart-1
    do k = 0,ipart-1
     do l = 1,natom
      nt = nt + 1
	select case(l)
	case(1)
    	 x(nt,1) = delta*i
     	 x(nt,2) = delta*j
     	 x(nt,3) = delta*k
        case(2)                                                                                       
         x(nt,1) = delta*i + bl
         x(nt,2) = delta*j
         x(nt,3) = delta*k
        case(3)                                                                                       
         x(nt,1) = delta*i + bl
         x(nt,2) = delta*j + bl
         x(nt,3) = delta*k
        case(4)                                                                                       
         x(nt,1) = delta*i
         x(nt,2) = delta*j + bl
         x(nt,3) = delta*k
        case(5)                                                                                       
         x(nt,1) = delta*i
         x(nt,2) = delta*j + bl
         x(nt,3) = delta*k + bl
        case(6)                                                                                       
         x(nt,1) = delta*i + bl
         x(nt,2) = delta*j + bl
         x(nt,3) = delta*k + bl
        case(7)                                                                                       
         x(nt,1) = delta*i + bl
         x(nt,2) = delta*j
         x(nt,3) = delta*k + bl
        case(8)
         x(nt,1) = delta*i
         x(nt,2) = delta*j
         x(nt,3) = delta*k + bl
	end select
     enddo
    enddo
   enddo
  enddo
  if(nt/=npart)then
     write(16,*)' stopping since nt is not equal to npart '
     write(16,*)' nt    = ',nt
     write(16,*)' npart = ',npart
     stop
  endif
 endif
! ------------------------------------------------------
  ecut = 4*(1/rc2**3-1)/(rc2**3) !cutoff value for the 6-12 potential
  write(16,*) 'ecut = ',ecut
  end subroutine init


!-----------------------------------------------------------------------------------
  subroutine potential   !calculate the potentials : consider LJ, harmonic
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r2i,r6i,addstr,addlj
  real(kind=8)               :: xr,r2
  integer                    :: i,j    !counters
!harmonic stretch pot:adjacent atoms in a chain.
!LJ : all nonbonded pairs, except atom pairs of harmonic stretch
!algorithm : molecule index = (atomindex+1)/natom + 1.
!ex) atom#1~8 : mol#1, 9~16:mol#2,...
!if diff of atomindex >=2, calc LJ. diff of atomindex=1-> if smaller one is multiple of 8, calc LJ
!otherwise: calc harmonic potential
  pot_lj = 0.d0
  pot_str = 0.d0
  do i = 1,npart-1
   do j = i+1,npart
       dr = x(i,:)-x(j,:)
       dr = dr-box*nint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
	if(j-i==1 .AND. mod(i,natom)/=0)then !harmonic
		addstr = 0.5d0*kvib*(dsqrt(r2)-bl)**2
		pot_str = pot_str + addstr
!		pot_str = pot_str + 0.5d0*kvib*(dsqrt(r2)-bl)**2
	else !LJ
	  if(r2<=rc2)then 
		r2i = 1/r2
		r6i=r2i**3
		addlj = 4*(r6i-1)*r6i-ecut
		pot_lj = pot_lj + addlj
!		pot_lj = pot_lj + 4*(r6i-1)*r6i-ecut
          endif
        endif
   enddo
  enddo
  end subroutine potential
!----------------------------------------------------------------------------------
  subroutine potential_fast(o,del_enolj,del_enostr)
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: del_enolj,del_enostr,r2i,r6i
  real(kind=8)               :: xr,r2
  integer                    :: i,j,o    !counters
!harmonic stretch pot:adjacent atoms in a chain.
!LJ : all nonbonded pairs, except atom pairs of harmonic stretch
!algorithm : molecule index = (atomindex+1)/natom + 1.
!ex) atom#1~8 : mol#1, 9~16:mol#2,...
!if diff of atomindex >=2, calc LJ. diff of atomindex=1-> if smaller one is multiple of 8, calc LJ
!otherwise: calc harmonic potential
  del_enolj=0.d0
  del_enostr=0.d0
  i = o
   do j = 1,npart
    if(j/=i)then
       dr = x(i,:)-x(j,:)
       dr = dr-box*nint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
        if(j-i==1 .AND. mod(i,natom)/=0)then !harmonic
                del_enostr = del_enostr + 0.5d0*kvib*(dsqrt(r2)-bl)**2
        else if(j-i==-1 .AND. mod(j,natom)/=0)then !another harmonic case
                del_enostr = del_enostr + 0.5d0*kvib*(dsqrt(r2)-bl)**2
	else!LJ
          if(r2<=rc2)then
                r2i = 1/r2
                r6i=r2i**3
                del_enolj = del_enolj + 4*(r6i-1)*r6i-ecut
          endif
        endif
    endif
   enddo
  end subroutine potential_fast
!-----------------------------------------------------------------------------------
  subroutine center   !keep particles in box
  implicit none
  integer                    :: i,j
  do i = 1,npart
   do j = 1,3
     if(x(i,j)>=box)then
       x(i,j) = x(i,j)-box
     endif
     if(x(i,j)<=0)then
       x(i,j) = x(i,j)+box
     endif
   enddo
  enddo
  end subroutine center

  subroutine pressure(p,PLRC)  !p is the pressure
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r2,r6i,p,rc3it,rc9it,PLRC
  integer                    :: i,j    !counters
  p = 0
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
  rc3it = (1.d0/rcut)**3/2.d0
  rc9it = (1.d0/rcut)**9/3.d0
  PLRC = (32*pi*density*npart*(rc9it-rc3it)/3.0D0)/volume
  p = (npart*tstar+16*p)/volume+PLRC	!P LRC
  end subroutine pressure

!-----------------------------------------------------------------------------------
  subroutine gr(iflag)
  implicit none
  integer, parameter                 :: num_bins=200000 
  integer                            :: iflag,jbin,i,j
  real*8, save,  dimension(num_bins) :: gn
  real*8,  dimension(3)              :: dr
  real*8, save                       :: rs,rmin,rmax,rdel,ridel,rb,const,rinner,router

! --------------------------------------------------------------------
  select case(iflag)
! ------------------------------------------------
  case(0)
   gn = 0   !initialization
   rmin = 0.8d0
   rmax = 0.5d0*box
   ridel = num_bins/(rmax-rmin)
   rdel = 1.d0/ridel
   const = 4.d0/3.d0*pi*density
   iflag = 1
! ------------------------------------------------
  case(1)
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
    open(16,file='gr.out')
    gn = gn/(ncalls*npart)
    rb = rmin+0.5*rdel
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

! ----------------------------------------------------------------
  subroutine mcmove(vpot,istep) !attempts to displace a particle
    implicit none
    integer                          :: j,o,istep
    real(kind=8), dimension(3)       :: xo
    real(kind=8)                     :: vpot,eno,enn,eno_lj,eno_str,enn_lj,enn_str,rnd

    call random_number(rnd)        !pick particle to be moved
    o=int(rnd*npart)+1
    call potential_fast(o,eno_lj,eno_str)  !calculate the potential energy
    xo = x(o,:)                         !current position
    eno=vpot
    do j = 1,3
      call random_number(rnd)
      x(o,j) = x(o,j) + (rnd-0.5)*delx  !new position
    enddo

    call potential_fast(o,enn_lj,enn_str)  !calculate the potential energy
    pot_lj = pot_lj + enn_lj - eno_lj
    pot_str = pot_str + enn_str - eno_str
    call random_number(rnd)
    enn = pot_lj + pot_str + e_lrc
    if(rnd<exp(-beta*(enn-eno)))then  !decide if you want to make a move
           istep = istep + 1  !count on accepted moves.
           vpot = enn
    else
           x(o,:)=xo
	   pot_lj = pot_lj - enn_lj + eno_lj
	   pot_str = pot_str - enn_str + eno_str
	   vpot = eno
    endif
  end subroutine mcmove

  subroutine com(zcm)
   integer  :: l
   real(kind=8), dimension(3)       :: zcm !center of mass
   do l = 1,3 !x,y,z
     zcm(l) =  sum(old_mol(:,l))/natom
   enddo
  end subroutine com

  subroutine mcvol(vpot,jstep,kstep) !attempts to displace a particle
    implicit none
    integer                          :: j,i,jstep,kstep
    real(kind=8)                     :: rnd       !random number
    real(kind=8)                     :: vpot,eno,enn,vo,lnvn,vn,boxo,arg
    real(kind=8)		     :: eno_lj,eno_str,e_lrcn
    real(kind=8), dimension(3)       :: zcm,dzcm !com and change in com
 !  -------------------------------------------------------------
 !  save the old values
    xm = x 
    boxo = box
    vo = box**3
 !  -------------------------------------------------------------
    call potential     !calculate the potential energy
    e_lrc = 8.0D0*pi*dfloat(npart)**2*((1/rcut)**9/3.0D0-(1/rcut)**3)/3.0D0/vo
    eno = pot_lj + pot_str + e_lrc
    eno_lj = pot_lj
    eno_str = pot_str
    call random_number(rnd)        !pick molecule to be moved
    lnvn=dlog(vo)+(rnd-0.5d0)*vmax
    vn = dexp(lnvn)
    box = vn**(1.d0/3.d0)         !new box length
!   --------------------------------------------------------------
! scale the coordinates : get com of each chain, scale coord of com and relocate each atom
!    x = box/boxo*xm               !scale the coordinates
 do i=1,nmol
   do j=1,natom
       old_mol(j,:) = x(natom*(i-1)+j,:)	!store old position of one molecule
   enddo
    call com(zcm)
    dzcm = zcm*(box/boxo-1.d0)                  !scale to new volume
   do j=1,natom
     x(natom*(i-1)+j,:) = old_mol(j,:)+dzcm    !new position
   enddo
 enddo
!   --------------------------------------------------------------
    call potential    !calculate the potential energy
    e_lrcn = 8.0D0*pi*dfloat(npart)**2*((1/rcut)**9/3.0D0-(1/rcut)**3)/3.0D0/vn !using constant rcut
    enn = pot_lj + pot_str + e_lrcn
!   --------------------------------------------------------------
    arg = -beta*((enn-eno)+pres*(vn-vo)-(npart+1)*dlog(vn/vo)/beta)
    call random_number(rnd)       
    kstep = kstep + 1
    if(rnd>dexp(arg))then
     box = boxo
     volume = vo
     x = xm                        !move is rejected
     pot_lj = eno_lj
     pot_str = eno_str
     vpot = eno
    else
     density = npart/vn
     volume = vn
     vpot = enn
     jstep = jstep + 1
    endif
  end subroutine mcvol

 END MODULE lj_met_chain_mod

