  MODULE water_spc_mod
      implicit none
      real*8, parameter ::clight=2.99792458D10,av=6.0221367D23,hb=1.05457266D-27 
      real*8, parameter :: pi=3.141592654,kbolt = 8.31447d0
      real*8, parameter :: coulomb=1.60219172D-19  !coulomb charge
      real*8, parameter :: KC=8.98755D9            !1.d0/(4 pi epsilon)
      integer           :: const_vol
      real*8            :: vmax      !used for max volume change attempt
      real(kind=8)      :: e_conv            !convert from program units to kJ/mol
      integer, parameter:: nbins=50
      real*8           :: box
      real(kind=8)          :: mass(20)
      real*8                                    :: density,tstar,pressure 
      real(kind=8)                              :: beta,delx
      real(kind=8)                              :: rcut,rc2,ecut,volume
      integer, dimension(2) :: seed
      integer               :: iread 
      integer               :: nsamp
      integer               :: npart,nmol,natom
      integer               :: first           !used in ewald_setup
      real(kind=8)                             :: pot_lj    !lj contribution to potential
   
      integer                                   :: move_option
      real(kind=8), allocatable, dimension(:,:)  :: x,xm,xx,v,vi,f,fewald,frwald
      real(kind=8), allocatable, dimension(:,:)  :: old_mol,xtemp,xatom,new_mol,c12,c06  !geom of atoms; 6-12 paramters
      real(kind=8), allocatable, dimension(:)   :: amass    !masses of the atoms in a molecule
      real(kind=8)                              :: cmass    !total mass of molecule
      real(kind=8), allocatable, dimension(:)   :: qi,zi      !charges on the atoms
      integer,      allocatable, dimension(:)   :: ni        !atomic numbers
      real(kind=8),  allocatable,  dimension(:)  :: atom_temp !temporary array used in rotation moves
      real(kind=8)                              :: e_lrc     !corrections due to cutoff

!   ------------------------------------------------------------------------------------------
!     used with ewald sums
      real(kind=8)                             :: nkappa,kappa   !width for charge distribution
      real(kind=8)                             :: twopi          !2*pi
      real(kind=8)                             :: vself          !the intramolecular self interation
      real(kind=8)                             :: pot_lr,pot_sr  !the long and short range part of ewald sum
      integer                                  :: kmax,KSQMAX    !used for ewald sums
      REAL(kind=8), allocatable, dimension(:)  :: kvec
      real(kind=8), allocatable, dimension(:)  :: z,rx,ry,rz,drx,dry,drz
      COMPLEX(kind=8), allocatable, dimension(:,:) ::    eikz,eiky,dikz,diky
      COMPLEX(kind=8), allocatable, dimension(:,:) ::    eikx,dikx
      COMPLEX(kind=8), allocatable, dimension(:)            ::    EIKR,eikt,fff,dikr,dikt
      COMPLEX(kind=8), allocatable, dimension(:)            ::    rhok_old,rhok_new   !used for updating ewald lr terms
!   ------------------------------------------------------------------------------------------

          
    contains

 subroutine masses
    implicit none
    mass(1)=1.008
    mass(2)=4.002602
    mass(3)=6.941	
    mass(4)=9.012182
    mass(5)=10.81 
    mass(6)=  12.011
    mass(7)=  14.007
    mass(8)=  15.999	
    mass(9)=  18.9984032
    mass(10)=	20.1797
    mass(11)=   22.98976928
    mass(12)=	24.3050
    mass(13)=	26.9815386
    mass(14)=	28.085	
    mass(15)=	30.973762
    mass(16)=	32.06
    mass(17)=	35.45
    mass(18)=	39.948
    mass(19)=	39.0983
    mass(20)=	40.078
  end subroutine masses

  subroutine input(ncycle)
    implicit none
    integer      :: ncycle,i
    real(kind=8) :: zcm,ri3
    real(kind=8) :: u_leng,u_time,u_mass  !conversions for mks to amu ps angrsroms
    open(15,file='water_spc.inp')
    call masses
!   ----------------------------------------------------------------
!   read in the data file
    read(15,*)iread        !if iread=1 read in old_configuration
    read(15,*)nmol         ! number of molecules
    read(15,*)natom        ! number of atoms per molecule 
    read(15,*)const_vol    ! set to zero for NPT otherwise NVT
    npart = natom*nmol
    allocate (qi(natom),zi(natom*2),ni(natom),xtemp(3,natom),xatom(3,natom),old_mol(3,natom),amass(3))
    allocate (new_mol(3,natom),atom_temp(natom))
    allocate (c06(natom,natom),c12(natom,natom))
    read(15,*)qi
    do i = 1,natom
     read(15,*)ni(i),xatom(:,i)
     amass(i) = mass(ni(i))
     write(16,'(a,f10.3)')' amass(i) = ',amass(i) 
    enddo
    cmass = sum(amass)
    read(15,*)c12(1,1),c06(1,1)  !lj interaction
    read(15,*)nkappa    !  ewald gaussian width
    read(15,*)ksqmax    ! used to setup k space
    read(15,*)ncycle    ! number of points in equilibration
    read(15,*)nsamp     ! sampling every nsamp points
    read(15,*)delx      ! max mc step size
    read(15,*)vmax      ! max mc volume step size
    read(15,*)move_option      ! controls types of moves.
    read(15,*)tstar     ! the desired temperature
    read(15,*)pressure  ! the desired pressure 
    read(15,*)density   ! the initial value of the density in grams / cm**3  
    read(15,*)seed      ! used for initializing random number generator
    close (15)
!   ----------------------------------------------------------------
!   convert to ps-A-amu units
    u_mass = 1.d-3/av
    u_leng = 1.d-10
    u_time = 1.d-12
    tstar = tstar*0.1d0*kbolt
    c12 = c12*418.4d0
    c06 = c06*418.4d0
    density = density*av/(1.d24*cmass)                   !conversion to number density
    pressure= pressure*1.01325d5*u_leng*u_time**2/u_mass  !convert the pressure from atm
    e_conv = 418.4*nmol                                  !energy conversion to kcal/mol
!   ----------------------------------------------------------------
    volume  = nmol/density
    npart = nmol*natom
    box = volume**(1.d0/3.d0)
    beta = 1.d0/tstar
!    rcut = 7.5d0                               !fixed value (can change is e_lrc is corrected.)
    rcut = box*0.5d0                           !this value is better, but don't use unless e_lrc is included.
    if(rcut>0.5d0*box)write(*,*)'box is too small for rcut value'
    if(rcut>0.5d0*box)stop
    kappa = nkappa/box
    rc2 = rcut*rcut
!   ----------------------------------------------------------------
!   need to calculate long range correction to 6-12
    ri3 = 1.d0/rcut**3
    e_lrc = 2*pi*density*(c12(1,1)*ri3/3.0D0-c06(1,1))*ri3/3.0D0   ! some function of r13 nmol c12 and c06  
    write(*,'(a,f10.5)')' e_lrc = ',e_lrc/e_conv/volume
    if(e_lrc==0) write(*,*)' You need to calcuate the long range correction.'
!   ----------------------------------------------------------------

    write(16,'(a,f8.4)')' the following are in A ps amu units'
    write(16,'(a,f8.4)')' the desired temp is    ',tstar
    write(16,'(a,i5)')' the move_option is    ',move_option
    write(16,'(a,f8.3)')' the box length is      ',box
    write(16,'(a,i8)')  ' number of particles is ',npart
    write(16,'(a,f8.5)')' the density is         ',density
    write(16,'(a,f8.5)')' the pressure is        ',pressure
    write(16,'(a,f8.5)')' the value of beta  is  ',beta
    write(16,'(a,f8.3)')' the volume is          ',volume
    write(16,'(a,f8.3)')' r cutoff is            ',rcut
    write(16,'(a,f8.3)')' molar mass is          ',cmass
!   ----------------------------------------------------------------
    allocate(x(npart,3),v(npart,3),vi(npart,3),xm(npart,3),xx(npart,3))
    allocate(f(npart,3),fewald(npart,3),frwald(npart,3))
    allocate(z(npart))
! ----------------------------------------------------------------
!   move atoms to center of mass frame
    write(16,*)' atoms in center of mass frame '
    do i = 1,3
     zcm = sum(xatom(i,:)*amass(:))/cmass
     xatom(i,:) = xatom(i,:)-zcm
    enddo
    write(16,'(3f10.5)')xatom
! ----------------------------------------------------------------
  end subroutine input

! ----------------------------------------------------------------
  subroutine init !set up an initial condition using fcc
  implicit none
  integer :: i,j,k,l,iat,ipart,nt
  real(kind=8)    :: delta  !spacing between atom(1)


! pick points on a  cubic grid
   
  ipart  = (1.00001d0*dfloat(nmol))**(1.0/3.0)
  delta = box/dfloat(ipart)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if(iread==1)then
   open(18,file='old_config.ini')
   read(18,*) box
   write(16,'(a,f10.4)')' the current box length is ',box
   volume = box**3
   read(18,*)x
   close(18)
 else
  nt = 0
  if(nmol==2)then ! just two molecules
     delta = 2
     do iat = 1,natom
       nt = nt + 1
       x(nt,1) = delta*i + xatom(1,iat)
       x(nt,2) = delta*j + xatom(2,iat)
       x(nt,3) = delta*k + xatom(3,iat)
     enddo
     do iat = 1,natom
       nt = nt + 1
       x(nt,1) = delta + xatom(1,iat)
       x(nt,2) = delta + xatom(2,iat)
       x(nt,3) = delta + xatom(3,iat)
     enddo
  else
   !put on a grid
    do i = 0,ipart-1
     do j = 0,ipart-1
      do k = 0,ipart-1
       do iat = 1,natom
         nt = nt + 1
         x(nt,1) = delta*i + xatom(1,iat)
         x(nt,2) = delta*j + xatom(2,iat)
         x(nt,3) = delta*k + xatom(3,iat)
       enddo
      enddo
     enddo
    enddo
  endif

  if(nt/=npart)then
     write(16,*)' stopping since nt is not equal to npart '
     write(16,*)' nt    = ',nt
     write(16,*)' npart = ',npart
     stop
  endif
 endif
! ------------------------------------------------------
  write(16,*)' here are the original charges'
  write(16,*)qi
  qi = qi*coulomb*dsqrt(0.1d0*av*kc*1.d10)    !convert to ps A amu units
  write(16,*)' here are the scaled charges'
  write(16,*)qi
  k = 0
  do i = 1,nmol
   do j = 1,natom
    k = k + 1
    if(k<=2*natom)zi(k) = qi(j)   !used in kwald_fast
    z(k) = qi(j)
   enddo
  enddo
!  ecut = (c12(1,1)/rc2**6-c06(1,1))/rc2**3  !cutoff value for the 6-12 potential     ! I think typo.
  ecut = (c12(1,1)/rc2**3-c06(1,1))/rc2**3
!  ecut = 0
  write(16,*)' ecut = ',ecut
  end subroutine init

!-----------------------------------------------------------------------------------
  subroutine potential    !calculate the potentials for the LJ potential
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: r6i
  real(kind=8)               :: r2
  integer                    :: i,j
  pot_lj = 0.d0
  do i = 1,npart-natom,natom  !just over the oxygen atoms here
   do j = i+natom,npart,natom  !just over the oxygen atoms here
       dr = x(i,:)-x(j,:)
       dr = dr-box*anint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
       if(r2<=rc2)then
         r6i=1.d0/r2**3
         pot_lj =  pot_lj + (c12(1,1)*r6i-c06(1,1))*r6i-ecut
       endif
    enddo
  enddo
  end subroutine potential

!-----------------------------------------------------------------------------------
  subroutine center   !keep particles in box
  implicit none
  integer                    :: i,j,ip
  do i = 1,npart,natom
   do j = 1,3
     if(x(i,j)>=box)then
       do ip = i,i+2
         x(ip,j) = x(ip,j)-box
       enddo
     endif
     if(x(i,j)<=0)then
       do ip = i,i+2
         x(ip,j) = x(ip,j)+box
       enddo
     endif
   enddo
  enddo
  end subroutine center

!*******************************************************************************
!* FICHE F.22.  ROUTINES TO PERFORM THE EWALD SUM                             **
!* This FORTRAN code is intended to illustrate points made in the text.       **
!* To our knowledge it works correctly.  However it is the responsibility of  **
!* the user to test it, if it is to be used in a research application.        **
!*******************************************************************************

!    *******************************************************************
!    ** REAL-SPACE AND RECIPROCAL-SPACE PARTS OF EWALD SUM FOR IONS.  **
!    **                                                               **
!    ** REFERENCES:                                                   **
!    **                                                               **
!    ** WOODCOCK AND SINGER, TRANS. FARADAY SOC. 67, 12, 1971.        **
!    ** DE LEEUW ET AL., PROC. ROY. SOC. A 373, 27, 1980.             **
!    ** HEYES, J. CHEM. PHYS. 74, 1924, 1981.                         **
!    ** SEE ALSO FINCHAM, MDIONS, CCP5 PROGRAM LIBRARY.               **
!    **                                                               **
!    ** ROUTINES SUPPLIED:                                            **
!    **                                                               **
!    ** SUBROUTINE SETUP ( KAPPA )                                    **
!    **    SETS UP THE WAVEVECTORS FOR USE IN THE EWALD SUM           **
!    ** SUBROUTINE RWALD ( KAPPA, VR )                                **
!    **    CALCULATES THE R-SPACE PART OF THE SUM                     **
!    ** SUBROUTINE KWALD ( KAPPA, VK )                                **
!    **    CALCULATES THE K-SPACE PART OF THE SUM                     **
!    ** REAL FUNCTION ERFC ( X )                                      **
!    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER  TOTK         THE TOTAL NUMBER OF K-VECTORS STORED    **
!    ** INTEGER  MAXK         MAXIMUM POSSIBLE NUMBER OF K-VECTORS    **
!    ** INTEGER  KMAX         MAX INTEGER COMPONENT OF THE K-VECTOR   **
!    ** INTEGER  KSQMAX       MAX SQUARE MOD OF THE K-VECTOR REQUIRED **
!    ** REAL     VR           ENERGY FROM R-SPACE SUM                 **
!    ** REAL     VK           ENERGY FROM K-SPACE SUM                 **
!    ** REAL     KVEC(MAXK)   ARRAY USED TO STORE K-VECTORS           **
!    ** REAL     KAPPA        WIDTH OF CANCELLING DISTRIBUTION        **
!    **                                                               **
!    ** USAGE:                                                        **
!    **                                                               **
!    ** SETUP IS CALLED ONCE AT THE BEGINNING OF THE SIMULATION       **
!    ** TO CALCULATE ALL THE K-VECTORS REQUIRED IN THE EWALD SUM.     **
!    ** THESE VECTORS ARE USED THROUGHOUT THE SIMULATION IN THE       **
!    ** SUBROUTINE KWALD TO CALCULATE THE K-SPACE CONTRIBUTION TO THE **
!    ** POTENTIAL ENERGY AT EACH CONFIGURATION. THE SELF TERM IS      **
!    ** SUBTRACTED FROM THE K-SPACE CONTRIBUTION IN KWALD.            **
!    ** THE SURFACE TERM FOR SIMULATIONS IN VACUUM IS NOT INCLUDED.   **
!    ** ROUTINE RWALD RETURNS THE R-SPACE CONTRIBUTION TO THE EWALD   **
!    ** SUM AND IS CALLED FOR EACH CONFIGURATION IN THE SIMULATION.   **
!    ** A CUBIC BOX AND UNIT BOX LENGTH ARE ASSUMED THROUGHOUT.       **
!    *******************************************************************



     SUBROUTINE ewald_setup
     implicit none
!    *******************************************************************
!    ** ROUTINE TO SET UP THE WAVE-VECTORS FOR THE EWALD SUM.         **
!    **                                                               **
!    ** THE WAVEVECTORS MUST FIT INTO A BOX OF UNIT LENGTH.           **
!    ** IN THIS EXAMPLE WE ALLOW A MAXIMUM OF 1000 WAVEVECTORS.       **
!    *******************************************************************
     INTEGER                        :: KSQ, KX, KY, KZ, TOTK
     REAL(kind=8)                   :: B, RKX, RKY, RKZ, RKSQ,twopibox

!    *******************************************************************
     if(first==1)then
      kmax = dsqrt(dfloat(ksqmax)) 
      allocate (eikz(npart,-kmax:kmax),eiky(npart,-kmax:kmax))
      allocate (eikx(npart,0:kmax))
      allocate (eikr(npart),rx(npart),ry(npart),rz(npart),eikt(npart))
      allocate (dikz(2*natom,-kmax:kmax),diky(2*natom,-kmax:kmax))
      allocate (dikx(2*natom,0:kmax))
      allocate (dikr(2*natom),drx(2*natom),dry(2*natom),drz(2*natom),dikt(2*natom))
      twopi = 2*pi
 
!     LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **
      TOTK = 0
      DO  KX = 0, KMAX
         DO KY = -KMAX, KMAX
           DO kz = -KMAX, KMAX
             KSQ = KX * KX + KY * KY + KZ * KZ
             IF ( ( KSQ < KSQMAX ) .AND. ( KSQ /= 0 ) )  TOTK = TOTK + 1
           enddo
         enddo
      enddo
      allocate(kvec(totk),rhok_old(totk),rhok_new(totk))  !here I assign the space
      first=2
     endif
 
     twopibox = twopi/box 
     B = 0.25d0/ (KAPPA*KAPPA)

        TOTK = 0
        DO  KX = 0, KMAX
           RKX = twopibox * REAL ( KX )
           DO KY = -KMAX, KMAX
              RKY = twopibox * REAL ( KY )
              DO kz = -KMAX, KMAX
                 RKZ = twopibox * REAL ( KZ )
                 KSQ = KX * KX + KY * KY + KZ * KZ
                 IF ( ( KSQ < KSQMAX ) .AND. ( KSQ /= 0 ) ) THEN
                    TOTK = TOTK + 1
                    RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
                    KVEC(TOTK) = twopi * EXP ( -B * RKSQ ) / RKSQ

                 ENDIF
              enddo
           enddo
        enddo

        
        end subroutine ewald_setup




        subroutine kwald

!    *******************************************************************
!    ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
!    **                                                               **
!    ** THE SELF TERM IS SUBTRACTED.                                  **
!    ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
!    ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
!    ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
!    ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER npart               NUMBER OF IONS                    **
!    ** REAL    RX(Npart),RY(Npart),RZ(Npart)   POSITIONS OF IONS     **
!    ** REAL    Z(N)                IONIC CHARGES                     **
!    ** REAL    VK                  K-SPACE POTENTIAL ENERGY          **
!    ** REAL    VKS                 SELF PART OF K-SPACE SUM          **
!    *******************************************************************
        implicit none

!       REAL(kind=8), dimension(npart)   ::  RX, RY, RZ, Z
!       REAL(kind=8), dimension(npart,3) ::  x,fewald 
        REAL(kind=8)                     ::  VK,t1,t2
        INTEGER                          ::  TOTK,o

        INTEGER                          ::  KX, KY, KZ, I, KSQ
        REAL(kind=8)                     ::  FACTOR, VD, VS
        real(kind=8),   PARAMETER        ::  RSQPI = 0.5641896 

!       COMPLEX(kind=8), dimension(npart,-kmax:kmax) ::    eikz,EIKY
!       COMPLEX(kind=8), dimension(npart,    0:kmax) ::    eikx
!       COMPLEX(kind=8), dimension(npart)            ::    EIKR,eikt,fff
        complex(kind=8)                              :: rhok    !rhok is complex conj of 12.1.18 Frankel
!    *******************************************************************
!    ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **
!    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **

       

        rx = x(:,1)/box
        ry = x(:,2)/box
        rz = x(:,3)/box

        eikx(:,0) = (1.0, 0.0)
        eiky(:,0) = (1.0, 0.0)
        eikz(:,0) = (1.0, 0.0)
        EIKX(:, 1) = dCMPLX ( dcOS ( twopi * rx) ,dSIN ( twopi * RX) )
        EIKY(:, 1) = dCMPLX ( dcOS ( twopi * RY) ,dSIN ( twopi * RY) )
        EIKZ(:, 1) = dCMPLX ( dCOS ( twopi * RZ) ,dSIN ( twopi * RZ) ) 
        EIKY(:, -1) = CONJG ( EIKY(:, 1) )
        EIKZ(:, -1) = CONJG ( EIKZ(:, 1) )


!    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **
        DO KX = 2, KMAX
           EIKX(:, KX) = EIKX(:, KX-1) * EIKX(:, 1)
        enddo
        DO KY = 2, KMAX
           EIKY(:,  KY) = EIKY(:, KY-1) * EIKY(:, 1)
           EIKY(:, -KY) = CONJG ( EIKY(:, KY) )
        enddo
        DO  KZ = 2, KMAX
           EIKZ(:,  KZ) = EIKZ(:, KZ-1) * EIKZ(:, 1)
           EIKZ(:, -KZ) = CONJG ( EIKZ(:, KZ) )
        enddo

!    ** SUM OVER ALL VECTORS **

        call cpu_time(t1)
        VD   = 0.0
        fewald = 0.d0
        TOTK = 0
        factor = 1.0
        DO KX =  0, KMAX
           DO KY =  -KMAX, KMAX
              EIKT(:) = EIKX(:, KX) * EIKY(:, KY) * z
              DO KZ =  -KMAX, KMAX
                 KSQ = KX * KX + KY * KY + Kz*kz
                 IF ( ( KSQ < KSQMAX ) .AND. ( KSQ /= 0 ) ) THEN
                    TOTK = TOTK + 1
                    EIKR(:) = EIKT(:) * EIKZ(:, KZ)
                    rhok = sum(eikr)
                    rhok_old(totk) = rhok

                    VD = VD + FACTOR * KVEC(TOTK) * cdabs(rhok)**2 

                    
                 ENDIF
             enddo
           enddo
           factor = 2.d0
        enddo
        vd = vd/volume
!       call cpu_time(t2)
!       print*,'t2-t1 = ',t2-t1

!    ** CALCULATES SELF PART OF K-SPACE SUM **
        VS = sum(z*z) 
        VS = RSQPI * KAPPA * vs
!       print*,' self part is ',vs
!    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **
        pot_lr = VD  - VS
        end SUBROUTINE KWALD 

        subroutine rwald
        implicit none
        REAL(kind=8), dimension(3)       ::  dr 
        REAL(kind=8)                     ::  RIJSQ, RIJ, KRIJ, VIJ
        real(kind=8)                     ::  phi_sr
        INTEGER                          ::  I, J  !counters
        INTEGER                          ::  o     !molecule being moved
!    *******************************************************************
!    ** CALCULATES R-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                     NUMBER OF IONS                  **
!    ** REAL    RX(N),RY(N),RZ(N)     POSITIONS OF IONS               **
!    ** REAL    Z(N)                  IONIC CHARGES                   **
!    ** REAL    pot_sr                R-SPACE POTENTIAL ENERGY        **
!    **                                                               **
!    ** ROUTINE REFERENCED:                                           **
!    **                                                               **
!    ** REAL FUNCTION ERFC ( X )                                      **
!    ** RETURNS THE COMPLEMENTARY ERROR FUNCTION                      **
!    *******************************************************************
!    *******************************************************************
        pot_sr = 0.0
        DO i = 1,npart-1
            dO j = i+1,npart
                dr = x(i,:)-x(j,:)
                dr = dr - box*anint(dr/box)
                RIJSQ = dot_product(dr,dr)
                RIJ   = SQRT ( RIJSQ )
                KRIJ  = KAPPA * RIJ
                phi_sr =  z(i) * erfc(krij)/rij
                pot_sr    = pot_sr + z(j)*phi_sr
           enddo
        enddo
        end subroutine rwald


     REAL(kind=8)  FUNCTION erfc(x)
!    *******************************************************************
!    ** APPROXIMATION TO THE COMPLEMENTARY ERROR FUNCTION             **
!    **                                                               **
!    ** REFERENCE:                                                    **
!    **                                                               **
!    ** ABRAMOWITZ AND STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,    **
!    **    NATIONAL BUREAU OF STANDARDS, FORMULA 7.1.26               **
!    *******************************************************************
        real(kind=8), PARAMETER ::  A1 = 0.254829592, A2 = -0.284496736, & 
                                    A3 = 1.421413741, A4 = -1.453152027, &
                                    A5 = 1.061405429, P  =  0.3275911   
        real(kind=8)            :: T, X, XSQ, TP
!    *******************************************************************
        T  = 1.0 / ( 1.0 + P * X )
        XSQ = X * X
        TP = T * ( A1 + T * ( A2 + T * ( A3 + T * ( A4 + T * A5 ) ) ) )
        eRFC = TP * EXP ( -XSQ )
        END function erfc


  subroutine self
  implicit none
  real(kind=8)                       :: rint    !related to H - O - H distances
  integer                            :: i,j 
  vself = 0.d0
  do i =1,natom-1
   do j = i+1,natom
    rint = dsqrt(sum((xatom(:,i)-xatom(:,j))**2))
    vself = qi(i)*qi(j)/rint + vself
   enddo
  enddo
  vself = vself*nmol
  write(16,*)' the self energy is ',vself/e_conv
  end subroutine self
  
  subroutine rotate
    implicit none
    integer                          :: i,j,l
    real(kind=8)                     :: rnd,z,sthe,psi,cpsi,spsi,phi,cphi,sphi
    real(kind=8), dimension(3,3)     :: a
    real(kind=8), dimension(3)       :: zcm !center of mass
     
  
    call random_number(rnd)        !pick angle to be moved
    phi = rnd*twopi
    call random_number(rnd)        !pick angle to be moved
    psi = rnd*twopi
    call random_number(rnd)        !pick angle to be moved
    z = -1.d0+2.d0*rnd
    cphi = dcos(phi)
    sphi = dsin(phi)
    cpsi = dcos(psi)
    spsi = dsin(psi)
    sthe = dsqrt(1.d0-z**2)
! Rotation matrix of Euler Angles, Z_1*X_2*Y_3. ref http://en.wikipedia.org/wiki/Euler_angles    
    a(1,1) =  cpsi*cphi-z*sphi*spsi  
    a(1,2) = -spsi*cphi-z*sphi*cpsi
    a(1,3) =  sthe*sphi

    a(2,1) =  cpsi*sphi+z*cphi*spsi
    a(2,2) = -spsi*sphi+z*cphi*cpsi
    a(2,3) = -sthe*cphi

    a(3,1) = sthe*spsi
    a(3,2) = sthe*cpsi
    a(3,3) = z


    xtemp   = matmul(a,xatom)


  end subroutine rotate

  subroutine rotate_trans
    implicit none
    integer                          :: i,j,l
    !somethin is missing
  end subroutine rotate_trans

  subroutine com(zcm)
   integer  :: l
   real(kind=8), dimension(3)       :: zcm !center of mass
   do l = 1,3 !x,y,z
     zcm(l) =  sum(old_mol(l,:)*amass)/cmass
   enddo
  end subroutine com


  subroutine mcmove(vpot,istep) !attempts to displace a particle
    implicit none
    integer                          :: j,i,o,istep
    real(kind=8), dimension(3)       :: zcm,dzcm !com and change in com
    real(kind=8)                     :: vpot,eno,enn,rnd,vr,vo
    real(kind=8)                     :: eno_sr,enn_sr,eno_lj,enn_lj,eno_lr,enn_lr
 
    eno_lr = pot_lr
!    eno_sr = pot_sr   ! I added

    call random_number(rnd)        !pick molecule to be moved
     o=int(rnd*nmol)+1
 !   ------------------------------------------------------------------
 !   calculate two of the potential changes for moving a molecule
     call potential_fast(o,eno_lj)        !calculate the potentials for the LJ potential
     call rwald_fast (o,eno_sr)           ! short range
 !   ------------------------------------------------------------------
 !   save the old total energy and select out the molecule being moved
     eno = vpot
     do j = 1,natom
       old_mol(:,j) = x(natom*(o-1)+j,:)                         !current position
     enddo

 !   ------------------------------------------------------------------
 !   calculate random translational displacement
     do j = 1,3
      call random_number(rnd)
      dzcm(j) =  (rnd-0.5)*delx  !randomly choose new com position
     enddo
 !   ------------------------------------------------------------------
     call com(zcm)    !find old com of old_mol
 !   ------------------------------------------------------------------
 !   make one of three kinds of moves
     select case(move_option)
     case(1)
       call rotate   !rotate to random orientation
       zcm = zcm+dzcm   !new com position
       do j = 1,natom
          new_mol(:,j) = xtemp(:,j) + zcm        !new position of atom
          x(natom*(o-1)+j,:) = new_mol(:,j)                  !update x vector position
       enddo
     case(2)
       call rotate_trans   !rotate to random orientation
       write(*,*)' this case does not exist'
     case(3)
       write(*,*)' this case does not exist'
       pause
     end select
!    ---------------------------------------------------------------
!    calculate the 3 contributions to the energies
     call rwald_fast (o,enn_sr)                ! short range
     pot_sr = pot_sr + enn_sr - eno_sr         ! new short range

     call potential_fast(o,enn_lj)        ! 6-12
     pot_lj = pot_lj + enn_lj - eno_lj         ! new 6-12

     call ewald_setup
     call kwald_fast                          ! new long range
     enn = pot_sr+pot_lr-vself+pot_lj         ! new potential
!     enn = pot_sr+pot_lr-vself+pot_lj+e_lrc/vo  ! I added
!    ---------------------------------------------------------------
!    decide to make a move based on the metropolis algorithm
     call random_number(rnd)
     if(rnd<exp(-beta*(enn-eno)))then     !decide if you want to make a move
       vpot = enn                         !accept move
       rhok_old = rhok_new
       istep = istep + 1                  !count on accepted moves.
     else                                 !reject move
       vpot = eno
       pot_lr = eno_lr
       pot_sr = pot_sr - enn_sr + eno_sr
       pot_lj = pot_lj - enn_lj + eno_lj
       do j = 1,natom
         x(natom*(o-1)+j,:)=old_mol(:,j)  !retain old position
       enddo
     endif
  end subroutine mcmove




        subroutine kwald_fast

!    *******************************************************************
!    ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
!    **                                                               **
!    ** THE SELF TERM IS SUBTRACTED.                                  **
!    ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
!    ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
!    ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
!    ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER npart               NUMBER OF IONS                    **
!    ** REAL    RX(Npart),RY(Npart),RZ(Npart)   POSITIONS OF IONS     **
!    ** REAL    Z(N)                IONIC CHARGES                     **
!    ** REAL    pot_lr                  K-SPACE POTENTIAL ENERGY          **
!    ** REAL    VKS                     SELF PART OF K-SPACE SUM          **
!    *******************************************************************
        implicit none

        REAL(kind=8)                     ::  VK,t1,t2
        INTEGER                          ::  TOTK,o

        INTEGER                          ::  KX, KY, KZ, I, KSQ
        REAL(kind=8)                     ::  FACTOR, VD, VS
        real(kind=8),   PARAMETER        ::  RSQPI = 0.5641896 

        complex(kind=8)                              ::    rhok    !rhok is form 12.1.18 Frankel
!    *******************************************************************
!    ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **
!    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **

       

        do i = 1,natom
         drx(i) = new_mol(1,i)/box
         dry(i) = new_mol(2,i)/box
         drz(i) = new_mol(3,i)/box
         drx(i+natom) = old_mol(1,i)/box
         dry(i+natom) = old_mol(2,i)/box
         drz(i+natom) = old_mol(3,i)/box
        enddo
         


        dikx(:,0) = (1.0, 0.0)
        diky(:,0) = (1.0, 0.0)
        dikz(:,0) = (1.0, 0.0)
        dikX(:, 1) = dCMPLX ( dcOS ( twopi * drx) ,dSIN ( twopi * dRX) )
        dikY(:, 1) = dCMPLX ( dcOS ( twopi * dRY) ,dSIN ( twopi * dRY) )
        dikZ(:, 1) = dCMPLX ( dCOS ( twopi * dRZ) ,dSIN ( twopi * dRZ) ) 
        dikY(:, -1) = CONJG ( dikY(:, 1) )
        dikZ(:, -1) = CONJG ( dikZ(:, 1) )


!    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **
        DO KX = 2, KMAX
           dikX(:, KX) = dikX(:, KX-1) * dikX(:, 1)
        enddo
        DO KY = 2, KMAX
           dikY(:,  KY) = dikY(:, KY-1) * dikY(:, 1)
           dikY(:, -KY) = CONJG ( dikY(:, KY) )
        enddo
        DO  KZ = 2, KMAX
           dikZ(:,  KZ) = dikZ(:, KZ-1) * dikZ(:, 1)
           dikZ(:, -KZ) = CONJG ( dikZ(:, KZ) )
        enddo

!    ** SUM OVER ALL VECTORS **

!       call cpu_time(t1)
        VD   = 0.0
        TOTK = 0
        factor = 1.0
        DO KX =  0, KMAX
           DO KY =  -KMAX, KMAX
              dikt(:) = dikx(:, KX) * diky(:, KY) * zi 
              DO KZ =  -KMAX, KMAX
                 KSQ = KX * KX + KY * KY + Kz*kz
                 IF ( ( KSQ < KSQMAX ) .AND. ( KSQ /= 0 ) ) THEN
                    TOTK = TOTK + 1
                    dikr(:) = dikt(:) * dikz(:, KZ)
                    rhok = rhok_old(totk)
                    do i = 1,natom
                     rhok = rhok + dikr(i) - dikr(i+natom)
                    enddo
                    rhok_new(totk) = rhok

                    VD = VD + FACTOR * KVEC(TOTK) * cdabs(rhok)**2 
                    
                 ENDIF
             enddo
           enddo
           factor = 2.d0
        enddo
        vd = vd/volume
!       call cpu_time(t2)
!       print*,'t2-t1 = ',t2-t1

!    ** CALCULATES SELF PART OF K-SPACE SUM **
        VS = sum(z*z) 
        VS = RSQPI * KAPPA * vs
!       print*,' self part is ',vs
!    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **
        pot_lr = VD  - VS
        end SUBROUTINE KWALD_fast

        SUBROUTINE RWALD_fast (o,del_pot_sr)
        implicit none
        REAL(kind=8), dimension(3)       ::  dr 
        REAL(kind=8)                     ::  VR, RIJSQ, RIJ, KRIJ, VIJ,del_pot_sr
        real(kind=8)                     ::  rinv, phi_sr
        INTEGER                          ::  I, J,init,ifin  !counters
        INTEGER                          ::  o               !molecule being moved


        VR = 0.0
        init = natom*(o-1)+1
        ifin = init+natom-1
        DO i = init,ifin 
           DO j = 1,npart
              if(j<init.or.j>ifin)then
                dr = x(i,:)-x(j,:)
                dr = dr - box*anint(dr/box)
                RIJSQ = dot_product(dr,dr)
                RIJ   = SQRT ( RIJSQ )
                KRIJ  = KAPPA * RIJ
                phi_sr =  z(i) * erfc(krij)/rij
                VR    = VR + z(j)*phi_sr
              endif
           enddo
         enddo
         del_pot_sr = vr
        end subroutine rwald_fast

!-----------------------------------------------------------------------------------
  subroutine potential_fast(o,del_eno)   !calculate the potentials for the LJ potential
  implicit none
  real(kind=8), dimension(3) :: dr
  real(kind=8)               :: del_eno,r2i,r6i
  real(kind=8)               :: r2
  integer                    :: i,j,o 
  del_eno = 0.d0
  i = natom*(o-1)+1
  do j = 1,npart,natom  !just over the oxygen atoms here
     if(j/=i)then
       dr = x(i,:)-x(j,:)
       dr = dr-box*nint(dr/box)  !periodic boundary conditions
       r2 = dot_product(dr,dr)
       if(r2<=rc2)then
         r2i = 1/r2
         r6i=r2i**3
         del_eno =  del_eno + (c12(1,1)*r6i-c06(1,1))*r6i-ecut
       endif
     endif
  enddo
  end subroutine potential_fast


  subroutine mcvol(vpot,jstep) !attempts to displace a particle
    implicit none
    integer                          :: j,i,jstep
    real(kind=8)                     :: rnd       !random number
    real(kind=8), dimension(3)       :: zcm,dzcm  !com and change in com
    real(kind=8)                     :: vpot,eno,enn,vo,lnvn,vn,boxo,boxn,arg
    real(kind=8)                     :: pot_lj_old,pot_sr_old,pot_lr_old,ri3,ri3n,e_lrcn
 !  --------------------------------------
 !  save the old values
    xm = x 
    pot_lj_old = pot_lj
    pot_lr_old = pot_lr
    pot_sr_old = pot_sr
    boxo = box
    vo = box**3
!    ---------------------------------------------------------------------------
!    ri3 = 8.d0/vo
    ri3 = 1.d0/(box/2.0D0)**3
    e_lrc = 2*pi*(density)*(c12(1,1)*ri3/3.0D0-c06(1,1))*ri3/3.0D0  !this needs to be fixed
!    print*, e_lrc/vo
    if(e_lrc==0)write(*,*)' fix the above expression '
    eno = pot_lr+pot_sr+pot_lj-vself+e_lrc
!   ------------------------------------------
    call random_number(rnd)            !pick a random number for volume move
    lnvn=dlog(vo)+(rnd-0.5d0)*vmax
    vn = dexp(lnvn)
    boxn = vn**(1.d0/3.d0)             !new box length
!   ------------------------------------------
!   determine the new configuration by scaling com
    do i = 1,nmol
      do j = 1,natom
        old_mol(:,j) = x(natom*(i-1)+j,:)          !old position
      enddo
      call com(zcm)                                !find com of old_mol
      dzcm = zcm*(boxn/boxo-1.d0)                  !scale to new volume
      do j = 1,natom
         x(natom*(i-1)+j,:) = old_mol(:,j)+dzcm    !new position
      enddo
    enddo
!   ------------------------------------------
!   calculate energies at new position
    box = boxn
    volume = vn
    call ewald_setup
    call kwald                      !long range ewald -  pot_lr
    call rwald                      !short range ewald - pot_sr
    call potential                  !6 12 potential    - pot_lj
!    enn = pot_lr+pot_sr+pot_lj-vself    !why e_lrc is missed??????
    ri3n = 1.d0/(box/2.0D0)**3 ! I added
    e_lrcn = 2*pi*(density*vo/vn)**2*(c12(1,1)*ri3n/3.0D0-c06(1,1)*ri3n)/3.0D0 ! I added
    enn = pot_lr+pot_sr+pot_lj-vself+e_lrcn    ! I added
!   ------------------------------------------
!   decide if you want to take a move
    arg = -beta*((enn-eno) + pressure*(vn-vo)-(nmol+1)*dlog(vn/vo)/beta)
    call random_number(rnd)        !pick a random number
    if(rnd>dexp(arg))then
!     vpot = eno     ! I added
     box = boxo
     volume = vo
     x = xm                        !move is rejected
     pot_lj = pot_lj_old           !restore the original values
     pot_lr = pot_lr_old 
     pot_sr = pot_sr_old 
    else
     vpot = enn     ! I added
     e_lrc = e_lrcn  ! I added
     jstep = jstep + 1             !keep track of the accepted moves
    endif
  end subroutine mcvol


  END MODULE water_spc_mod

