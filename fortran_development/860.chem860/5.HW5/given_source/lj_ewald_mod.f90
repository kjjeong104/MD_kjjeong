  MODULE lj_ewald_mod
      implicit none
      real*8, parameter     :: clight=2.99792458D10,av=6.0221367D23,hb=1.05457266D-27 
      real*8, parameter     :: mass_h=1.007825d0
      real*8, parameter     :: pi=3.141592654
      integer, parameter    :: npart=128,nbins=50
      integer, parameter    :: wp = 8

      integer, parameter             :: MAXK = 16000 
!     integer, PARAMETER             :: kmax = 5, KSQMAX = 27   !used for ewald sums
      integer, PARAMETER             :: kmax = 9, KSQMAX = 84   !used for ewald sums
      REAL(kind=8), dimension(maxk)  :: kvec    
      real(kind=8)                   :: kappa                   !width for charge distribution

      real*8                :: box,twopi,twopibox
      real*8 :: density,tstar 
      real(kind=8)                              :: trun,time_test,tcor,tsample,tstu,tequil,dt
      real(kind=8)                              :: rcut,rc2,volume
      integer, dimension(2) :: seed
      integer               :: nsamp,ncalls,mcalls
   

!   ------------------------------------------------------------------------------------------

          
  contains




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
     REAL(kind=8)                   :: B, RKX, RKY, RKZ, RKSQ

!    *******************************************************************
     twopi = 2*pi
     twopibox = 2*pi/box                !code was updated here
     B = 0.25d0/ (KAPPA*KAPPA)          !code was updated here

!    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **

        TOTK = 0
        DO  KX = 0, KMAX
           RKX = twopibox * REAL ( KX )          !code was updated here
           DO KY = -KMAX, KMAX  
              RKY = twopibox * REAL ( KY )       !code was updated here
              DO kz = -KMAX, KMAX
                 RKZ = twopibox * REAL ( KZ )    !code was updated here
                 KSQ = KX * KX + KY * KY + KZ * KZ
                 IF ( ( KSQ < KSQMAX ) .AND. ( KSQ /= 0 ) ) THEN
                    TOTK = TOTK + 1
                    IF ( TOTK > MAXK ) STOP 'KVEC IS TOO SMALL'
                    RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
                    KVEC(TOTK) = twopi * EXP ( -B * RKSQ ) / RKSQ

!                   write(*,'(3i4,f12.7)')kx,ky,kz,kvec(totk)
                 ENDIF
              enddo
           enddo
        enddo
        pause

        WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
        WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) TOTK
        
        end subroutine ewald_setup




        subroutine kwald (x,z,VK,fewald)

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

        REAL(kind=8), dimension(npart)   ::  RX, RY, RZ, Z
        REAL(kind=8), dimension(npart,3) ::  x,fewald 
        REAL(kind=8)                     ::  VK,t1,t2
        INTEGER                          ::  TOTK

        INTEGER                          ::  KX, KY, KZ, I, KSQ
        REAL(kind=8)                     ::  FACTOR, VD, VS
        real(kind=8),   PARAMETER        ::  RSQPI = 0.5641896 

        COMPLEX(kind=8), dimension(npart,-kmax:kmax) ::    eikz,EIKY
        COMPLEX(kind=8), dimension(npart,    0:kmax) ::    eikx
        COMPLEX(kind=8), dimension(npart)              ::    EIKR,eikt,fff
        complex(kind=8)                                :: rhok    !rhok is complex conj of 12.1.18 Frankel
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
!                   EIKR(:) = EIKX(:, KX) * EIKY(:, KY) * EIKZ(:, KZ)
                    EIKR(:) = EIKT(:) * EIKZ(:, KZ)
                    rhok = dconjg(sum(eikr))

                    VD = VD + FACTOR * KVEC(TOTK) * cdabs(rhok)**2 

                    fff = factor*kvec(totk)*dimag(eikr(:)*rhok)
                    fewald(:,1) = fewald(:,1) + kx*fff   !for problem 5, you need to fix if incorrect
                    fewald(:,2) = fewald(:,2) + ky*fff   !for problem 5, you need to fix if incorrect
                    fewald(:,3) = fewald(:,3) + kz*fff   !for problem 5, you need to fix if incorrect
                    
                 ENDIF
             enddo
           enddo
           factor = 2.d0
        enddo
         vd = vd/volume
        fewald = -fewald*4.d0*pi/(volume*box)
!       call cpu_time(t2)
!       print*,'t2-t1 = ',t2-t1

!    ** CALCULATES SELF PART OF K-SPACE SUM **
        VS = sum(z*z) 
        VS = RSQPI * KAPPA * vs
!       print*,' self part is ',vs
!    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **
        VK = VD  - VS
        end SUBROUTINE KWALD 

        SUBROUTINE RWALD (x,z, VR,frwald )
        implicit none
        REAL(kind=8), dimension(npart)   ::  z
        REAL(kind=8), dimension(npart,3) ::  x,frwald 
        REAL(kind=8), dimension(3)       ::  dr , derv
        REAL(kind=8)                     ::  VR, RIJSQ, RIJ, KRIJ, VIJ
        real(kind=8)                     ::  rinv, phi_sr,pre,fr
        INTEGER                          ::  I, J  !counters
!    *******************************************************************
!    ** CALCULATES R-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                     NUMBER OF IONS                  **
!    ** REAL    RX(N),RY(N),RZ(N)     POSITIONS OF IONS               **
!    ** REAL    Z(N)                  IONIC CHARGES                   **
!    ** REAL    VR                    R-SPACE POTENTIAL ENERGY        **
!    **                                                               **
!    ** ROUTINE REFERENCED:                                           **
!    **                                                               **
!    ** REAL FUNCTION ERFC ( X )                                      **
!    ** RETURNS THE COMPLEMENTARY ERROR FUNCTION                      **
!    *******************************************************************
!    *******************************************************************
        VR = 0.0
        frwald = 0
        pre = 2.d0*kappa/dsqrt(pi)
        DO i = 1, npart-1               !code was updated here
           DO j = i+1,npart             !code was updated here
              dr = x(i,:)-x(j,:)
              dr = dr - box*anint(dr/box)
              RIJSQ = dot_product(dr,dr)
              RIJ   = SQRT ( RIJSQ )
              rinv = 1.d0/rij
              KRIJ  = KAPPA * RIJ
              phi_sr =  z(i) * erfc(krij)*rinv
              fr     =  z(i) * pre* dexp (-KRIJ*krij ) 
              VR    = VR + z(j)*phi_sr
              derv = -z(j)*(fr+phi_sr) * dr*rinv**2
              frwald(i,:) = frwald(i,:) + derv
              frwald(j,:) = frwald(j,:) - derv
           enddo
         enddo
        end subroutine rwald

!*******************************************************************************
! Complementary error-function  erfc(x) = 1 - erf(x)
! Reference: W.J. Cody, Mathematics of Computation 22 (1969), 631-637
! Taken from CERN library.
!*******************************************************************************

      real(wp) function erfc(x)

      implicit none
      real(wp), intent(in)  :: x
      integer               :: j
      real(wp)              :: a,b,v,y
      real(wp), parameter   :: isqrtpi = 0.56418958354775629_wp   ! 1/sqrt(pi)
      real(wp), parameter   :: p1(0:3) = (/  2.426679552305318e+2_wp, &
                                             2.197926161829415e+1_wp, &
                                             6.996383488619136e+0_wp, &
                                            -3.560984370181538e-2_wp  /)
      real(wp), parameter   :: q1(0:3) = (/  2.150588758698612e+2_wp, &
                                             9.116490540451490e+1_wp, &
                                             1.508279763040779e+1_wp, &
                                             1.000000000000000e+0_wp  /)
      real(wp), parameter   :: p2(0:7) = (/  3.004592610201616e+2_wp, &
                                             4.519189537118729e+2_wp, &
                                             3.393208167343437e+2_wp, &
                                             1.529892850469404e+2_wp, &
                                             4.316222722205674e+1_wp, &
                                             7.211758250883094e+0_wp, &
                                             5.641955174789740e-1_wp, &
                                            -1.368648573827167e-7_wp  /)
      real(wp), parameter   :: q2(0:7) = (/  3.004592609569833e+2_wp, &
                                             7.909509253278980e+2_wp, &
                                             9.313540948506096e+2_wp, &
                                             6.389802644656312e+2_wp, &
                                             2.775854447439876e+2_wp, &
                                             7.700015293522947e+1_wp, &
                                             1.278272731962942e+1_wp, &
                                             1.000000000000000e+0_wp  /)
      real(wp), parameter   :: p3(0:4) = (/ -2.996107077035422e-3_wp, &
                                            -4.947309106232507e-2_wp, &
                                            -2.269565935396869e-1_wp, &
                                            -2.786613086096478e-1_wp, &
                                            -2.231924597341847e-2_wp  /)
      real(wp), parameter   :: q3(0:4) = (/  1.062092305284679e-2_wp, &
                                             1.913089261078298e-1_wp, &
                                             1.051675107067932e+0_wp, &
                                             1.987332018171353e+0_wp, &
                                             1.000000000000000e+0_wp  /)


      v = abs(x)
      if (v <= 0.46875_wp) then
          y = v**2
          a = p1(3)
          b = q1(3)
          do j=2,0,-1
              a = a*y + p1(j)
              b = b*y + q1(j)
          enddo
          erfc = 1.0_wp - x*a/b
          return
      elseif (v <= 4.0_wp) then
          a = p2(7)
          b = q2(7)
          do j=6,0,-1
              a = a*v + p2(j)
              b = b*v + q2(j)
          enddo
          erfc = exp(-v**2)*a/b
      elseif (v <= 10.0_wp) then
          y = 1.0_wp/(v**2)
          a = p3(4)
          b = q3(4)
          do j=3,0,-1
              a = a*y + p3(j)
              b = b*y + q3(j)
          enddo
          erfc = exp(-v**2)*(isqrtpi+y*a/b)/v
      else
          erfc = 0.0_wp
      endif
      if (x <= 0.0_wp) erfc = 2.0_wp - erfc

      end function erfc  


     REAL(kind=8)  FUNCTION erfcnn(x)
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
        ERFCnn = TP * EXP ( -XSQ )
        END function erfcnn


    END MODULE lj_ewald_mod




