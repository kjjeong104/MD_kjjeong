     module fourier_mod
     implicit none
     real*8, parameter :: pi = 3.141592654

     contains

     subroutine FFTET(ft,split,hold,isign)
!	this program is the driver program for taking the fast fourier 
!	transform (time to energy) of the provided
!	function,  the fourier transform itself
!	is done by the subroutine four1 from numerical recipes
!	timestep is in femtoseconds when read in and is then
!	converted to atomic time units 
!	final energies are in wavenumbes
!	
!	-----------------------------------------------------------
	implicit double precision (a-h,o-z)
	real(kind=8), dimension(:)  :: split, hold
	real(kind=8), dimension(:)       :: ft 
        integer                                       :: i,k,ntot,isign

        ntot = ubound(ft,1)
   
!	-------------------------------------------------------
!	-------------------------------------------------------
!	pack the array to be fed to the subroutine.  The real parts of
!	the function must be in odd array boxes, the imaginary in even boxes
!	also they must go in increasing order of the independent variable
!	i.e. for f(x) must go in increasing value of x
	k=1
	do i=1,ntot
	  split(k)=ft(i)
	  split(k+1)=0.d0
	  k=k+2
	enddo
!	---------------------------------------------------------------
!	actually do the subroutine call
	call four1(split,isign)
!	---------------------------------------------------------------
!	unpack split  Split is originally presented in 
!	order of increasing + frequencies, and then decreasing negative
!	frequencies.  The multiplication of 1/n is to renormalize.  Only
!	needed for inverse f.t.
        k=1
        do i=ntot+1, 2*ntot
          hold(k)=split(i)
          k=k+1
        enddo
      
        do i=1,ntot
          hold(k)=split(i)
          k=k+1
        enddo
!	-------------------------------------------------------------
      
      end subroutine FFTET
!     ----------------------------------------------------------------------
      SUBROUTINE four1(DATA,ISIGN)
	implicit none
!	this program performs the fourier transform of data
!	nn is the number of complex data points.  nn must be a 
!	power of two.  this is not checked for!!  isign controls
!	whether this is a forward (1) or reverse (-1) transform
!	data must be packed with the real portions of the function
!	in odd array elements and the imaginary portion in the even
!	array elements
      integer       ::  i,n,nn,j,mmax,m,istep,isign
      real(kind=8)  ::  WR,WI,WPR,WPI,WTEMP,THETA,tempr,tempi
      real(kind=8), dimension(:)  :: data

      n = ubound(data,1)


      J=1
      DO 11 I=1,N,2
        IF(J.GT.I)THEN
          TEMPR=DATA(J)
          TEMPI=DATA(J+1)
          DATA(J)=DATA(I)
          DATA(J+1)=DATA(I+1)
          DATA(I)=TEMPR
          DATA(I+1)=TEMPI
        ENDIF
        M=N/2
1       IF ((M.GE.2).AND.(J.GT.M)) THEN
          J=J-M
          M=M/2
        GO TO 1
        ENDIF
        J=J+M
11    CONTINUE
      MMAX=2
2     IF (N.GT.MMAX) THEN
        ISTEP=2*MMAX
        THETA=6.28318530717959D0/(ISIGN*MMAX)
        WPR=-2.D0*DSIN(0.5D0*THETA)**2
        WPI=DSIN(THETA)
        WR=1.D0
        WI=0.D0
        DO 13 M=1,MMAX,2
          DO 12 I=M,N,ISTEP
            J=I+MMAX
            TEMPR=(WR)*DATA(J)-(WI)*DATA(J+1)
            TEMPI=(WR)*DATA(J+1)+(WI)*DATA(J)
            DATA(J)=DATA(I)-TEMPR
            DATA(J+1)=DATA(I+1)-TEMPI
            DATA(I)=DATA(I)+TEMPR
            DATA(I+1)=DATA(I+1)+TEMPI
12        CONTINUE
          WTEMP=WR
          WR=WR*WPR-WI*WPI+WR
          WI=WI*WPR+WTEMP*WPI+WI
13      CONTINUE
        MMAX=ISTEP
      GO TO 2
      ENDIF

      end subroutine four1 
 
     end module fourier_mod
