  program fpl
! ----------------------------------------------------------------SECTION A
! declare all the variables that will be used.
  implicit none
  integer, parameter                 :: npart=6,nint=5  !number of atoms and number internal dof
  integer                            :: ncycle,nsamp,istart=0
  integer                            :: i,j             !counter
  real*8                             :: fii,req,delt,dt,etot,time_tot,vpot
  real*8,  dimension(npart)          :: x,xm,xx,v,f
  real*8,  dimension(nint)           :: r,fr
  real*8,  dimension(npart,nint)     :: drx
  real*8,  dimension(npart,npart)    :: umatrix          !used for normal modes
  real*8,  dimension(npart)          :: root,bk,q,pq,enm    !used for normal modes

! ----------------------------------------------------------------SECTION A'
! this file is new.  You can you the output file init.xyz to make a jmol movie.
  open(25,file='init.xyz')   !jmol movie file
! ----------------------------------------------------------------SECTION B

! read in the date file
  open(15,file='fpl_simp.inp')
  read(15,*)time_tot  !total integration time
  read(15,*)dt        !time step
  read(15,*)delt      !sampling time
  read(15,*)fii       !force constant
  read(15,*)req       !equilibrium bond length
  close (15)
! ----------------------------------------------------------------SECTION C
! set up derivatives of r wrt x
  drx = 0
  do i = 1,npart-1
   drx(i,i) = -1
   drx(i+1,i) = 1
  enddo
! ----------------------------------------------------------------SECTION D
  ncycle = time_tot/dt   !number of integration steps
  nsamp =  delt/dt       !sampling rate
! ----------------------------------------------------------------SECTION E
! calculate the normal modes of the linear chain
! root - are the eigenvalues
! umatrix - hessian on input, returned as transformation between normal modes and cartesian on return
  call normal_modes(npart,fii,umatrix,root,bk)
! ---------------------------------------------------------------SECTION F
  call init(npart,x,xm,v,dt,req) !set up an initial condition
! ---------------------------------------------------------------SECTION G
  do i = 1,ncycle
   call force(npart,nint,x,f,r,fr,drx,fii,vpot,req) 
   call integrate(npart,dt,f,x,xx,xm,v ,etot,vpot)   !integrate equations of motion
   if(mod(i,10)==1)call nm_time(npart,req,x,v,q,pq,i,dt,root,bk,enm,umatrix,istart)
   if(mod(i,ncycle/40)==1)call jmol(npart,x)    !make a movie
  enddo
    
  end program fpl


  subroutine jmol(npart,x)
  implicit none
  integer                          :: npart,i
  character(len=4)                 :: atom=' C  '
  real(kind=8), dimension(npart)   :: x
  write(25,*)npart
  write(25,*)
  do i = 1,npart
    write(25,'(a4,3f10.3)')atom,x(i),0.d0,0.d0
  enddo
  end subroutine jmol

! ----------------------------------------------------------------
  subroutine init(npart,x,xm,v,dt,req) !set up an initial condition
  implicit none
  real(kind=8), dimension(npart) :: x,xm,v
  integer                        :: i,npart
  real(kind=8)                   :: req,dt
  v = 0
  x(1) = 0.d0
  do i = 2,npart
   x(i) = x(i-1)+req
  enddo
  x(1)=x(1)-1.d0                   !stretch one bond
  xm = x - v*dt
! ------------------------------------------------------
  end subroutine init


!  ---------------------------------------------------------------------------------
!  Subroutine to calculate the forces.
!  Variables are: x cartesion coordinates; f forces wrt x; 
!                 r internals coordinates; f forces wrt r; 
!                 req equilbrium bond length; vpot potential energy
!  ------------------------------------------------------------------------------
  subroutine force(npart,nint,x,f,r,fr,drx,fii,vpot,req)
  implicit none
  real(kind=8), dimension(npart)      :: x,f
  real(kind=8), dimension(nint)       :: r,fr
  real(kind=8), dimension(npart,nint) :: drx 
  integer                             :: i,npart,nint   
  real(kind=8)         :: req,fii,vpot
   
  do i = 1,npart-1  !calculate internals
    r(i) = x(i+1)-x(i)-req
  enddo

  fr = -fii*r                        !calculate internal forces
  vpot = 0.5d0*fii*dot_product(r,r)  !calculate potential
  f = matmul(drx,fr)                 !calculate cartesian forces
  
  end subroutine force

!-----------------------------------------------------------------------------------
  subroutine integrate(npart,dt,f,x,xx,xm,v,etot,vpot)   !integrate equations of motion
  implicit none
  real(kind=8), dimension(npart)   :: xx,x,f,xm,v
  real(kind=8)                     ::  dt,dt2, dit, etot, vpot
  integer                          :: npart,mdim
    dt2 = dt*dt
    dit = 1/(2*dt)
    xx = 2*x-xm+dt2*f
    v = (xx-xm)*dit
    xm = x
    x = xx
    etot=vpot+0.5*sum(v*v)   !total energy
  end subroutine integrate

!-----------------------------------------------------------------------------------
  subroutine normal_modes(n,fii,u,root,bk)
  real(kind=8), dimension(n,n)     :: u
  real(kind=8), dimension(n)       :: bk,root
  real(kind=8)                     :: fii 
  integer                          :: n,i
! ---------------------------------------------------------
! set up the hessian
  u = 0
  do i = 1,n-1
   u(i,i) = 2.d0*fii
   u(i+1,i) = -fii
   u(i,i+1) = -fii
  enddo
  u(1,1) = fii
  u(n,n) = fii
! ---------------------------------------------------------
! print out the hessian
  write(*,*)' if the number of atoms is less that 15 the Hessian is now printed.'
  if(n<=15)then
    do i = 1,n
      write(*,'(15f7.3)')u(i,:)
    enddo
  endif
! ---------------------------------------------------------
! find the normal modes by diagonalizing the Hessian
  call house(u,n,root,bk) !on return u is the transformation matrix to normal modes.

! print out the roots
  write(*,*)' here are the eigenvalues of the hessian.'
  write(*,*)root

  write(*,*)' If the number of atoms is less that 15 the normal modes are now printed.'
  write(*,*)' Each column corresponds to a normal mode vector.'
  if(n<=15)then
    do i = 1,n
      write(*,'(15f8.4)')u(i,:)
    enddo
  endif

  end subroutine normal_modes

!-----------------------------------------------------------------------------------
! subroutine converts from cartesian to normal coordinates
  subroutine nm_time(npart,req,x,v,q,pq,ic,dt,root,bk,enm,u,istart)
  implicit none
  real(kind=8), dimension(npart,npart)    :: u
  real(kind=8), dimension(npart)          :: bk,root,x,v,q,pq,enm
  real(kind=8)                            :: dt,req
  integer                                 :: i,ic,istart,npart
  if(istart==0)then
   open(55,file='nm_time.dat')
   write(55,*)' time	energy(1)  energy(2)  energy(3)  energy(4)  energy(5)  energy(6)'
   print*,' root(3) = ',root(:)
   istart=1
  endif
  do i = 1,npart
    bk(i) = x(i) - req*(i-1)
  enddo
  q = matmul(transpose(u),bk)
  pq = matmul(transpose(u),v)
  enm = 0.5d0*(pq*pq+root*q*q)
  write(55,'(7f10.5)')ic*dt,enm(1),enm(2),enm(3),enm(4),enm(5),enm(6)  !write results for 6 modes
  end subroutine nm_time


!
!    A - The (N,N) matrix to be diagonalized.
!    D - Array of dimension N containing eigenvalues on output.
!    
     subroutine house(a,n,d,e)
     implicit none
     real(kind=8), dimension(n,n)  :: a
     real(kind=8), dimension(n)    :: d,e
     real(kind=8)                  :: b,dd,h,f,g,hh,scale
     real(kind=8)                  :: p,r,s,c
     integer                       :: j,k,i,l,n,iter,m

     if(n>=1)then
      do  i=n,2,-1  
          l=i-1
          h=0.d0
          scale=0.d0
          if(l.gt.1)then
            do k=1,l
              scale=scale+dabs(a(i,k))
            enddo
            if(scale.eq.0.d0)then
              e(i)=a(i,l)
            else
              do k=1,l
                A(I,K)=A(I,K)/SCALE
                H=H+A(I,K)**2
              enddo
              F=A(I,L)
              G=-SIGN(dSQRT(H),F)
              E(I)=SCALE*G
              H=H-F*G
              A(I,L)=F-G
              F=0.d0
              DO J=1,L
                A(J,I)=A(I,J)/H
                G=0.d0
                DO  K=1,J
                  G=G+A(J,K)*A(I,K)
                enddo
                IF(L.GT.J)THEN
                  DO K=J+1,L
                    G=G+A(K,J)*A(I,K)
                  enddo
                ENDIF
                E(J)=G/H
                F=F+E(J)*A(I,J)
              enddo
              HH=F/(H+H)
              DO J=1,L
                F=A(I,J)
                G=E(J)-HH*F
                E(J)=G
                DO K=1,J
                  A(J,K)=A(J,K)-F*E(K)-G*A(I,K)
                enddo
              enddo
            ENDIF
          ELSE
            E(I)=A(I,L)
          ENDIF
          D(I)=H
    enddo
   ENDIF
      D(1)=0.d0
      E(1)=0.d0
      DO 23 I=1,N
        L=I-1
        IF(D(I).NE.0.d0)THEN
          DO 21 J=1,L
            G=0.d0
            DO 19 K=1,L
              G=G+A(I,K)*A(K,J)
19          CONTINUE
            DO 20 K=1,L
              A(K,J)=A(K,J)-G*A(K,I)
20          CONTINUE
21        CONTINUE
        ENDIF
        D(I)=A(I,I)
        A(I,I)=1.d0
        IF(L.GE.1)THEN
          DO 22 J=1,L
            A(I,J)=0.d0
            A(J,I)=0.d0
22        CONTINUE
        ENDIF
23    CONTINUE
!
!  Now diagonalize the triadiagonal matrix produced above.
!
      IF (N.GT.1) THEN
        DO I=2,N
          E(I-1)=E(I)
        enddo
        E(N)=0.d0
        DO 45 L=1,N
          ITER=0
1         DO 42 M=L,N-1
            DD=dABS(D(M))+dABS(D(M+1))
            IF (dABS(E(M))+DD==DD) GO TO 2
42        CONTINUE
          M=N
2         IF(M.NE.L)THEN
            IF(ITER.EQ.30)PAUSE 'too many iterations'
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.d0*E(L))
            R=SQRT(G**2+1.d0)
            G=D(M)-D(L)+E(L)/(G+SIGN(R,G))
            S=1.d0
            C=1.d0
            P=0.d0
            DO 44 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(dABS(F).GE.dABS(G))THEN
                C=G/F
                R=dSQRT(C**2+1.d0)
                E(I+1)=F*R
                S=1.d0/R
                C=C*S
              ELSE
                S=F/G
                R=dSQRT(S**2+1.d0)
                E(I+1)=G*R
                C=1.d0/R  
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 43 K=1,N
                F=A(K,I+1)
                A(K,I+1)=S*A(K,I)+C*F
                A(K,I)=C*A(K,I)-S*F
43            CONTINUE
44          CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.d0
            GO TO 1
          ENDIF
45      CONTINUE
      ENDIF
!
!  Now sort the eigenvalues and eigenvectors increasing order.
!
      DO I=1,N-1
        K=I
        P=D(I)
        DO J=I+1,N
          IF(D(J)<=P)THEN
            K=J
            P=D(J)
          ENDIF
        enddo
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO J=1,N
            P=A(J,I)
            A(J,I)=A(J,K)
            A(J,K)=P
          enddo
        ENDIF
      enddo

   end subroutine house

