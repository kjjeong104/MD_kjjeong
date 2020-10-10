! ----------------------------------------------------------------
  subroutine mcmove(x,vpot,istep) !attempts to displace a particle
    use lj_metro_mod
    implicit none
    integer                          :: j,o,istep
    real(kind=8), dimension(3)       :: xn
    real(kind=8)                     :: vpot,eno,enn,rnd
    real(kind=8), dimension(npart,3) :: x

    call random_number(rnd)        !pick particle to be moved
    o=int(rnd*npart)+1
    call potential(o,x,eno)  !calculate the potential energy
   
    xn = x(o,:)                         !current position

    do j = 1,3
      call random_number(rnd)
      x(o,j) = x(o,j) + (rnd-0.5)*delx  !new position
    enddo

    call potential(o,x,enn)  !calculate the potential energy
    call random_number(rnd)
    if(rnd<exp(-beta*(enn-eno)))then  !decide if you want to make a move
           vpot = enn
           istep = istep + 1  !count on accepted moves.
    else
           x(o,:)=xn
    endif
  end subroutine mcmove
