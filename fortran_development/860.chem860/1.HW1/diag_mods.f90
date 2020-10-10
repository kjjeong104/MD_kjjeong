   module diag_mods 
!
!  Reduction of a real symmetric matrix, A, to tridagonal
!  form by the Householder method. This is followed by the
!  evaluation of the eigenvalues and eigenvectors.
!
!    A - The (N,N) matrix to be diagonalized.
!    D - Array of dimension N containing eigenvalues on output.
!    
     contains

     subroutine house(a,d,e)
     implicit none
     real(kind=8), dimension(:,:)  :: a
     real(kind=8), dimension(:)    :: d,e
     real(kind=8)                  :: b,dd,h,f,g,hh,scale
     real(kind=8)                  :: p,r,s,c
     integer                       :: j,k,i,l,n,iter,m
     n = ubound(a,1)

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

   end module diag_mods
