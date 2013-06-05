SUBROUTINE diaghHF( H,S,N,EIGENVAL,EIGENVECT)
          IMPLICIT NONE
          EXTERNAL :: DSYGV
          INTEGER, INTENT(IN) :: N
          DOUBLE PRECISION, INTENT(IN) :: H(N,N),S(N,N)
          DOUBLE PRECISION, INTENT(OUT) :: EIGENVAL(N),EIGENVECT(N,N)
          DOUBLE PRECISION :: SS(N,N),worktemp(4)
          INTEGER :: lwork, info
          DOUBLE PRECISION, ALLOCATABLE :: work(:)
          
          
          EIGENVECT = H
          SS = S
          
          if( size(H,1)/=size(H,2) ) then
             print*,'H is not a square'
             stop
          endif
          
          lwork =  - 1
          
          ! Work space query:
          
          call DSYGV( 1, 'V', 'U', N, EIGENVECT, N, SS, N, EIGENVAL, worktemp, lwork, INFO )
          
          lwork = INT(worktemp(1))
          allocate(work(lwork))
          
          ! Calculation of eigenvalues and eigenvectors:
          
          call DSYGV( 1, 'V', 'U', N, EIGENVECT, N, SS, N, EIGENVAL, work, lwork, INFO )
          
          if( info /= 0 ) then
             print*,'Something wrong in diagh',info
             stop
          endif
          
          deallocate( work )
          
END SUBROUTINE diaghHF
