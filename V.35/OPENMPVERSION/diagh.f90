SUBROUTINE diagh( H,N,EIGENVAL,EIGENVECT)
          IMPLICIT NONE
          EXTERNAL :: dsyev
          INTEGER, INTENT(IN) :: N
          DOUBLE PRECISION, INTENT(IN) :: H(N,N)
          DOUBLE PRECISION, INTENT(OUT) :: EIGENVAL(N),EIGENVECT(N,N)
          INTEGER :: lwork, info
          DOUBLE PRECISION, ALLOCATABLE :: work(:)
          
          EIGENVECT = H
          if( size(H,1)/=size(H,2) ) then
             print*,'H is not a square'
             stop
          endif
          
          lwork = 3*N - 1
          allocate( work(lwork) )
          
          work = 0.0d0

          call dsyev( 'V', 'U', N, EIGENVECT, N, EIGENVAL, work, lwork, info )
          
          if( info /= 0 ) then
             !print*,'Something wrong in diagh',info
             !stop
          endif
          
          deallocate( work )

END SUBROUTINE diagh
