SUBROUTINE diaghc( H,N,EIGENVAL,EIGENVECT,info)
          IMPLICIT NONE
          EXTERNAL :: zheev
          INTEGER, INTENT(IN) :: N
          COMPLEX*16, INTENT(IN) :: H(N,N)
          DOUBLE PRECISION, INTENT(OUT) :: EIGENVAL(N)
          COMPLEX*16, INTENT(OUT) :: EIGENVECT(N,N)
          INTEGER, INTENT(OUT) :: info
          INTEGER :: lwork,M
          COMPLEX*16, ALLOCATABLE :: work(:)
          DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
          
          EIGENVECT = H
          if( size(H,1)/=size(H,2) ) then
             print*,'H is not a square'
             stop
          endif
          
          lwork = 6*N - 1
          M = 3*N-2
          allocate( work(lwork),rwork(M) )
          
          work = (0.0d0,0.0d0)
          rwork = 0.0d0

          !call dsyev( 'V', 'U', N, EIGENVECT, N, EIGENVAL, work, lwork, info )
          
          call zheev( 'V', 'U', N, EIGENVECT, N, EIGENVAL, work, lwork, rwork, info )
          
          !if( info /= 0 ) then
          !   print*,'Something wrong in diagh',info
          !   stop
          !endif
          
          deallocate( work,rwork )

END SUBROUTINE diaghc
