SUBROUTINE diaghcc( H,N,EIGENVAL,EIGENVECT,info)
          IMPLICIT NONE
          EXTERNAL :: zheev
          INTEGER, INTENT(IN) :: N
          COMPLEX*16, INTENT(IN) :: H(N,N)
          COMPLEX*16, INTENT(OUT) :: EIGENVAL(N)
          COMPLEX*16, INTENT(OUT) :: EIGENVECT(N,N)
          COMPLEX*16  :: A(N,N),VL(N,N)
          INTEGER, INTENT(OUT) :: info
          INTEGER :: lwork,M
          COMPLEX*16, ALLOCATABLE :: work(:)
          DOUBLE PRECISION, ALLOCATABLE :: rwork(:)
          
          A = H
          if( size(H,1)/=size(H,2) ) then
             print*,'H is not a square'
             stop
          endif
          
          lwork = 2*N
          M = 2*N
          allocate( work(lwork),rwork(M) )
          
          work = (0.0d0,0.0d0)
          rwork = 0.0d0

          !call dsyev( 'V', 'U', N, EIGENVECT, N, EIGENVAL, work, lwork, info )
          !call zheev( 'V', 'U', N, EIGENVECT, N, EIGENVAL, work, lwork, rwork, info )
          
          !call ZGEEV( 'N', 'V', N, A, N, EIGENVAL, EIGENVECT, N, work, lwork, rwork, info )
          call ZGEEV( 'N', 'V', N, A, N, EIGENVAL, VL, N, EIGENVECT, N, work, lwork, rwork, info )
          !if( info /= 0 ) then
          !   print*,'Something wrong in diagh',info
          !   stop
          !endif
          
          deallocate( work,rwork )

END SUBROUTINE diaghcc
