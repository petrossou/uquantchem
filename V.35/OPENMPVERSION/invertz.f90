!============================================================================
! Returns the inverse of a matrix calculated by finding the LU decomposition.  
!============================================================================
FUNCTION invertz(A)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: A
  DOUBLE PRECISION, DIMENSION(SIZE(A,1),SIZE(A,2)) :: invertz
  DOUBLE PRECISION, DIMENSION(SIZE(A,1)) :: WORK  ! work array for LAPACK
  INTEGER, DIMENSION(SIZE(A,1)) :: IPIV   ! pivot indices
  INTEGER :: N, INFO
 
  ! Remeber on needs to use an interface to use array valued functions
 
  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  invertz = A
  N = SIZE(A,1)
  
  ! DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  call DGETRF(N, N, invertz, N, IPIV, INFO)

  IF (INFO .NE. 0) THEN
     STOP 'Matrix is numerically singular!'
  ENDIF

  ! DGETRI computes the inverse of a matrix using the LU factorization
  ! computed by DGETRF.
  call DGETRI(N, invertz, N, IPIV, WORK, N, INFO)

  IF (INFO .NE. 0) THEN
     STOP 'Matrix inversion failed!'
  ENDIF
END FUNCTION invertz
