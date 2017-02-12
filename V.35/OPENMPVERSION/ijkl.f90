FUNCTION ijkl(I,J,K,L)
      ! This function maps the four indexes i,j,k,l onto the symmetry reduced
      ! intex P = ijkl, ijkl:(i,j,k,l) ---> P. This function can then be used
      ! together with the one dimensional array Intsv(:) produced as output to
      ! the call of the subroutine  eeints. Thus any tensor element
      ! Ints(I,J,K,L) can be accessed through the composite "function"
      ! Ints(I,J,K,L) = Intsv(ijkl(I,J,K,L))
      IMPLICIT NONE
      INTEGER :: ijkl
      INTEGER :: I,J,K,L
      INTEGER :: M,N

      IF ( I .LE. J ) THEN
              M = (J*(J-1))/2 + I
      ELSE
              M = (I*(I-1))/2 + J
      ENDIF

      IF ( K .LE. L ) THEN
              N = (L*(L-1))/2 + K
      ELSE
              N = (K*(K-1))/2 + L
      ENDIF

     IF ( M .LE. N ) THEN
        ijkl = (N*(N-1))/2 + M
     ELSE
        ijkl = (M*(M-1))/2 + N
     ENDIF
END FUNCTION ijkl
