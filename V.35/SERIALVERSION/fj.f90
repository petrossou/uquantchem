FUNCTION fj(J,L,M,A,B) 
      ! The function f_j(l,m,a,b) on page 237 in the Cook book
      IMPLICIT NONE
      DOUBLE PRECISION :: fj
      DOUBLE PRECISION :: A,B
      INTEGER :: J,L,M
      INTEGER :: K
      DOUBLE PRECISION, EXTERNAL :: binomfac
      fj = 0.0d0
      DO K=MAX(0,J-M),MIN(J,L)
        fj = fj + binomfac(L,K)*binomfac(M,J-K)*(A**(L-K))*(B**(M+K-J))
      ENDDO
END FUNCTION fj
