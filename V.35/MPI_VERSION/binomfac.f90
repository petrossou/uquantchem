FUNCTION binomfac(N,M)
      ! calculates the binomial factor
      ! N over M = N!/(M!*(N-M)!)
      IMPLICIT NONE
      DOUBLE PRECISION :: binomfac
      INTEGER :: N,M
      DOUBLE PRECISION, EXTERNAL :: fac

      IF ( M .GT. N ) THEN
              print*,'WARNING M > N IN binomfac.f90 ABORTING!'
              STOP
      ENDIF
      binomfac = fac(N)/(fac(M)*fac(N-M))
      
END FUNCTION binomfac
