FUNCTION fac(N)
      ! simple factorial function: n!
      IMPLICIT NONE
      DOUBLE PRECISION :: fac
      INTEGER :: N,I
      fac = 1.0d0
      IF ( N .GT. 0 ) THEN
        DO I=1,N
                fac = fac*I
        ENDDO
      ENDIF
END FUNCTION fac
