FUNCTION dfac(N)
      ! This function calculates the 
      ! double factorial n!!
      IMPLICIT NONE
      DOUBLE PRECISION :: dfac
      INTEGER :: N, I
      
      dfac = 1.0d0
      IF ( N .GT. 0 ) THEN
        I = N
        DO WHILE ( I .GT. 0  )
                dfac = I*dfac
                I = I - 2
        ENDDO
      ENDIF
      IF ( N .LT. -1 ) THEN
              I = -N-2
              DO WHILE ( I .GT. 0  )
                dfac = I*dfac
                I = I-2
              ENDDO
              dfac = (-1)**MOD((-N-1)/2,2)/dfac
      ENDIF

END FUNCTION dfac
