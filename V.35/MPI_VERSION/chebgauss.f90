SUBROUTINE chebgauss(CGORDER,CGQ)
      ! This subroutine calculates the 
      ! Abscissae and weights of a 
      ! CGORDER Chebyshev-Gauss quadrature.
      ! The only imput is the integer CGORDER
      ! which secifies the order of the quadrature.
      ! The output is organized as follows:
      ! CGQ(I,1) = Abscissae
      ! CGQ(I,2) = Weight.
      ! I = index of weights and abscissae.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: CGORDER
      DOUBLE PRECISION, INTENT(OUT) :: CGQ(CGORDER,2)
      INTEGER :: I
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0

      DO I=1,CGORDER
                CGQ(I,1) = cos(PI*I/(CGORDER*1.0d0 + 1.0d0))
                CGQ(I,2) = (PI/(CGORDER*1.0d0 + 1.0d0))*( 1.0d0 - CGQ(I,1)**2 )
      ENDDO
END SUBROUTINE chebgauss
