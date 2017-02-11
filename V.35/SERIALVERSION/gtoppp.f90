FUNCTION gtoppp(PSI,r)
! Calculates the second radial derivative, at argument r, of a spherically
! symmetric gto basis function consisting of contracted primitive gaussians
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: gtoppp
      TYPE(BASFUNCT) :: PSI
      DOUBLE PRECISION :: r
      INTEGER :: I
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

      
      gtoppp = 0.0d0
      DO I=1,PSI%NPRIM
                gtoppp = PSI%CONTRCOEFF(I)*((2.0d0*PSI%EXPON(I)/pi)**0.750d0)*((2.0d0*PSI%EXPON(I))**2)*r*( 3 - 2.0d0*PSI%EXPON(I)*r**2 )*(PSI%EXPON(I)**(0.50d0*SUM(PSI%L)))*exp(-PSI%EXPON(I)*r**2)
      ENDDO

END FUNCTION gtoppp
