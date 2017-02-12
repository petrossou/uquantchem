FUNCTION gtopp(PSI,r)
! Calculates the second radial derivative, at argument r, of a spherically
! symmetric gto basis function consisting of contracted primitive gaussians
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: gtopp
      TYPE(BASFUNCT) :: PSI
      DOUBLE PRECISION :: r
      INTEGER :: I
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

      
      gtopp = 0.0d0
      DO I=1,PSI%NPRIM
                gtopp = gtopp -2.0d0*PSI%EXPON(I)*PSI%CONTRCOEFF(I)*((2.0d0*PSI%EXPON(I)/pi)**0.750d0)*( PSI%EXPON(I)**(0.50d0*SUM(PSI%L)) )*exp(-PSI%EXPON(I)*r**2)
                gtopp = gtopp + ((2.0d0*r*PSI%EXPON(I))**2)*PSI%CONTRCOEFF(I)*((2.0d0*PSI%EXPON(I)/pi)**0.750d0)*( PSI%EXPON(I)**(0.50d0*SUM(PSI%L)) )*exp(-PSI%EXPON(I)*r**2)
      ENDDO

END FUNCTION gtopp
