FUNCTION gto(PSI,r)
! Calculates the value at argument for a spherically symmetric gto
! basis function consisting of contracted primitive gaussians
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: gto
      TYPE(BASFUNCT) :: PSI
      DOUBLE PRECISION :: r
      INTEGER :: I
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      gto = 0.0d0

      DO I=1,PSI%NPRIM
                gto = gto + PSI%CONTRCOEFF(I)*((2.0d0*PSI%EXPON(I)/pi)**0.750d0)*( PSI%EXPON(I)**(0.50d0*SUM(PSI%L)) )*exp(-PSI%EXPON(I)*r**2)
      ENDDO

END FUNCTION gto
