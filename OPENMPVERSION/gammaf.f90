FUNCTION gammaf(I)
      ! calculates the value og the GAMMA-function G(I+0.5)
      ! gfortan has the gamma function as an intrinsic function 
      ! whereas ifort on osx does not seem to have it.
      IMPLICIT NONE
      DOUBLE PRECISION :: gammaf
      INTEGER :: I 
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, EXTERNAL :: dfac
    
      gammaf = sqrt(pi)*dfac(2*I-1)/2.0d0**I

END FUNCTION  gammaf
