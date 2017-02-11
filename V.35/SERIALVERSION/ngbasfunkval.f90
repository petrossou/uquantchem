SUBROUTINE ngbasfunkval(PSI,R,nucgradient)
      ! Calculates the value of the gradient of a real basis-function at the spatial
      ! point R, with respect to the nuclear coordinates at which the basis
      ! function is centered.
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(IN) :: R(3)
      DOUBLE PRECISION, INTENT(OUT) :: nucgradient(3)
      DOUBLE PRECISION :: grad(3)

      CALL gradbasfunkval(PSI,R,grad)
        
      nucgradient = -grad

END SUBROUTINE ngbasfunkval

