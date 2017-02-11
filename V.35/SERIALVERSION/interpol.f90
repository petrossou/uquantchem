SUBROUTINE interpol(rc,deltar,PSI,Acoeff,LAMDA,POLY)
      ! This routine calculates the interpolating 
      ! cubic polynomial linking together the 
      ! second derivative of an sto-type basis function:
      ! sto''(1) = Acoeff*(LAMDA**2)*EXP(-LAMDA*r) at the 
      ! point r = rc, and a the second derivative of a 
      ! gto-type basis function: gto''(r) at r = rc + deltar
      ! so that both second and third derivative of the polunomial 
      ! matches tha sto at r=rc and the gto at r=rc+deltar. This is 
      ! done according to J. Chem. Phys. 115, 5362 (2001).
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: rc,deltar,Acoeff,LAMDA
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(OUT) :: POLY(4)
      DOUBLE PRECISION :: A(4,4), B(4)
      INTEGER :: INFO, IPIV(4)
      DOUBLE PRECISION, EXTERNAL :: gtopp,gtoppp
      A(:,:) = 0.0d0
      
      A(1,1) = 1.0d0
      A(1,2) = rc
      A(1,3) = rc**2
      A(1,4) = rc**3

      A(2,1) = 1.0d0
      A(2,2) = rc+deltar
      A(2,3) = (rc+deltar)**2
      A(2,4) = (rc+deltar)**3

      A(3,1) = 0.0d0
      A(3,2) = 1.0d0
      A(3,3) = 2.0d0*rc
      A(3,4) = 3.0d0*rc**2
      
      A(4,1) = 0.0d0
      A(4,2) = 1.0d0
      A(4,3) = 2.0d0*(rc+deltar)
      A(4,4) = 3.0d0*(rc+deltar)**2

      B(1) = Acoeff*(LAMDA**2)*EXP(-LAMDA*rc)
      B(2) = gtopp(PSI,rc+deltar)

      B(3) = -Acoeff*(LAMDA**3)*EXP(-LAMDA*rc)
      B(4) = gtoppp(PSI,rc+deltar)

      CALL DGESV( 4, 1, A, 4, IPIV, B, 4, INFO )
        
      POLY = B
      IF ( INFO .NE. 0 ) THEN
              WRITE(*,*)'FAILURE IN POLYNOMIAL INTERPOLATION'
              WRITE(*,*)'WHEN CONSTRUCTING CUSP CORRECTION'
              STOP
      ENDIF
END SUBROUTINE interpol
