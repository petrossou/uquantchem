RECURSIVE FUNCTION pvoronoi(I,NATOMS,ATOMS,r)
      ! This function is given by Eqn. (13)
      ! in A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: pvoronoi
      INTEGER :: I,NATOMS
      TYPE(ATOM) :: ATOMS(NATOMS)
      DOUBLE PRECISION :: r(3)
      DOUBLE PRECISION, EXTERNAL :: sofmu,atomicradii
      DOUBLE PRECISION :: pp(NATOMS),mu,nu,u12,a12,x,r1,r2,R12,sumpp
      INTEGER :: J,kk
      
      DO kk=1,NATOMS
        pp(kk) = 1.0d0
        DO J=1,NATOMS
                IF ( J .NE. kk ) THEN
                        r1 = sqrt(DOT_PRODUCT(r-ATOMS(kk)%R,r-ATOMS(kk)%R))
                        r2 = sqrt(DOT_PRODUCT(r-ATOMS(J)%R,r-ATOMS(J)%R))
                        R12 =sqrt(DOT_PRODUCT(ATOMS(kk)%R-ATOMS(J)%R,ATOMS(kk)%R-ATOMS(J)%R))
                        mu = (r1-r2)/R12
                        ! Here we adjust for the atomic size according to the appendix 
                        ! in A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
                        x = atomicradii(ATOMS(kk)%Z)/atomicradii(ATOMS(J)%Z)    ! (A4)
                        u12 = (x - 1.0d0)/(x + 1.0d0 )                          ! (A6)
                        a12 = u12/( (u12**2) - 1.0d0 )                          ! (A5)
                        IF ( dabs(a12) .GT. 0.50d0 ) a12 = 0.50d0*a12/dabs(a12) ! (A3)
                        nu = mu + a12*(1.0d0 - mu**2 )                          ! (A2)
                        pp(kk) = pp(kk)*sofmu(nu)                                         ! (A1)
                ENDIF
        ENDDO
      ENDDO
      
      sumpp = SUM(pp)
      IF ( sumpp .NE. 0.0d0 ) THEN
                pvoronoi = pp(I)/SUM(pp)
      ELSE
                pvoronoi = 0.0d0
      ENDIF
END FUNCTION pvoronoi
