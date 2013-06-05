SUBROUTINE gradpvoronoi(I,NATOMS,ATOMS,r,grad)
      ! This subroutine calculates all the nuclear gradients of the Becke weight 
      ! function centered at atom I. See Appendix B in J. Chem. Phys. 98, 5612 (1993)
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION,INTENT(OUT) :: grad(NATOMS,3)
      INTEGER, INTENT(IN) :: I,NATOMS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: r(3)
      DOUBLE PRECISION, EXTERNAL :: sofmu,atomicradii,pnonorm,tmu
      DOUBLE PRECISION :: p(NATOMS),mu,nu,u12,a12,x,r1,r2,R12,sumpp
      DOUBLE PRECISION :: gradp(NATOMS,NATOMS,3),Z,gnu(3),gmu(3)
      INTEGER :: J,kk,LL,B,I1,I2

      gradp = 0.0d0
      DO J=1,NATOMS
        p(J) = pnonorm(J,NATOMS,ATOMS,r)
      ENDDO

      Z = SUM(p)
      
      DO kk=1,NATOMS                    ! loop over nuclear gradients
                DO J=1,NATOMS           ! loop over non-normalized Becke-weight fuctions centered at J
                        IF ( kk .EQ. J .OR. kk .EQ. I ) THEN
                                I1 = 1
                                I2 = NATOMS
                        ENDIF
                        IF ( kk .NE. J .AND. kk .NE. I ) THEN
                                I1 = kk
                                I2 = kk
                        ENDIF
                        DO B=I1,I2      ! loop over atomic cites not equal to J
                                IF ( B .NE. J ) THEN
                                        r1 = sqrt(DOT_PRODUCT(ATOMS(J)%R-r,ATOMS(J)%R-r))
                                        r2 = sqrt(DOT_PRODUCT(ATOMS(B)%R-r,ATOMS(B)%R-r))
                                        R12 =sqrt(DOT_PRODUCT(ATOMS(J)%R-ATOMS(B)%R,ATOMS(J)%R-ATOMS(B)%R))
                                        mu = (r1-r2)/R12
                                        ! Here we adjust for the atomic size according to the appendix 
                                        ! in A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
                                
                                        x = atomicradii(ATOMS(J)%Z)/atomicradii(ATOMS(B)%Z)     ! (A4)
                                        u12 = (x - 1.0d0)/(x + 1.0d0 )                          ! (A6)
                                        a12 = u12/( (u12**2) - 1.0d0 )                          ! (A5)
                                        IF ( dabs(a12) .GT. 0.50d0 ) a12 = 0.50d0*a12/dabs(a12) ! (A3)
                                        nu = mu + a12*(1.0d0 - mu**2 )                          ! (A2)
                                        
                                        gmu = 0.0d0
                                        ! Calculating the gradient of mu
                                        IF ( kk .EQ. J .AND. J .EQ. I ) THEN
                                                IF ( r2 .NE. 0.0d0 ) gmu = (1.0d0/R12)*(ATOMS(B)%R-r)/r2 
                                                gmu = gmu - (r1-r2)*(1.0d0/R12**3)*(ATOMS(J)%R-ATOMS(B)%R)
                                        ENDIF
                                        IF ( kk .EQ. J .AND. J .NE. I ) THEN
                                                IF ( r1 .NE. 0.0d0 ) gmu = (1.0d0/R12)*(ATOMS(J)%R-r)/r1 
                                                gmu = gmu - (r1-r2)*(1.0d0/R12**3)*(ATOMS(J)%R-ATOMS(B)%R)
                                        ENDIF
                                        IF ( kk .NE. J .AND. J .NE. I ) THEN
                                                IF ( r2 .NE. 0.0d0 ) gmu = -(1.0d0/R12)*(ATOMS(B)%R-r)/r2 
                                                gmu = gmu + (r1-r2)*(1.0d0/R12**3)*(ATOMS(J)%R-ATOMS(B)%R)
                                        ENDIF
                                        IF ( kk .NE. J .AND. J .EQ. I ) THEN
                                                IF ( r2 .NE. 0.0d0 )  gmu = -(1.0d0/R12)*(ATOMS(B)%R-r)/r2 
                                                gmu = gmu + (r1-r2)*(1.0d0/R12**3)*(ATOMS(J)%R-ATOMS(B)%R)
                                        ENDIF
                                        IF ( kk .EQ. I .AND. J .NE. I ) THEN
                                                IF ( r2 .NE. 0.0d0 )  gmu = gmu +(1.0d0/R12)*(ATOMS(B)%R-r)/r2
                                                IF ( r1 .NE. 0.0d0 )  gmu = gmu -(1.0d0/R12)*(ATOMS(J)%R-r)/r1
                                        ENDIF
                                        ! Transforming the gradient og mu to the  gradient of nu
                                        gnu = (1.0d0 - 2.0d0*a12*mu)*gmu
                                        
                                        gradp(kk,J,:) = gradp(kk,J,:) + p(J)*tmu(nu)*gnu
                                
                                ENDIF
                        ENDDO
                ENDDO
      ENDDO

      ! Finally calculating all the nuclear gradients of the Becke weight
      ! funcion centered at atom I:
      
      grad = 0.0d0
      DO kk=1,NATOMS
                grad(kk,:) = gradp(kk,I,:)/Z
                DO J=1,NATOMS
                        grad(kk,:) = grad(kk,:) - p(I)*gradp(kk,J,:)/Z**2
                ENDDO
      ENDDO     
      
END SUBROUTINE gradpvoronoi
