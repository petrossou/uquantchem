FUNCTION laplacejastrow(I,N,N3,r,b,c)
        ! This function calculates the laplacian of the Jastrow factor,
        ! with respect to the electron coordinate I.
        ! N = total number of electrons
        ! r = coordinates of all electrons, r is a Nx3-array
        ! b,c = jastrow parameters. See Eqn.(31,32), p.11 in my QDMC notes
        IMPLICIT NONE
        DOUBLE PRECISION :: laplacejastrow
        INTEGER :: I,N,N3
        DOUBLE PRECISION :: r(N3),b,c
        DOUBLE PRECISION :: rij,V(3)
        INTEGER :: J

        laplacejastrow = 0.0d0

        IF ( I .LE. ( N - MOD(N,2) )/2 ) THEN
                DO J=1,( N - MOD(N,2) )/2
                        V(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        V(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        V(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        rij = sqrt(DOT_PRODUCT(V,V))
                        IF ( J .NE. I ) laplacejastrow = laplacejastrow + 2*b*(1.0d0/4.0d0)/(rij*(1 + c*rij)**3)
                ENDDO
                DO J=( N - MOD(N,2) )/2+1,N
                        V(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        V(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        V(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        rij = sqrt(DOT_PRODUCT(V,V))
                        IF ( J .NE. I ) laplacejastrow = laplacejastrow + 2*b*(1.0d0/2.0d0)/(rij*(1 + c*rij)**3)
                ENDDO
        ELSE
                DO J=( N - MOD(N,2) )/2+1,N
                        V(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        V(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        V(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        rij = sqrt(DOT_PRODUCT(V,V))
                        IF ( J .NE. I ) laplacejastrow = laplacejastrow + 2*b*(1.0d0/4.0d0)/(rij*(1 + c*rij)**3)
                ENDDO
                DO J=1,( N - MOD(N,2) )/2
                        V(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        V(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        V(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        rij = sqrt(DOT_PRODUCT(V,V))
                        IF ( J .NE. I ) laplacejastrow = laplacejastrow + 2*b*(1.0d0/2.0d0)/(rij*(1 + c*rij)**3)
                ENDDO
        ENDIF

END FUNCTION laplacejastrow
