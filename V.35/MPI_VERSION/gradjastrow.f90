SUBROUTINE gradjastrow(I,N,N3,r,b,c,grad)
        ! This subroutine calculates the gradient of the Jastrow factor,
        ! with respect to the electron coordinate I.
        ! N = total number of electrons
        ! r = coordinates of all electrons, r is a 3*N-long array
        ! b,c = jastrow parameters. See Eqn.(20,21), p.8 in my QDMC notes
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
        INTEGER, INTENT(IN) :: I,N,N3
        DOUBLE PRECISION, INTENT(IN) :: r(N3),b,c
        DOUBLE PRECISION :: rij(3),leng
        INTEGER :: J

        grad(:) = 0.0d0

        IF ( I .LE. ( N - MOD(N,2) )/2 ) THEN
                DO J=1,( N - MOD(N,2) )/2
                        rij(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        rij(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        rij(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        leng = sqrt(DOT_PRODUCT(rij,rij))
                        IF ( leng .NE. 0.0d0 ) THEN
                                grad = grad + b*(1.0d0/4.0d0)*(rij/leng)/((1 + c*leng)**2)
                        ENDIF
                ENDDO
                DO J=( N - MOD(N,2) )/2+1,N
                        rij(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        rij(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        rij(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        leng = sqrt(DOT_PRODUCT(rij,rij))
                        IF ( leng .NE. 0.0d0 ) THEN
                                grad = grad + b*(1.0d0/2.0d0)*(rij/leng)/((1 + c*leng)**2)
                        ENDIF
                ENDDO
        ELSE
                DO J=( N - MOD(N,2) )/2+1,N
                        rij(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        rij(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        rij(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        leng = sqrt(DOT_PRODUCT(rij,rij))
                        IF ( leng .NE. 0.0d0 ) THEN
                                 grad = grad + b*(1.0d0/4.0d0)*(rij/leng)/((1 + c*leng)**2)
                        ENDIF
                ENDDO
                DO J=1,( N - MOD(N,2) )/2
                        rij(1) = r(3*(I-1) + 1 ) - r(3*(J-1) + 1)
                        rij(2) = r(3*(I-1) + 2 ) - r(3*(J-1) + 2)
                        rij(3) = r(3*(I-1) + 3 ) - r(3*(J-1) + 3)
                        leng = sqrt(DOT_PRODUCT(rij,rij))
                        IF ( leng .NE. 0.0d0 ) THEN 
                                grad = grad + b*(1.0d0/2.0d0)*(rij/leng)/((1 + c*leng)**2)
                        ENDIF
                                       
                ENDDO
        ENDIF

END SUBROUTINE gradjastrow
