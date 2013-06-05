SUBROUTINE rhoo(BAS,P,R,rho)
        ! Calculates the charge density from the basis functions 
        ! and the density matrix P. Has to be normalized with 1/sqrt(Ne!)
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: rho
        DOUBLE PRECISION, INTENT(IN) :: R(3)
        INTEGER :: NBAS
        TYPE(BASIS),INTENT(IN) :: BAS
        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS) 
        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS)
        INTEGER :: I
        DOUBLE PRECISION :: basfunkval

        DO I=1,BAS%NBAS
                V1(I) = basfunkval(BAS%PSI(I),R)
        ENDDO
        V2 = V1
        rho = DOT_PRODUCT(V2,MATMUL(P,V1))
END SUBROUTINE rhoo

