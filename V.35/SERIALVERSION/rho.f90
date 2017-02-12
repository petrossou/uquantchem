FUNCTION rho(BAS,P,R)
        ! Calculates the charge density from the basis functions 
        ! and the density matrix P. Has to be normalized with 1/sqrt(Ne!)
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: rho
        DOUBLE PRECISION :: R(3)
        INTEGER :: NBAS
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: P(BAS%NBAS,BAS%NBAS) 
        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS)
        INTEGER :: I
        DOUBLE PRECISION, EXTERNAL :: basfunkval

        DO I=1,BAS%NBAS
                V1(I) = basfunkval(BAS%PSI(I),R)
        ENDDO
        V2 = V1
        rho = DOT_PRODUCT(V2,MATMUL(P,V1))
END FUNCTION rho

