SUBROUTINE gradprimkinetic(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,gradient)
        ! This function returns the integrals
        ! gradient< prim1 | -0.5*NABLA**2 | prim2 > 
        ! As given by The cook Book on pages 234-239
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(OUT) :: gradient(2,3)
        INTEGER, INTENT(IN) :: L1,M1,N1,L2,M2,N2
        DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,A(3),B(3)
        DOUBLE PRECISION :: grad(2,3),grad1(2,3),grad2(2,3),grad3(2,3)
        EXTERNAL :: gradprimoverlap
        
        gradient = 0.0d0
        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,grad)
        
        gradient = alpha2*( 2*(L2+M2+N2) + 3 )*grad
        
        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2+2,M2,N2,B,alpha2,grad1)
        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2+2,N2,B,alpha2,grad2)
        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+2,B,alpha2,grad3)
        
        gradient = gradient - 2*(alpha2**2)*( grad1 + grad2 + grad3 )

        IF ( L2 .GE. 2 ) THEN
                CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2-2,M2,N2,B,alpha2,grad1) 
                gradient = gradient - 0.50d0*L2*(L2-1)*grad1
        ENDIF
        IF ( M2 .GE. 2 ) THEN
                CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2-2,N2,B,alpha2,grad2)
                gradient = gradient - 0.50d0*M2*(M2-1)*grad2
        ENDIF
        IF ( N2 .GE. 2 ) THEN
                CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-2,B,alpha2,grad3)
                gradient = gradient - 0.50d0*N2*(N2-1)*grad3
        ENDIF
END SUBROUTINE gradprimkinetic
