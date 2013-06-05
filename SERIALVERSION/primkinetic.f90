FUNCTION primkinetic(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)
        ! This function returns the integral
        ! < prim1 | -0.5*NABLA**2 | prim2 > 
        ! As given by The cook Book on pages 234-239
        IMPLICIT NONE
        DOUBLE PRECISION :: primkinetic
        INTEGER :: L1,M1,N1,L2,M2,N2
        DOUBLE PRECISION :: alpha1,alpha2,A(3),B(3)
        DOUBLE PRECISION, EXTERNAL :: primoverlap
        
        primkinetic = alpha2*( 2*(L2+M2+N2) + 3 )*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)

        primkinetic = primkinetic - 2*(alpha2**2)*( primoverlap(L1,M1,N1,A,alpha1,L2+2,M2,N2,B,alpha2) )
        primkinetic = primkinetic - 2*(alpha2**2)*( primoverlap(L1,M1,N1,A,alpha1,L2,M2+2,N2,B,alpha2) + primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+2,B,alpha2))

        IF ( L2 .GE. 2 ) primkinetic = primkinetic - 0.50d0*L2*(L2-1)*primoverlap(L1,M1,N1,A,alpha1,L2-2,M2,N2,B,alpha2)
        IF ( M2 .GE. 2 ) primkinetic = primkinetic - 0.50d0*M2*(M2-1)*primoverlap(L1,M1,N1,A,alpha1,L2,M2-2,N2,B,alpha2)
        IF ( N2 .GE. 2 ) primkinetic = primkinetic - 0.50d0*N2*(N2-1)*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-2,B,alpha2)
END FUNCTION primkinetic
