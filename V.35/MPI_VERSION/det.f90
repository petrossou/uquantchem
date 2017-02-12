FUNCTION det(MATRIX,N)
        ! This function calculates the determinant of a NxN matrix A,
        ! by first transforming it into upper diagonal form
        IMPLICIT NONE
        DOUBLE PRECISION :: DET
        INTEGER :: N
        DOUBLE PRECISION :: MATRIX(N,N)
        INTEGER :: I,J,K
        DOUBLE PRECISION :: A(N,N),B(N,N),C(N,N),TECKEN

        A = MATRIX

        TECKEN = 1.0d0
        DET = 1.0d0
        IF ( N .GT. 1 ) THEN
                DO I=1,N-1
                        DO J=I+1,N
                                IF ( A(I,I) .EQ. 0.0d0 .AND. J .EQ. I+1 ) THEN
                                        K = I
                                        B(I,:) = A(I,:)
                                        DO WHILE ( A(K,I) .EQ. 0.0d0 .AND. K .LT. N )
                                                K=K+1
                                                C(K,:) = A(K,:)
                                        ENDDO
                                        A(I,:) = C(K,:)
                                        A(K,:) = B(I,:)
                                        TECKEN = -1.0d0*TECKEN
                                ENDIF
                                IF ( A(I,I) .NE. 0.0d0  ) THEN
                                        A(J,:) = A(J,:) - (A(J,I)/A(I,I))*A(I,:)
                                ELSE
                                        DET = 0.0d0
                                ENDIF
                        ENDDO
                        IF ( I .EQ. 1 ) THEN 
                                DET = A(1,1)
                        ELSE
                                DET = DET*A(I,I)
                        ENDIF
                ENDDO
                DET = DET*A(N,N)*TECKEN
        ELSE
                DET = MATRIX(1,1)
        ENDIF
END FUNCTION DET
