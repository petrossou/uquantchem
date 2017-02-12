SUBROUTINE gradhforbitalval(BAS,C,N,R,grad)
        ! Calculates the value of the gradient of N:th Hartree-Fock orbital PHI_N(R)
        ! at the spatial point R.  C(:,N) = the expansion coefficients of the 
        ! N:th HF-orbital. See QDMC-notes p.9 Eqn.26
        USE datatypemodule
        IMPLICIT NONE
        TYPE(BASIS), INTENT(IN) :: BAS
        DOUBLE PRECISION, INTENT(IN) :: C(BAS%NBAS,BAS%NBAS) 
        INTEGER, INTENT(IN) :: N
        DOUBLE PRECISION, INTENT(IN) :: R(3)
        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
        INTEGER :: I,L1,L2,L3
        DOUBLE PRECISION :: gradbas(3),RR(3),RB,gx,gy,gz,LAMDA,STONORM,BETA,GAMA,factor1,factor2,factor3,P(6)
        EXTERNAL :: gradbasfunkval
        
        grad(:) = 0.0d0

        DO I=1,BAS%NBAS
                RR = R-BAS%PSI(I)%R
                RB = sqrt(DOT_PRODUCT(RR,RR))
                IF ( RB .GT. BAS%PSI(I)%RC + BAS%PSI(I)%DR .OR. .not.  BAS%PSI(I)%HYBRID ) THEN
                        CALL gradbasfunkval(BAS%PSI(I),R,gradbas)
                        grad  = grad + C(I,N)*gradbas
                ENDIF
                IF ( RB .LE. BAS%PSI(I)%RC .AND. BAS%PSI(I)%HYBRID ) THEN
                        L1 = BAS%PSI(I)%L(1)
                        L2 = BAS%PSI(I)%L(2)
                        L3 = BAS%PSI(I)%L(3)
                        LAMDA = BAS%PSI(I)%STOEXPON
                        STONORM = BAS%PSI(I)%A
                        BETA = BAS%PSI(I)%BETA
                        GAMA = BAS%PSI(I)%GAMA
                        IF ( RB .NE. 0.0d0 ) THEN
                                gx = C(I,N)*STONORM*BAS%PSI(I)%NORM*( L1*RR(1)**(L1-1) - (LAMDA/RB)*RR(1)**(L1+1) )*(RR(2)**L2)*(RR(3)**L3)*EXP(-LAMDA*RB)
                                grad(1) = grad(1) + gx
                                gx = C(I,N)*STONORM*BAS%PSI(I)%NORM*BETA*( L1*RB*RR(1)**(L1-1) + (1/RB)*RR(1)**(L1+1) )*(RR(2)**L2)*(RR(3)**L3)
                                grad(1) = grad(1) + gx
                                gx = C(I,N)*STONORM*BAS%PSI(I)%NORM*GAMA*L1*(RR(1)**(L1-1))*(RR(2)**L2)*(RR(3)**L3)
                                grad(1) = grad(1) + gx

                                gy = C(I,N)*STONORM*BAS%PSI(I)%NORM*( L2*RR(2)**(L2-1) - (LAMDA/RB)*RR(2)**(L2+1) )*(RR(1)**L1)*(RR(3)**L3)*EXP(-LAMDA*RB)
                                grad(2) = grad(2) + gy
                                gy = C(I,N)*STONORM*BAS%PSI(I)%NORM*BETA*( L2*RB*RR(2)**(L2-1) + (1/RB)*RR(2)**(L2+1) )*(RR(1)**L1)*(RR(3)**L3)
                                grad(2) = grad(2) + gy
                                gy = C(I,N)*STONORM*BAS%PSI(I)%NORM*GAMA*L2*(RR(1)**L1)*(RR(2)**(L2-1))*(RR(3)**L3)
                                grad(2) = grad(1) + gy
                                
                                gz = C(I,N)*STONORM*BAS%PSI(I)%NORM*( L3*RR(3)**(L3-1) - (LAMDA/RB)*RR(3)**(L3+1) )*(RR(1)**L1)*(RR(2)**L2)*EXP(-LAMDA*RB)
                                grad(3) = grad(3) + gz
                                gz = C(I,N)*STONORM*BAS%PSI(I)%NORM*BETA*( L3*RB*RR(3)**(L3-1) + (1/RB)*RR(3)**(L3+1) )*(RR(1)**L1)*(RR(2)**L2)
                                grad(3) = grad(3) + gz
                                gz = C(I,N)*STONORM*BAS%PSI(I)%NORM*GAMA*L3*(RR(1)**L1)*(RR(2)**L2)*(RR(3)**(L3-1))
                                grad(3) = grad(3) + gz
                        ENDIF
                ENDIF
                IF ( RB .GT. BAS%PSI(I)%RC .AND. RB .LE. BAS%PSI(I)%RC + BAS%PSI(I)%DR .AND. BAS%PSI(I)%HYBRID ) THEN
                        P = BAS%PSI(I)%P
                        L1 = BAS%PSI(I)%L(1)
                        L2 = BAS%PSI(I)%L(2)
                        L3 = BAS%PSI(I)%L(3)
                        factor1 = ( P(1) + P(2)*RB + (1.0d0/2.0d0)*P(3)*RB**2 + (1.0d0/6.0d0)*P(4)*RB**3 + (1.0d0/12.0d0)*P(5)*RB**4 + (1.0d0/20.0d0)*P(6)*RB**5 )
                        factor2 = ( RR(1)**L1 )*( RR(2)**L2 )*( RR(3)**L3 )
                        factor3 = (P(2)/RB + P(3) + (1.0d0/2.0d0)*P(4)*RB + (1.0d0/3.0d0)*P(5)*RB**2 + (1.0d0/4.0d0)*P(6)*RB**3 )
                        
                        gx = C(I,N)*BAS%PSI(I)%NORM*factor1*L1*(RR(1)**(L1-1))*(RR(2)**L2)*(RR(3)**L3)
                        gx = gx + C(I,N)*BAS%PSI(I)%NORM*factor2*factor3*RR(1)
                        grad(1) = grad(1) + gx

                        gy = C(I,N)*BAS%PSI(I)%NORM*factor1*L2*(RR(1)**L1)*(RR(2)**(L2-1))*(RR(3)**L3)
                        gy = gy + C(I,N)*BAS%PSI(I)%NORM*factor2*factor3*RR(2)
                        grad(2) = grad(2) + gy

                        gz = C(I,N)*BAS%PSI(I)%NORM*factor1*L3*(RR(1)**L1)*(RR(2)**L2)*(RR(3)**(L3-1))
                        gz = gz + C(I,N)*BAS%PSI(I)%NORM*factor2*factor3*RR(3)
                        grad(3) = grad(3) + gz
                ENDIF
        ENDDO
END SUBROUTINE gradhforbitalval


