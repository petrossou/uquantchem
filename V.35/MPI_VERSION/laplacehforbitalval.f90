FUNCTION laplacehforbitalval(BAS,C,N,R)
        ! Calculates the value of the laplacian of N:th Hartree-Fock orbital PHI_N(R)
        ! at the spatial point R.  C(:,N) = the expansion coefficients of the 
        ! N:th HF-orbital. See QDMC-notes p.12 Eqn. (34.5)
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: laplacehforbitalval
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: C(BAS%NBAS,BAS%NBAS) 
        INTEGER :: N
        DOUBLE PRECISION :: R(3)
        INTEGER :: I,L1,L2,L3
        DOUBLE PRECISION :: gradbas(3),lapx,lapy,lapz,LAMDA,STONORM,RR(3),RB,BETA,GAMA,factor1,factor2,factor3,factor4,P(6)
        DOUBLE PRECISION, EXTERNAL :: laplacebasfunkval
        
        laplacehforbitalval = 0.0d0

        DO I=1,BAS%NBAS
                RR = R-BAS%PSI(I)%R
                RB = sqrt(DOT_PRODUCT(RR,RR))
                IF ( RB .GT. BAS%PSI(I)%RC + BAS%PSI(I)%DR .OR. .not.  BAS%PSI(I)%HYBRID ) THEN
                        laplacehforbitalval  = laplacehforbitalval + C(I,N)*laplacebasfunkval(BAS%PSI(I),R)
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
                           lapx = ( L1*(L1-1)*RR(1)**(L1-2) - (LAMDA/RB)*(2*L1 + 1)*RR(1)**L1 + ( (LAMDA/RB)**2 + LAMDA/(RB**3) )*RR(1)**(L1+2)  )*(RR(2)**L2)*(RR(3)**L3)*EXP(-LAMDA*RB) 
                           lapy = ( L2*(L2-1)*RR(2)**(L2-2) - (LAMDA/RB)*(2*L2 + 1)*RR(2)**L2 + ( (LAMDA/RB)**2 + LAMDA/(RB**3) )*RR(2)**(L2+2)  )*(RR(1)**L1)*(RR(3)**L3)*EXP(-LAMDA*RB) 
                           lapz = ( L3*(L3-1)*RR(3)**(L3-2) - (LAMDA/RB)*(2*L3 + 1)*RR(3)**L3 + ( (LAMDA/RB)**2 + LAMDA/(RB**3) )*RR(3)**(L3+2)  )*(RR(1)**L1)*(RR(2)**L2)*EXP(-LAMDA*RB)
                           laplacehforbitalval = laplacehforbitalval + C(I,N)*BAS%PSI(I)%NORM*STONORM*(lapx + lapy + lapz )
                           
                           lapx = ( RB*L1*(L1-1)*RR(1)**(L1-2) + (1.0d0/RB)*(2*L1 + 1)*RR(1)**L1 - ( 1.0d0/(RB**3) )*RR(1)**(L1+2)  )*(RR(2)**L2)*(RR(3)**L3)
                           lapy = ( RB*L2*(L2-1)*RR(2)**(L2-2) + (1.0d0/RB)*(2*L2 + 1)*RR(2)**L2 - ( 1.0d0/(RB**3) )*RR(2)**(L2+2)  )*(RR(1)**L1)*(RR(3)**L3)
                           lapz = ( RB*L3*(L3-1)*RR(3)**(L3-2) + (1.0d0/RB)*(2*L3 + 1)*RR(3)**L3 - ( 1.0d0/(RB**3) )*RR(3)**(L3+2)  )*(RR(1)**L1)*(RR(2)**L2)
                           laplacehforbitalval = laplacehforbitalval + C(I,N)*BAS%PSI(I)%NORM*BETA*(lapx + lapy + lapz )

                           lapx = L1*(L1-1)*( RR(1)**(L1-2) )*( RR(2)**L2 )*( RR(3)**L3 )
                           lapy = L2*(L2-1)*( RR(1)**L1 )*( RR(2)**(L2-2) )*( RR(3)**L3 )
                           lapz = L3*(L3-1)*( RR(1)**L1 )*( RR(2)**L2 )*( RR(3)**(L3-2) )
                           laplacehforbitalval = laplacehforbitalval + C(I,N)*BAS%PSI(I)%NORM*GAMA*(lapx + lapy + lapz )
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
                        factor4 = (-P(2)/(RB**3) + (1.0d0/(2.0d0*RB))*P(4) + (2.0d0/3.0d0)*P(5) + (3.0d0/4.0d0)*P(6)*RB )
                    
                        lapx = factor1*L1*(L1-1)*(RR(1)**(L1-2))*(RR(2)**L2)*(RR(3)**L3)
                        lapx = lapx + factor3*(2*L1+1)*factor2 
                        lapx = lapx + factor4*(RR(1)**(L1+2))*(RR(2)**L2)*(RR(3)**L3)

                        lapy = factor1*L2*(L2-1)*(RR(1)**L1)*(RR(2)**(L2-2))*(RR(3)**L3)
                        lapy = lapy + factor3*(2*L2+1)*factor2 
                        lapy = lapy + factor4*(RR(1)**L1)*(RR(2)**(L2+2))*(RR(3)**L3)

                        lapz = factor1*L3*(L3-1)*(RR(1)**L1)*(RR(2)**L2)*(RR(3)**(L3-2))
                        lapz = lapz + factor3*(2*L3+1)*factor2 
                        lapz = lapz + factor4*(RR(1)**L1)*(RR(2)**L2)*(RR(3)**(L3+2))
                        
                        laplacehforbitalval = laplacehforbitalval + C(I,N)*BAS%PSI(I)%NORM*(lapx + lapy + lapz )
                ENDIF
        ENDDO
END FUNCTION laplacehforbitalval


