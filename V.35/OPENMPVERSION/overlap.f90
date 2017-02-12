SUBROUTINE overlap(NATOMS,BAS,S,gradS)
      ! This subroutine calculates the overlap matrix
      USE datatypemodule
      IMPLICIT NONE
      INTEGER :: NATOMS
      DOUBLE PRECISION, EXTERNAL :: primoverlap
      EXTERNAL :: gradprimoverlap
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(OUT) :: S(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,L1,M1,N1,L2,M2,N2,M
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(2,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      
      !print*,'S_ij='

      S(:,:) = 0.0d0
      gradS(:,:,:,:) = 0.0d0
      
      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=I,BAS%NBAS
                        
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R
                        NO2 = BAS%PSI(J)%NORM
                        
                        DO N=1,BAS%PSI(I)%NPRIM
                                alpha1 = BAS%PSI(I)%EXPON(N)
                                coeff1 = BAS%PSI(I)%CONTRCOEFF(N)
                                NP1 = BAS%PSI(I)%PRIMNORM(N)
                                DO K=1,BAS%PSI(J)%NPRIM
                                        alpha2 = BAS%PSI(J)%EXPON(K)
                                        coeff2 = BAS%PSI(J)%CONTRCOEFF(K)
                                        NP2 = BAS%PSI(J)%PRIMNORM(K)
                                        !----------------------------------
                                        ! Here calculate the overlap matrix
                                        !----------------------------------
                                        
                                        S(I,J) = S(I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)
                                        
                                        !-----------------------------------------------------
                                        ! Here we calculate the gradient of the overlap matrix
                                        !-----------------------------------------------------
                                        
                                        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,gradient)
                                        
                                        gradS(BAS%PSI(I)%ATYPE,:,I,J) = gradS(BAS%PSI(I)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                        gradS(BAS%PSI(J)%ATYPE,:,I,J) = gradS(BAS%PSI(J)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                S(J,I) = S(I,J)
                                DO M=1,NATOMS
                                        gradS(M,:,J,I) = gradS(M,:,I,J)
                                ENDDO
                                !gradS(BAS%PSI(I)%ATYPE,:,J,I) = gradS(BAS%PSI(I)%ATYPE,:,I,J)
                                !gradS(BAS%PSI(J)%ATYPE,:,J,I) = gradS(BAS%PSI(J)%ATYPE,:,I,J)
                        ENDIF
                        !print*,I,J,S(I,J)
                ENDDO

      ENDDO
END SUBROUTINE overlap
