SUBROUTINE kinetic(NATOMS,BAS,T,gradT,NSI,CFORCE)
      ! This subroutine calculates the kinetic energy matrix
      USE datatypemodule
      IMPLICIT NONE
      INTEGER :: NATOMS,NSI
      LOGICAL, INTENT(IN) :: CFORCE
      DOUBLE PRECISION, EXTERNAL :: primkinetic
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(OUT) :: T(BAS%NBAS,BAS%NBAS),gradT(NATOMS,3,NSI,NSI)
      INTEGER :: I,J,N,K,L1,M1,N1,L2,M2,N2,M
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(2,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      
      !print*,'T_ij='
      gradT(:,:,:,:) = 0.0d0

      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=I,BAS%NBAS
                        
                        T(I,J) = 0.0d0
                        
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

                                        !-------------------------------------
                                        ! Calculation of Kinetic energy matrix
                                        !-------------------------------------

                                        T(I,J) = T(I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*primkinetic(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)

                                        !-------------------------------------------------
                                        ! Calculation of gradient of Kinetic energy matrix
                                        !--------------------------------------------------
                                        IF ( CFORCE ) THEN
                                         CALL gradprimkinetic(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,gradient)

                                          gradT(BAS%PSI(I)%ATYPE,:,I,J) = gradT(BAS%PSI(I)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                          gradT(BAS%PSI(J)%ATYPE,:,I,J) = gradT(BAS%PSI(J)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                        ENDIF
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                T(J,I) = T(I,J)
                                IF ( CFORCE ) THEN
                                  DO M=1,NATOMS
                                        gradT(M,:,J,I) = gradT(M,:,I,J)
                                  ENDDO
                                ENDIF
                                !gradT(BAS%PSI(I)%ATYPE,:,J,I) = gradT(BAS%PSI(I)%ATYPE,:,I,J)
                                !gradT(BAS%PSI(J)%ATYPE,:,J,I) = gradT(BAS%PSI(J)%ATYPE,:,I,J)
                        ENDIF
                        !print*,I,J,T(I,J)
                ENDDO

      ENDDO
END SUBROUTINE kinetic
