SUBROUTINE potential(BAS,NATOMS,ATOMS,V,gradV)
      ! This subroutine calculates the potential energy matrix
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primpotential
      INTEGER, INTENT(IN) :: NATOMS
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: V(BAS%NBAS,BAS%NBAS),gradV(NATOMS,3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,M,L1,M1,N1,L2,M2,N2
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(3,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      
      !print*,'V_ij='

      gradV(:,:,:,:) = 0.0d0
      
      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=I,BAS%NBAS
                        
                        V(I,J) = 0.0d0
                        
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
                                        DO M=1,NATOMS
                                                !----------------------------------------------
                                                ! Here we calculate the potential energy matrix
                                                !----------------------------------------------
                                                
                                                V(I,J) = V(I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*primpotential(L1,M1,N1,A,alpha1,ATOMS(M)%R,L2,M2,N2,B,alpha2)
                                                
                                                !---------------------------------------------------------------
                                                ! Here we calculate the gradient  of the potential energy matrix
                                                !---------------------------------------------------------------

                                                CALL gradprimpotential(L1,M1,N1,A,alpha1,ATOMS(M)%R,L2,M2,N2,B,alpha2,gradient)

                                                gradV(BAS%PSI(I)%ATYPE,:,I,J) = gradV(BAS%PSI(I)%ATYPE,:,I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                                gradV(BAS%PSI(J)%ATYPE,:,I,J) = gradV(BAS%PSI(J)%ATYPE,:,I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                                gradV(M,:,I,J)                = gradV(M,:,I,J)                - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(3,:)
                                        ENDDO
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                V(J,I) = V(I,J)
                                DO M=1,NATOMS
                                        gradV(M,:,J,I) = gradV(M,:,I,J) 
                                ENDDO
                        ENDIF
                        !print*,I,J,V(I,J)
                ENDDO

      ENDDO
END SUBROUTINE potential
