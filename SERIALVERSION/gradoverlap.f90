SUBROUTINE gradoverlap(NATOMS,ATOMS,BAS,SHIFT,gradS)
      ! This subroutine calculates the overlap matrix
      USE datatypemodule
      IMPLICIT NONE
      INTEGER :: NATOMS,SHIFT
      DOUBLE PRECISION, EXTERNAL :: primoverlap
      EXTERNAL :: gradprimoverlap
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,L1,M1,N1,L2,M2,N2,M
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(2,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      
      !print*,'S_ij='

      gradS(:,:,:,:) = 0.0d0
      
      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R-ATOMS(SHIFT)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=I,BAS%NBAS
                        
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R-ATOMS(SHIFT)%R
                        NO2 = BAS%PSI(J)%NORM
                        
                        DO N=1,BAS%PSI(I)%NPRIM
                                alpha1 = BAS%PSI(I)%EXPON(N)
                                coeff1 = BAS%PSI(I)%CONTRCOEFF(N)
                                NP1 = BAS%PSI(I)%PRIMNORM(N)
                                DO K=1,BAS%PSI(J)%NPRIM
                                        alpha2 = BAS%PSI(J)%EXPON(K)
                                        coeff2 = BAS%PSI(J)%CONTRCOEFF(K)
                                        NP2 = BAS%PSI(J)%PRIMNORM(K)
                                        !---------------------------------------------------------------------
                                        ! Here we calculate the gradient of the overlap matrix
                                        ! when all the basis functions have been
                                        ! shifted with the vector ATOMS(SHIFT)%R  = R, i.e PHI(r) --> PHI(r+R)
                                        !---------------------------------------------------------------------
                                        
                                        CALL gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,gradient)
                                        
                                        IF ( BAS%PSI(I)%ATYPE .NE. SHIFT ) THEN
                                                gradS(SHIFT,:,I,J) = gradS(SHIFT,:,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                                gradS(BAS%PSI(I)%ATYPE,:,I,J) = gradS(BAS%PSI(I)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                        ENDIF

                                        IF ( BAS%PSI(J)%ATYPE .NE. SHIFT ) THEN
                                                gradS(SHIFT,:,I,J) = gradS(SHIFT,:,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                                gradS(BAS%PSI(J)%ATYPE,:,I,J) = gradS(BAS%PSI(J)%ATYPE,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                        ENDIF
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                DO M=1,NATOMS
                                        gradS(M,:,J,I) = gradS(M,:,I,J)
                                ENDDO
                                !gradS(BAS%PSI(I)%ATYPE,:,J,I) = gradS(BAS%PSI(I)%ATYPE,:,I,J)
                                !gradS(BAS%PSI(J)%ATYPE,:,J,I) = gradS(BAS%PSI(J)%ATYPE,:,I,J)
                        ENDIF
                        !print*,I,J,S(I,J)
                ENDDO

      ENDDO
END SUBROUTINE gradoverlap
