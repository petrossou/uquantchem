SUBROUTINE normalize(BAS)
      ! This subroutine calculates the normalization constants 
      ! for all of the basis functions stored in BAS
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primoverlap
      TYPE(BASIS), INTENT(INOUT) :: BAS
      INTEGER :: I,J,K,L,M,N
      DOUBLE PRECISION :: NORMC
      DOUBLE PRECISION :: A(3),alpha1,alpha2,coeff1,coeff2,NP1,NP2
        
      DO I=1,BAS%NBAS
                NORMC = 0.0d0
                L= BAS%PSI(I)%L(1)
                M= BAS%PSI(I)%L(2)
                N= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                DO J=1,BAS%PSI(I)%NPRIM
                        alpha1 = BAS%PSI(I)%EXPON(J)
                        coeff1 = BAS%PSI(I)%CONTRCOEFF(J)
                        NP1 = BAS%PSI(I)%PRIMNORM(J)
                        DO K=1,BAS%PSI(I)%NPRIM
                                alpha2 = BAS%PSI(I)%EXPON(K)
                                coeff2 = BAS%PSI(I)%CONTRCOEFF(K)
                                NP2 = BAS%PSI(I)%PRIMNORM(K)
                                NORMC = NORMC + coeff1*coeff2*NP1*NP2*primoverlap(L,M,N,A,alpha1,L,M,N,A,alpha2)
                        ENDDO
                ENDDO

                BAS%PSI(I)%NORM = 1.0d0/SQRT(NORMC)
      ENDDO
END SUBROUTINE normalize
