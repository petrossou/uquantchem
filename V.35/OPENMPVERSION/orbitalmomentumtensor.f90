SUBROUTINE orbitalmomentumtensor(BAS,LTENSOR)
      ! This subroutine calculates the orbital momentum 
      ! matrices Lx,Ly,Lz  (see p. 56-57 in dark blue notebook)
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primoverlap
      EXTERNAL :: gradprimoverlap
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(OUT) :: LTENSOR(3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,L1,M1,N1,L2,M2,N2,M
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(2,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      

      LTENSOR(:,:,:) = 0.0d0
      
      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                
                DO J=1,BAS%NBAS
                        
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
                                        ! Here calculate the Lx-tensor
                                        !----------------------------------
                                        IF ( N2 .GE. 1 ) THEN
                                           LTENSOR(1,I,J) = LTENSOR(1,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( N2*primoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2-1,B,alpha2) )
                                           LTENSOR(1,I,J) = LTENSOR(1,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(2)*N2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2) )
                                        ENDIF
                                        LTENSOR(1,I,J) = LTENSOR(1,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(2)*(-2.0d0*alpha2)*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2) )
                                        
                                        IF ( M2 .GE. 1 ) THEN
                                           LTENSOR(1,I,J) = LTENSOR(1,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( M2*primoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2+1,B,alpha2) )
                                           LTENSOR(1,I,J) = LTENSOR(1,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( B(3)*M2*primoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2) )
                                        ENDIF
                                        LTENSOR(1,I,J) = LTENSOR(1,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(3)*2.0d0*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2) )
                                        !----------------------------------
                                        ! Here calculate the Ly-tensor
                                        !----------------------------------
                                        IF ( L2 .GE. 1 ) THEN
                                           LTENSOR(2,I,J) = LTENSOR(2,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( L2*primoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2+1,B,alpha2) )
                                           LTENSOR(2,I,J) = LTENSOR(2,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(3)*L2*primoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2) )
                                        ENDIF
                                        LTENSOR(2,I,J) = LTENSOR(2,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(3)*(-2.0d0*alpha2)*primoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2) )
                                        
                                        IF ( N2 .GE. 1 ) THEN
                                           LTENSOR(2,I,J) = LTENSOR(2,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( N2*primoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2-1,B,alpha2) )
                                           LTENSOR(2,I,J) = LTENSOR(2,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( B(1)*N2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2) )
                                        ENDIF
                                        LTENSOR(2,I,J) = LTENSOR(2,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(1)*2.0d0*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2) )
                                        !----------------------------------
                                        ! Here calculate the Lz-tensor
                                        !----------------------------------
                                        IF ( M2 .GE. 1 ) THEN
                                           LTENSOR(3,I,J) = LTENSOR(3,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( M2*primoverlap(L1,M1,N1,A,alpha1,L2+1,M2-1,N2,B,alpha2) )
                                           LTENSOR(3,I,J) = LTENSOR(3,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(1)*M2*primoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2) )
                                        ENDIF
                                        LTENSOR(3,I,J) = LTENSOR(3,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(1)*(-2.0d0*alpha2)*primoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2) )
                                        
                                        IF ( L2 .GE. 1 ) THEN
                                           LTENSOR(3,I,J) = LTENSOR(3,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( L2*primoverlap(L1,M1,N1,A,alpha1,L2-1,M2+1,N2,B,alpha2) )
                                           LTENSOR(3,I,J) = LTENSOR(3,I,J) - NO1*NO2*coeff1*coeff2*NP1*NP2*( B(2)*L2*primoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2) )
                                        ENDIF
                                        LTENSOR(3,I,J) = LTENSOR(3,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*( B(2)*2.0d0*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2) )
                                ENDDO
                        ENDDO
                ENDDO

      ENDDO
END SUBROUTINE orbitalmomentumtensor
