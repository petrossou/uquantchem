SUBROUTINE graddipoletensorinhom(BAS,NATOMS,d,omega,t,gDPTENSOR)
      ! This subroutine calculates the nuclear gradient of the 
      ! dipole-tensor d_mn = <m| W-field |n>
      ! where W-field = e*sin( omega*t + (c/omega)*r ).
      ! Here:  if d = 1, then r = x, e = y 
      ! Here:  if d = 2, then r = y, e = z 
      ! Here:  if d = 3, then r = z, e = x 
      ! The result is stored in the array: DPTENSOR
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primoverlap
      EXTERNAL :: gradprimoverlap
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: NATOMS,d
      DOUBLE PRECISION :: omega,t
      DOUBLE PRECISION, INTENT(OUT) :: gDPTENSOR(NATOMS,3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,N,K,L1,M1,N1,L2,M2,N2,M,BA,AA
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(2,3),Dpole,gDpole(2,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      

      gDPTENSOR(:,:,:,:) = 0.0d0
      
      DO I=1,BAS%NBAS
                
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
                AA = BAS%PSI(I)%ATYPE

                DO J=1,BAS%NBAS
                        
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R
                        NO2 = BAS%PSI(J)%NORM
                        BA = BAS%PSI(J)%ATYPE

                        DO N=1,BAS%PSI(I)%NPRIM
                                alpha1 = BAS%PSI(I)%EXPON(N)
                                coeff1 = BAS%PSI(I)%CONTRCOEFF(N)
                                NP1 = BAS%PSI(I)%PRIMNORM(N)
                                DO K=1,BAS%PSI(J)%NPRIM
                                        alpha2 = BAS%PSI(J)%EXPON(K)
                                        coeff2 = BAS%PSI(J)%CONTRCOEFF(K)
                                        NP2 = BAS%PSI(J)%PRIMNORM(K)
                                        !------------------------------------------------------
                                        ! Here calculate the grad<psi(N)| W-field |psi(K)> -tensor elements
                                        !-------------------------------------------------------
                                        CALL gradpriminhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,d,OMEGA,t,gDpole)
                                        gDPTENSOR(AA,:,I,J) = gDPTENSOR(AA,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gDpole(1,:)
                                        gDPTENSOR(BA,:,I,J) = gDPTENSOR(BA,:,I,J) + NO1*NO2*coeff1*coeff2*NP1*NP2*gDpole(2,:)
                                ENDDO
                        ENDDO
                ENDDO

      ENDDO
END SUBROUTINE graddipoletensorinhom
