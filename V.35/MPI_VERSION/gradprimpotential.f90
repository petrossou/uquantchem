SUBROUTINE gradprimpotential(L1,M1,N1,A,alpha1,C,L2,M2,N2,B,alpha2,gradient)
      ! This function calculates the gradient of the potential energy integral given on the
      ! bottom of p. 244 in the Cook Book.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L1,M1,N1,L2,M2,N2
      DOUBLE PRECISION, INTENT(OUT) :: gradient(3,3)
      DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,A(3),B(3),C(3)
      DOUBLE PRECISION :: P(3),PA(3),PB(3),PC(3)
      DOUBLE PRECISION :: gama, X,TERM0,TERM1,TERM2,TERM3,AB(3)
      DOUBLE PRECISION, ALLOCATABLE :: FNUVEC(:)
      INTEGER :: L,R,I,M,S,J,N,T,K
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, EXTERNAL :: Afunc,dAfunc,primpotential
      EXTERNAL :: Fnurec

      !===========================================================================
      ! First we calculate the gradient with respect to the nuclear coordinate  C.
      ! This is stored in gradient(3,:)
      !===========================================================================

      gama = alpha1+alpha2
      P = (A*alpha1+B*alpha2)/gama
      PA = P-A
      PB = P-B
      PC = P-C
      AB = A-B

      X = gama*DOT_PRODUCT(PC,PC)

      ALLOCATE(FNUVEC(L1+L2+M1+M2+N1+N2+2))
      
      CALL Fnurec(L1+L2+M1+M2+N1+N2+1,FNUVEC,X)
      
      gradient(:,:) = 0.0d0
      
      DO L=0,L1+L2
         DO R=0,(L-MOD(L,2))/2
            DO I=0,(L-2*R-MOD(L-2*R,2))/2
               DO M=0,M1+M2
                  DO S=0,(M-MOD(M,2))/2
                     DO J=0,(M-2*S-MOD(M-2*S,2))/2
                        DO N=0,N1+N2
                           DO T=0,(N-MOD(N,2))/2
                              DO K=0,(N-2*T-MOD(N-2*T,2))/2

                                        TERM0 = Afunc(L,R,I,L1,L2,PA(1),PB(1),PC(1),gama)*Afunc(M,S,J,M1,M2,PA(2),PB(2),PC(2),gama)*Afunc(N,T,K,N1,N2,PA(3),PB(3),PC(3),gama) 
                                        TERM1 = dAfunc(L,R,I,L1,L2,PA(1),PB(1),PC(1),gama)*Afunc(M,S,J,M1,M2,PA(2),PB(2),PC(2),gama)*Afunc(N,T,K,N1,N2,PA(3),PB(3),PC(3),gama) 
                                        TERM2 = Afunc(L,R,I,L1,L2,PA(1),PB(1),PC(1),gama)*dAfunc(M,S,J,M1,M2,PA(2),PB(2),PC(2),gama)*Afunc(N,T,K,N1,N2,PA(3),PB(3),PC(3),gama) 
                                        TERM3 = Afunc(L,R,I,L1,L2,PA(1),PB(1),PC(1),gama)*Afunc(M,S,J,M1,M2,PA(2),PB(2),PC(2),gama)*dAfunc(N,T,K,N1,N2,PA(3),PB(3),PC(3),gama) 
                                
                                        ! Derivative with respect to Cx
                                        gradient(3,1) = gradient(3,1) - TERM1*FNUVEC(L+M+N-2*(R+S+T)-(I+J+K)+1)

                                        ! Derivative with respect to Cy
                                        gradient(3,2) = gradient(3,2) - TERM2*FNUVEC(L+M+N-2*(R+S+T)-(I+J+K)+1)
                                
                                        ! Derivative with respect to Cz
                                        gradient(3,3) = gradient(3,3) - TERM3*FNUVEC(L+M+N-2*(R+S+T)-(I+J+K)+1)

                                        ! Adding the gradient of the  F_nu-function:
                                        gradient(3,:) = gradient(3,:) + 2*gama*PC*TERM0*FNUVEC(L+M+N-2*(R+S+T)-(I+J+K)+2)

                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      
      X = (alpha1*alpha2/gama)*DOT_PRODUCT(AB,AB)
      gradient(3,:) = gradient(3,:)*(2.0d0*pi/gama)*EXP(-X)

      !=========================================================================
      ! Secondly we calculate the gradients with respect to the nuclear
      ! coordinates A and B. They are stored in gradient(1,:) and geadient(2,:)
      !=========================================================================

      ! Derivative with respect to A:
      IF ( L1 .GT. 0 ) gradient(1,1) = -L1*primpotential(L1-1,M1,N1,A,alpha1,C,L2,M2,N2,B,alpha2) 
      gradient(1,1) =  gradient(1,1) + 2*alpha1*primpotential(L1+1,M1,N1,A,alpha1,C,L2,M2,N2,B,alpha2)
      
      IF ( M1 .GT. 0 ) gradient(1,2) = -M1*primpotential(L1,M1-1,N1,A,alpha1,C,L2,M2,N2,B,alpha2) 
      gradient(1,2) = gradient(1,2) + 2*alpha1*primpotential(L1,M1+1,N1,A,alpha1,C,L2,M2,N2,B,alpha2)
      
      IF ( N1 .GT. 0 ) gradient(1,3) = -N1*primpotential(L1,M1,N1-1,A,alpha1,C,L2,M2,N2,B,alpha2) 
      gradient(1,3) = gradient(1,3) + 2*alpha1*primpotential(L1,M1,N1+1,A,alpha1,C,L2,M2,N2,B,alpha2)

      ! Derivative with respect to B:
      IF ( L2 .GT. 0 ) gradient(2,1) = -L2*primpotential(L1,M1,N1,A,alpha1,C,L2-1,M2,N2,B,alpha2) 
      gradient(2,1) = gradient(2,1) + 2*alpha2*primpotential(L1,M1,N1,A,alpha1,C,L2+1,M2,N2,B,alpha2)
      
      IF ( M2 .GT. 0 ) gradient(2,2) = -M2*primpotential(L1,M1,N1,A,alpha1,C,L2,M2-1,N2,B,alpha2) 
      gradient(2,2) = gradient(2,2) + 2*alpha2*primpotential(L1,M1,N1,A,alpha1,C,L2,M2+1,N2,B,alpha2)
      
      IF ( N2 .GT. 0 ) gradient(2,3) = -N2*primpotential(L1,M1,N1,A,alpha1,C,L2,M2,N2-1,B,alpha2) 
      gradient(2,3) = gradient(2,3) + 2*alpha2*primpotential(L1,M1,N1,A,alpha1,C,L2,M2,N2+1,B,alpha2)
      
      END SUBROUTINE gradprimpotential
