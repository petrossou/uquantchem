FUNCTION primpotential(L1,M1,N1,A,alpha1,C,L2,M2,N2,B,alpha2)
      ! This function calculates the potential energy integral given on the
      ! bottom of p. 244 in the Cook Book.
      IMPLICIT NONE
      INTEGER :: L1,M1,N1,L2,M2,N2
      DOUBLE PRECISION :: primpotential
      DOUBLE PRECISION :: alpha1,alpha2,A(3),B(3),C(3),P(3),PA(3),PB(3),PC(3)
      DOUBLE PRECISION :: gama, X,TERM,AB(3)
      DOUBLE PRECISION, ALLOCATABLE :: FNUVEC(:)
      INTEGER :: L,R,I,M,S,J,N,T,K
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, EXTERNAL :: Afunc
      EXTERNAL :: Fnurec

      gama = alpha1+alpha2
      P = (A*alpha1+B*alpha2)/gama
      PA = P-A
      PB = P-B
      PC = P-C
      AB = A-B

      X = gama*DOT_PRODUCT(PC,PC)

      ALLOCATE(FNUVEC(L1+L2+M1+M2+N1+N2+1))
      
      CALL Fnurec(L1+L2+M1+M2+N1+N2,FNUVEC,X)
      
      primpotential = 0.0d0
      
      DO L=0,L1+L2
         DO R=0,(L-MOD(L,2))/2
            DO I=0,(L-2*R-MOD(L-2*R,2))/2
               DO M=0,M1+M2
                  DO S=0,(M-MOD(M,2))/2
                     DO J=0,(M-2*S-MOD(M-2*S,2))/2
                        DO N=0,N1+N2
                           DO T=0,(N-MOD(N,2))/2
                              DO K=0,(N-2*T-MOD(N-2*T,2))/2
                                TERM = Afunc(L,R,I,L1,L2,PA(1),PB(1),PC(1),gama)*Afunc(M,S,J,M1,M2,PA(2),PB(2),PC(2),gama)*Afunc(N,T,K,N1,N2,PA(3),PB(3),PC(3),gama)
                                primpotential = primpotential + TERM*FNUVEC(L+M+N-2*(R+S+T)-(I+J+K)+1)
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
      primpotential = primpotential*(2.0d0*pi/gama)*EXP(-X)

      END FUNCTION primpotential
