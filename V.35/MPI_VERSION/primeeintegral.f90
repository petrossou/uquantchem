FUNCTION primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      ! This function calculates the two-electron integral for primitive gaussians 
      ! primeeintegral = (i j | k l ) = ( 1 2 | 3 4 ) = ( phi_1(r)*phi_2(r) | phi_3(r')*phi_4(r') ) 
      ! through Rys quadrature, J. Comp. Chem. 11, 972-977 (1990)
      IMPLICIT NONE
      DOUBLE PRECISION :: primeeintegral
      INTEGER :: L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4
      LOGICAL :: SAMECITE,PECHNIDI
      DOUBLE PRECISION :: alpha1,alpha2,alpha3,alpha4,PRYSR(25,25),PRYSW(25,25)
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3)
      DOUBLE PRECISION, EXTERNAL :: oneint
      EXTERNAL :: RysQuad
      INTEGER :: LAMBDA,N,I
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION :: Rho,al1,al2,ROOTS(25),WEIGHT(25),X,P(3),Q(3),TERM
        
      al1 = alpha1 + alpha2
      al2 = alpha3 + alpha4

      P = (A*alpha1 + B*alpha2)/al1
      Q = (C*alpha3 + D*alpha4)/al2

      Rho = al1*al2/(al1 + al2 )

      X = Rho*DOT_PRODUCT(P-Q,P-Q)

      LAMBDA = L1 + M1 + N1 + L2 + M2 + N2 + L3 + M3 + N3 + L4 + M4 + N4

      N = (LAMBDA+MOD(LAMBDA,2))/2 + 1 - MOD(LAMBDA,2)
      
      IF ( .not. SAMECITE ) THEN
                ! Here the rys roots and weights are calculated from 
                ! scratch, a very slow process
                IF ( N .GT. 5 ) THEN 
                        call RysQuad(N,X,ROOTS,WEIGHT)
                        PECHNIDI = .FALSE.
                ENDIF
                !-----------------------------------------------------
                ! Here the rys roots and weights are calculated from 
                ! polynomial fits of pree-calculated roots and weights.
                ! The fits are taken from the GAMES-package ( Same as 
                ! in Pyquante), a very fast process.
                !------------------------------------------------------
                IF ( N .LE. 3 ) THEN
                        call ROOTWEIGHTMAX3(N,X,ROOTS,WEIGHT)
                        PECHNIDI = .TRUE.
                ENDIF
                IF ( N .EQ. 4 ) THEN
                        call ROOTWEIGHT4(N,X,ROOTS,WEIGHT)
                        PECHNIDI = .TRUE.
                ENDIF
                IF ( N .EQ. 5 ) THEN
                        call ROOTWEIGHT5(N,X,ROOTS,WEIGHT)
                        PECHNIDI = .TRUE.
                ENDIF
                
                primeeintegral = 0.0d0
                DO I=1,N
                        TERM =      oneint(1,L1,L2,L3,L4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,ROOTS(I))
                        TERM = TERM*oneint(2,M1,M2,M3,M4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,ROOTS(I))
                        TERM = TERM*oneint(3,N1,N2,N3,N4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,ROOTS(I))*WEIGHT(I)
                        primeeintegral = primeeintegral + TERM
                ENDDO
      ELSE
                ! Here we use the pre-calculated roots and weights of the rys
                ! polynomials in the case when A=B=C=D:
                
                primeeintegral = 0.0d0
                DO I=1,N
                        PECHNIDI = .FALSE.
                        TERM =      oneint(1,L1,L2,L3,L4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,PRYSR(N,I))
                        TERM = TERM*oneint(2,M1,M2,M3,M4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,PRYSR(N,I))
                        TERM = TERM*oneint(3,N1,N2,N3,N4,A,B,C,D,alpha1,alpha2,alpha3,alpha4,PECHNIDI,PRYSR(N,I))*PRYSW(N,I)
                        primeeintegral = primeeintegral + TERM
                ENDDO
        ENDIF

      primeeintegral = 2.0d0*sqrt(Rho/pi)*primeeintegral

      END FUNCTION primeeintegral
