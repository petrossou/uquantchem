SUBROUTINE gradprimeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW,grad)
      ! This subroutine calculates the gradient of the two-electron integral for primitive gaussians 
      ! primeeintegral = (i j | k l ) = ( 1 2 | 3 4 ) = ( phi_1(r)*phi_2(r) | phi_3(r')*phi_4(r') ) 
      ! through Rys quadrature, J. Comp. Chem. 11, 972-977 (1990)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4
      LOGICAL,INTENT(IN) :: SAMECITE
      DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,alpha3,alpha4,PRYSR(25,25),PRYSW(25,25)
      DOUBLE PRECISION, INTENT(IN) :: A(3),B(3),C(3),D(3)
      DOUBLE PRECISION, INTENT(OUT) :: grad(4,3)
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      
      grad(:,:)  = 0.0d0
      ! Calculation of gradient with respect to A(=B):
      IF ( A(1) .EQ. B(1) .AND. A(2) .EQ. B(2) .AND. A(3) .EQ. B(3) ) THEN
              IF ( L1 .GT. 0 ) THEN
                      grad(1,1) = -(L1+L2)*primeeintegral(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                      IF ( L2 .GT. 0 ) grad(1,1) = -(L1+L2)*primeeintegral(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(1,1) = grad(1,1) +2*(alpha1+alpha2)*primeeintegral(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)


              IF ( M1 .GT. 0 ) THEN
                      grad(1,2) = -(M1+M2)*primeeintegral(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                      IF ( M2 .GT. 0 ) grad(1,2) = -(M1+M2)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(1,2) = grad(1,2) +2*(alpha1+alpha2)*primeeintegral(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)


              IF ( N1 .GT. 0 ) THEN
                      grad(1,3) = -(N1+N2)*primeeintegral(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                      IF ( N2 .GT. 0 ) grad(1,3) = -(N1+N2)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(1,3) = grad(1,3) +2*(alpha1+alpha2)*primeeintegral(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)


              grad(2,:) = 0.0d0

      ELSE
              ! Calculation of gradient with respect to A
                IF ( L1 .GT. 0 ) grad(1,1) = -L1*primeeintegral(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(1,1) = grad(1,1) +2*alpha1*primeeintegral(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

                IF ( M1 .GT. 0 ) grad(1,2) = -M1*primeeintegral(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(1,2) = grad(1,2) +2*alpha1*primeeintegral(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                IF ( N1 .GT. 0 ) grad(1,3) = -N1*primeeintegral(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(1,3) = grad(1,3) +2*alpha1*primeeintegral(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                ! Calculation of gradient with respect to B:
                IF ( L2 .GT. 0 ) grad(2,1) = -L2*primeeintegral(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(2,1) = grad(2,1) +2*alpha2*primeeintegral(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

                IF ( M2 .GT. 0 ) grad(2,2) = -M2*primeeintegral(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(2,2) = grad(2,2) +2*alpha2*primeeintegral(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                IF ( N2 .GT. 0 ) grad(2,3) = -N2*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(2,3) = grad(2,3) +2*alpha2*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      ENDIF
      
      ! Calculation of gradient with respect to C(=D):
      IF ( C(1) .EQ. D(1) .AND. C(2) .EQ. D(2) .AND. C(3) .EQ. D(3) ) THEN
              IF ( L3 .GT. 0 ) THEN
                      grad(3,1) = -(L3+L4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3-1,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                      IF ( L4 .GT. 0 ) grad(3,1) = -(L3+L4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4-1,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(3,1) = grad(3,1) +2*(alpha3+alpha4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3+1,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

              
              IF ( M3 .GT. 0 ) THEN
                        grad(3,2) = -(M3+M4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3-1,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                        IF ( M4 .GT. 0 ) grad(3,2) = -(M3+M4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4-1,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(3,2) = grad(3,2) +2*(alpha3+alpha4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3+1,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

              
              IF ( N3 .GT. 0 ) THEN 
                      grad(3,3) = -(N3+N4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3-1,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ELSE
                      IF ( N4 .GT. 0 )grad(3,3) = -(N3+N4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4-1,D,alpha4,SAMECITE,PRYSR,PRYSW)
              ENDIF
              grad(3,3) = grad(3,3) +2*(alpha3+alpha4)*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3+1,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

              
              grad(4,:) = 0.0d0

      ELSE
                ! Calculation of gradient with respect to C:
                IF ( L3 .GT. 0 ) grad(3,1) = -L3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3-1,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(3,1) = grad(3,1) +2*alpha3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3+1,M3,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

                IF ( M3 .GT. 0 ) grad(3,2) = -M3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3-1,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(3,2) = grad(3,2) +2*alpha3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3+1,N3,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                IF ( N3 .GT. 0 ) grad(3,3) = -N3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3-1,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(3,3) = grad(3,3) +2*alpha3*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3+1,C,alpha3,L4,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                ! Calculation of gradient with respect to D:
                IF ( L4 .GT. 0 ) grad(4,1) = -L4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4-1,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(4,1) = grad(4,1) +2*alpha4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4+1,M4,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)

                IF ( M4 .GT. 0 ) grad(4,2) = -M4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4-1,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(4,2) = grad(4,2) +2*alpha4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4+1,N4,D,alpha4,SAMECITE,PRYSR,PRYSW)
      
                IF ( N4 .GT. 0 ) grad(4,3) = -N4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4-1,D,alpha4,SAMECITE,PRYSR,PRYSW)
                grad(4,3) = grad(4,3) +2*alpha4*primeeintegral(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4+1,D,alpha4,SAMECITE,PRYSR,PRYSW)
    ENDIF
      END SUBROUTINE gradprimeeintegral
