SUBROUTINE linesearchmin(POLYORDER,ENERGY,DRF,NPOINTS,DR,DRM)
     ! This subroutine fits a number of NPOINTS energies and moves 
     ! to 4:th order polynomial in the intervall [0,DR] and calcuates
     ! the position of the minimum of this polynomial in the intervall 
     ! [0,DR]
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: NPOINTS,POLYORDER
     DOUBLE PRECISION, INTENT(IN) :: DR,ENERGY(NPOINTS),DRF(NPOINTS)
     DOUBLE PRECISION, INTENT(OUT) :: DRM
     DOUBLE PRECISION :: A(NPOINTS,POLYORDER),AA(POLYORDER,POLYORDER),B(POLYORDER),EMIN,ETEST,R,V(POLYORDER)
     INTEGER :: I,J,IPIV(5),NSEARCH,INFO
     EXTERNAL :: DGESV
     
     NSEARCH = 10000

     !=======================================
     ! Setting up the least squares matrices:
     !=======================================
     DO I=1,NPOINTS
        DO J=1,POLYORDER
                A(I,J) = DRF(I)**(J-1)
        ENDDO
     ENDDO

     AA = MATMUL(TRANSPOSE(A),A)
      B = MATMUL(TRANSPOSE(A),ENERGY)
     
     !======================================


     !======================================
     ! Performing the least squares fit
     !======================================

     CALL DGESV( POLYORDER, 1, AA, POLYORDER, IPIV, B, POLYORDER, INFO )
     
     !======================================

     !==================================================
     ! Searching for the minimum in the intervall [0,DR]
     !==================================================
     IF ( INFO .EQ. 0 ) THEN
             EMIN = ENERGY(1)
             DRM = DRF(1)
             DO I=2,NSEARCH
                R = DR*(I-1)/(NSEARCH-1)
                DO J=1,POLYORDER
                        V(J) = R**(J-1)
                ENDDO
                ! Estimating the energy at the point R
                ! using the Fourth order fit

                ETEST = DOT_PRODUCT(V,B)
                IF ( ETEST .LT. EMIN ) THEN
                        EMIN = ETEST
                        DRM  = R
                ENDIF
             ENDDO
     ELSE
             DRM = 0.0010d0
     ENDIF

END SUBROUTINE linesearchmin
