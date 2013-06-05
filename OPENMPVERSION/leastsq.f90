SUBROUTINE leastsq(rc,N,PSI,Acoeff,LAMDA,C)
      ! Performs a least square fit of the function 
      ! f(r) = Acoeff*EXP(-LAMDA*r) + C to the N data points
      ! r(i) = (rc/(N-1))*(i-1), y(i) = GTO(r(i))
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: rc
      INTEGER, INTENT(IN) :: N
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(OUT) :: Acoeff,LAMDA,C
      DOUBLE PRECISION :: r(N), A1(N,2),A2(N,3),b(N),ATA1(2,2),ATA2(3,3),x(2),y(2),dx(3),dy(3)
      INTEGER :: I,J,INFO,IPIV1(2),IPIV2(3),NITERMAX
      DOUBLE PRECISION :: tol,MIX,x2old
      DOUBLE PRECISION, EXTERNAL :: gto
       
      NITERMAX = 200

      !==================================================
      ! Performing a simple least square fit to the 
      ! function: f(r) = Acoeff*EXP(-LAMDA*r)
      !==================================================

      DO I=1,N
        r(I) = rc*(I-1)/(1.0d0*(N-1))
        A1(I,1) = 1.0d0
        A1(I,2) = -r(I)
        b(I) = LOG(gto(PSI,r(I)))
       ENDDO

       ATA1 = MATMUL(TRANSPOSE(A1),A1)
       y = MATMUL(TRANSPOSE(A1),b)
       CALL DGESV( 2, 1, ATA1, 2, IPIV1, y, 2, INFO )

       Acoeff = EXP(y(1))
       LAMDA = y(2)
       C = 0.00d0

       !====================================================
       ! Performing a non-linear least square fit to the 
       ! function: f(r) = Acoeff*EXP(-LAMDA*r) + C 
       !====================================================

       tol = 1.0d0
       J = 0
       
       DO WHILE ( tol .GT. 1.0E-10 .AND. J .LT. NITERMAX )
                dy = 0.0d0
                DO I=1,N
                        r(I) = rc*(I-1)/(1.0d0*(N-1))
                        A2(I,1) = EXP(-LAMDA*r(I))
                        A2(I,2) = -Acoeff*R(I)*EXP(-LAMDA*r(I))
                        A2(I,3) = 1.0d0
                        dy(1) =  dy(1) + 2.0d0*( Acoeff*EXP(-LAMDA*r(I)) + C - gto(PSI,r(I)) )*EXP(-LAMDA*r(I))
                        dy(2) =  dy(2) - 2.0d0*( Acoeff*EXP(-LAMDA*r(I)) + C - gto(PSI,r(I)) )*Acoeff*r(I)*EXP(-LAMDA*r(I))
                        dy(3) =  dy(3) + 2.0d0*( Acoeff*EXP(-LAMDA*r(I)) + C - gto(PSI,r(I)) )
                        b(I) =  Acoeff*EXP(-LAMDA*r(I)) + C - gto(PSI,r(I))
                ENDDO

                IF ( J .EQ. 0 ) THEN
                        MIX = 0.001d0
                ELSE
                       IF ( DOT_PRODUCT(b,b) .GE. x2old ) THEN
                                MIX = MIX*10.0d0
                       ELSE
                                MIX = MIX/10.0d0
                       ENDIF
                ENDIF

                x2old = DOT_PRODUCT(b,b)
                ATA2 = MATMUL(TRANSPOSE(A2),A2)

                ATA2(1,1) = ATA2(1,1)*(1.0d0 + MIX )
                ATA2(2,2) = ATA2(2,2)*(1.0d0 + MIX )
                ATA2(3,3) = ATA2(3,3)*(1.0d0 + MIX )

                CALL DGESV( 3, 1, ATA2, 3, IPIV2, dy, 3, INFO )
                
                Acoeff = Acoeff + dy(1)
                LAMDA  = LAMDA +  dy(2)
                C      =  C     + dy(3)

                tol = sqrt(DOT_PRODUCT(dy,dy))
                J = J+1
        ENDDO
        IF ( J .GE. NITERMAX ) THEN
                WRITE(*,*)'CUSP CORRECTION FAILED!'
        ENDIF
END SUBROUTINE leastsq
