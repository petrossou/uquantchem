FUNCTION Fnu(N,X)
      ! The function F_nu(x) described in the Cook Book 
      ! on pages 244 and 280, is here calculated from 
      ! the Kummer confluent hypergeometric function.
      ! See WIREs Computational Molecular Science, 2012, 2: 290-303
      IMPLICIT NONE
      DOUBLE PRECISION :: Fnu
      DOUBLE PRECISION :: X
      INTEGER :: N
      DOUBLE PRECISION :: DF
      DOUBLE PRECISION, EXTERNAL :: gammaf
      DOUBLE PRECISION, PARAMETER :: Tol= 1.0E-16
      INTEGER, PARAMETER :: MAXNUMBEROFTERMS = 10000000
      DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION :: ak,bk,k,M
      INTEGER :: I

      IF ( N .NE. 0 .AND. X .NE. 0.0d0 .AND. X .LE. 5.0d0 ) THEN
                ak = N + 1.0d0/2.0d0 
                bk = N + 3.0d0/2.0d0 
                DF = -X*ak/bk             ! k = 1
                M = 1.0d0 + DF            ! k = 1
                k = 1.0d0                 ! k = 1
                !=========================================================
                ! Calculating the Kummer confluent hypergeometric function
                ! M(n+1/2,n+3/2,-X). Here the sum from 2 to "infinity"
                !=========================================================
               DO WHILE ( DABS(DF) .GT. Tol .AND. INT(k) .LE. MAXNUMBEROFTERMS .OR. k .EQ. 1 )
                        ak = ak + 1.0d0
                        bk = bk + 1.0d0
                        k =  k  + 1.0d0
                        DF =  DF*(-ak*X)/(k*bk)
                        M = M + DF
               ENDDO
               
               Fnu = M/(2.0d0*N + 1.0d0)
               
               IF ( INT(k) .GT. MAXNUMBEROFTERMS ) THEN
                        print*,'WARNING F_nu(x) function did not reach converge to'
                        print*,'target tolerance=',Tol,'ABORTING IN Fnu.f90'
                        STOP
               ENDIF
       ENDIF
       IF ( N .NE. 0 .AND. X .GT. 5.0d0 ) THEN
               Fnu = sqrt(PI/(4.0d0*X))*DERF(sqrt(X)) ! The Boys function of order 0
               ! Using the forward recursion formula (9.8.13), p. 367 in T. Helgaker et al
               I = 0
               DO WHILE ( I .LT. N ) 
                        Fnu = (( 2*I + 1 )*Fnu - DEXP(-X))/(2.0d0*X)
                        I = I + 1
               ENDDO
       ENDIF
    
       IF ( N .EQ. 0 .AND. X .NE. 0.0d0 ) Fnu = sqrt(PI/(4.0d0*X))*DERF(sqrt(X))
       IF ( X .EQ. 0.0d0 ) Fnu = 1.0d0/(2.0d0*N + 1.0d0)
     
END FUNCTION Fnu
