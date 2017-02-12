FUNCTION oneint(coord,I,J,K,L,RI,RJ,RK,RL,alphaI,alphaJ,alphaK,alphaL,PECHNIDI,ta)
      ! This function calculates the Ix(ta) term on p. 974 in J. Comp. Chem. 11, 972-977 (1990)
      ! which is used to calculate the two electron integrals through
      ! Rys-quadrature.
      IMPLICIT NONE
      DOUBLE PRECISION :: oneint
      INTEGER :: I,J,K,L,coord,N,M
      LOGICAL :: PECHNIDI
      DOUBLE PRECISION :: RI(3),RJ(3),RK(3),RL(3),alphaI,alphaJ,alphaK,alphaL,ta
      DOUBLE PRECISION :: A,B,Rho,P(3),Q(3),taa
      DOUBLE PRECISION :: B00,B10,C00,B10p,C00p
      DOUBLE PRECISION, ALLOCATABLE :: Icoord(:,:),ixjxm(:)
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, EXTERNAL :: binomfac
      
      ALLOCATE(Icoord(I+J+1,K+L+1),ixjxm(K+L+1))


      A = alphaI + alphaJ
      B = alphaK + alphaL

      P = (alphaI*RI+alphaJ*RJ)/A
      Q = (alphaK*RK+alphaL*RL)/B

      Rho = A*B/(A+B)
       
      IF ( .not. PECHNIDI ) THEN
                !-----------------------------------------------------------
                ! Recursion parameters to be used with the roots and weights 
                ! calculated from scratch with rysquad.f
                !-----------------------------------------------------------
                ! Equations in (12-14) J. Comp. Chem. 11, 972-977 (1990):
                !-----------------------------------------------------------
                B00  = (1.0d0/(2.0d0*(A+B)))*ta**2

                B10  = 1.0d0/(2.0d0*A) - (B/(2.0d0*A*(A+B)))*ta**2
      
                B10p = 1.0d0/(2.0d0*B) - (A/(2.0d0*B*(A+B)))*ta**2

                C00  = P(coord)-RI(coord) + (B*(Q(coord)-P(coord))/(A+B))*ta**2
      
                C00p = Q(coord)-RK(coord) + (A*(P(coord)-Q(coord))/(A+B))*ta**2
      ELSE
              !----------------------------------------------------------------
              ! Recursion parameters to be used with the roots and weights 
              ! obtained from the polynomial fits taken from GAMESS, i.e 
              ! together with the routines ROOTWEIGHTMAX3.f90, ROOTWEIGHT4.f90 
              ! and ROOTWEIGHT5.f90.
              !----------------------------------------------------------------
                taa = ta/((A+B)*(ta+1.0d0))

                B00 = 0.50d0*taa 

                B10 = 1.0d0/(2.0d0*A*(ta+1.0d0)) + 0.50d0*taa

                B10p = 1.0d0/(2.0d0*B*(ta+1.0d0)) + 0.50d0*taa

                C00  = (P(coord)-RI(coord))/(1.0d0 + ta ) + ( B*(Q(coord)-RI(coord)) + A*(P(coord)-RI(coord)) )*taa
                
                C00p  = (Q(coord)-RK(coord))/(1.0d0 + ta ) + ( B*(Q(coord)-RK(coord)) + A*(P(coord)-RK(coord)) )*taa
      ENDIF


      ! Equation (11) in J. Comp. Chem. 11, 972-977 (1990):

      Icoord(1,1) = (pi/sqrt(A*B))*exp( -(alphaI*alphaJ/A)*(RI(coord)-RJ(coord))**2 - (alphaK*alphaL/B)*(RK(coord)-RL(coord))**2 )

      ! Reccursion according to Eqn (15) in J. Comp. Chem. 11, 972-977 (1990),
      ! in order to calculate Ix(1,0,0,0):
      
      IF ( I+J .GT. 0 ) THEN
        Icoord(2,1) = C00*Icoord(1,1)
      ENDIF

      ! Reccursion according to Eqn (15) in J. Comp. Chem. 11, 972-977 (1990),
      ! in order to calculate Ix(0,1,0,0):
      
      IF ( K+L .GT. 0 ) THEN
              Icoord(1,2) = C00p*Icoord(1,1)
      ENDIF
      
      ! Reccursion according to Eqn (15) in J. Comp. Chem. 11, 972-977
      ! (1990), ! in order to calculate Ix(ix+jx,0,0,0):
      IF ( I+J .GT. 1 ) THEN
                DO N=2,I+J
                        Icoord(1+N,1) = (N-1)*B10*Icoord(1+N-2,1) + C00*Icoord(1+N-1,1)
                ENDDO
      ENDIF

      ! Reccursion according to Eqn (16) in J. Comp. Chem. 11, 972-977
      ! (1990), ! in order to calculate Ix(0,kx+lx,0,0):
      IF ( K+L .GT. 0 ) THEN
              DO M=2,K+L
                       Icoord(1,1+M) = (M-1)*B10p*Icoord(1,1+M-2) + C00p*Icoord(1,1+M-1)
              ENDDO
      ENDIF
      
      ! Reccursion according to Eqn (16) in J. Comp. Chem. 11, 972-977 (1990),
      ! in order to calculate Ix(i,0,k,0) for 0 <= i <= ix+jx and 0 <= k <=  kx+lx
      IF ( I+J .NE. 0 .AND. K+L .NE. 0 ) THEN
                DO N=1,I+J
                        Icoord(1+N,2) = N*B00*Icoord(1+N-1,1) + C00p*Icoord(1+N,1)
                        DO M=2,K+L
                                Icoord(1+N,1+M) = (M-1)*B10p*Icoord(1+N,1+M-2) + N*B00*Icoord(1+N-1,1+M-1) + C00p*Icoord(1+N,1+M-1)
                        ENDDO
                ENDDO
      ENDIF

      ! Here we calculate Ix(ix,jx,kx+M,0), for M =0,..,lx according to Eqn. (19) in J. Comp. Chem. 11, 972-977 (1990)
      ! and from Ix(ix,jx,kx,0), ! Ix(ix,jx,kx+1,0),Ix(ix,jx,kx+2,0),...,Ix(ix,jx,kx+lx,0) we calculate with a
      ! similar recursion formula Ix(ix,jx,kx,lx)
      oneint = 0.0d0

      DO M=0,L
        ixjxm(K+M+1) = 0.0d0
         DO N=0,J
                !Here we calculate Ix(ix,jx,kx+M,0):
                ixjxm(K+M+1) = ixjxm(K+M+1) + binomfac(J,N)*Icoord(1+I+N,1+K+M)*(RI(coord)-RJ(coord))**(J-N)
         ENDDO
         ! Here we calculate Ix(ix,jx,kx,lx):
         oneint = oneint + binomfac(L,M)*ixjxm(K+M+1)*(RK(coord)-RL(coord))**(L-M)
      ENDDO
        
      END FUNCTION oneint
