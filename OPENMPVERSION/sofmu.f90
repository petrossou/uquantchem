FUNCTION sofmu(mu)
      ! This is the function s(mu) = (1/2)*(1 - f(mu) )
      ! described by Eqn (3) 
      ! in A. D. Becke, J. Chem. Phys. 88, 2547 (1988).
      ! Here the recommended third order (k=3) is used 
      IMPLICIT NONE
      DOUBLE PRECISION :: sofmu
      DOUBLE PRECISION :: mu
      DOUBLE PRECISION :: p1,p2,p3

      p1 = 1.50d0*mu - 0.50d0*mu**3
      p2 = 1.50d0*p1 - 0.50d0*p1**3
      p3 = 1.50d0*p2 - 0.50d0*p2**3

      sofmu = 0.50d0*(1.0d0 - p3 )

END FUNCTION sofmu
