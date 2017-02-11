FUNCTION tmu(mu)
      ! This function is the function t(mu)
      ! described in appendix B Eqn (B9)
      ! In J. Chem. Phys. 98, 5612, (1993)
      IMPLICIT NONE
      DOUBLE PRECISION :: tmu
      DOUBLE PRECISION :: mu,muu
      DOUBLE PRECISION :: p1,p2,p3,sofmu

      
      p1 = 1.50d0*mu - 0.50d0*mu**3
      p2 = 1.50d0*p1 - 0.50d0*p1**3
      p3 = 1.50d0*p2 - 0.50d0*p2**3

      sofmu = 0.50d0*(1.0d0 - p3 )
      
      IF ( sofmu .NE. 0.0d0 ) THEN
        tmu = -(1.0d0/sofmu)*(27.0d0/16.0d0)*(1.0d0 - p2**2 )*( 1.0d0 - p1**2 )*( 1.0d0 - mu**2 )
      ELSE
             tmu = 0.0d0
      ENDIF

END FUNCTION tmu
