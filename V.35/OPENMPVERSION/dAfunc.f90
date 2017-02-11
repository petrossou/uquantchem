FUNCTION dAfunc(L,R,I,L1,L2,PA,PB,PC,gama)
      ! This is the function A_lri(l1,l2,Ax,Bx,Cx,gamma) defined
      ! at the top of page 245 in the Cook Book.
      IMPLICIT NONE
      DOUBLE PRECISION :: dAfunc
      INTEGER :: L,R,I,L1,L2
      DOUBLE PRECISION :: PA,PB,PC,gama
      DOUBLE PRECISION, EXTERNAL :: fj,fac
      DOUBLE PRECISION :: eps
      eps = 1.0d0/(4.0d0*gama)
      
      IF ( L-2*R-2*I .NE. 0 ) THEN
                dAfunc = (((-1.0d0)**L)*fj(L,L1,L2,PA,PB)*((-1.0d0)**I)*fac(L)*(PC**(L-2*R-2*I-1))*eps**(R+I))/(fac(R)*fac(I)*fac(L-2*R-2*I-1))
      ELSE
              dAfunc = 0.0d0
      ENDIF

END FUNCTION dAfunc
