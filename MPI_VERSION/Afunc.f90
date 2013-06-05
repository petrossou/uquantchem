FUNCTION Afunc(L,R,I,L1,L2,PA,PB,PC,gama)
      ! This is the function A_lri(l1,l2,Ax,Bx,Cx,gamma) defined
      ! at the top of page 245 in the Cook Book.
      IMPLICIT NONE
      DOUBLE PRECISION :: Afunc
      INTEGER :: L,R,I,L1,L2
      DOUBLE PRECISION :: PA,PB,PC,gama
      DOUBLE PRECISION, EXTERNAL :: fj,fac
      DOUBLE PRECISION :: eps
      eps = 1.0d0/(4.0d0*gama)

      Afunc = (((-1.0d0)**L)*fj(L,L1,L2,PA,PB)*((-1.0d0)**I)*fac(L)*(PC**(L-2*R-2*I))*eps**(R+I))/(fac(R)*fac(I)*fac(L-2*R-2*I))

END FUNCTION Afunc
