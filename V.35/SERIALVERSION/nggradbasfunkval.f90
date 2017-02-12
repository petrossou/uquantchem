SUBROUTINE nggradbasfunkval(PSI,R,nggrad)
      ! Calculates the value of nuclear gradient of the gradient of a real basis-function at the spatial
      ! point R. I.e, The function = GRAD_nuc[GRAD(Psi(R))]
      ! The output is stored in the 3x3 matrix nggrad, according to the
      ! following ordering:
      ! nggrad(1,1) = d^2(Psi)/(dRxdx), nggrad(1,2) = d^2(Psi)/(dRydx), nggrad(1,3) = d^2(Psi)/(dRzdx)
      ! nggrad(2,1) = d^2(Psi)/(dRxdy), nggrad(2,2) = d^2(Psi)/(dRydy), nggrad(2,3) = d^2(Psi)/(dRzdy)
      ! nggrad(3,1) = d^2(Psi)/(dRxdz), nggrad(3,2) = d^2(Psi)/(dRydz), nggrad(3,3) = d^2(Psi)/(dRzdz)
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(IN) :: R(3)
      DOUBLE PRECISION, INTENT(OUT) :: nggrad(3,3)
      DOUBLE PRECISION :: hh,xv(3),yv(3),zv(3),gp1(3),gm1(3),gp2(3),gm2(3)

      xv(:) = 0.0d0
      yv(:) = 0.0d0
      zv(:) = 0.0d0
      xv(1) = 1.0d0
      yv(2) = 1.0d0
      zv(3) = 1.0d0

      hh = 0.00010d0
      
      ! Calculating the derivative of the gradient with respect to the
      ! x-component of the nuclear coordinate, Rx.
      CALL gradbasfunkval(PSI,R+xv*hh,gp1)
      CALL gradbasfunkval(PSI,R-xv*hh,gm1)
      CALL gradbasfunkval(PSI,R+2.0d0*xv*hh,gp2)
      CALL gradbasfunkval(PSI,R-2.0d0*xv*hh,gm2)

      nggrad(:,1) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )
     
      ! Calculating the derivative of the gradient with respect to the
      ! y-component of the nuclear coordinate, Ry.
      CALL gradbasfunkval(PSI,R+yv*hh,gp1)
      CALL gradbasfunkval(PSI,R-yv*hh,gm1)
      CALL gradbasfunkval(PSI,R+2.0d0*yv*hh,gp2)
      CALL gradbasfunkval(PSI,R-2.0d0*yv*hh,gm2)

      nggrad(:,2) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )
      
      ! Calculating the derivative of the gradient with respect to the
      ! z-component of the nuclear coordinate, Rz.
      CALL gradbasfunkval(PSI,R+zv*hh,gp1)
      CALL gradbasfunkval(PSI,R-zv*hh,gm1)
      CALL gradbasfunkval(PSI,R+2.0d0*zv*hh,gp2)
      CALL gradbasfunkval(PSI,R-2.0d0*zv*hh,gm2)

      nggrad(:,3) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )

END SUBROUTINE nggradbasfunkval

