SUBROUTINE nglaplacebasfunkval(PSI,R,nglapl)
      ! Calculates the value of nuclear gradient of the laplacian of a real basis-function at the spatial
      ! point R. I.e, The function = GRAD_nuc[Laplace(Psi(R))]
      ! The output is stored in the array, nglapl, according to the
      ! following ordering:
      !nglapl(1) = dLaplace(Psi(R))]/dRx
      !nglapl(2) = dLaplace(Psi(R))]/dRy
      !nglapl(3) = dLaplace(Psi(R))]/dRz
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(IN) :: R(3)
      DOUBLE PRECISION, INTENT(OUT) :: nglapl(3)
      DOUBLE PRECISION :: hh,xv(3),yv(3),zv(3),gp1,gm1,gp2,gm2
      DOUBLE PRECISION, EXTERNAL  :: laplacebasfunkval
      xv(:) = 0.0d0
      yv(:) = 0.0d0
      zv(:) = 0.0d0
      xv(1) = 1.0d0
      yv(2) = 1.0d0
      zv(3) = 1.0d0

      hh = 0.00010d0
      
      ! Calculating the derivative of the laplacian with respect to the
      ! x-component of the nuclear coordinate, Rx.
      gp1 =  laplacebasfunkval(PSI,R+xv*hh)
      gm1 =  laplacebasfunkval(PSI,R-xv*hh)
      gp2 =  laplacebasfunkval(PSI,R+2.0d0*xv*hh)
      gm2 =  laplacebasfunkval(PSI,R-2.0d0*xv*hh)

      nglapl(1) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )
     
      ! Calculating the derivative of the laplacian with respect to the
      ! y-component of the nuclear coordinate, Ry.
      gp1 =  laplacebasfunkval(PSI,R+yv*hh)
      gm1 =  laplacebasfunkval(PSI,R-yv*hh)
      gp2 =  laplacebasfunkval(PSI,R+2.0d0*yv*hh)
      gm2 =  laplacebasfunkval(PSI,R-2.0d0*yv*hh)
      nglapl(2) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )
      
      ! Calculating the derivative of the laplacian with respect to the
      ! z-component of the nuclear coordinate, Rz.
      gp1 =  laplacebasfunkval(PSI,R+zv*hh)
      gm1 =  laplacebasfunkval(PSI,R-zv*hh)
      gp2 =  laplacebasfunkval(PSI,R+2.0d0*zv*hh)
      gm2 =  laplacebasfunkval(PSI,R-2.0d0*zv*hh)

      nglapl(3) = -(1.0d0/(12.0d0*hh))*( gm2 - 8.0d0*gm1 + 8.0d0*gp1 - gp2 )

END SUBROUTINE nglaplacebasfunkval

