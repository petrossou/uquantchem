SUBROUTINE gradbasfunkval(PSI,R,grad)
      ! Calculates the value of the gradient of a real basis-function at the spatial
      ! point R. I.e, The function = GRAD(Psi(R))
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(IN) :: R(3)
      DOUBLE PRECISION, INTENT(OUT) :: grad(3)
      DOUBLE PRECISION :: RP(3),RP2,dp,dpp,dppp
      INTEGER :: I,Lp,Lpp,Lppp,Mp,Mpp,Mppp,Np,Npp,Nppp
     
      grad(:) = 0.0d0

      RP = R - PSI%R
      RP2 = DOT_PRODUCT(RP,RP)
      DO  I=1,PSI%NPRIM
                ! The derivative with respect to the x-coordinate (See  Eqn. 27, P.9 on QDMC-notes )
                dp   = PSI%L(1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                dpp  =       -2*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*PSI%EXPON(I)
                Lp   = PSI%L(1) - 1
                Lpp  = PSI%L(1) + 1
                IF ( Lp .GE. 0 ) grad(1) = grad(1) +  dp*( RP(1)**Lp  )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                grad(1) =                  grad(1) + dpp*( RP(1)**Lpp )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                
                ! The derivative with respect to the y-coordinate (See  Eqn. 27, P.9 on QDMC-notes )
                dp   = PSI%L(2)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                Mp   = PSI%L(2) - 1
                Mpp  = PSI%L(2) + 1
                IF ( Mp .GE. 0 ) grad(2) = grad(2) +  dp*( RP(1)**PSI%L(1) )*( RP(2)**Mp  )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                grad(2) =                  grad(2) + dpp*( RP(1)**PSI%L(1) )*( RP(2)**Mpp )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)

                ! The derivative with respect to the z-coordinate (See  Eqn. 27, P.9 on QDMC-notes )
                dp   = PSI%L(3)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                Np   = PSI%L(3) - 1
                Npp  = PSI%L(3) + 1
                IF ( Np .GE. 0 ) grad(3) = grad(3) +  dp*( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**Np  )*EXP(-PSI%EXPON(I)*RP2)
                grad(3) =                  grad(3) + dpp*( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**Npp )*EXP(-PSI%EXPON(I)*RP2)
      ENDDO

      END SUBROUTINE gradbasfunkval

