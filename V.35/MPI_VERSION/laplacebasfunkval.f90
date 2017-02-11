FUNCTION laplacebasfunkval(PSI,R)
      ! Calculates the value of the Laplacian of a real basis-function at the spatial
      ! point R. I.e, The function = NABLA^2(Psi(R))
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: laplacebasfunkval
      TYPE(BASFUNCT) :: PSI
      DOUBLE PRECISION :: R(3),RP(3),RP2,dp,dpp,dppp
      INTEGER :: I,Lp,Lpp,Lppp,Mp,Mpp,Mppp,Np,Npp,Nppp
     
      laplacebasfunkval = 0.0d0

      RP = R - PSI%R
      RP2 = DOT_PRODUCT(RP,RP)
      DO  I=1,PSI%NPRIM
                ! The second derivative with respect to the x-coordinate (See  Eqn. 36, P.13 on QDMC-notes )
                dp   = PSI%L(1)*(PSI%L(1) - 1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                dpp  =       -2*(2*PSI%L(1)+1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*PSI%EXPON(I)
                dppp =                       4*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*PSI%EXPON(I)**2
                Lp   = PSI%L(1) - 2
                Lpp  = PSI%L(1)
                Lppp = PSI%L(1) + 2
                IF ( Lp   .GE. 0 ) laplacebasfunkval = laplacebasfunkval +   dp*( RP(1)**Lp   )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Lpp  .GE. 0 ) laplacebasfunkval = laplacebasfunkval +  dpp*( RP(1)**Lpp  )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Lppp .GE. 0 ) laplacebasfunkval = laplacebasfunkval + dppp*( RP(1)**Lppp )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                
                ! The second derivative with respect to the y-coordinate (See  Eqn. 36, P.13 on QDMC-notes )
                dp   = PSI%L(2)*(PSI%L(2) - 1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                dpp  =       -2*(2*PSI%L(2)+1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*PSI%EXPON(I)
                Mp   = PSI%L(2) - 2
                Mpp  = PSI%L(2)
                Mppp = PSI%L(2) + 2
                IF ( Mp   .GE. 0 ) laplacebasfunkval = laplacebasfunkval +   dp*( RP(1)**PSI%L(1) )*( RP(2)**Mp   )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Mpp  .GE. 0 ) laplacebasfunkval = laplacebasfunkval +  dpp*( RP(1)**PSI%L(1) )*( RP(2)**Mpp  )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Mppp .GE. 0 ) laplacebasfunkval = laplacebasfunkval + dppp*( RP(1)**PSI%L(1) )*( RP(2)**Mppp )*( RP(3)**PSI%L(3) )*EXP(-PSI%EXPON(I)*RP2)

                ! The second derivative with respect to the z-coordinate (See  Eqn. 36, P.13 on QDMC-notes )
                dp   = PSI%L(3)*(PSI%L(3) - 1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM
                dpp  =       -2*(2*PSI%L(3)+1)*PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*PSI%EXPON(I)
                Np   = PSI%L(3) - 2
                Npp  = PSI%L(3)
                Nppp = PSI%L(3) + 2
                IF ( Np   .GE. 0 ) laplacebasfunkval = laplacebasfunkval +   dp*( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**Np   )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Npp  .GE. 0 ) laplacebasfunkval = laplacebasfunkval +  dpp*( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**Npp  )*EXP(-PSI%EXPON(I)*RP2)
                IF ( Nppp .GE. 0 ) laplacebasfunkval = laplacebasfunkval + dppp*( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**Nppp )*EXP(-PSI%EXPON(I)*RP2)
      ENDDO
      END FUNCTION laplacebasfunkval

