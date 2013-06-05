SUBROUTINE hessianbasfunk(PSI,R,hb)
      ! Calculates the hessian of a real basis-function at the spatial
      ! point R. I.e, The the matrix hb_ij = d^2 PHI/(dx_i*dx_j)
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASFUNCT), INTENT(IN) :: PSI
      DOUBLE PRECISION, INTENT(IN) :: R(3)
      DOUBLE PRECISION, INTENT(OUT) :: hb(3,3)
      DOUBLE PRECISION :: RP(3),RP2,A,B,C,D
      INTEGER :: I
     
      hb(:,:) = 0.0d0
      RP = R - PSI%R
      RP2 = DOT_PRODUCT(RP,RP)
      DO  I=1,PSI%NPRIM
                !-------
                ! hb_xx:
                !-------
                B = -2*PSI%EXPON(I)*(2*PSI%L(1)+1)*RP(1)**PSI%L(1)
                
                A = 0.0d0
                IF ( PSI%L(1) .GE. 2 ) A = PSI%L(1)*(PSI%L(1)-1)*RP(1)**(PSI%L(1)-2)
                
                C = (2*PSI%EXPON(I))**2*RP(1)**(PSI%L(1)+2)

                hb(1,1) = hb(1,1) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B + C )*RP(2)**PSI%L(2)*RP(3)**PSI%L(3)*EXP(-PSI%EXPON(I)*RP2)

                !--------------------------------------------------------------------------------------
                ! hb_yy:
                !-------
                B = -2*PSI%EXPON(I)*(2*PSI%L(2)+1)*RP(2)**PSI%L(2)
                
                A = 0.0d0
                IF ( PSI%L(2) .GE. 2 ) A = PSI%L(2)*(PSI%L(2)-1)*RP(2)**(PSI%L(2)-2)
                
                C = (2*PSI%EXPON(I))**2*RP(2)**(PSI%L(2)+2)

                hb(2,2) = hb(2,2) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B + C )*RP(1)**PSI%L(1)*RP(3)**PSI%L(3)*EXP(-PSI%EXPON(I)*RP2)
               
                !--------------------------------------------------------------------------------------
                ! hb_zz:
                !-------
                B = -2*PSI%EXPON(I)*(2*PSI%L(3)+1)*RP(3)**PSI%L(3)
               
                A = 0.0d0
                IF ( PSI%L(3) .GE. 2 ) A = PSI%L(3)*(PSI%L(3)-1)*RP(3)**(PSI%L(3)-2)
                
                C = (2*PSI%EXPON(I))**2*RP(3)**(PSI%L(3)+2)

                hb(3,3) = hb(3,3) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B + C )*RP(1)**PSI%L(1)*RP(2)**PSI%L(2)*EXP(-PSI%EXPON(I)*RP2)
              
                !--------------------------------------------------------------------------------------
                ! hb_xy:
                !------
                A = 0.0d0
                IF ( PSI%L(1) .GE. 1 ) A = PSI%L(1)*RP(1)**(PSI%L(1)-1)
                B = -2*PSI%EXPON(I)*RP(1)**(PSI%L(1)+1)

                C = 0.0d0
                IF ( PSI%L(2) .GE. 1 ) C = PSI%L(2)*RP(2)**(PSI%L(2)-1)
                D = -2*PSI%EXPON(I)*RP(2)**(PSI%L(2)+1)

                hb(1,2) = hb(1,2) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B )*( C + D )*RP(3)**PSI%L(3)*EXP(-PSI%EXPON(I)*RP2)
                hb(2,1) = hb(1,2)

                !--------------------------------------------------------------------------------------
                ! hb_xz:
                !------
                A = 0.0d0
                IF ( PSI%L(1) .GE. 1 ) A = PSI%L(1)*RP(1)**(PSI%L(1)-1)
                B = -2*PSI%EXPON(I)*RP(1)**(PSI%L(1)+1)

                C = 0.0d0
                IF ( PSI%L(3) .GE. 1 ) C = PSI%L(3)*RP(3)**(PSI%L(3)-1)
                D = -2*PSI%EXPON(I)*RP(3)**(PSI%L(3)+1)

                hb(1,3) = hb(1,3) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B )*( C + D )*RP(2)**PSI%L(2)*EXP(-PSI%EXPON(I)*RP2)
                hb(3,1) = hb(1,3)

                !--------------------------------------------------------------------------------------
                ! hb_yz:
                !------
                A = 0.0d0
                IF ( PSI%L(2) .GE. 1 ) A = PSI%L(2)*RP(2)**(PSI%L(2)-1)
                B = -2*PSI%EXPON(I)*RP(2)**(PSI%L(2)+1)

                C = 0.0d0
                IF ( PSI%L(3) .GE. 1 ) C = PSI%L(3)*RP(3)**(PSI%L(3)-1)
                D = -2*PSI%EXPON(I)*RP(3)**(PSI%L(3)+1)

                hb(2,3) = hb(2,3) + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*( A + B )*( C + D )*RP(1)**PSI%L(1)*EXP(-PSI%EXPON(I)*RP2)
                hb(3,2) = hb(2,3)
                
      ENDDO
      END SUBROUTINE hessianbasfunk
