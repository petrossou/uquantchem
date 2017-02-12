MODULE exchcorrmodule
      IMPLICIT NONE
      
      CONTAINS

!==============================================================================
! In the first part of this module all the exchange correlation potentials are
! listed. In the second part the auxilliary functions are listed
!==============================================================================
        SUBROUTINE VWN(densu,densd,Vc)
             !------------------------------------------------------------!
             ! This subroutine calculates the LSDA correlation            !
             ! potential of Vosko, Wilk, and Nusair (VWM) at the point r  !
             ! See for example  Phys. Rev. B. 42, 3205 (1990).            !
             ! densu = spin up density, densd = spin down density         !
             !------------------------------------------------------------!
             USE datatypemodule
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN) :: densu,densd
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION :: rs,zeta,dens
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             
             dens = densu + densd 
             ! calculating zeta and rs
             IF ( dens .NE. 0.0d0 ) THEN
                     zeta = (densu-densd)/dens
                     rs = (3.0d0/(4.0d0*PI*dens))**(1.0d0/3.0d0)
             ELSE
                     zeta = 0.0d0
                     rs = 1.0E20
             ENDIF

             ! The spin up component of the correlation potential
             Vc(1) = muP(rs) + deps(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betaf(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdr(rs)*zeta**(4.0d0) ) )

             Vc(1) = Vc(1) + (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betaf(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betaf(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 - zeta )
             

             ! The spin down component of the correlation potential
             Vc(2) = muP(rs) + deps(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betaf(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdr(rs)*zeta**(4.0d0) ) )

             Vc(2) = Vc(2) - (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betaf(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betaf(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 + zeta )
                

        END SUBROUTINE VWN
        
        SUBROUTINE VLYP(densu,densd,gdensu,gdensd,ldensu,ldensd,Vc)
             !------------------------------------------------------------!
             ! This subroutine calculates Gradient corrected correlation  !
             ! potential of Lee, Yang and Parr (LYP) at the point r       !
             ! For details see   Phys. Rev. B. 37, 785 (1988).            !
             ! densu = spin up density, densd = spin down density         !
             ! gdensu = spin up gradient, gdensd = spin down gradient     !
             ! ldensu = spin up laplacian, ldensd = spin down laplacian   !
             !------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN)  ::  densu,densd,gdensu(3),gdensd(3),ldensu,ldensd
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION :: g,rh,grh(3),lrh
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             DOUBLE PRECISION, PARAMETER :: a = 0.0490d0
             DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
             DOUBLE PRECISION :: CF,df2,dg21,dg22,dg23,dg24,gg2,ff2,grf(3),grg(3),lf,lg
             DOUBLE PRECISION :: rha,rhb,grha(3),grhb(3),lrha,lrhb
            
             rha = densu
             rhb = densd
             grha = gdensu
             grhb = gdensd
             lrha = ldensu
             lrhb = ldensd
             
             CF = (3.0d0/10.0d0)*(3.0d0*PI**2)**(2/3)

             ! the total density
             rh = rha + rhb
             ! rhe gamma function of eqn (13) in Phys. Rev. B. 37, 785 (1988)
             IF ( rh .GT. 1.0E-20 ) THEN
                     g = 2.0d0*(1.0d0 - (rha**2 + rhb**2)/rh**2 )
             ELSE
                     g = 2.0d0
             ENDIF
             
             !=============================================================================
             ! In what follows we use the expression (25) in Phys. Rev. B. 37, 785 (1988)
             ! to calculate the spinpolarized correlation potential
             !==============================================================================
             !----------------------------------
             ! rhe spin up part of the potential
             !----------------------------------
             
             df2 = F2pa(rh,rha,rhb,1.0d0) 
             ff2 = F2(g,rh)
            
             gg2 = G2(g,rh)

             dg21 = G2pa( rha**(8.0d0/3.0d0) + rhb**(8.0d0/3.0d0),rha,rhb, (8.0d0/3.0d0)*rha**(5.0d0/3.0d0) ) 

             ! Calculating the gradient and laplacian
             ! of the total charge-density
             lrh = lrha + lrhb
             grh = grha + grhb

             dg22 = G2pa( rh*lrh - DOT_PRODUCT(grh,grh),rha,rhb, lrh )
        
             dg23 = G2pa( rha*lrha + rhb*lrhb,rha,rhb,lrha)

             dg24 = G2pa( DOT_PRODUCT(grha,grha) + DOT_PRODUCT(grhb,grhb),rha,rhb,0.0d0)

             CALL gradG2(g,rh,grh,grg)
             lg = lapG2(g,rh,grh,lrh)

             ! here the expression (25) in Phys. Rev. B. 37, 785 (1988) is used
             ! for the spin up density:
               

             IF ( rh .GT. 1.0E-20 ) THEN
                Vc(1) = -a*( df2*rh + ff2 ) - (2.0d0**(5.0d0/3.0d0))*a*b*CF*( dg21 + (8.0d0/3.0d0)*gg2*rha**(5.0d0/3.0d0) )
                Vc(1) = Vc(1) - (a*b/4.0d0)*( rh*lg + 4.0d0*(DOT_PRODUCT(grg,grh) + gg2*lrh ) + dg22 )
                Vc(1) = Vc(1) - (a*b/36.0d0)*( 3.0d0*rha*lg + 4.0d0*( DOT_PRODUCT(grha,grg) + gg2*lrha ) + 3.0d0*dg23 + dg24 )
             ELSE
                Vc(1) = 0.0d0
             ENDIF
             !----------------------------------
             ! the spin down part of the potential
             !----------------------------------
             
             !----------------------------------------------------------
             ! switching up and down densities, gradients and laplacians
             !----------------------------------------------------------
             rha = densd
             rhb = densu
             !---------------
             grha = gdensd
             grhb = gdensu
             !---------------
             lrha = ldensd
             lrhb = ldensu
             !---------------


             df2 = F2pa(rh,rha,rhb,1.0d0) 

             dg21 = G2pa( rha**(8.0d0/3.0d0) + rhb**(8.0d0/3.0d0),rha,rhb, (8.0d0/3.0d0)*rha**(5.0d0/3.0d0) ) 

             dg22 = G2pa( rh*lrh - DOT_PRODUCT(grh,grh),rha,rhb, lrh )

             dg23 = G2pa( rha*lrha + rhb*lrhb,rha,rhb,lrha)

             dg24 = G2pa( DOT_PRODUCT(grha,grha) + DOT_PRODUCT(grhb,grhb),rha,rhb,0.0d0)

             ! here the expression (25) in Phys. Rev. B. 37, 785 (1988) is used
             ! for the spin down density:
              
             IF ( rh .GT. 1.0E-20 ) THEN
                Vc(2) = -a*( df2*rh + ff2 ) - (2.0d0**(5.0d0/3.0d0))*a*b*CF*( dg21 + (8.0d0/3.0d0)*gg2*rha**(5.0d0/3.0d0) )
                Vc(2) = Vc(2) - (a*b/4.0d0)*( rh*lg + 4.0d0*(DOT_PRODUCT(grg,grh) + gg2*lrh ) + dg22 )
                Vc(2) = Vc(2) - (a*b/36.0d0)*( 3.0d0*rha*lg + 4.0d0*( DOT_PRODUCT(grha,grg) + gg2*lrha ) + 3.0d0*dg23 + dg24 )
             ELSE
                     Vc(2) = 0.0d0
             ENDIF
             
        END SUBROUTINE VLYP

        SUBROUTINE VB3LYP(densu,densd,gdensu,gdensd,ldensu,ldensd,Vc)
             !-----------------------------------------------------------------------!
             ! This subroutine calculates Gradient corrected echange correlation     !
             ! potential B3LYP without 20 % of the non-local hartree fock exchange   !
             ! potential, 0.20*VXHF, which has to be added later, i.e:               !
             !  VB3LYP ---> VB3LYP + 0.20*VXHF.                                      !
             ! See Eqn(2) in J.Phys. Chem., Vol. 98, No. 45, 1994                    !
             ! densu = spin up density, densd = spin down density                    !
             ! gdensu = spin up gradient, gdensd = spin down gradient                !
             ! ldensu = spin up laplacian, ldensd = spin down laplacian              !
             !-----------------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN)  ::  densu,densd,gdensu(3),gdensd(3),ldensu,ldensd
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION :: VVLYP(2), VVXB88(2), VVWN(2),VXLSDA(2)
             DOUBLE PRECISION, PARAMETER :: a = 0.20d0
             DOUBLE PRECISION, PARAMETER :: b = 0.720d0
             DOUBLE PRECISION, PARAMETER :: c = 0.810d0

             CALL VWN(densu,densd,VVWN)

             CALL VXB88(densu,densd,gdensu,gdensd,VVXB88)

             CALL VLYP(densu,densd,gdensu,gdensd,ldensu,ldensd,VVLYP)

             CALL Vx(densu,densd,VXLSDA)
                
             Vc = (1.0d0 - a )*VXLSDA + b*VVXB88 + c*VVLYP + (1.0d0 - c)*VVWN

       END SUBROUTINE VB3LYP

        FUNCTION EB3LYP(densu,densd,gdensu,gdensd,ldensu,ldensd)
                ! Calculationg the B3LYP energy density acording to 
                ! Eqn (22) in Phys. Rev. B. 37, 785 (1988) + the rest, i.e Becke
                ! exchange ...., except 20 % of the exact exchange that has to
                ! be added later.
                IMPLICIT NONE
                DOUBLE PRECISION :: densu,densd,gdensu(3),gdensd(3),ldensu,ldensd
                DOUBLE PRECISION :: EB3LYP
                DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
                DOUBLE PRECISION, PARAMETER :: a = 0.0490d0
                DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
                DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                DOUBLE PRECISION :: CF,gam,twu,twd,tw,dens,LYPE,VWNE,LSDAE,XB88E
                DOUBLE PRECISION :: rs,zeta
                DOUBLE PRECISION, PARAMETER :: aa = 0.20d0
                DOUBLE PRECISION, PARAMETER :: bb = 0.720d0
                DOUBLE PRECISION, PARAMETER :: cc = 0.810d0
               
                dens = densu+densd
                IF ( dens .LT. 1.0E-20 ) THEN
                        EB3LYP = 0.0d0
                        RETURN
                ENDIF
                ! calculating zeta and rs
                IF ( dens .NE. 0.0d0 ) THEN
                     zeta = (densu-densd)/dens
                     rs = (3.0d0/(4.0d0*PI*dens))**(1.0d0/3.0d0)
                ELSE
                     zeta = 0.0d0
                     rs = 1.0E20
                ENDIF

                CF = (3.0d0/10.0d0)*(3.0d0*PI**2)**(2/3)

                gam = 2.0d0*( 1.0d0 - ( densu**2 + densd**2)/( (densu+densd)**2) )

                twu = (1.0d0/8.0d0)*( DOT_PRODUCT(gdensu,gdensu)/densu - ldensu)
                twd = (1.0d0/8.0d0)*( DOT_PRODUCT(gdensd,gdensd)/densd - ldensd)

                tw = twu + twd

                ! The LYP correlation energy density acording to Eqn (22) in Phys. Rev. B. 37, 785 (1988)
                ! is calculated here
                
                LYPE = -a*(gam/(1.0d0 + d*dens**(-1.0d0/3.0d0) ) )*( dens + 2.0d0*b*dens**(-5.0d0/3.0d0)*( CF*2**(2.0d0/3.0d0)*( densu**(8.0d0/3.0d0) + densd**(8.0d0/3.0d0) ) &
                & - dens*tw + (1.0d0/9.0d0)*(densu*twu + densd*twd) + (1.0d0/18.0d0)*(densu*ldensu+densd*ldensd) )*EXP(-c*dens**(-1.0d0/3.0d0) ) )

                ! The LSDA exchange energy density is calculated here:
                LSDAE = -(3.0d0/2.0d0)*( densd*((3.0d0/(8.0d0*PI))*densd)**(1.0d0/3.0d0) + densu*((3.0d0/(8.0d0*PI))*densu)**(1.0d0/3.0d0) )

                ! The correlation energy density of of Vosko, Wilk, and Nusair is calculated here
                VWNE = ( epsP(rs) + deps(rs,zeta) )*(densu+densd)

                ! The Becke exchange energy density is calculated here:
                XB88E = EXB88(densu,densd,gdensu,gdensd)

                ! The B3LYB energy density excluding 20% exact exchange is calculated here
                EB3LYP = (1.0d0 - aa )*LSDAE + bb*XB88E + cc*LYPE + (1.0d0 - cc)*VWNE
                
        END FUNCTION EB3LYP

        SUBROUTINE VXB88(rha,rhb,grha,grhb,Vxb)
             !-----------------------------------------------------------!
             ! This subroutine calculates the Becke exchange potential   !
             ! as given by the functional derivative with respect to     !
             ! the density of ( Ex - Ex(LDA) ) given by Eqn. (8) in      !
             ! Phys. Rev. A 38, 3098 (1988)                              !
             ! rha = spin up density, rhb = spin down density            !
             ! grha = spin up gradient, grhb = spin down gradient        !
             !-----------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN) :: rha,rhb,grha(3),grhb(3)
             DOUBLE PRECISION, INTENT(OUT) :: Vxb(2)
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: rh,xa,xb,asha,dasha,ashb,dashb
             DOUBLE PRECISION, PARAMETER :: bb = 0.00420d0
             
             IF ( rha .GT. 1.0E-20 ) THEN
                xa = sqrt(DOT_PRODUCT(grha,grha))/rha
             ELSE
                xa = 0.0d0
             ENDIF
             IF ( rhb .GT. 1.0E-20 ) THEN
                xb = sqrt(DOT_PRODUCT(grhb,grhb))/rhb
             ELSE
                     xb = 0.0d0
             ENDIF

             ! calculating of the inverse hyperbolic functions:
             ! asinh(xa) and asinh(xb)
             asha = LOG( xa + sqrt(xa**2 + 1) )
             ashb = LOG( xb + sqrt(xb**2 + 1) )

             ! calculating the derivatives of the 
             ! inverse hyperbolic functions 
             ! asinh(xa) and asinh(xb)
             dasha = 1.0d0/sqrt(xa**2 + 1)
             dashb = 1.0d0/sqrt(xb**2 + 1)

             Vxb(1) =( ((4.0d0/3.0d0)*bb*(rha**(1.0d0/3.0d0))*xa**2)/(1.0d0 + 6.0d0*bb*xa*asha ) )*(6.0d0*bb*xa*(asha + xa*dasha)/(1.0d0 + 6.0d0*bb*xa*asha ) - 1.0d0 )
             
             Vxb(2) =( ((4.0d0/3.0d0)*bb*(rhb**(1.0d0/3.0d0))*xb**2)/(1.0d0 + 6.0d0*bb*xb*ashb ) )*(6.0d0*bb*xb*(ashb + xb*dashb)/(1.0d0 + 6.0d0*bb*xb*ashb ) - 1.0d0 )

     END SUBROUTINE VXB88
      
     FUNCTION EXB88(rha,rhb,grha,grhb)
             !---------------------------------------------------------------!
             ! This subroutine calculates the Becke exchange energy density  !
             ! of ( Ex - Ex(LDA) ) given by Eqn. (8) in                      !
             ! Phys. Rev. A 38, 3098 (1988)                                  !
             ! rha = spin up density, rhb = spin down density                !
             ! grha = spin up gradient, grhb = spin down gradient            !
             !---------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION  :: rha,rhb,grha(3),grhb(3)
             DOUBLE PRECISION  :: EXB88
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: rh,xa,xb,asha,dasha,ashb,dashb
             DOUBLE PRECISION, PARAMETER :: bb = 0.00420d0
             
             xa = sqrt(DOT_PRODUCT(grha,grha))/rha
             xb = sqrt(DOT_PRODUCT(grhb,grhb))/rhb

             ! calculating of the inverse hyperbolic functions:
             ! asinh(xa) and asinh(xb)
             asha = LOG( xa + sqrt(xa**2 + 1) )
             ashb = LOG( xb + sqrt(xb**2 + 1) )

             EXB88 = -bb*( rha**(4.0d0/3.0d0)*xa**2/(1.0d0 + 6.0d0*bb*xa*asha ) + rhb**(4.0d0/3.0d0)*xb**2/(1.0d0 + 6.0d0*bb*xb*ashb ) )

     END FUNCTION EXB88

     SUBROUTINE VPBE(rha,rhb,lgrh,Vc)
             !----------------------------------------------------------------------------------!
             ! This subroutine calculates Gradient corrected exchange correlation               !
             ! potential of Perdew, Burke and Ernzerhof (PBE), Phys. Rev. Lett. 77, 3865 (1996) !
             ! rha = spin up density, rhb = spin down density                                   !
             ! lgrh = length of the total density gradient. The potential is calculated using   !
             ! finite differences rather than analytically in order to avoid  bugs.             !
             !----------------------------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN) :: rha,rhb,lgrh
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: g,rh,grh(3),lh
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             DOUBLE PRECISION ::  hh,rs,zeta
             DOUBLE PRECISION :: vp2,vm2,vp1,vm1,v0
            
             !print*,'(2)',lgrh

             IF ( rha + rhb .EQ. 0.0d0 ) THEN
                     Vc = 0.0d0
                     RETURN
             ENDIF

             IF ( rha .GT. 1.0E-20 ) THEN
                     hh = 0.0010d0*rha
             ELSE
                     hh = 0.0010d0
             ENDIF

             !f(0):
             rh = rha + rhb
             IF ( rh .GT. 1.0E-20 ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = (rha - rhb)/rh
             ELSE
                        rs = 1.0E20
                        zeta = 0.0d0
             ENDIF
             
             v0 = rh*( HPBE(rha,rhb,lgrh) + EUEG(rs,zeta) + EXPBE(rha,rhb,lgrh) )

             !f(-h):
             rh = (rha-hh) + rhb
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( (rha-hh) - rhb)/rh
             ELSE
                        rs = 1.0E20
                        zeta = 0.0d0
             ENDIF
          
             vm1 = rh*( HPBE(rha-hh,rhb,lgrh) + EUEG(rs,zeta)  + EXPBE(rha-hh,rhb,lgrh) )

             !f(-2*h):
             rh = (rha-2*hh) + rhb
             IF ( rh .GT. 1.0E-20  ) THEN
                rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                zeta = ( (rha-2*hh) - rhb)/rh
             ELSE
                     rs = 1.0E20
                     zeta = 0.0d0
             ENDIF
          
             vm2 = rh*( HPBE(rha-2*hh,rhb,lgrh) + EUEG(rs,zeta)  + EXPBE(rha-2*hh,rhb,lgrh) )

             !f(+h):
             rh = (rha+hh) + rhb
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( (rha+hh) - rhb)/rh
             ELSE
                        rs = 1.0E20
                        zeta = 0.0d0
             ENDIF
          
             vp1 = rh*( HPBE(rha+hh,rhb,lgrh) + EUEG(rs,zeta)  + EXPBE(rha+hh,rhb,lgrh) )

             !f(+2*h):
             rh = (rha+2*hh) + rhb
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( (rha+2*hh) - rhb)/rh
             ELSE
                        rs = 1.0E20
                        zeta = 0.0d0
             ENDIF
          
             vp2 = rh*( HPBE(rha+2*hh,rhb,lgrh) + EUEG(rs,zeta)  + EXPBE(rha+2*hh,rhb,lgrh) )

             ! Here we calculate the derivative of the exchange correlation  integrand
             ! with respect to the spin up density rha, using the finite
             ! difference approximation of O(h**4) see Koonin, p.6 Table 1.2

             Vc(1) = (1.0d0/(12.0d0*hh))*( vm2 - 8*Vm1 + 8*Vp1 - vp2 )

             !----------------------------
             !the spin down component:
             !----------------------------
             IF ( rhb .LT. 1.0d0 ) THEN
                        hh = 0.0010d0*rhb
             ELSE
                        hh = 0.0010d0
             ENDIF

             !f(-h):
             rh = rha + (rhb-hh)
             IF ( rh .GT. 1.0E-20  ) THEN
                rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                zeta = ( rha - (rhb-hh))/rh
             ELSE
                     rs = 1.0E20
                     zeta = 0.0d0
             ENDIF
          
             vm1 = rh*( HPBE(rha,rhb-hh,lgrh) + EUEG(rs,zeta)  + EXPBE(rha,rhb-hh,lgrh) )

             !f(-2*h):
             rh = rha + (rhb-2*hh)
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( rha - (rhb-2*hh) )/rh
             ELSE
                     rs = 1.0E20
                     zeta = 0.0d0
             ENDIF
             vm2 = rh*( HPBE(rha,rhb-2*hh,lgrh) + EUEG(rs,zeta)  + EXPBE(rha,rhb-2*hh,lgrh) )

             !f(+h):
             rh = rha + (rhb+hh)
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( rha - (rhb+hh) )/rh
             ELSE
                     rs = 1.0E20
                     zeta = 0.0d0
             ENDIF
          
             vp1 = rh*( HPBE(rha,rhb+hh,lgrh) + EUEG(rs,zeta)  + EXPBE(rha,rhb+hh,lgrh) )

             !f(+2*h):
             rh = rha + (rhb+2*hh)
             IF ( rh .GT. 1.0E-20  ) THEN
                        rs = (3.0d0/(4.0d0*rh))**(1.0d0/3.0d0)
                        zeta = ( rha - (rhb+2*hh))/rh
             ELSE
                     rs = 1.0E20
                     zeta = 0.0d0
             ENDIF
          
             vp2 = rh*( HPBE(rha,rhb+2*hh,lgrh) + EUEG(rs,zeta)  + EXPBE(rha,rhb+2*hh,lgrh) )

             ! Here we calculate the derivative of the exchange correlation  integrand
             ! with respect to the spin down density rhb, using the finite
             ! difference approximation of O(h**4) see Koonin, p.6 Table 1.2

             Vc(2) = (1.0d0/(12.0d0*hh))*( vm2 - 8*Vm1 + 8*Vp1 - vp2 )
             
      END SUBROUTINE VPBE


         SUBROUTINE Vx(densu,densd,Vxx)
             !-----------------------------------------------------------!
             ! This subroutine calculates the LSDA exchange  potential   !
             ! densu = spin up density, densd = spin down density        !
             !-----------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN) :: densu,densd
             DOUBLE PRECISION, INTENT(OUT) :: Vxx(2)
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: dens
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             
             dens = densu + densd 
             ! The exchange potential for spin up and down is calculated
             Vxx(1) = -2.0d0*((3.0d0/(8.0d0*PI))*densu)**(1.0d0/3.0d0)
             Vxx(2) = -2.0d0*((3.0d0/(8.0d0*PI))*densd)**(1.0d0/3.0d0)
         END SUBROUTINE Vx
!======================================================================================
! Here follows the auxilliary functions used to define the exchange correlation
! potential of Vosko, Wilk, and Nusair (VWM).
! potentials listed above. See pages 3206-3207 in Phys. Rev. B (R) 42, 3205 ! (1990)
!======================================================================================
        FUNCTION ffunk(zeta)
                !--------------------------------------------!
                ! Eqn (6),in PRB 42, 3205 (1990)             !
                !--------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: ffunk
                DOUBLE PRECISION :: zeta
                ffunk = ( (1.0d0 + zeta)**(4.0d0/3.0d0) + (1.0d0 - zeta)**(4.0d0/3.0d0) - 2.0d0 )/( 2.0d0**(4.0d0/3.0d0) - 2.0d0 ) 
        END FUNCTION ffunk

        FUNCTION ffunkp(zeta)
                !--------------------------------------------!
                ! the first derivtive of ffunk               !
                !--------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: ffunkp
                DOUBLE PRECISION :: zeta
                ffunkp = (4.0d0/3.0d0)*( (1.0d0 + zeta)**(1.0d0/3.0d0) - (1.0d0 - zeta)**(1.0d0/3.0d0) )/( 2.0d0**(4.0d0/3.0d0) - 2.0d0 ) 
        END FUNCTION ffunkp
        
        FUNCTION ffunkbis(zeta)
                ! the second derivative of ffunk
                IMPLICIT NONE
                DOUBLE PRECISION :: ffunkbis
                DOUBLE PRECISION :: zeta
                ffunkbis = (4.0d0/9.0d0)*( (1.0d0 + zeta)**(-2.0d0/3.0d0) + (1.0d0 - zeta)**(-2.0d0/3.0d0)  )/( 2.0d0**(4.0d0/3.0d0) - 2.0d0 ) 
        END FUNCTION ffunkbis

        FUNCTION epsP(rs)
                !--------------------------------------------!
                ! Eqn (11),     i=P,  in PRB 42, 3205 (1990) !
                !--------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: epsP
                DOUBLE PRECISION :: rs
                DOUBLE PRECISION, PARAMETER :: A =  0.03110d0
                DOUBLE PRECISION, PARAMETER :: B = -0.0480d0
                DOUBLE PRECISION, PARAMETER :: C =  0.00200d0
                DOUBLE PRECISION, PARAMETER :: D = -0.01160d0
                DOUBLE PRECISION, PARAMETER :: beta1 =  1.05290d0
                DOUBLE PRECISION, PARAMETER :: beta2 =  0.33340d0
                DOUBLE PRECISION, PARAMETER :: gama  = -0.14230d0

                IF ( rs .LT. 1.0d0 .AND. rs .NE. 0.0d0 ) THEN 
                        epsP = A*log(rs) + B + C*rs*log(rs) + D*rs
                ELSE
                        epsP = gama/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                ENDIF

                IF ( rs .EQ. 0.0d0 )  epsP  = -1.0E20
        END FUNCTION epsP
        
        FUNCTION epsF(rs)
                !--------------------------------------------!
                ! Eqn (11),     i=F,  in PRB 42, 3205 (1990) !
                !--------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: epsF
                DOUBLE PRECISION :: rs
                DOUBLE PRECISION, PARAMETER :: A =  0.015550d0
                DOUBLE PRECISION, PARAMETER :: B = -0.02690d0
                DOUBLE PRECISION, PARAMETER :: C =  0.00070d0
                DOUBLE PRECISION, PARAMETER :: D = -0.0048
                DOUBLE PRECISION, PARAMETER :: beta1 =  1.39810d0
                DOUBLE PRECISION, PARAMETER :: beta2 =  0.26110d0
                DOUBLE PRECISION, PARAMETER :: gama  = -0.08430d0

                IF ( rs .LT. 1.0d0 .AND. rs .NE. 0.0d0 ) THEN 
                        epsF = A*log(rs) + B + C*rs*log(rs) + D*rs
                ELSE
                        epsF = gama/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                ENDIF

                IF ( rs .EQ. 0.0d0 )  epsF  = -1.0E20

        END FUNCTION epsF
       
        FUNCTION alphac(rs)
                !---------------------------------!
                ! Eqn (7), in PRB 42, 3205 (1990) !
                !---------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: alphac
                DOUBLE PRECISION :: rs
                DOUBLE PRECISION, PARAMETER :: A = -0.0337737278807790d0
                DOUBLE PRECISION, PARAMETER :: b =  1.131070d0
                DOUBLE PRECISION, PARAMETER :: c =  13.00450d0
                DOUBLE PRECISION, PARAMETER :: x0 = -0.004758400d0
                DOUBLE PRECISION :: Q,Xrs,Xx0

                Q = sqrt(4.0d0*c - b**2.0d0 )
                Xrs = rs + b*sqrt(rs)+c
                Xx0 = x0**2 + b*dabs(x0)+c

                alphac = log(rs/Xrs) + (2.0d0*b/Q)*atan(Q/(2*sqrt(rs) + b ))
                alphac = alphac - (b*x0/Xx0)*( log( ( (sqrt(rs) - x0)**2.0d0 )/Xrs ) + (2.0d0*(b+2.0d0*x0)/Q)*atan(Q/(2*sqrt(rs) + b )) )
                alphac = alphac*A
       END FUNCTION alphac

                FUNCTION betaf(rs)
                        !---------------------------------!
                        ! Eqn (4), in PRB 42, 3205 (1990) !
                        !---------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: betaf
                        DOUBLE PRECISION :: rs
                        betaf = ffunkbis(0.0d0)*( epsF(rs) - epsP(rs) )/alphac(rs) - 1.0d0
                END FUNCTION betaf
                
                FUNCTION deps(rs,zeta)
                        !---------------------------------!
                        ! Eqn (2), in PRB 42, 3205 (1990) !
                        !---------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: deps
                        DOUBLE PRECISION :: rs,zeta
                        deps = alphac(rs)*(ffunk(zeta)/ffunkbis(0.0d0))*( 1.0d0 + betaf(rs)*zeta**4.0d0 )
                END FUNCTION deps

                FUNCTION muP(rs)
                        !--------------------------------------------!
                        ! Eqn (12,a,b), i=P,  in PRB 42, 3205 (1990) !
                        !--------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: muP
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION, PARAMETER :: A =  0.03110d0
                        DOUBLE PRECISION, PARAMETER :: B = -0.0480d0
                        DOUBLE PRECISION, PARAMETER :: C =  0.00200d0
                        DOUBLE PRECISION, PARAMETER :: D = -0.01160d0
                        DOUBLE PRECISION, PARAMETER :: beta1 =  1.05290d0
                        DOUBLE PRECISION, PARAMETER :: beta2 =  0.33340d0
                        DOUBLE PRECISION, PARAMETER :: gama  = -0.14230d0

                        IF ( rs .LT. 1.0d0 ) THEN 
                                muP = A*log(rs) + (B - (1.0d0/3.0d0)*A ) + (2.0d0/3.0d0)*C*rs*log(rs) + (1.0d0/3.0d0)*(2.0d0*D -C)*rs
                        ELSE
                                mup = epsP(1.0d0 + (7.0d0/6.0d0)*beta1*sqrt(rs) + (4.0d0/3.0d0)*beta2*rs )/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                        ENDIF
                END FUNCTION muP
                
                FUNCTION muF(rs)
                        !--------------------------------------------!
                        ! Eqn (12,a,b), i=F,  in PRB 42, 3205 (1990) !
                        !--------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: muF
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION, PARAMETER :: A =  0.015550d0
                        DOUBLE PRECISION, PARAMETER :: B = -0.02690d0
                        DOUBLE PRECISION, PARAMETER :: C =  0.00070d0
                        DOUBLE PRECISION, PARAMETER :: D = -0.0048
                        DOUBLE PRECISION, PARAMETER :: beta1 =  1.39810d0
                        DOUBLE PRECISION, PARAMETER :: beta2 =  0.26110d0
                        DOUBLE PRECISION, PARAMETER :: gama  = -0.08430d0
                        IF ( rs .LT. 1.0d0 ) THEN 
                                muF = A*log(rs) + (B - (1.0d0/3.0d0)*A ) + (2.0d0/3.0d0)*C*rs*log(rs) + (1.0d0/3.0d0)*(2.0d0*D -C)*rs
                        ELSE
                                muF = epsF(1.0d0 + (7.0d0/6.0d0)*beta1*sqrt(rs) + (4.0d0/3.0d0)*beta2*rs )/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                        ENDIF
                END FUNCTION muF

                FUNCTION rdadr(rs)
                        !-------------------------------------------------------------------!
                        ! Eqn (9),in PRB 42, 3205 (1990)                                    !
                        ! the derivative of alphac times rc, i.e rdadr(rs) = rc*alphac'(rc) !
                        !-------------------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: rdadr
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION, PARAMETER :: A = -0.0337737278807790d0
                        DOUBLE PRECISION, PARAMETER :: b =  1.131070d0
                        DOUBLE PRECISION, PARAMETER :: c =  13.00450d0
                        DOUBLE PRECISION, PARAMETER :: x0 = -0.004758400d0
                        DOUBLE PRECISION :: b1,b2,b3

                        b1 = (b*x0 - c)/(c*x0)
                        b2 = (x0 - b)/(c*x0)
                        b3 = - 1.0d0/(c*x0)

                        rdadr = A*( 1.0d0 + b1*sqrt(rs) )/( 1.0d0 + b1*sqrt(rs) + b2*rs + b3*rs**(1.50d0) )

                END FUNCTION rdadr

                FUNCTION rdbdr(rs)
                        !-------------------------------------------------------------------!
                        ! Eqn (10),in PRB 42, 3205 (1990)                                   !
                        ! rdbdr(rs) = rs*betaf'(rs)                                         !
                        !-------------------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: rdbdr
                        DOUBLE PRECISION :: rs
                        rdbdr = (ffunkbis(0.0d0)/alphac(rs))*3.0d0*( epsF(rs) - muF(rs) - ( epsP(rs) - muP(rs) ) )
                        rdbdr = rdbdr - ( ffunkbis(0.0d0)*(epsF(rs)-epsP(rs))/(alphac(rs)**2.0d0) )*rdadr(rs)
                END FUNCTION rdbdr
                !=================================================================================================
                ! Here follows the auxilliary functions used to calculate the
                ! LYP correlation potential Phys. Rev. B 37, 785 (1988)  ! a ! =0.04918 b=0.132
                !================================================================================================

                FUNCTION G2(g,rh)
                        !------------------------------------!
                        ! Eqn (26), in PRB 37, 785 (1988)    !
                        !------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: G2
                        DOUBLE PRECISION :: g,rh
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        IF ( rh .GT. 0.0d0 ) THEN
                                G2 = F2(g,rh)*(rh**(-5.0d0/3.0d0))*EXP(-c*rh**(-1.0d0/3.0d0))
                        ELSE
                                G2 = 0.0d0
                        ENDIF
                END FUNCTION G2
                
                FUNCTION F2(g,rh)
                        !------------------------------------!
                        ! Eqn (26), in PRB 37, 785 (1988)    !
                        !------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: F2
                        DOUBLE PRECISION :: g,rh
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        IF ( rh .GT. 0.0d0 ) THEN
                                F2 = g/(1.0d0 + d*rh**(-1.0d0/3.0d0) )
                        ELSE
                                F2 = 0.0d0
                        ENDIF
                        
                END FUNCTION F2

                SUBROUTINE gradrho(BAS,P,r,grho)
                        ! Calculates the gradient of the charge density with
                        ! respect to the electronic coordinate, r.
                        USE datatypemodule
                        IMPLICIT NONE 
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS,3),grad(3)
                        INTEGER :: I
                        
                        DO I=1,BAS%NBAS
                                V1(I) = basfunkval(BAS%PSI(I),r)
                                CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                V2(I,:) = grad
                         ENDDO

                         grho(1) = 2.0d0*DOT_PRODUCT(V2(:,1),MATMUL(P,V1))
                         grho(2) = 2.0d0*DOT_PRODUCT(V2(:,2),MATMUL(P,V1))
                         grho(3) = 2.0d0*DOT_PRODUCT(V2(:,3),MATMUL(P,V1))
                 END SUBROUTINE gradrho
                                 
                SUBROUTINE nucgradrho(NATOMS,BAS,P,gradS,r,grho)
                        ! Calculates the nuclear gradient of the charge density.
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(NATOMS,BAS%NBAS,3),grad(3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J
                        
                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) V1(I) = basfunkval(BAS%PSI(I),r)
                                        IF ( BAS%PSI(I)%ATYPE .EQ. J ) THEN
                                                CALL ngbasfunkval(BAS%PSI(I),r,grad)
                                                V2(J,I,:) = grad
                                        ELSE
                                                V2(J,I,:) = 0.0d0
                                        ENDIF
                              ENDDO
                              ! Calculating the nuclear gradient of the density
                              ! matrix using Eqn (9) , page 125 in black note
                              ! book, i.e grad(P) = -P*grad(S)*P
                              gradP(J,1,:,:) = -MATMUL(P,MATMUL(gradS(J,1,:,:),P))
                              gradP(J,2,:,:) = -MATMUL(P,MATMUL(gradS(J,2,:,:),P))
                              gradP(J,3,:,:) = -MATMUL(P,MATMUL(gradS(J,3,:,:),P))
                        ENDDO
                        
                        ! See derivation of Eqn (6) and Eqn(6) p. 124 in the black note book.
                        DO J=1,NATOMS
                                grho(J,1) = DOT_PRODUCT(V1,MATMUL(gradP(J,1,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,1),MATMUL(P,V1))
                                grho(J,2) = DOT_PRODUCT(V1,MATMUL(gradP(J,2,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,2),MATMUL(P,V1))
                                grho(J,3) = DOT_PRODUCT(V1,MATMUL(gradP(J,3,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,3),MATMUL(P,V1))
                        ENDDO
                 END SUBROUTINE nucgradrho
                
                SUBROUTINE nucgradgradrho(NATOMS,BAS,P,gradS,r,grho)
                        ! Calculates the nuclear gradient of the gradient of the charge density.
                        ! the output corresponding to the nuclear gradient with
                        ! respect to atom J, is organized as follows:
                        ! grho(J,1,1) = d^2(dens)/(dRxdx), grho(J,1,2) = d^2(dens)/(dRydx), grho(J,1,3) = d^2(dens)/(dRzdx)
                        ! grho(J,2,1) = d^2(dens)/(dRxdy), grho(J,2,2) = d^2(dens)/(dRydy), grho(J,2,3) = d^2(dens)/(dRzdy)
                        ! grho(J,3,1) = d^2(dens)/(dRxdz), grho(J,3,2) = d^2(dens)/(dRydz), grho(J,3,3) = d^2(dens)/(dRzdz)
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS,3),V3(NATOMS,BAS%NBAS,3),V4(NATOMS,BAS%NBAS,3,3)
                        DOUBLE PRECISION :: grad(3),ggrad(3,3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J,n,m
                        
                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) THEN
                                                V1(I) = basfunkval(BAS%PSI(I),r)
                                                CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                V2(BAS%NBAS,:) = grad
                                        ENDIF
                                        IF ( BAS%PSI(I)%ATYPE .EQ. J ) THEN
                                                CALL ngbasfunkval(BAS%PSI(I),r,grad)
                                                V3(J,I,:) = grad
                                                CALL nggradbasfunkval(BAS%PSI(I),r,ggrad)
                                                V4(J,I,:,:) = ggrad
                                        ELSE
                                                V3(J,I,:) = 0.0d0
                                                V4(J,I,:,:) = 0.0d0
                                        ENDIF
                              ENDDO
                              ! Calculating the nuclear gradient of the density
                              ! matrix using Eqn (9) , page 125 in black note
                              ! book, i.e grad(P) = -P*grad(S)*P
                              gradP(J,1,:,:) = -MATMUL(P,MATMUL(gradS(J,1,:,:),P))
                              gradP(J,2,:,:) = -MATMUL(P,MATMUL(gradS(J,2,:,:),P))
                              gradP(J,3,:,:) = -MATMUL(P,MATMUL(gradS(J,3,:,:),P))
                        ENDDO
                        
                        ! See derivation of Eqn (7) and Eqn(7) p. 124 in the black note book.
                        DO J=1,NATOMS
                                DO n=1,3
                                        DO m=1,3
                                                grho(J,n,m) = 2.0d0*DOT_PRODUCT(V2(:,n),MATMUL(gradP(J,m,:,:),V1)) 
                                                grho(J,n,m) = grho(J,n,m)+ 2.0d0*( DOT_PRODUCT(V4(J,:,n,m),MATMUL(P,V1)) + DOT_PRODUCT(V2(:,n),MATMUL(P,V3(J,:,m)) ) )
                                        ENDDO
                                ENDDO
                        ENDDO
                 END SUBROUTINE nucgradgradrho
                
                SUBROUTINE nucgradlaprho(NATOMS,BAS,P,gradS,r,grho)
                        ! Calculates the nuclear gradient of the Laplacian of the charge density.
                        ! the output corresponding to the nuclear gradient with
                        ! respect to atom J, is organized as follows:
                        ! grho(J,1) = d(laplacedens)/dRx, grho(J,2) = d(laplacedens)/dRy, grho(J,3) = d(laplacedens)/dRz
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval,laplacebasfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS,3),V3(NATOMS,BAS%NBAS,3),V4(NATOMS,BAS%NBAS,3,3),V5(BAS%NBAS),V6(NATOMS,BAS%NBAS,3)
                        DOUBLE PRECISION :: grad(3),ggrad(3,3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J,n,m
                        
                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) THEN
                                                ! The value of the basis function 
                                                V1(I) = basfunkval(BAS%PSI(I),r)
                                                ! the electron coordinate gradien of the basis function
                                                CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                V2(BAS%NBAS,:) = grad
                                                ! the electron coordinate laplacian of the basis function
                                                V5(I) = laplacebasfunkval(BAS%PSI(I),r)
                                        ENDIF
                                        IF ( BAS%PSI(I)%ATYPE .EQ. J ) THEN
                                                ! the nuclear coordinate gradient of the basis function
                                                CALL ngbasfunkval(BAS%PSI(I),r,grad)
                                                V3(J,I,:) = grad
                                                ! the nuclear gradient of the electron gradient of the basis function
                                                CALL nggradbasfunkval(BAS%PSI(I),r,ggrad)
                                                V4(J,I,:,:) = ggrad
                                                ! the nuclear gradient of the laplacian of the basis function
                                                CALL nglaplacebasfunkval(BAS%PSI(I),r,grad)
                                                V6(J,I,:) = grad
                                        ELSE
                                                V3(J,I,:) = 0.0d0
                                                V4(J,I,:,:) = 0.0d0
                                                V6(J,I,:) = 0.0d0
                                        ENDIF
                              ENDDO
                              ! Calculating the nuclear gradient of the density
                              ! matrix using Eqn (9) , page 125 in black note
                              ! book, i.e grad(P) = -P*grad(S)*P
                              gradP(J,1,:,:) = -MATMUL(P,MATMUL(gradS(J,1,:,:),P))
                              gradP(J,2,:,:) = -MATMUL(P,MATMUL(gradS(J,2,:,:),P))
                              gradP(J,3,:,:) = -MATMUL(P,MATMUL(gradS(J,3,:,:),P))
                        ENDDO
                        
                        ! Here everything is put together according to eqn
                        ! (8-10), p. 125 in the black note book.
                        DO J=1,NATOMS
                              DO m=1,3
                                grho(J,m) = 2.0d0*DOT_PRODUCT(V5,MATMUL(gradP(J,m,:,:),V1)) 
                                DO n=1,3
                                    grho(J,m) = grho(J,m) + 2.0d0*DOT_PRODUCT(V2(:,n),MATMUL(gradP(J,m,:,:),V2(:,n))) 
                                ENDDO
                                DO n=1,3
                                    grho(J,m) = grho(J,m) + 4.0d0*DOT_PRODUCT(V4(J,:,n,m),MATMUL(P,V2(:,n)))
                                ENDDO
                                grho(J,m) = grho(J,m) + 2.0d0*( DOT_PRODUCT(V6(J,:,m),MATMUL(P,V1)) + DOT_PRODUCT(V5,MATMUL(P,V3(J,:,m)) ) )
                             ENDDO
                        ENDDO
                 END SUBROUTINE nucgradlaprho
                
                 FUNCTION laprho(BAS,P,r)
                        ! Calculates the laplacian of the charge density with
                        ! respect to the electronic coordinate, r.
                        USE datatypemodule
                        IMPLICIT NONE
                        DOUBLE PRECISION :: laprho 
                        TYPE(BASIS) :: BAS
                        DOUBLE PRECISION :: P(BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION :: r(3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval,laplacebasfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(BAS%NBAS,3),V3(BAS%NBAS),grad(3)
                        INTEGER :: I
                        
                        DO I=1,BAS%NBAS
                                V1(I) = basfunkval(BAS%PSI(I),r)
                                V3(I) = laplacebasfunkval(BAS%PSI(I),r)
                                CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                V2(I,:) = grad
                         ENDDO
                         laprho = 0.0d0
                         DO I=1,3
                                laprho = laprho + DOT_PRODUCT(V2(:,I),MATMUL(P,V2(:,I))) 
                         ENDDO
                         laprho = laprho + DOT_PRODUCT(V3,MATMUL(P,V1))
                         laprho = 2.0d0*laprho
                 END FUNCTION laprho

                 FUNCTION F2pa(rh,rha,rhb,drhdrha)
                        ! Calculating the derivative of F2 with 
                        ! respect to the spin-up density, rha.
                        ! drhdrha = derivative of total density with respect to 
                        ! spin-up density. rh = total density
                        ! rha = spin-up density, rhb = spin-down density
                        IMPLICIT NONE
                        DOUBLE PRECISION :: F2pa
                        DOUBLE PRECISION :: rh,rha,rhb,drhdrha
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        DOUBLE PRECISION :: rt,g
                        rt = rha + rhb
                        IF ( rt .GT. 1.0E-20 ) THEN
                                g = 2.0d0*( 1.0d0 - ( rha**2 + rhb**2 )/rt )
                        ELSE
                                g = 2.0d0
                        ENDIF
                        IF ( rh .GT. 0.0d0 ) THEN
                                F2pa = (1.0d0/(1.0d0 + d*rh**(-1.0d0/3.0d0) ))*( (4.0d0*rhb/(rt**3))*(rhb-rha) + (d/3.0d0)*F2(g,rh)*drhdrha*rh**(-4.0d0/3.0d0) )
                        ELSE
                                F2pa = 0.0d0
                        ENDIF
                        !print*,F2pa
                END FUNCTION F2pa

                 FUNCTION G2pa(rh,rha,rhb,drhdrha)
                        ! Calculating the derivative of G2 with 
                        ! respect to the spin-up density, rha.
                        ! drhdrha = derivative of total density with respect to 
                        ! spin-up density. rh = total density
                        ! rha = spin-up density, rhb = spin-down density
                        IMPLICIT NONE
                        DOUBLE PRECISION :: G2pa
                        DOUBLE PRECISION :: rh,rha,rhb,drhdrha
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        DOUBLE PRECISION :: rt,g
                        rt = rha + rhb
                        IF ( rt .GT. 1.0E-20 ) THEN
                                g = 2.0d0*( 1.0d0 - ( rha**2 + rhb**2 )/rt )
                        ELSE
                                g = 2.0d0
                        ENDIF
                       
                        IF ( rh .GT. 1.0E-20  ) THEN
                                G2pa = ( F2pa(rh,rha,rhb,drhdrha)/F2(g,rh) + ( (c/3.0d0)*rh**(-4.0d0/3.0d0) - (5.0d0/(3.0d0*rh)) )*drhdrha )*G2(g,rh)
                        ELSE
                                G2pa = 0.0d0
                        ENDIF
                
                 END FUNCTION G2pa

                SUBROUTINE gradF2(g,rh,graddens,grad)
                        ! calculates the gradient of the F2 function
                        ! graddens = gradient of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION, INTENT(IN) :: g,graddens(3),rh
                        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        grad = graddens*( (2.0d0/rh)*(1.0d0-(g/2.0d0) - 1.0d0/rh ) + F2(g,rh)*(d/3.0d0)**(-4.0d0/3.0d0) )*(1.0d0/(1.0d0 + d*rh**(-1.0d0/3.0d0) ) )
                END SUBROUTINE gradF2

                FUNCTION lapF2(g,rh,graddens,lapdens)
                        ! Calculates the laplacian of the function F2. 
                        ! lapdens = laplacian of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION :: lapF2
                        DOUBLE PRECISION :: g,rh,graddens(3),lapdens
                        DOUBLE PRECISION :: grad(3),DF22,Drh2
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        DOUBLE PRECISION :: fa
                        
                        fa = (1.0d0 + d*rh**(-1.0d0/3.0d0) )
                        Drh2 = DOT_PRODUCT(graddens,graddens)
                        CALL gradF2(g,rh,graddens,grad)
                        DF22 =  DOT_PRODUCT(grad,grad)

                        lapF2 = DF22*( (d/3.0d0)*rh**(-4.0d0/3.0d0) + (lapdens/Drh2)*fa )
                        lapF2 = lapF2 + (1.0d0/fa)*Drh2*( 8.0d0/rh**3 - (6.0d0/rh**2)*(1.0d0-g/2.0d0) - (4.0d0/9.0d0)*F2(g,rh)*d*rh**(-7.0d0/3.0d0) )
                        lapF2 = lapF2 + (1.0d0/fa)*(d/3.0d0)*DOT_PRODUCT(graddens,grad)*rh**(-4.0d0/3.0d0) 
                END FUNCTION lapF2

                SUBROUTINE gradG2(g,rh,graddens,grad)
                        ! calculates the gradient of the G2 function
                        ! graddens = gradient of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION, INTENT(IN) :: g,graddens(3),rh
                        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        DOUBLE PRECISION :: gf(3),gg

                        CALL gradF2(g,rh,graddens,gf)

                        gg = G2(g,rh)
                        
                        grad = (gf/F2(g,rh))*gg + ( (c/3.0d0)*rh**(-4.0d0/3.0d0) - (5.0d0/(3.0d0*rh)) )*gg*graddens
                END SUBROUTINE gradG2

                FUNCTION lapG2(g,rh,graddens,lapdens)
                        ! Calculates the laplacian of the function G2. 
                        ! lapdens = laplacian of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION :: lapG2
                        DOUBLE PRECISION :: g,rh,graddens(3),lapdens
                        DOUBLE PRECISION :: gradf(3),gradg(3),DF22,Drh2
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        DOUBLE PRECISION :: fa,fb,fc,ff,gg
                        
                        fa =  5.0d0/(3.0d0*rh**2) -(4.0d0*c/9.0d0)*rh**(-7.0d0/3.0d0)
                        fb = (c/3.0d0)*rh**(-4.0d0/3.0d0) - 5.0d0/(3.0d0*rh)

                        ff = F2(g,rh)
                        gg = G2(g,rh)

                        Drh2 = DOT_PRODUCT(graddens,graddens)
                        
                        CALL gradF2(g,rh,graddens,gradf)
                        CALL gradG2(g,rh,graddens,gradg)

                        DF22 =  DOT_PRODUCT(gradf,gradf)
                        
                        lapG2 = (lapF2(g,rh,graddens,lapdens)/ff)*gg - (DF22/(ff**2))*gg + DOT_PRODUCT(gradf,gradg)/ff

                        lapG2 = lapG2 + fa*gg*Drh2 + fb*( DOT_PRODUCT(gradg,graddens) + gg*lapdens )

                END FUNCTION lapG2
                !===================================================================
                ! Here we define som auxilliary functions used to define the
                ! VPBE exchange correlation potential 
                !==================================================================

                FUNCTION EUEG(rs,zeta)
                        !---------------------------------------------------------!
                        ! The LDA correlation energy of the uniform electron gas  !
                        ! as given by Eqn (1) in Phys. Rev B, 42, 3205 (1990)     !
                        !---------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: EUEG
                        DOUBLE PRECISION :: rs,zeta
                        EUEG = epsP(rs) + deps(rs,zeta)
                END FUNCTION EUEG

                FUNCTION HPBE(rha,rhb,lgrh)
                        !-----------------------------------------------------------!
                        ! the function H as given by Eqn (7) and (8) in             !
                        ! Phys. Rev. Lett. 77, 3865 (1996)                          !
                        ! rha = spin up density, rhb = spin down density            !
                        ! lgrh = the length of the gradient of the total            !
                        ! density, i,e lgrh = sqrt(DOT_PRODUCT(grad,grad))          ! 
                        !-----------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: HPBE
                        DOUBLE PRECISION :: rha,rh,rhb,lgrh
                        DOUBLE PRECISION :: rs,zeta,A,PHH,ks,kf,t,n
                        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
                        DOUBLE PRECISION, PARAMETER :: bb = 0.0667250d0
                        DOUBLE PRECISION, PARAMETER :: gg = 0.0310910d0
                      
                        rh = rha+rhb
                        
                        IF ( rh .GT. 1.0E-20 ) THEN
                                zeta = (rha - rhb)/(rha+rhb)
                                rs = (3.0d0/(4.0d0*(rha+rhb)))**(1.0d0/3.0d0)
                        ELSE
                                zeta = 0.0d0
                                rs = 1.0E20
                        ENDIF

                        PHH = 0.50d0*( (1+zeta)**(2.0d0/3.0d0) + (1-zeta)**(2.0d0/3.0d0) )

                        kf = (1.0d0/rs)*(9.0d0*PI/4.0d0)**(1.0d0/3.0d0)

                        ks = sqrt(4.0d0*kf/PI)

                        n = 3.0d0/(4.0d0*PI*rs**3)

                        t = lgrh/(2.0d0*PHH*ks*n)
                        
                        IF (  EXP(-EUEG(rs,zeta)/(gg*PHH**3)) .NE. 1.0d0   ) THEN
                                A = (bb/gg)/( EXP(-EUEG(rs,zeta)/(gg*PHH**3)) - 1.0d0 )
                        ELSE
                                A = 1.0E20
                        ENDIF

                        HPBE = gg*(PHH**3)*LOG( 1.0d0 + (bb/gg)*(t**2)*(1.0d0 + A*t**2 )/( 1.0d0 + A*t**2 + (A*t**2)**2 ) )
                        
                END FUNCTION HPBE

                FUNCTION EXPBE(rha,rhb,lgrh)
                        !-----------------------------------------------------------!
                        ! the integrand/n   as given by Eqn (10) in                 !
                        ! Phys. Rev. Lett. 77, 3865 (1996), i.e                     !
                        ! epsilon_x(n)*F_X(s)                                       !
                        ! rha = spin up density, rhb = spin down density            !
                        ! lgrh = the length of the gradient of the total            !
                        ! density, i,e lgrh = sqrt(DOT_PRODUCT(grad,grad))          ! 
                        !-----------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: EXPBE
                        DOUBLE PRECISION :: rha,rhb,lgrh
                        DOUBLE PRECISION :: EXUEG,kf,s,rs,n,FXX
                        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
                        DOUBLE PRECISION, PARAMETER :: mm = 0.219510d0
                        DOUBLE PRECISION, PARAMETER :: kk = 1.245 ! Here we actually choose the kappa value of the revPBE ( PRL 80, 890 (1998) )
                       
                        IF ( rha+rhb .GT. 1.0E-20 ) THEN
                                rs = (3.0d0/(4.0d0*(rha+rhb)))**(1.0d0/3.0d0)
                        ELSE
                                rs = 1.0E20
                        ENDIF

                        kf = (1.0d0/rs)*(9.0d0*PI/4.0d0)**(1.0d0/3.0d0)

                        n = 3.0d0/(4.0d0*PI*rs**3)

                        s = lgrh/(2.0d0*kf*n)
                        
                        EXPBE = -3.0d0*kf/(4.0d0*PI)
                        
                        FXX = 1.0d0 + kk - kk/(1.0d0 + mm*(s**2)/kk )
                        
                        EXPBE = EXPBE*FXX
                        
               END FUNCTION EXPBE

END MODULE exchcorrmodule
