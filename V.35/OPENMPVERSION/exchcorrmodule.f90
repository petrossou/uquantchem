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

             IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
             IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0

             ! The spin up component of the correlation potential
             Vc(1) = muP(rs) + deps(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betaf(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdr(rs)*zeta**(4.0d0) ) )

             Vc(1) = Vc(1) + (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betaf(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betaf(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 - zeta )
             

             ! The spin down component of the correlation potential
             Vc(2) = muP(rs) + deps(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betaf(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdr(rs)*zeta**(4.0d0) ) )

             Vc(2) = Vc(2) - (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betaf(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betaf(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 + zeta )
                

        END SUBROUTINE VWN
        
        SUBROUTINE VWNC(densu,densd,Vc)
             !-------------------------------------------------------------------!
             ! This subroutine calculates the LSDA correlation                   !
             ! potential of Vosko, Wilk, and Nusair (VWM) at the point r         !
             ! using the interpolation refered to as VWN3 for the ferro magnetic !
             ! corr energy epsF and the paramagnetic corr energy epsP this fitt  ! 
             ! is different from the one used in the VWN subroutine above, which !
             ! uses the form given by eqn (11) in Phys. Rev. B. 42, 3205 (1990). ! 
             ! and the parameters given by table II in the same paper. Here      !
             ! instead we use the parametrization given by Eqn. 4.4 in           !
             ! Can. J. Phys 58, 1200  (1980), or Eqn. (39) Comp. Phys. Comm. 66, !
             ! 383 (1991)                                                        !
             ! densu = spin up density, densd = spin down density                !
             !-------------------------------------------------------------------!
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

             IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
             IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
             
             ! The spin up component of the correlation potential
             Vc(1) = muPC(rs) + depsC(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betafC(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdrC(rs)*zeta**(4.0d0) ) )

             Vc(1) = Vc(1) + (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betafC(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betafC(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 - zeta )
             

             ! The spin down component of the correlation potential
             Vc(2) = muPC(rs) + depsC(rs,zeta) + (ffunk(zeta)/ffunkbis(0.0d0))*( -(1.0d0/3.0d0)*rdadr(rs)*(1.0d0 + betafC(rs)*zeta**(4.0d0) ) + &
             & alphac(rs)*( -(1.0d0/3.0d0)*rdbdrC(rs)*zeta**(4.0d0) ) )

             Vc(2) = Vc(2) - (alphac(rs)/ffunkbis(0.0d0))*( 4.0d0*betafC(rs)*ffunk(zeta)*zeta**(3.0d0) + (1.0d0 + betafC(rs)*zeta**(4.0d0) )*ffunkp(zeta) )*(1.0d0 + zeta )
                

        END SUBROUTINE VWNC
        
        SUBROUTINE VLYP(rha,rhb,grha,grhb,lrha,lrhb,Vc)
             !------------------------------------------------------------!
             ! This subroutine calculates Gradient corrected correlation  !
             ! potential of Lee, Yang and Parr (LYP) at the point r       !
             ! For details see   Phys. Rev. B. 37, 785 (1988).            !
             ! densu = spin up density, densd = spin down density         !
             ! gdensu = spin up gradient, gdensd = spin down gradient     !
             ! ldensu = spin up laplacian, ldensd = spin down laplacian   !
             !------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN)  :: rha,rhb,grha(3),grhb(3),lrha,lrhb
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION :: g,rh,grh(3),lrh
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             DOUBLE PRECISION, PARAMETER :: a = 0.04918d0
             DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
             DOUBLE PRECISION :: CF,df2,dg21,dg22,dg23,dg24,gg2,ff2,grf(3),grg(3),lf,lg,GPRIME,FPRIME
            
             
             CF = (3.0d0/10.0d0)*(3.0d0*PI**2)**(2.0d0/3.0d0)

             ! the total density
             rh = rha + rhb
             IF ( rh .LT. 1.0E-20 ) THEN
                     Vc = 0.0d0
                     RETURN
             ENDIF
             ! rhe gamma function of eqn (13) in Phys. Rev. B. 37, 785 (1988)
             IF ( rh .GT. 1.0E-20 ) THEN
                     g = 2.0d0*(1.0d0 - (rha**2 + rhb**2)/rh**2 )
             ELSE
                     g = 1.0d0
             ENDIF
             
             !=============================================================================
             ! In what follows we use the expression (25) in Phys. Rev. B. 37, 785 (1988)
             ! to calculate the spinpolarized correlation potential
             !==============================================================================
             !----------------------------------
             ! rhe spin up part of the potential
             !----------------------------------
             
             ff2 = F2(g,rh)
             gg2 = G2(g,rh)
             

             ! Calculating the gradient and laplacian
             ! of the total charge-density
             lrh = lrha + lrhb
             grh = grha + grhb

             FPRIME = F2pa(rh,rha,rhb,1.0d0) 
             GPRIME = G2pa( rh,rha,rhb, 1.0d0 )
            
             dg21 = GPRIME*(rha**(8.0d0/3.0d0) + rhb**(8.0d0/3.0d0))
             dg22 = GPRIME*(rh*lrh - DOT_PRODUCT(grh,grh))
             dg23 = GPRIME*(rha*lrha + rhb*lrhb)
             dg24 = GPRIME*(DOT_PRODUCT(grha,grha) + DOT_PRODUCT(grhb,grhb))

             CALL gradG2(g,rha,rhb,grha,grhb,grg)
             lg = lapG2(g,rha,rhb,grha,grhb,lrha,lrhb)
             
             ! here the expression (25) in Phys. Rev. B. 37, 785 (1988) is used
             ! for the spin up density:

             IF ( rh .GT. 1.0E-20 ) THEN
                Vc(1) = -a*( FPRIME*rh + ff2 ) - (2.0d0**(5.0d0/3.0d0))*a*b*CF*( dg21 + (8.0d0/3.0d0)*gg2*rha**(5.0d0/3.0d0) )
                Vc(1) = Vc(1) - (a*b/4.0d0)*( rh*lg + 4.0d0*(DOT_PRODUCT(grg,grh) + gg2*lrh ) + dg22 )
                Vc(1) = Vc(1) - (a*b/36.0d0)*( 3.0d0*rha*lg + 4.0d0*( DOT_PRODUCT(grha,grg) + gg2*lrha ) + 3.0d0*dg23 + dg24 )
             ELSE
                Vc(1) = 0.0d0
             ENDIF
             !----------------------------------
             ! the spin down part of the potential
             !----------------------------------
             FPRIME = F2pa(rh,rhb,rha,1.0d0)
             GPRIME = G2pa( rh,rhb,rha, 1.0d0 )
             
             dg21 = GPRIME*(rha**(8.0d0/3.0d0) + rhb**(8.0d0/3.0d0))
             dg22 = GPRIME*(rh*lrh - DOT_PRODUCT(grh,grh))
             dg23 = GPRIME*(rha*lrha + rhb*lrhb)
             dg24 = GPRIME*(DOT_PRODUCT(grha,grha) + DOT_PRODUCT(grhb,grhb))

             IF ( rh .GT. 1.0E-20 ) THEN
                Vc(2) = -a*( FPRIME*rh + ff2 ) - (2.0d0**(5.0d0/3.0d0))*a*b*CF*( dg21 + (8.0d0/3.0d0)*gg2*rhb**(5.0d0/3.0d0) )
                Vc(2) = Vc(2) - (a*b/4.0d0)*( rh*lg + 4.0d0*(DOT_PRODUCT(grg,grh) + gg2*lrh ) + dg22 )
                Vc(2) = Vc(2) - (a*b/36.0d0)*( 3.0d0*rhb*lg + 4.0d0*( DOT_PRODUCT(grhb,grg) + gg2*lrhb ) + 3.0d0*dg23 + dg24 )
             ELSE
                     Vc(2) = 0.0d0
             ENDIF
             
        END SUBROUTINE VLYP

SUBROUTINE VLYPA(rha,rhb,grha,grhb,Vc)
             !------------------------------------------------------------!
             ! This subroutine calculates Gradient corrected correlation  !
             ! potential of Lee, Yang and Parr (LYP) at the point r       !
             ! using the alternative expression in CHEM. PHYS. LETT. 157, !
             ! 200 (1989), that is independent of the laplacian           !
             ! densu = spin up density, densd = spin down density         !
             ! gdensu = spin up gradient, gdensd = spin down gradient     !
             !------------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN)  :: rha,rhb,grha(3),grhb(3)
             DOUBLE PRECISION, INTENT(OUT) :: Vc(2)
             DOUBLE PRECISION :: g,rh,grh(3),lrh
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             DOUBLE PRECISION, PARAMETER :: a = 0.04918d0
             DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
             DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
             DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
             DOUBLE PRECISION :: CF,df2,dg21,dg22,dg23,dg24,gg2,ff2,grf(3),grg(3),lf,lg,GPRIME,FPRIME
             DOUBLE PRECISION :: ga,gb,gz,ep,zi,it,th,ka,la
             DOUBLE PRECISION :: omeg,delta,domeg,ddelta,ep1,ep2,ep0,dedrh,dep1drha,dep1drhb,dep2drha,dep2drhb
             DOUBLE PRECISION :: dep0drhb,dep0drha

             CF = (3.0d0/10.0d0)*(3.0d0*PI**2)**(2.0d0/3.0d0)

             ! the total density
             rh = rha + rhb
             IF ( rh .LT. 1.0E-20 ) THEN
                     Vc = 0.0d0
                     RETURN
             ENDIF
             
             gz = DOT_PRODUCT(grha+grhb,grha+grhb)
             ga = DOT_PRODUCT(grha,grha)
             gb = DOT_PRODUCT(grhb,grhb)
                
             omeg  = rh**(-11.0d0/3.0d0)*EXP(-c*rh**(-1.0d0/3.0d0))/(1.0d0 + d*rh**(-1.0d0/3.0d0) )
             delta = c*rh**(-1.0d0/3.0d0)    + d*rh**(-1.0d0/3.0d0)/(1.0d0 + d*rh**(-1.0d0/3.0d0) )
             
             domeg = omeg*( (1.0d0/3.0d0)*rh**(-4.0d0/3.0d0)*( c  + d/ (1.0d0 + d*rh**(-1.0d0/3.0d0) ) ) - 11.0d0/(3.0d0*rh) )
             ddelta = -(1.0d0/3.0d0)*(delta/rh) + (d**2/3.0d0)*rh**(-5.0d0/3.0d0)/( (1.0d0 + d*rh**(-1.0d0/3.0d0) )**2 )
        
             ep0 = -4.0d0*a*rha*rhb/(rh + d*rh**(2.0d0/3.0d0) )

             ep1 = 2**(11.0d0/3.0d0)*CF*(rha**(8.0d0/3.0d0) + rhb**(8.0d0/3.0d0)) + ( (47.0d0/18.0d0) - 7.0d0*delta/18.0d0 )*gz - ( 2.50d0 - delta/18.0d0 )*( ga + gb ) &
             & - ((delta-11.0d0)/9.0d0 )*( (rha/rh)*ga + (rhb/rh)*gb )

             ep2 = -(2.0d0/3.0d0)*rh**2*gz + ( (2.0d0/3.0d0)*rh**2 - rha**2 )*gb + ( (2.0d0/3.0d0)*rh**2 - rhb**2 )*ga
        
             dep0drha = ep0*( 1.0d0/rha - ( 1.0d0 + (2.0d0*d/3.0d0)*rh**(-1.0d0/3.0d0) )/(rh + d*rh**(2.0d0/3.0d0) ) )
             
             dep0drhb = ep0*( 1.0d0/rhb - ( 1.0d0 + (2.0d0*d/3.0d0)*rh**(-1.0d0/3.0d0) )/(rh + d*rh**(2.0d0/3.0d0) ) )

             dep1drha = (8.0d0/3.0d0)*2**(11.0d0/3.0d0)*CF*rha**(5.0d0/3.0d0) + (1.0d0/18.0d0)*ddelta*( ga + gb - 7.0d0*gz ) + &
             & ( (delta-11.0d0)/(9.0d0*rh) - (1.0d0/9.0d0)*ddelta )*( (rha/rh)*ga + (rhb/rh)*gb ) - ((delta-11)/(9.0d0*rh))*ga

             dep1drhb = (8.0d0/3.0d0)*2**(11.0d0/3.0d0)*CF*rhb**(5.0d0/3.0d0) + (1.0d0/18.0d0)*ddelta*( ga + gb - 7.0d0*gz ) + &
             & ( (delta-11.0d0)/(9.0d0*rh) - (1.0d0/9.0d0)*ddelta )*( (rha/rh)*ga + (rhb/rh)*gb ) - ((delta-11)/(9.0d0*rh))*gb
             
             dep2drha = (4.0d0/3.0d0)*rh*( ga - gz ) + 2.0d0*( (2.0d0/3.0d0)*rh - rha )*gb
      
             dep2drhb = (4.0d0/3.0d0)*rh*( gb - gz ) + 2.0d0*( (2.0d0/3.0d0)*rh - rhb )*ga
        
             Vc(1) = dep0drha - a*b*( domeg*( rha*rhb*ep1 + ep2 ) + omeg*( rhb*ep1 + rha*rhb*dep1drha + dep2drha ) )
             Vc(2) = dep0drhb - a*b*( domeg*( rha*rhb*ep1 + ep2 ) + omeg*( rha*ep1 + rha*rhb*dep1drhb + dep2drhb ) )

END SUBROUTINE VLYPA

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
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             IF ( densu+densd .LT. 1.0E-20 ) THEN
                     Vc = 0.0d0
                     RETURN
             ENDIF
             
             !CALL VWN(densu,densd,VVWN)
             CALL VWNC(densu,densd,VVWN)

             CALL VXB88(densu,densd,gdensu,gdensd,VVXB88)

             !CALL VLYP(densu,densd,gdensu,gdensd,ldensu,ldensd,VVLYP)
             CALL VLYPA(densu,densd,gdensu,gdensd,VVLYP)
             
             CALL Vx(densu,densd,VXLSDA)
             
             !VXLSDA(1) = -2.0d0*(3.0d0/(4.0d0*PI))**(1.0d0/3.0d0)*densu**(1.0d0/3.0d0)
             !VXLSDA(2) = -2.0d0*(3.0d0/(4.0d0*PI))**(1.0d0/3.0d0)*densd**(1.0d0/3.0d0)
                
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
                DOUBLE PRECISION, PARAMETER :: a = 0.04918d0
                DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
                DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                DOUBLE PRECISION :: CF,gam,twu,twd,tw,dens,LYPE,VWNE,LSDAE,XB88E
                DOUBLE PRECISION :: rs,zeta,fz
                DOUBLE PRECISION, PARAMETER :: aa = 0.20d0
                DOUBLE PRECISION, PARAMETER :: bb = 0.720d0
                DOUBLE PRECISION, PARAMETER :: cc = 0.810d0
                DOUBLE PRECISION :: om,del,dpa,dpb,dpt
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

                IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
                IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
                
                CF = (3.0d0/10.0d0)*(3.0d0*PI**2)**(2.0d0/3.0d0)

                !gam = 2.0d0*( 1.0d0 - ( densu**2 + densd**2)/( (densu+densd)**2) )

                !twu = (1.0d0/8.0d0)*( DOT_PRODUCT(gdensu,gdensu)/densu - ldensu )
                !twd = (1.0d0/8.0d0)*( DOT_PRODUCT(gdensd,gdensd)/densd - ldensd )

                !tw = (1.0d0/8.0d0)*( DOT_PRODUCT(gdensu+gdensd,gdensu+gdensd)/dens - (ldensu+ldensd))

                !---------------------------------------------------------------------------------------------
                !ALT1: The LYP correlation energy density acording to Eqn (22) in Phys. Rev. B. 37, 785 (1988):
                !----------------------------------------------------------------------------------------------
                
                !LYPE = -a*(gam/(1.0d0 + d*dens**(-1.0d0/3.0d0) ) )*( dens + 2.0d0*b*dens**(-5.0d0/3.0d0)*( CF*2.0d0**(2.0d0/3.0d0)*( densu**(8.0d0/3.0d0) + densd**(8.0d0/3.0d0) ) &
                !& - dens*tw + (1.0d0/9.0d0)*(densu*twu + densd*twd) + (1.0d0/18.0d0)*(densu*ldensu+densd*ldensd) )*EXP(-c*dens**(-1.0d0/3.0d0) ) )
              
                !-----------------------------------------------------
                !ALT2: (Eqn (2) page 201 in CHEM. PHYS. LETT. 157, 200)
                !------------------------------------------------------

                om = dens**(-11.0d0/3.0d0)*EXP(-c*dens**(-1.0d0/3.0d0))/(1.0d0 + d*dens**(-1.0d0/3.0d0)) 
                del = c*dens**(-1.0d0/3.0d0) + d*dens**(-1.0d0/3.0d0)/(1.0d0 + d*dens**(-1.0d0/3.0d0)) 
                dpa = DOT_PRODUCT(gdensu,gdensu)
                dpb = DOT_PRODUCT(gdensd,gdensd)
                dpt = DOT_PRODUCT(gdensu+gdensd,gdensu+gdensd)

                LYPE = -4*a*(densu*densd/dens)/(1.0d0 + d*dens**(-1.0d0/3.0d0))
                LYPE = LYPE -a*b*om*( densu*densd*( 2.0d0**(11.0d0/3.0d0)*CF*( densu**(8.0d0/3.0d0) + densd**(8.0d0/3.0d0) ) + ( 47.0d0/18.0d0 - 7*del/18.0d0)*dpt &
                & - (5.0d0/2.0d0-del/18.0d0)*(dpa+dpb) - ((del-11)/9.0d0)*( (densu/dens)*dpa + (densd/dens)*dpb ) ) &
                & - (2.0d0/3.0d0)*dens**2*dpt + ( (2.0d0/3.0d0)*dens**2 - densu**2)*dpb + ( (2.0d0/3.0d0)*dens**2 - densd**2)*dpa )
                
                ! The LSDA exchange energy density is calculated here:
                ! according to Eqn (26) in J.P Perdue and Y. Wang PRB, 45, 13244 (1992)
                fz = 0.50d0*( (1.0d0 + zeta)**(4.0d0/3.0d0) + (1.0d0 - zeta)**(4.0d0/3.0d0) )
                LSDAE =  - (3.0d0/4.0d0)*(3*dens/PI)**(1.0d0/3.0d0)*fz*dens

                !LSDAE = -(3.0d0/2.0d0)*(3.0d0/(4.0d0*PI))**(1.0d0/3.0d0)*(densu**(4.0d0/3.0d0) + densd**(4.0d0/3.0d0) )

                ! The correlation energy density of of Vosko, Wilk, and Nusair is calculated here
                !VWNE = EUEG(rs,zeta)*dens
                VWNE = EUEGC(rs,zeta)*dens
                
                ! The Becke exchange energy density is calculated here:
                XB88E = EXB88(densu,densd,gdensu,gdensd)

                ! The B3LYB energy density excluding 20% exact exchange is calculated here
                EB3LYP = (1.0d0 - aa )*LSDAE + bb*XB88E + cc*LYPE + (1.0d0 - cc)*VWNE
                
        END FUNCTION EB3LYP

        SUBROUTINE gVB3LYP(densu,densd,gdensu,gdensd,gV)
                ! The derivative with respect to the total gradient length the gradient length of the spin down and up densities 
                ! of the B3LYB energy density excluding 20% exact exchange is calculated here
                IMPLICIT NONE
                DOUBLE PRECISION,INTENT(IN) :: densu,densd,gdensu(3),gdensd(3)
                DOUBLE PRECISION,INTENT(OUT) :: gV(3)
                DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
                DOUBLE PRECISION, PARAMETER :: a = 0.04918d0
                DOUBLE PRECISION, PARAMETER :: b = 0.1320d0
                DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                DOUBLE PRECISION :: CF,gam,twu,twd,tw,dens,gVLYPE(3),gVXB88E(2)
                DOUBLE PRECISION :: rs,zeta,fz
                DOUBLE PRECISION, PARAMETER :: aa = 0.20d0
                DOUBLE PRECISION, PARAMETER :: bb = 0.720d0
                DOUBLE PRECISION, PARAMETER :: cc = 0.810d0
                DOUBLE PRECISION :: om,del,dpa,dpb,dpt
                dens = densu+densd
                
                IF ( dens .LT. 1.0E-20 ) THEN
                        gVLYPE = 0.0d0
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

                IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
                IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
                
                om = dens**(-11.0d0/3.0d0)*EXP(-c*dens**(-1.0d0/3.0d0))/(1.0d0 + d*dens**(-1.0d0/3.0d0)) 
                del = c*dens**(-1.0d0/3.0d0) + d*dens**(-1.0d0/3.0d0)/(1.0d0 + d*dens**(-1.0d0/3.0d0)) 
                dpa = sqrt(DOT_PRODUCT(gdensu,gdensu))
                dpb = sqrt(DOT_PRODUCT(gdensd,gdensd))
                dpt = sqrt(DOT_PRODUCT(gdensu+gdensd,gdensu+gdensd))

                gVLYPE(1) = -2.0d0*a*b*om*dpa*( densu*densd*(  - 5.0d0/2.0d0 + del/18.0d0 - ((del-11.0d0)/9.0d0)*(densu/dens) ) +  (2.0d0/3.0d0)*dens**2 - densd**2 )*cc
                
                gVLYPE(2) = -2.0d0*a*b*om*dpb*( densu*densd*(  - 5.0d0/2.0d0 + del/18.0d0 - ((del-11.0d0)/9.0d0)*(densd/dens) ) +  (2.0d0/3.0d0)*dens**2 - densu**2 )*cc
                
                gVLYPE(3) = -2.0d0*a*b*om*dpt*( densu*densd*( 47.0d0/18.0d0 - 7.0d0*del/18.0d0 )  - (2.0d0/3.0d0)*dens**2 )*cc

                CALL gVXB88(densu,densd,gdensu,gdensd,gVXB88E)

                gVLYPE(1) = gVLYPE(1) + bb*gVXB88E(1)
                gVLYPE(2) = gVLYPE(2) + bb*gVXB88E(2)
                
                gv = gVLYPE

        END SUBROUTINE gVB3LYP

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
             DOUBLE PRECISION :: rh,xa,xb,asha,dasha,ashb,dashb,VXLSDA(2),nna,nnb
             DOUBLE PRECISION, PARAMETER :: bb = 0.00420d0
             
             IF ( rha .GT. 1.0E-20 ) THEN
                xa = sqrt(DOT_PRODUCT(grha,grha))/(rha**(4.0d0/3.0d0))
             ELSE
                xa = 0.0d0
             ENDIF
             IF ( rhb .GT. 1.0E-20 ) THEN
                xb = sqrt(DOT_PRODUCT(grhb,grhb))/(rhb**(4.0d0/3.0d0))
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

             nna = (1.0d0 + 6.0d0*bb*xa*asha )
             nnb = (1.0d0 + 6.0d0*bb*xb*ashb )
                
             Vxb(1) =  -( ((4.0d0/3.0d0)*bb*(rha**(1.0d0/3.0d0))*xa**2)/nna )*(6.0d0*bb*xa*(asha + xa*dasha)/nna - 1.0d0 )
             
             Vxb(2) =  -( ((4.0d0/3.0d0)*bb*(rhb**(1.0d0/3.0d0))*xb**2)/nnb )*(6.0d0*bb*xb*(ashb + xb*dashb)/nnb - 1.0d0 )
             
     END SUBROUTINE VXB88
     
     SUBROUTINE gVXB88(rha,rhb,grha,grhb,gV)
                ! This subroutine calculates the derivatives of the Becke exchange
                ! correction with respect to the lenght of the density gradients
                ! for the spin up and down chanells. See Eqn (15) p. "7" in the 
                ! other black note book.
                IMPLICIT NONE
                DOUBLE PRECISION, INTENT(IN) :: rha,rhb,grha(3),grhb(3)
                DOUBLE PRECISION, INTENT(OUT) :: gV(2)
                DOUBLE PRECISION :: rh,xa,xb,asha,dasha,ashb,dashb,VXLSDA(2),nna,nnb
                DOUBLE PRECISION, PARAMETER :: bb = 0.00420d0

                IF ( rha .GT. 1.0E-20 ) THEN
                        xa = sqrt(DOT_PRODUCT(grha,grha))/(rha**(4.0d0/3.0d0))
                ELSE
                        xa = 0.0d0
                ENDIF
                IF ( rhb .GT. 1.0E-20 ) THEN
                        xb = sqrt(DOT_PRODUCT(grhb,grhb))/(rhb**(4.0d0/3.0d0))
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

                nna = (1.0d0 + 6.0d0*bb*xa*asha )
                nnb = (1.0d0 + 6.0d0*bb*xb*ashb )

                gV(1) = -bb*xa*( 2.0d0 + 6.0d0*bb*xa*( asha - xa*dasha ) )/(nna**2)
                gV(2) = -bb*xb*( 2.0d0 + 6.0d0*bb*xb*( ashb - xa*dashb ) )/(nnb**2)

     END  SUBROUTINE gVXB88

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
             DOUBLE PRECISION :: rh,rs,zeta,xa,xb,asha,dasha,ashb,dashb,LSDAEX,fz
             DOUBLE PRECISION, PARAMETER :: bb = 0.00420d0
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             
             xa = sqrt(DOT_PRODUCT(grha,grha))/(rha**(4.0d0/3.0d0))
             xb = sqrt(DOT_PRODUCT(grhb,grhb))/(rhb**(4.0d0/3.0d0))
             
             
             ! calculating of the inverse hyperbolic functions:
             ! asinh(xa) and asinh(xb)
             asha = LOG( xa + sqrt(xa**2 + 1) )
             ashb = LOG( xb + sqrt(xb**2 + 1) )


             EXB88 =  -bb*( rha**(4.0d0/3.0d0)*xa**2/(1.0d0 + 6.0d0*bb*xa*asha ) + rhb**(4.0d0/3.0d0)*xb**2/(1.0d0 + 6.0d0*bb*xb*ashb ) )

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
             DOUBLE PRECISION :: A,PHH,PHHP,ks,kf,t,s,n,FXX,FXXP
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: g,rh,grh(3),lh,Vhom(2)
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             DOUBLE PRECISION, PARAMETER :: bb = 0.0667250d0
             DOUBLE PRECISION, PARAMETER :: gg = 0.0310910d0
             DOUBLE PRECISION, PARAMETER :: mm = 0.219510d0
             DOUBLE PRECISION, PARAMETER :: kk = 1.245 ! Here we actually choose the kappa value of the revPBE ( PRL 80, 890 (1998) )
                                                       ! instead of kk = 0.804, as is used in the original PBE (Phys. Rev. Lett. 77, 3865 (1996) )
             DOUBLE PRECISION ::  hh,rs,zeta,EXHOM,VXHOM(2),GFAKTOR,FAKTOR,HHPBE
             DOUBLE PRECISION ::  dtdn(2),dAdn(2),dHdn(2),EU,fzeta,fzetap
            
             n = rha + rhb

             IF ( n .LT. 1.0E-20 ) THEN
                     Vc = 0.0d0
                     RETURN
             ENDIF

             zeta = (rha - rhb)/n
             
             IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
             IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
             
             rs = (3.0d0/(4.0d0*PI*n))**(1.0d0/3.0d0)

             PHH = 0.50d0*( (1+zeta)**(2.0d0/3.0d0) + (1-zeta)**(2.0d0/3.0d0) )
             PHHP = (1.0d0/3.0d0)*( (1+zeta)**(-1.0d0/3.0d0) - (1-zeta)**(-1.0d0/3.0d0) )

             kf = (3*PI**2*n)**(1.0d0/3.0d0)
             
             ks = sqrt(4.0d0*kf/PI)
             
             t = lgrh/(2.0d0*PHH*ks*n)
             
             s = lgrh/(2.0d0*kf*n)
        
             !EU = EUEG(rs,zeta)
             EU = EUEGC(rs,zeta)

             IF (  EXP(-EU/(gg*PHH**3)) .NE. 1.0d0   ) THEN
                     A = (bb/gg)/( EXP(-EU/(gg*PHH**3)) - 1.0d0 )
             ELSE
                     A = 1.0E20
             ENDIF
             
             !=======================================================
             ! The nonpolarized exchange energy density is used here
             !=======================================================

             EXHOM = - 3.0d0*kf/(4*PI)      !non-pol
             
             FXX  = 1.0d0 + kk - kk/(1.0d0 + (mm/kk)*s**2 )
             FXXP =  2.0d0*mm*s/((1.0d0 + (mm/kk)*s**2 )**2)
             
             VXHOM = (4.0d0/3.0d0)*EXHOM*( FXX - FXXP*s ) !non-pol
             
             !===================================================
             ! Here we use the polarized exchange energy density:
             !===================================================
             !fzeta = 0.50d0*( (1+zeta)**(4.0d0/3.0d0) + (1-zeta)**(4.0d0/3.0d0) )   ! polarized
             !fzetap = (2.0d0/3.0d0)*( (1+zeta)**(1.0d0/3.0d0) - (1-zeta)**(1.0d0/3.0d0) ) ! polrized 
             
             !VXHOM(1) = -(3.0d0*n/PI)**(1.0d0/3.0d0)*( fzeta*( FXX - FXXP*s ) + (3.0d0/4.0d0)*fzetap*(-zeta+1)*FXX ) ! polrized  
             !VXHOM(2) = -(3.0d0*n/PI)**(1.0d0/3.0d0)*( fzeta*( FXX - FXXP*s ) + (3.0d0/4.0d0)*fzetap*(-zeta-1)*FXX ) ! polrized 

             ! Calculating the derivatives dA/dn and dt/dn for the up and down
             ! densities:
                
             !CALL VWN(rha,rhb,Vhom)
             CALL VWNC(rha,rhb,Vhom)
             dAdn(1) = (A**2/(bb*n*PHH**3))*EXP(-EU/(gg*PHH**3))*( Vhom(1) - ((PHHP/PHH)*(-zeta+1.0d0)+1.0d0)*EU )
             dAdn(2) = (A**2/(bb*n*PHH**3))*EXP(-EU/(gg*PHH**3))*( Vhom(2) - ((PHHP/PHH)*(-zeta-1.0d0)+1.0d0)*EU )

             dtdn(1) = -(t/n)*( 0.0d0*(PHHP/PHH)*(-zeta+1) + 7.0d0/6.0d0 )
             dtdn(2) = -(t/n)*( 0.0d0*(PHHP/PHH)*(-zeta-1) + 7.0d0/6.0d0 )
             
             dHdn(1) = (PHHP/n)*(-zeta+1)
             dHdn(2) = (PHHP/n)*(-zeta-1)
                
             HHPBE = HPBE(rha,rhb,lgrh)
             dHdn  = dHdn*(3*HHPBE/PHH)
             
             GFAKTOR = 1.0d0 + (bb/gg)*(t**2)*(1.0d0 + A*t**2)/(1.0d0 + A*t**2 + (A*t**2)**2 )

             GFAKTOR = GFAKTOR*(1.0d0 + A*t**2 + (A*t**2)**2 )**2

             FAKTOR = (bb*PHH**3/GFAKTOR)

             dHdn = dHdn + FAKTOR*( 2*t*(1.0d0+2*A*t**2)*dtdn - A*t**6*(2.0d0+A*t**2)*dAdn )

             !===================================================================
             ! Here we put the entire PBE exchange correlation potential together
             !===================================================================
              Vc(:) =  VXHOM + HHPBE
              Vc = Vc + Vhom +  n*dHdn

      END SUBROUTINE VPBE

      FUNCTION gVPBE(rha,rhb,lgrh)
               ! This function calculates the derivative of: rho*( HPBE + E_c(uniform) + F_x*E_x(uniform) )
               ! with respect to the length of the total density gradient,|grad*rho|. See my notes page "4"
               ! equation (7) in the other black note book
               IMPLICIT NONE
               DOUBLE PRECISION :: gVPBE
               DOUBLE PRECISION :: rha,rhb,lgrh
               DOUBLE PRECISION :: A,PHH,PHHP,ks,kf,t,s,n,EXHOM,EU,term,rs,zeta
               DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
               DOUBLE PRECISION, PARAMETER :: bb = 0.0667250d0
               DOUBLE PRECISION, PARAMETER :: gg = 0.0310910d0
               DOUBLE PRECISION, PARAMETER :: mm = 0.219510d0
               DOUBLE PRECISION, PARAMETER :: kk = 1.245 ! Here we actually choose the kappa value of the revPBE ( PRL 80, 890 (1998) )
                                                         ! instead of kk = ! 0.804, as is used in  the original PBE (Phys. Rev. Lett. 77, 3865 (1996) )
               n = rha + rhb

               IF ( n .LT. 1.0E-20 ) THEN
                       gVPBE = 0.0d0
                       RETURN
               ENDIF

               zeta = (rha - rhb)/n

               IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
               IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
                
               rs = (3.0d0/(4.0d0*PI*n))**(1.0d0/3.0d0)

               PHH = 0.50d0*( (1+zeta)**(2.0d0/3.0d0) + (1-zeta)**(2.0d0/3.0d0))

               kf = (3*PI**2*n)**(1.0d0/3.0d0)

               ks = sqrt(4.0d0*kf/PI)

               t = lgrh/(2.0d0*PHH*ks*n)

               s = lgrh/(2.0d0*kf*n)

               EXHOM = - 3.0d0*kf/(4*PI)

               !EU = EUEG(rs,zeta)
               EU = EUEGC(rs,zeta)

               IF (  EXP(-EU/(gg*PHH**3)) .NE. 1.0d0   ) THEN
                       A = (bb/gg)/( EXP(-EU/(gg*PHH**3)) - 1.0d0 )
               ELSE
                       A = 1.0E20
               ENDIF

               term = 1.0d0 + A*t**2 + A**2*t**4

               gVPBE = bb*PHH**2*t*( 1.0d0 + 2*A*t**2)/( ks*( 1.0d0 + (bb/gg)*t**2*( 1 + A*t**2 )/term )*term**2 )

               gVPBE = gVPBE + mm*s*EXHOM/( kf*( 1.0d0 + (mm/kk)*s**2 )**2 )

       END FUNCTION gVPBE

         SUBROUTINE Vx(densu,densd,Vxx)
             !-----------------------------------------------------------!
             ! This subroutine calculates the LSDA exchange  potential   !
             ! densu = spin up density, densd = spin down density        !
             !-----------------------------------------------------------!
             IMPLICIT NONE
             DOUBLE PRECISION, INTENT(IN) :: densu,densd
             DOUBLE PRECISION, INTENT(OUT) :: Vxx(2)
             DOUBLE PRECISION, EXTERNAL :: rho
             DOUBLE PRECISION :: dens,funkz,pfunkz,zeta
             DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
             
             dens = densu + densd
             zeta = (densu-densd)/dens
             
             IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
             IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
             
             ! Here we have written down the polarized version of the Exchange 
             ! potential according to Eqn (26) in J.P Perdue and Y. Wang PRB, 45, 13244 (1992)
             funkz = 0.50d0*( (1+zeta)**(4.0d0/3.0d0) + (1-zeta)**(4.0d0/3.0d0))
             pfunkz = (2.0d0/3.0d0)*( (1+zeta)**(1.0d0/3.0d0) - (1-zeta)**(1.0d0/3.0d0))
             Vxx(1) = -(3*dens/PI)**(1.0d0/3.0d0)*( funkz + (3.0d0/4.0d0)*pfunkz*(-zeta+1.0d0) )
             Vxx(2) = -(3*dens/PI)**(1.0d0/3.0d0)*( funkz + (3.0d0/4.0d0)*pfunkz*(-zeta-1.0d0) )
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
        
        FUNCTION epsPC(rs)
                !-------------------------------------------------------------------!
                ! Eqn. (39) i= P,  in Comp. Phys. Comm. 66, 383 (1991)              !
                !-------------------------------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: epsPC
                DOUBLE PRECISION :: rs
                DOUBLE PRECISION, PARAMETER :: A = 0.03109070d0   
                DOUBLE PRECISION, PARAMETER :: b = 3.727440d0    
                DOUBLE PRECISION, PARAMETER :: c = 12.93520d0    
                DOUBLE PRECISION, PARAMETER :: x0 = -0.104980d0  
                DOUBLE PRECISION :: Q,Xrs,Xx0

                Q = sqrt(4.0d0*c - b**2.0d0 )
                Xrs = rs + b*sqrt(rs)+c
                Xx0 = x0**2 + b*dabs(x0)+c

                epsPC = log(rs/Xrs) + (2.0d0*b/Q)*atan(Q/(2*sqrt(rs) + b ))
                epsPC = epsPC - (b*x0/Xx0)*( log( ( (sqrt(rs) - x0)**2.0d0 )/Xrs ) + (2.0d0*(b+2.0d0*x0)/Q)*atan(Q/(2*sqrt(rs) + b )) )
                epsPC = epsPC*A
       END FUNCTION epsPC

       FUNCTION epsFC(rs)
                !-------------------------------------------------------------------!
                ! Eqn. (39) i= F,  in Comp. Phys. Comm. 66, 383 (1991)              !
                !-------------------------------------------------------------------!
                IMPLICIT NONE
                DOUBLE PRECISION :: epsFC
                DOUBLE PRECISION :: rs
                DOUBLE PRECISION, PARAMETER :: A = 0.015545350d0
                DOUBLE PRECISION, PARAMETER :: b = 7.060420d0
                DOUBLE PRECISION, PARAMETER :: c = 18.05780d0
                DOUBLE PRECISION, PARAMETER :: x0 = -0.32500d0
                DOUBLE PRECISION :: Q,Xrs,Xx0

                Q = sqrt(4.0d0*c - b**2.0d0 )
                Xrs = rs + b*sqrt(rs)+c
                Xx0 = x0**2 + b*dabs(x0)+c

                epsFC = log(rs/Xrs) + (2.0d0*b/Q)*atan(Q/(2*sqrt(rs) + b ))
                epsFC = epsFC - (b*x0/Xx0)*( log( ( (sqrt(rs) - x0)**2.0d0 )/Xrs ) + (2.0d0*(b+2.0d0*x0)/Q)*atan(Q/(2*sqrt(rs) + b )) )
                epsFC = epsFC*A
       END FUNCTION epsFC

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
                
                FUNCTION betafC(rs)
                        !---------------------------------!
                        ! Eqn (4), in PRB 42, 3205 (1990) !
                        !---------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: betafC
                        DOUBLE PRECISION :: rs
                        betafC = ffunkbis(0.0d0)*( epsFC(rs) - epsPC(rs) )/alphac(rs) - 1.0d0
                END FUNCTION betafC
                
                FUNCTION deps(rs,zeta)
                        !---------------------------------!
                        ! Eqn (2), in PRB 42, 3205 (1990) !
                        !---------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: deps
                        DOUBLE PRECISION :: rs,zeta
                        deps = alphac(rs)*(ffunk(zeta)/ffunkbis(0.0d0))*( 1.0d0 + betaf(rs)*zeta**4.0d0 )
                END FUNCTION deps

                FUNCTION depsC(rs,zeta)
                        !---------------------------------!
                        ! Eqn (2), in PRB 42, 3205 (1990) !
                        !---------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: depsC
                        DOUBLE PRECISION :: rs,zeta
                        depsC = alphac(rs)*(ffunk(zeta)/ffunkbis(0.0d0))*( 1.0d0 + betafC(rs)*zeta**4.0d0 )
                END FUNCTION depsC
                
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
                                mup = epsP(rs)*(1.0d0 + (7.0d0/6.0d0)*beta1*sqrt(rs) + (4.0d0/3.0d0)*beta2*rs )/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                        ENDIF
                END FUNCTION muP
                
                FUNCTION muF(rs)
                        !-----------------------------------------------!
                        ! Eqn (12,a,b), i=F,  in PRB 42, 3205 (1990)    !
                        !-----------------------------------------------!
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
                                muF = epsF(rs)*(1.0d0 + (7.0d0/6.0d0)*beta1*sqrt(rs) + (4.0d0/3.0d0)*beta2*rs )/( 1.0d0 + beta1*sqrt(rs) + beta2*rs )
                        ENDIF
                END FUNCTION muF

                FUNCTION muPC(rs)
                        !-------------------------------------------------------------------!
                        ! Eqn. (40) i= P,  in Comp. Phys. Comm. 66, 383 (1991)              !
                        !-------------------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: mupC
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION, PARAMETER :: A = 0.03109070d0   
                        DOUBLE PRECISION, PARAMETER :: b = 3.727440d0    
                        DOUBLE PRECISION, PARAMETER :: c = 12.93520d0    
                        DOUBLE PRECISION, PARAMETER :: x0 = -0.104980d0  
                        DOUBLE PRECISION :: b1,b2,b3

                        b1 = (b*x0 - c)/(c*x0)
                        b2 = (x0 - b)/(c*x0)
                        b3 = - 1.0d0/(c*x0)

                        muPC = epsPC(rs) - (A/3.0d0)*( 1.0d0 + b1*sqrt(rs) )/( 1.0d0 + b1*sqrt(rs) + b2*rs + b3*rs**(1.50d0) )

                END FUNCTION muPC
                
                FUNCTION muFC(rs)
                        !-------------------------------------------------------------------!
                        ! Eqn. (40) i= P,  in Comp. Phys. Comm. 66, 383 (1991)              !
                        !-------------------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: muFC
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION, PARAMETER :: A = 0.01554535000
                        DOUBLE PRECISION, PARAMETER :: b = 7.060420d0
                        DOUBLE PRECISION, PARAMETER :: c = 18.05780d0
                        DOUBLE PRECISION, PARAMETER :: x0 = -0.32500d0
                        DOUBLE PRECISION :: b1,b2,b3

                        b1 = (b*x0 - c)/(c*x0)
                        b2 = (x0 - b)/(c*x0)
                        b3 = - 1.0d0/(c*x0)

                        muFC = epsFC(rs) - (A/3.0d0)*( 1.0d0 + b1*sqrt(rs) )/( 1.0d0 + b1*sqrt(rs) + b2*rs + b3*rs**(1.50d0) )

                END FUNCTION muFC
                
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
                        DOUBLE PRECISION :: depp,deff
                        DOUBLE PRECISION, PARAMETER :: Ap =  0.03110d0
                        DOUBLE PRECISION, PARAMETER :: Bp = -0.0480d0
                        DOUBLE PRECISION, PARAMETER :: Cp =  0.00200d0
                        DOUBLE PRECISION, PARAMETER :: Dp = -0.01160d0
                        DOUBLE PRECISION, PARAMETER :: beta1p =  1.05290d0
                        DOUBLE PRECISION, PARAMETER :: beta2p =  0.33340d0
                        DOUBLE PRECISION, PARAMETER :: gamap  = -0.14230d0
                        DOUBLE PRECISION, PARAMETER :: Af =  0.015550d0
                        DOUBLE PRECISION, PARAMETER :: Bf = -0.02690d0
                        DOUBLE PRECISION, PARAMETER :: Cf =  0.00070d0
                        DOUBLE PRECISION, PARAMETER :: Df = -0.0048
                        DOUBLE PRECISION, PARAMETER :: beta1f =  1.39810d0
                        DOUBLE PRECISION, PARAMETER :: beta2f =  0.26110d0
                        DOUBLE PRECISION, PARAMETER :: gamaf  = -0.08430d0

                        IF ( rs .LT. 1.0d0 ) THEN
                                depp = -gamap*( 0.050d0*beta1p/sqrt(rs) + beta2p )/(( 1.0d0 + beta1p*sqrt(rs) + beta2p*rs )**2)
                        ELSE
                                depp = Ap/rs + Cp*( log(rs) + 1.0d0 ) + Dp
                        ENDIF
                        
                        IF ( rs .LT. 1.0d0 ) THEN
                                deff = -gamaf*( 0.050d0*beta1f/sqrt(rs) + beta2f )/(( 1.0d0 + beta1f*sqrt(rs) + beta2f*rs )**2)
                        ELSE
                                deff = Af/rs + Cf*( log(rs) + 1.0d0 ) + Df
                        ENDIF
                        
                        rdbdr = (ffunkbis(0.0d0)/alphac(rs))*rs*( deff - depp )
                        rdbdr = rdbdr - ( ffunkbis(0.0d0)*(epsF(rs)-epsP(rs))/(alphac(rs)**2.0d0) )*rdadr(rs)
                END FUNCTION rdbdr
                
                FUNCTION rdbdrC(rs)
                        !-------------------------------------------------------------------!
                        ! Eqn (10),in PRB 42, 3205 (1990)                                   !
                        ! rdbdr(rs) = rs*betaf'(rs)                                         !
                        !-------------------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: rdbdrC
                        DOUBLE PRECISION :: rs
                        DOUBLE PRECISION :: depp,deff
                        DOUBLE PRECISION, PARAMETER :: Ap = 0.06218140d0
                        DOUBLE PRECISION, PARAMETER :: bp = 3.727440d0
                        DOUBLE PRECISION, PARAMETER :: cp = 12.93520d0
                        DOUBLE PRECISION, PARAMETER :: x0p = -0.104980d0
                        DOUBLE PRECISION, PARAMETER :: Af =  0.03109070d0
                        DOUBLE PRECISION, PARAMETER :: bf = 7.060420d0
                        DOUBLE PRECISION, PARAMETER :: cf = 18.05780d0
                        DOUBLE PRECISION, PARAMETER :: x0f = -0.325000d0
                        DOUBLE PRECISION :: b1,b2,b3

                        
                        b1 = (bp*x0p - cp)/(cp*x0p)
                        b2 = (x0p - bp)/(cp*x0p)
                        b3 = - 1.0d0/(cp*x0p)

                        depp = Ap*( 1.0d0 + b1*sqrt(rs) )/( 1.0d0 + b1*sqrt(rs) + b2*rs + b3*rs**(1.50d0) )
                        
                        b1 = (bf*x0f - cf)/(cf*x0f)
                        b2 = (x0f - bf)/(cf*x0f)
                        b3 = - 1.0d0/(cf*x0f)

                        deff = Af*( 1.0d0 + b1*sqrt(rs) )/( 1.0d0 + b1*sqrt(rs) + b2*rs + b3*rs**(1.50d0) )
                        
                        rdbdrC = (ffunkbis(0.0d0)/alphac(rs))*( deff - depp )
                        rdbdrC = rdbdrC - ( ffunkbis(0.0d0)*(epsFC(rs)-epsPC(rs))/(alphac(rs)**2.0d0) )*rdadr(rs)
                END FUNCTION rdbdrC
                
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
                        
                        F2 = g/(1.0d0 + d*rh**(-1.0d0/3.0d0) )
                        
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
                
                SUBROUTINE gradrho2(BAS,P,r,grho,V1,V2)
                        ! Calculates the gradient of the charge density with
                        ! respect to the electronic coordinate, r.
                        USE datatypemodule
                        IMPLICIT NONE 
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(3),V1(BAS%NBAS),V2(BAS%NBAS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: grad(3)
                        INTEGER :: I
                        
                        DO I=1,BAS%NBAS
                                V1(I) = basfunkval(BAS%PSI(I),r)
                                CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                V2(I,:) = grad
                         ENDDO

                         grho(1) = 2.0d0*DOT_PRODUCT(V2(:,1),MATMUL(P,V1))
                         grho(2) = 2.0d0*DOT_PRODUCT(V2(:,2),MATMUL(P,V1))
                         grho(3) = 2.0d0*DOT_PRODUCT(V2(:,3),MATMUL(P,V1))
                 END SUBROUTINE gradrho2
                                 
                SUBROUTINE nucgradrho(NATOMS,BAS,P,gradS,r,BECKECENTER,grho)
                        ! Calculates the nuclear gradient of the charge density.
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS,BECKECENTER
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(NATOMS,BAS%NBAS,3),grad(3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J
                        
                        grho = 0.0d0
                        V2 = 0.0d0

                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) V1(I) = basfunkval(BAS%PSI(I),r)
                                        IF ( BAS%PSI(I)%ATYPE .EQ. BECKECENTER ) THEN
                                                V2(J,I,:) = 0.0d0
                                        ELSE
                                                IF ( J .EQ. BECKECENTER ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = grad
                                                ENDIF
                                                IF ( J .EQ. BAS%PSI(I)%ATYPE ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = -grad
                                                ENDIF
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
                
                SUBROUTINE nucgradrhohole(NATOMS,BAS,P,PHOLE,gradS,r,BECKECENTER,grho)
                        ! Calculates the nuclear gradient of the charge density.
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS,BECKECENTER
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),PHOLE(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(NATOMS,BAS%NBAS,3),grad(3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J
                        
                        grho = 0.0d0
                        V2 = 0.0d0

                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) V1(I) = basfunkval(BAS%PSI(I),r)
                                        IF ( BAS%PSI(I)%ATYPE .EQ. BECKECENTER ) THEN
                                                V2(J,I,:) = 0.0d0
                                        ELSE
                                                IF ( J .EQ. BECKECENTER ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = grad
                                                ENDIF
                                                IF ( J .EQ. BAS%PSI(I)%ATYPE ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = -grad
                                                ENDIF
                                        ENDIF
                              ENDDO
                              ! Calculating the nuclear gradient of the density
                              ! matrix using Eqn (9) , page 125 in black note
                              ! book, i.e grad(P) = -P*grad(S)*P
                              gradP(J,1,:,:) = -MATMUL(P,MATMUL(gradS(J,1,:,:),P)) + MATMUL(PHOLE,MATMUL(gradS(J,1,:,:),PHOLE))
                              gradP(J,2,:,:) = -MATMUL(P,MATMUL(gradS(J,2,:,:),P)) + MATMUL(PHOLE,MATMUL(gradS(J,2,:,:),PHOLE))
                              gradP(J,3,:,:) = -MATMUL(P,MATMUL(gradS(J,3,:,:),P)) + MATMUL(PHOLE,MATMUL(gradS(J,3,:,:),PHOLE))
                        ENDDO
                        
                        ! See derivation of Eqn (6) and Eqn(6) p. 124 in the black note book.
                        DO J=1,NATOMS
                                grho(J,1) = DOT_PRODUCT(V1,MATMUL(gradP(J,1,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,1),MATMUL(P-PHOLE,V1)) 
                                grho(J,2) = DOT_PRODUCT(V1,MATMUL(gradP(J,2,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,2),MATMUL(P-PHOLE,V1)) 
                                grho(J,3) = DOT_PRODUCT(V1,MATMUL(gradP(J,3,:,:),V1)) + 2.0d0*DOT_PRODUCT(V2(J,:,3),MATMUL(P-PHOLE,V1)) 
                        ENDDO
                 END SUBROUTINE nucgradrhohole

                SUBROUTINE nucgradrhonip(NATOMS,BAS,P,r,BECKECENTER,grho)
                        ! Calculates the nuclear gradient of the charge density.
                        USE datatypemodule
                        IMPLICIT NONE
                        INTEGER, INTENT(IN) :: NATOMS,BECKECENTER
                        TYPE(BASIS), INTENT(IN) :: BAS
                        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS)
                        DOUBLE PRECISION, INTENT(IN) :: r(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grho(NATOMS,3)
                        DOUBLE PRECISION, EXTERNAL :: basfunkval
                        EXTERNAL :: gradbasfunkval
                        DOUBLE PRECISION :: V1(BAS%NBAS),V2(NATOMS,BAS%NBAS,3),grad(3),gradP(NATOMS,3,BAS%NBAS,BAS%NBAS)
                        INTEGER :: I,J
                        
                        grho = 0.0d0
                        V2 = 0.0d0

                        DO J=1,NATOMS
                              DO I=1,BAS%NBAS
                                        IF ( J .EQ. 1 ) V1(I) = basfunkval(BAS%PSI(I),r)
                                        IF ( BAS%PSI(I)%ATYPE .EQ. BECKECENTER ) THEN
                                                V2(J,I,:) = 0.0d0
                                        ELSE
                                                IF ( J .EQ. BECKECENTER ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = grad
                                                ENDIF
                                                IF ( J .EQ. BAS%PSI(I)%ATYPE ) THEN
                                                        CALL gradbasfunkval(BAS%PSI(I),r,grad)
                                                        V2(J,I,:) = -grad
                                                ENDIF
                                        ENDIF
                              ENDDO
                              ! Calculating the nuclear gradient of the density
                              ! matrix using Eqn (9) , page 125 in black note
                              ! book, i.e grad(P) = -P*grad(S)*P
                              !gradP(J,1,:,:) = -MATMUL(P,MATMUL(gradS(J,1,:,:),P))
                              !gradP(J,2,:,:) = -MATMUL(P,MATMUL(gradS(J,2,:,:),P))
                              !gradP(J,3,:,:) = -MATMUL(P,MATMUL(gradS(J,3,:,:),P))
                        ENDDO
                        
                        ! See derivation of Eqn (6) and Eqn(6) p. 124 in the black note book.
                        DO J=1,NATOMS
                                grho(J,1) =  2.0d0*DOT_PRODUCT(V2(J,:,1),MATMUL(P,V1)) 
                                grho(J,2) =  2.0d0*DOT_PRODUCT(V2(J,:,2),MATMUL(P,V1)) 
                                grho(J,3) =  2.0d0*DOT_PRODUCT(V2(J,:,3),MATMUL(P,V1)) 
                        ENDDO
                 END SUBROUTINE nucgradrhonip

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
                                F2pa = (1.0d0/(1.0d0 + d*rh**(-1.0d0/3.0d0) ))*( 4.0d0*rhb*(rhb-rha)/rt**3 + (d/3.0d0)*F2(g,rh)*drhdrha*rh**(-4.0d0/3.0d0) )
                        ELSE
                                F2pa = 0.0d0
                        ENDIF
                        
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
                                g = 1.0d0
                        ENDIF
                       
                        IF ( rh .GT. 1.0E-20  ) THEN
                                G2pa = ( F2pa(rh,rha,rhb,drhdrha)/F2(g,rh) + ( (c/3.0d0)*rh**(-4.0d0/3.0d0) - (5.0d0/(3.0d0*rh)) )*drhdrha )*G2(g,rh)
                        ELSE
                                G2pa = 0.0d0
                        ENDIF
                
                 END FUNCTION G2pa

                SUBROUTINE gradF2(g,rha,rhb,grha,grhb,grad)
                        ! calculates the gradient of the F2 function
                        ! graddens = gradient of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION, INTENT(IN) :: g,grha(3),grhb(3),rha,rhb
                        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        DOUBLE PRECISION :: rh,graddens(3),ggam(3),nnn
                        
                        CALL gradgammus(rha,rhb,grha,grhb,ggam)
                        
                        rh = rha + rhb
                        graddens = grha + grhb
                        nnn = 1.0d0 + d*rh**(-1.0d0/3.0d0)

                        grad = ( ggam/g + ( (d/3.0d0)*rh**(-1.0d0/3.0d0)/nnn )*graddens/rh )*F2(g,rh)

                END SUBROUTINE gradF2

                FUNCTION lapF2(g,rha,rhb,grha,grhb,lrha,lrhb)
                        ! Calculates the laplacian of the function F2. 
                        ! lapdens = laplacian of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION :: lapF2
                        DOUBLE PRECISION :: g,rha,rhb,grha(3),grhb(3),lrha,lrhb
                        DOUBLE PRECISION :: gf(3),DF22,Drh2,ggam(3)
                        DOUBLE PRECISION, PARAMETER :: d = 0.3490d0
                        DOUBLE PRECISION :: xx,ff,rh,graddens(3),lapdens
                        
                        graddens = grha + grhb
                        rh = rha + rhb
                        lapdens = lrha + lrhb
                        ff = F2(g,rh)
                        xx = (d/3.0d0)*rh**(-1.0d0/3.0d0)/(1.0d0 + d*rh**(-1.0d0/3.0d0) )
                        Drh2 = DOT_PRODUCT(graddens,graddens)/(rh**2)
                        CALL gradF2(g,rha,rhb,grha,grhb,gf)
                        CALL gradgammus(rha,rhb,grha,grhb,ggam)
                        ggam = ggam/g
                        lapF2 = DOT_PRODUCT(gf,gf)/ff + ( laplacegammus(rha,rhb,grha,grhb,lrha,lrhb)/g - DOT_PRODUCT(ggam,ggam) + (xx**2 -(4.0d0/3.0d0)*xx)*Drh2 + xx*lapdens/rh )*ff
                        
                END FUNCTION lapF2

                SUBROUTINE gradG2(g,rha,rhb,grha,grhb,grad)
                        ! calculates the gradient of the G2 function
                        ! graddens = gradient of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION, INTENT(IN) :: g,rha,rhb,grha(3),grhb(3)
                        DOUBLE PRECISION, INTENT(OUT) :: grad(3)
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        DOUBLE PRECISION :: gf(3),gg,GdividedbyF
                        DOUBLE PRECISION :: graddens(3),rh
                        
                        rh = rha + rhb
                        CALL gradF2(g,rha,rhb,grha,grhb,gf)
                        gg = G2(g,rh)
                        graddens = grha + grhb
                        
                        !ALT2:
                        grad = ( gf*rh**(-5.0d0/3.0d0) + F2(g,rh)*( (c/3.0d0)/rh**3.0d0 - (5.0d0/3.0d0)/rh**(8.0d0/3.0d0) )*graddens )*EXP(-c*rh**(-1.0d0/3.0d0))
                        
                        !ALT1: 
                        !grad = (gf/F2(g,rh))*gg + ( (c/3.0d0)/rh**(4.0d0/3.0d0) - (5.0d0/3.0d0)/rh )*gg*graddens

                END SUBROUTINE gradG2

                FUNCTION lapG2(g,rha,rhb,grha,grhb,lrha,lrhb)
                        ! Calculates the laplacian of the function G2. 
                        ! lapdens = laplacian of the chargedensity
                        IMPLICIT NONE
                        DOUBLE PRECISION :: lapG2
                        DOUBLE PRECISION :: g,rha,rhb,grha(3),grhb(3),lrha,lrhb
                        DOUBLE PRECISION :: gradf(3),gradg(3),DF22,Drh2
                        DOUBLE PRECISION, PARAMETER :: c = 0.25330d0
                        DOUBLE PRECISION :: fa,fb,fc,ff,gg,lf2,rh,graddens(3),lapdens
                        
                        rh = rha + rhb
                        fa =  5.0d0/(3.0d0*rh**2) -(4.0d0*c/9.0d0)*rh**(-7.0d0/3.0d0)
                        fb = (c/3.0d0)*rh**(-4.0d0/3.0d0) - 5.0d0/(3.0d0*rh)

                        ff = F2(g,rh)
                        gg = G2(g,rh)
                        lapdens = lrha + lrhb
                        graddens = grha + grhb
                        lf2 = lapF2(g,rha,rhb,grha,grhb,lrha,lrhb)

                        Drh2 = DOT_PRODUCT(graddens,graddens)
                        
                        CALL gradF2(g,rha,rhb,grha,grhb,gradf)
                        CALL gradG2(g,rha,rhb,grha,grhb,gradg)

                        DF22 =  DOT_PRODUCT(gradf,gradf)
                        
                        !ALT2:
                        lapG2 = lf2*rh**(-5.0d0/3.0d0) + (1.0d0/3.0d0)*rh**(-3)*( c - 5.0d0*rh**(1.0d0/3.0d0) )*( 2.0d0*DOT_PRODUCT(gradf,graddens) + ff*lapdens )
                        lapG2 = lapG2 + (ff/9.0d0)*rh**(-4)*( rh**(1.0d0/3.0d0)*( 40.0d0 + c**2*rh*(-2.0d0/3.0d0) ) - 14.0d0*c )*Drh2
                        lapG2 = lapG2*EXP(-c*rh**(-1.0d0/3.0d0))
                        
                        !ALT1:
                        !lapG2 = (lf2/ff)*gg - (DF22/(ff**2))*gg + DOT_PRODUCT(gradf,gradg)/ff
                        !lapG2 = lapG2 + fa*gg*Drh2 + fb*( DOT_PRODUCT(gradg,graddens) + gg*lapdens )
                        
                END FUNCTION lapG2
                
                FUNCTION gammus(rha,rhb)
                        ! calculates the function gamma as defined by 
                        ! Eqn (13) in Phys. Rev. B. 37, 785 (1988)
                        IMPLICIT NONE
                        DOUBLE PRECISION  :: gammus 
                        DOUBLE PRECISION  :: rha,rhb
                        gammus = 2.0d0*( 1.0d0 - (rha**2 + rhb**2)/((rha+rhb)**2) )
                END FUNCTION gammus

                SUBROUTINE gradgammus(rha,rhb,grha,grhb,grad)
                        ! calculates the gradient of the function gamma as defined by
                        ! Eqn (13) in Phys. Rev. B. 37, 785 (1988)
                        IMPLICIT NONE
                        DOUBLE PRECISION, INTENT(IN)  :: rha,rhb,grha(3),grhb(3)
                        DOUBLE PRECISION, INTENT(OUT)  :: grad(3) 
                        DOUBLE PRECISION :: rh,grh(3)

                        grh = grha + grhb
                        rh  = rha + rhb

                        grad = 4.0d0*( (1.0d0 - gammus(rha,rhb)/2.0d0 )*(grh/rh) - ( rha*grha + rhb*grhb )/(rh**2) )

                END SUBROUTINE gradgammus

                FUNCTION laplacegammus(rha,rhb,grha,grhb,lrha,lrhb)
                        ! calculates the laplacian of the function gamma as defined by
                        ! Eqn (13) in Phys. Rev. B. 37, 785 (1988)
                        IMPLICIT NONE
                        DOUBLE PRECISION  :: laplacegammus
                        DOUBLE PRECISION :: rha,rhb,grha(3),grhb(3),lrha,lrhb
                        DOUBLE PRECISION :: rh,grh(3),d2,d2a,d2b,lrh,gradg(3)

                        grh = grha + grhb
                        rh  = rha + rhb
                        lrh = lrha + lrhb

                        d2 = (1.0d0/(rh**2))*DOT_PRODUCT(grh,grh)
                        d2a = DOT_PRODUCT(grha,grha)
                        d2b = DOT_PRODUCT(grhb,grhb)
                        CALL gradgammus(rha,rhb,grha,grhb,gradg)

                        laplacegammus = 4.0d0*( (1.0d0 - gammus(rha,rhb)/2.0d0)*( d2 + lrh/rh ) - (1.0d0/rh)**2*(d2a + d2b + rha*lrha + rhb*lrhb ) - DOT_PRODUCT(gradg,grh)/rh  )

                END FUNCTION laplacegammus
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
                
                FUNCTION EUEGC(rs,zeta)
                        !---------------------------------------------------------!
                        ! The LDA correlation energy of the uniform electron gas  !
                        ! as given by Eqn (1) in Phys. Rev B, 42, 3205 (1990)     !
                        ! using Vosko, Wilks and Nusair's original                !
                        ! parametrizatin/interpolation as given in Can. J. Phys.  !
                        ! 58, 1200 (1980)                                         ! 
                        !---------------------------------------------------------!
                        IMPLICIT NONE
                        DOUBLE PRECISION :: EUEGC
                        DOUBLE PRECISION :: rs,zeta
                        
                        EUEGC = epsPC(rs) + depsC(rs,zeta)
                
                END FUNCTION EUEGC

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
                        
                        n = rha+rhb
                        
                        IF ( n .GT. 1.0E-20 ) THEN
                                zeta = (rha - rhb)/n
                                rs = (3.0d0/(4.0d0*PI*n))**(1.0d0/3.0d0)
                        ELSE
                                zeta = 0.0d0
                                rs = 1.0E20
                        ENDIF

                        IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
                        IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
                        
                        PHH = 0.50d0*( (1+zeta)**(2.0d0/3.0d0) + (1-zeta)**(2.0d0/3.0d0) )

                        kf = (3.0d0*PI**2*n)**(1.0d0/3.0d0)
                        
                        ks = sqrt(4.0d0*kf/PI)
                        
                        IF ( n .GT. 1.0E-20 ) THEN
                                t = lgrh/(2.0d0*PHH*ks*n)
                        ELSE
                                t = 1.0E20
                        ENDIF
                        
                        IF (  EXP(-EUEGC(rs,zeta)/(gg*PHH**3)) .NE. 1.0d0   ) THEN
                                A = (bb/gg)/( EXP(-EUEGC(rs,zeta)/(gg*PHH**3)) - 1.0d0 )
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
                        DOUBLE PRECISION :: kf,s,rs,n,FXX,funkz,zeta
                        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
                        DOUBLE PRECISION, PARAMETER :: mm = 0.219510d0
                        DOUBLE PRECISION, PARAMETER :: kk = 1.245 ! Here we actually choose the kappa value of the revPBE ( PRL 80, 890 (1998) )
                                                                  ! instead of kk = 0.804, as is used in the original PBE (Phys. Rev. Lett. 77, 3865 (1996) )
                        
                        n = rha + rhb

                        IF ( n .GT. 1.0E-20 ) THEN
                                rs = (3.0d0/(4.0d0*PI*n))**(1.0d0/3.0d0)
                                zeta = (rha-rhb)/n
                        ELSE
                                rs = 1.0E20
                                zeta = 0.0d0
                        ENDIF

                        IF ( zeta .LT. 0.0d0 ) zeta = 0.0d0
                        IF ( zeta .GT. 1.0d0 ) zeta = 1.0d0
                        
                        kf = (1.0d0/rs)*(9.0d0*PI/4.0d0)**(1.0d0/3.0d0)
                        
                        IF ( kf*n .GT. 1.0E-20 ) THEN
                                s = lgrh/(2.0d0*kf*n)
                        ELSE
                                s = 1.0E20
                        ENDIF
                        
                        !======================================
                        ! The non-polarized exchange density/n
                        !======================================
                        EXPBE = -3.0d0*kf/(4.0d0*PI)
                        
                        !=================================
                        ! The polarized exchange density/n
                        !=================================
                        !funkz  = 0.50d0*( (1+zeta)**(4.0d0/3.0d0) + (1-zeta)**(4.0d0/3.0d0) )
                        !EXPBE  = - (3.0d0/4.0d0)*(3.0d0*n/PI)**(1.0d0/3.0d0)*funkz
                        !==================================
                        ! End of plarized energy expression
                        !==================================

                        !==========================================
                        ! This part is independent of polarization:
                        !==========================================
                        FXX = 1.0d0 + kk - kk/(1.0d0 + (mm/kk)*s**2 )
                        
                        EXPBE = EXPBE*FXX
                        
               END FUNCTION EXPBE

END MODULE exchcorrmodule
