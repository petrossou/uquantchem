SUBROUTINE dftforcedensnip(CORRLEVEL,NATOMS,BAS,Pup,Pdown,r,BECKECENTER,excforcedens)
        !===================================================================!
        ! This subroutine calculates the exchange correlation force density !
        ! at the point r in space.                                          !
        !===================================================================!
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NATOMS,BECKECENTER
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION, INTENT(OUT) :: excforcedens(NATOMS,3)
        DOUBLE PRECISION :: dens,densu,densd,gdensu(3),gdensd(3),gdens(3),lengthgdu,lengthgdd,lengthgd,ngrho(NATOMS,3),ggb3lyp(3)
        DOUBLE PRECISION :: ngrhou(NATOMS,3),ngrhod(NATOMS,3),V(2),Vcc(2),Vex(2),hru(3,3,NATOMS),hrd(3,3,NATOMS),hr(3,3,NATOMS),gg
        DOUBLE PRECISION, EXTERNAL :: rho
        INTEGER :: I,J,n,m

        densu = rho(BAS,Pup,r)
        densd = rho(BAS,Pdown,r)

        dens = densd + densu
        excforcedens = 0.0d0

        ! Calculation of the nuclear gradient of the charge densities
        CALL nucgradrhonip(NATOMS,BAS,Pup,r,BECKECENTER,ngrhou)
        CALL nucgradrhonip(NATOMS,BAS,Pdown,r,BECKECENTER,ngrhod)
        ngrho = ngrhou + ngrhod
        
        IF ( dens .GT. 1.0E-20 )  THEN
                SELECT CASE (CORRLEVEL)
                        CASE ('LDA')
                                CALL VWNC(densu,densd,Vcc)
                                CALL Vx(densu,densd,Vex)
                                V = Vcc + Vex
                                
                                excforcedens =  -V(1)*ngrhou - V(2)*ngrhod

                        CASE('PBE')
                                CALL gradrho(BAS,Pup+Pdown,r,gdens)
                                lengthgd = sqrt(DOT_PRODUCT(gdens,gdens))
                                
                                CALL VPBE(densu,densd,lengthgd,Vcc)
                                V = Vcc
                                
                                ! Calculating the hessian of the charge density
                                CALL hessianrhonip(NATOMS,BAS,Pup,r,BECKECENTER,hru)
                                CALL hessianrhonip(NATOMS,BAS,Pdown,r,BECKECENTER,hrd)

                                hr = hru + hrd
                                
                                ! Calculating the derivative of the PBE exchange
                                ! correlation energy density  with respect to 
                                ! the lenght of the total density gradient
                                gg = gVPBE(densu,densd,lengthgd)
                                
                                excforcedens = -V(1)*ngrhou -V(2)*ngrhod

                                DO I=1,NATOMS
                                        excforcedens(I,:) = excforcedens(I,:) - gg*MATMUL(hr(:,:,I),gdens)/lengthgd
                                ENDDO

                        CASE('B3LYP')
                                CALL gradrho(BAS,Pup,r,gdensu)
                                CALL gradrho(BAS,Pdown,r,gdensd)
                                CALL gradrho(BAS,Pup+Pdown,r,gdens)
                                lengthgdu = sqrt(DOT_PRODUCT(gdensu,gdensu))
                                lengthgdd = sqrt(DOT_PRODUCT(gdensd,gdensd))
                                lengthgd  = sqrt(DOT_PRODUCT(gdens,gdens))
                                
                                CALL VB3LYP(densu,densd,gdensu,gdensd,0.0d0,0.0d0,Vcc)
                                V = Vcc
                                
                                ! Calculating the hessian of the charge density
                                CALL hessianrhonip(NATOMS,BAS,Pup,r,BECKECENTER,hru)
                                CALL hessianrhonip(NATOMS,BAS,Pdown,r,BECKECENTER,hrd)
                                
                                hr = hru + hrd

                                ! Calculating the derivative of the PBE exchange
                                ! correlation energy density  with respect to 
                                ! the lenght of the total density gradient
                                CALL gVB3LYP(densu,densd,gdensu,gdensd,ggb3lyp)
                                
                                excforcedens = - V(1)*ngrhou - V(2)*ngrhod
                                
                                
                                DO I=1,NATOMS
                                        excforcedens(I,:) = excforcedens(I,:) - ggb3lyp(1)*MATMUL(hru(:,:,I),gdensu)/lengthgdu
                                        excforcedens(I,:) = excforcedens(I,:) - ggb3lyp(2)*MATMUL(hrd(:,:,I),gdensd)/lengthgdd
                                        excforcedens(I,:) = excforcedens(I,:) - ggb3lyp(3)*MATMUL(hr(:,:,I),gdens)/lengthgd
                                ENDDO

                        END SELECT
                ELSE
                        excforcedens = 0.0d0
                ENDIF

END SUBROUTINE dftforcedensnip

