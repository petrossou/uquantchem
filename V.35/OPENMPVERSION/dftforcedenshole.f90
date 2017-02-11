SUBROUTINE dftforcedenshole(CORRLEVEL,NATOMS,BAS,Pup,Pdown,gradS,r,BECKECENTER,excforcedens,IORBNR,PHOLE)
        !===================================================================!
        ! This subroutine calculates the exchange correlation force density !
        ! at the point r in space.                                          !
        !===================================================================!
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NATOMS,BECKECENTER,IORBNR
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION, INTENT(OUT) :: excforcedens(NATOMS,3),PHOLE(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: dens,densu,densd,gdensu(3),gdensd(3),gdens(3),lengthgdu,lengthgdd,lengthgd,ngrho(NATOMS,3),ggb3lyp(3)
        DOUBLE PRECISION :: ngrhou(NATOMS,3),ngrhod(NATOMS,3),V(2),Vcc(2),Vex(2),hru(3,3,NATOMS),hrd(3,3,NATOMS),hr(3,3,NATOMS),gg,ngrhole(NATOMS,3)
        DOUBLE PRECISION, EXTERNAL :: rho
        INTEGER :: I,J,n,m

        IF ( IORBNR .GT. 0 ) THEN
                        densu = rho(BAS,Pup-PHOLE,r)
                        densd = rho(BAS,Pdown,r)
        ELSE      
                        densu = rho(BAS,Pup,r)
                        densd = rho(BAS,Pdown-PHOLE,r)
        ENDIF

        dens = densd + densu
        excforcedens = 0.0d0

        ! Calculation of the nuclear gradient of the charge densities
        IF ( IORBNR .GT. 0 ) THEN
                        CALL nucgradrhohole(NATOMS,BAS,Pup,PHOLE,gradS,r,BECKECENTER,ngrhou)
                        CALL nucgradrho(NATOMS,BAS,Pdown,gradS,r,BECKECENTER,ngrhod)
        ELSE
                        CALL nucgradrho(NATOMS,BAS,Pup,gradS,r,BECKECENTER,ngrhou)
                        CALL nucgradrhohole(NATOMS,BAS,Pdown,PHOLE,gradS,r,BECKECENTER,ngrhod)
        ENDIF

        ngrho = ngrhou + ngrhod
        
        IF ( dens .GT. 1.0E-20 )  THEN
                SELECT CASE (CORRLEVEL)
                        CASE ('LDA')
                                CALL VWNC(densu,densd,Vcc)
                                CALL Vx(densu,densd,Vex)
                                V = Vcc + Vex
                                
                                excforcedens =  -V(1)*ngrhou - V(2)*ngrhod

                        CASE('PBE')
                                CALL gradrho(BAS,Pup+Pdown-PHOLE,r,gdens)
                                
                                lengthgd = sqrt(DOT_PRODUCT(gdens,gdens))
                                
                                CALL VPBE(densu,densd,lengthgd,Vcc)
                                V = Vcc
                                
                                ! Calculating the hessian of the charge density
                                IF ( IORBNR .GT. 0 ) THEN
                                                CALL hessianrhohole(NATOMS,BAS,Pup,PHOLE,gradS,r,BECKECENTER,hru)
                                                CALL hessianrho(NATOMS,BAS,Pdown,gradS,r,BECKECENTER,hrd)
                                ELSE
                                                CALL hessianrho(NATOMS,BAS,Pup,gradS,r,BECKECENTER,hru)
                                                CALL hessianrhohole(NATOMS,BAS,Pdown,PHOLE,gradS,r,BECKECENTER,hrd)
                                ENDIF

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
                                IF ( IORBNR .GT. 0  ) THEN
                                                CALL gradrho(BAS,Pup-PHOLE,r,gdensu)
                                                CALL gradrho(BAS,Pdown,r,gdensd)
                                ELSE
                                                CALL gradrho(BAS,Pup,r,gdensu)
                                                CALL gradrho(BAS,Pdown-PHOLE,r,gdensd)
                                ENDIF
                                
                                CALL gradrho(BAS,Pup+Pdown-PHOLE,r,gdens)
                       
                                lengthgdu = sqrt(DOT_PRODUCT(gdensu,gdensu))
                                lengthgdd = sqrt(DOT_PRODUCT(gdensd,gdensd))
                                lengthgd  = sqrt(DOT_PRODUCT(gdens,gdens))
                                
                                CALL VB3LYP(densu,densd,gdensu,gdensd,0.0d0,0.0d0,Vcc)
                                V = Vcc
                                
                                ! Calculating the hessian of the charge density
                                IF ( IORBNR .GT. 0 ) THEN
                                                CALL hessianrhohole(NATOMS,BAS,Pup,PHOLE,gradS,r,BECKECENTER,hru)
                                                CALL hessianrho(NATOMS,BAS,Pdown,gradS,r,BECKECENTER,hrd)
                                ELSE
                                                CALL hessianrho(NATOMS,BAS,Pup,gradS,r,BECKECENTER,hru)
                                                CALL hessianrhohole(NATOMS,BAS,Pdown,PHOLE,gradS,r,BECKECENTER,hrd)
                                ENDIF
                                
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

END SUBROUTINE dftforcedenshole

