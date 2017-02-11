SUBROUTINE vxc(CORRLEVEL,BAS,Pup,Pdown,r,V)
        ! this subroutine calculates the exchange correlation 
        ! potential.
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION, INTENT(OUT) :: V(2)
        DOUBLE PRECISION, EXTERNAL :: rho
        DOUBLE PRECISION :: dens,gdens(3),ldensu,ldensd,densu,densd,gdensu(3),gdensd(3)
        DOUBLE PRECISION :: Vcc(2),Vex(2),P(BAS%NBAS,BAS%NBAS),lgrh
        
        densu = rho(BAS,Pup,r)
        densd = rho(BAS,Pdown,r)
        
        dens = densd + densu
        V = 0.0d0

        IF ( dens .GT. 1.0E-20 )  THEN
                SELECT CASE (CORRLEVEL)
                        CASE ('LDA')
                                CALL VWNC(densu,densd,Vcc)
                                CALL Vx(densu,densd,Vex)
                                V = Vcc + Vex

                        CASE('PBE')
                                P = Pup + Pdown
                                CALL gradrho(BAS,P,r,gdens)
                                lgrh = sqrt(DOT_PRODUCT(gdens,gdens))
                                !print*,'(1)',lgrh
                                CALL VPBE(densu,densd,lgrh,Vcc)
                                V = Vcc

                        CASE('B3LYP')
                                CALL gradrho(BAS,Pup,r,gdensu)
                                CALL gradrho(BAS,Pdown,r,gdensd)
                                !----------------------------------
                                ! The Laplacian is no longer needed
                                ! since we use the expression from 
                                ! Chem. Phys. Lett. 157, 200 (1989)
                                !-----------------------------------
                                !ldensu = laprho(BAS,Pup,r)
                                !ldensd = laprho(BAS,Pdown,r)
                                ldensu = 0.0d0
                                ldensd = 0.0d0
                                !------------------------------------------------------------
                                ! remember that 20% of the non-local exact exchange potential
                                ! is still missing from the B3LYP exchange-correlation.
                                !------------------------------------------------------------
                                CALL VB3LYP(densu,densd,gdensu,gdensd,ldensu,ldensd,Vcc)
                                V = Vcc

                        CASE DEFAULT
                                P = Pup + Pdown
                                CALL gradrho(BAS,P,r,gdens)
                                lgrh = sqrt(DOT_PRODUCT(gdens,gdens))
                                CALL VPBE(densu,densd,lgrh,Vcc)
                                V = Vcc
                END SELECT
         ELSE
                        V = 0.0d0
         ENDIF

END SUBROUTINE vxc

