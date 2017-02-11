SUBROUTINE vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd,V)
        ! this subroutine calculates the exchange correlation 
        ! potential. This is the alternative version to be used 
        ! to calculate the gradient of the exchange correlation potential 
        ! with respect to nuclear coordinates.
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(IN) :: densu,densd,gdensu(3),gdensd(3),lgrh,ldensu,ldensd
        DOUBLE PRECISION, INTENT(OUT) :: V(2)
        DOUBLE PRECISION, EXTERNAL :: rho
        DOUBLE PRECISION :: gdens(3)
        DOUBLE PRECISION :: Vcc(2),Vex(2)
        
        SELECT CASE (CORRLEVEL)
                CASE ('LDA')
                        CALL VWN(densu,densd,Vcc)
                        CALL Vx(densu,densd,Vex)
                        V = Vcc + Vex

                CASE('PBE')
                        CALL VPBE(densu,densd,lgrh,Vcc)
                        V = Vcc

                CASE('B3LYP')
                        !------------------------------------------------------------
                        ! remember that 20% of the non-local exact exchange potential
                        ! is still missing from the B3LYP exchange-correlation.
                        !------------------------------------------------------------
                        CALL VB3LYP(densu,densd,gdensu,gdensd,ldensu,ldensd,Vcc)
                        V = Vcc

                CASE DEFAULT
                        CALL VPBE(densu,densd,lgrh,Vcc)
                        V = Vcc
        END SELECT

END SUBROUTINE vxcalt

