FUNCTION excdr(CORRLEVEL,BAS,Pup,Pdown,r)
        ! this subroutine calculates the exchange correlation energy density
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        DOUBLE PRECISION :: excdr
        CHARACTER(LEN=20) :: CORRLEVEL
        TYPE(BASIS)  :: BAS
        DOUBLE PRECISION :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION :: densu,densd,gdensu(3),gdensd(3),gdens(3),ldensu,ldensd
        DOUBLE PRECISION :: P(BAS%NBAS,BAS%NBAS),lgrh
        DOUBLE PRECISION, EXTERNAL :: rho,excdens
        
        densu = rho(BAS,Pup,r)
        densd = rho(BAS,Pdown,r)
        
        SELECT CASE (CORRLEVEL)
                CASE ('LDA')
                        excdr = excdens(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd)

                CASE('PBE')
                        P = Pup + Pdown
                        CALL gradrho(BAS,P,r,gdens)
                        lgrh = sqrt(DOT_PRODUCT(gdens,gdens))
                        excdr = excdens(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd)
                        

                CASE('B3LYP')
                        CALL gradrho(BAS,Pup,r,gdensu)
                        CALL gradrho(BAS,Pdown,r,gdensd)
                        
                        !---------------------------------------------------
                        ! Laplacian is no longer needed
                        ! since we use the Chem. Phys. Lett 157, 200 (1989)
                        !---------------------------------------------------
                        !ldensu = laprho(BAS,Pup,r)
                        !ldensd = laprho(BAS,Pdown,r)
                        ldensu = 0.0d0
                        ldensd = 0.0d0
                        !------------------------------------------------------------
                        ! remember that 20% of the non-local exact exchange potential
                        ! is still missing from the B3LYP exchange-correlation.
                        !------------------------------------------------------------
                        excdr = excdens(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd)
                        

                CASE DEFAULT
                        P = Pup + Pdown
                        CALL gradrho(BAS,P,r,gdens)
                        lgrh = sqrt(DOT_PRODUCT(gdens,gdens))
                        excdr = excdens(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd)
        END SELECT
        IF ( excdr .NE. excdr ) excdr = 0.0d0
END FUNCTION excdr

