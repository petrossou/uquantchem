FUNCTION excdens(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd)
        ! this function calculates the exchange correlation energy density. 
        USE exchcorrmodule
        IMPLICIT NONE
        DOUBLE PRECISION :: excdens
        CHARACTER(LEN=20) :: CORRLEVEL
        DOUBLE PRECISION :: densu,densd,gdensu(3),gdensd(3),lgrh,ldensu,ldensd,dens
        DOUBLE PRECISION :: exchalt
        DOUBLE PRECISION :: gdens(3)
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        DOUBLE PRECISION :: rs,zeta,fz
        
        dens = densu+densd

        IF ( dens .NE. 0.0d0 ) THEN
                rs = (3.0d0/(4.0d0*PI*dens))**(1.0d0/3.0d0)
                zeta = (densu-densd)/dens
        ELSE
                zeta = 0.0d0
                rs = 1.0E20
        ENDIF
        
        IF ( dens .GT. 1.0E-20 ) THEN
                SELECT CASE (CORRLEVEL)
                        CASE ('LDA')
                                ! The correlation energy:
                                excdens = EUEGC(rs,zeta)*dens
                                ! Adding the LSDA exchange contribution:
                                ! Here we have written down the polarized version of the Exchange energy density
                                ! according to Eqn (26) in J.P Perdue and Y. Wang PRB, 45, 13244 (1992)
                                fz = 0.50d0*( (1.0d0 + zeta)**(4.0d0/3.0d0) + (1.0d0 - zeta)**(4.0d0/3.0d0) )
                                excdens = excdens - (3.0d0/4.0d0)*(3*dens/PI)**(1.0d0/3.0d0)*fz*dens
                        
                        CASE('PBE')
                                excdens = dens*( HPBE(densu,densd,lgrh) + EUEGC(rs,zeta) +  EXPBE(densu,densd,lgrh) )
                        

                        CASE('B3LYP')
                                !------------------------------------------------------------
                                ! remember that 20% of the non-local exact exchange potential
                                ! is still missing from the B3LYP exchange-correlation.
                                !------------------------------------------------------------
                                excdens = EB3LYP(densu,densd,gdensu,gdensd,ldensu,ldensd)

                        CASE DEFAULT
                                excdens = dens*( HPBE(densu,densd,lgrh) + EUEGC(rs,zeta) +  EXPBE(densu,densd,lgrh) )
                END SELECT
         ELSE
                 excdens = 0.0d0
         ENDIF

END FUNCTION excdens

