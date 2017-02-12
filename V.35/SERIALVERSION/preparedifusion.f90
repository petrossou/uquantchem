SUBROUTINE preparedifusion(atomnr,rnuc,DT,vlength,vvect,zvect,z,zetaa,vmeanvect,aa,vmeanz,vmeanrho,rhovect,zbis,rhobis,dr,qtilde)
        !--------------------------------------------------------------
        ! All the numbers in parenthesis to the right are in  reference 
        ! to the numbering of the steps in Umrigars algorithm, 
        ! of the so called improved algorithm on pages 2886-2887 in
        ! J. Chem. Phys. 99, 2865-2890, (1993)
        !--------------------------------------------------------------
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: atomnr
        DOUBLE PRECISION, INTENT(IN) :: rnuc(3),DT,vlength,vvect(3),zvect(3),z
        DOUBLE PRECISION, INTENT(OUT) :: aa,vmeanvect(3),vmeanz,vmeanrho,rhovect(3),zbis,rhobis,dr(3),zetaa
        REAL, INTENT(OUT) :: qtilde 

        zetaa = sqrt( atomnr**2 + 1.0d0/DT )     ! (12)

        aa = 0.50d0*(1.0d0 + (1.0d0/vlength)*DOT_PRODUCT(vvect,zvect) ) + (atomnr*z)**2/(10.0d0*( 4.0d0 + (atomnr*z)**2 ) ) ! (14)

        ! Calculating the mean drift velocity (row below 14)
        vmeanvect = ( -1.0d0 + sqrt(1.0d0+2.0d0*aa*DT*vlength**2) )*vvect/(aa*DT*vlength**2)

        ! Projecting out the components <v>z znd <v>rho (15)
        vmeanz = DOT_PRODUCT(vmeanvect,zvect)
        vmeanrho = sqrt(DOT_PRODUCT(vmeanvect - vmeanz*zvect,vmeanvect - vmeanz*zvect))
        
        ! Calculating the rho unit vector

        rhovect = 0.0d0

        IF ( vmeanrho .NE. 0.0d0 ) rhovect = (vmeanvect - vmeanz*zvect)/vmeanrho
        
        ! Calculating z''  (16)
        zbis = z + vmeanz*DT
        IF ( zbis .LT. 0.0d0 ) zbis = 0.0d0

        ! calculating rhobis  (17)
        rhobis = 2.0d0*vmeanrho*DT*zbis/( z + zbis )

        ! calculating d(r) (18)
        dr = rnuc + rhobis*rhovect + zbis*zvect

        ! Calculating qtilde (19)
        qtilde = 0.50*ERFC(REAL((z+vmeanz*DT)/sqrt(2.0d0*DT)))
        
        END SUBROUTINE preparedifusion
