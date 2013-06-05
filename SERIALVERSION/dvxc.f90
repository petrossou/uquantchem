SUBROUTINE dvxc(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd,dV)
        ! this subroutine calculates the derivative of the exchange correlation 
        ! potential of the spin up and spin down chanel with respect to:
        ! (1) (a) the spin up density for the spin up chanel     V(1,1,1) =  dV(1)/d(densu)
        !     (b) the spin down density for the spin up chanel   V(1,2,1) =  dV(1)/d(densd)
        !     (c) the spin up density for the spin down chanel   V(2,1,1) =  dV(2)/d(densu)
        !     (d) the spin down density for the spin down chanel V(2,2,1) =  dV(2)/d(densd)
        ! (2) (a) the x-derivative of the spin up density for the spin up chanel     V(1,1,2) =  dV(1)/d(gdensu(1))
        !     (b) the x-derivative of the spin down density for the spin up chanel   V(1,2,2) =  dV(1)/d(gdensd(1))
        !     (c) the x-derivative of the spin up density for the spin down chanel   V(2,1,2) =  dV(2)/d(gdensu(1))
        !     (d) the x-derivative of the spin down density for the spin down chanel V(2,2,2) =  dV(2)/d(gdensd(1))
        ! (3) (a) the y-derivative of the spin up density for the spin up chanel     V(1,1,3) =  dV(1)/d(gdensu(2))
        !     (b) the y-derivative of the spin down density for the spin up chanel   V(1,2,3) =  dV(1)/d(gdensd(2))
        !     (c) the y-derivative of the spin up density for the spin down chanel   V(2,1,3) =  dV(2)/d(gdensu(2))
        !     (d) the y-derivative of the spin down density for the spin down chanel V(2,2,3) =  dV(2)/d(gdensd(2))
        ! (4) (a) the z-derivative of the spin up density for the spin up chanel     V(1,1,4) =  dV(1)/d(gdensu(3))
        !     (b) the z-derivative of the spin down density for the spin up chanel   V(1,2,4) =  dV(1)/d(gdensd(3))
        !     (c) the z-derivative of the spin up density for the spin down chanel   V(2,1,4) =  dV(2)/d(gdensu(3))
        !     (d) the z-derivative of the spin down density for the spin down chanel V(2,2,4) =  dV(2)/d(gdensd(3))
        ! (5) (a) the length of the total gradient for the spin up chanel            V(1,1,5) =  V(1,2,5) = dV(1)/d(lgrh)
        !     (b) the length of the total gradient for the spin down chanel          V(2,1,5) =  V(2,2,5) = dV(2)/d(lgrh)
        ! (6) (a) laplacian of the spin up density for the spin up chanel            V(1,1,6) =  dV(1)/d(ldensu)
        !     (b) laplacian of the spin down density for the spin up chanel          V(1,2,6) =  dV(1)/d(ldensd)
        !     (c) laplacian of the spin up density for the spin down chanel          V(2,1,6) =  dV(2)/d(ldensu)
        !     (d) laplacian of the spin down density for the spin down chanel        V(2,2,6) =  dV(2)/d(ldensd)
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(IN) :: densu,densd,gdensu(3),gdensd(3),lgrh,ldensu,ldensd
        DOUBLE PRECISION, INTENT(OUT) :: dV(2,2,6)
        DOUBLE PRECISION :: gdens(3)
        DOUBLE PRECISION :: Vcc(2),Vex(2)
        DOUBLE PRECISION :: V0(2),Vp1(2),Vm1(2),Vp2(2),Vm2(2),hh,xv(3),yv(3),zv(3)
        
        xv(:) = 0.0d0
        yv(:) = 0.0d0
        zv(:) = 0.0d0
        
        xv(1) = 1.0d0
        yv(2) = 1.0d0
        zv(3) = 1.0d0
        
        dV(:,:,:) = 0.0d0

        SELECT CASE (CORRLEVEL)
                CASE ('LDA')
                       ! derivative with respect to the spin up density
                       IF ( densu .LT. 1.0d0 ) THEN
                               hh = densu*0.00010d0
                       ELSE
                               hh = 0.00010d0
                       ENDIF
                       CALL vxcalt(CORRLEVEL,densu+hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                       CALL vxcalt(CORRLEVEL,densu-hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                       CALL vxcalt(CORRLEVEL,densu+2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                       CALL vxcalt(CORRLEVEL,densu-2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                       dV(1,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                       dV(2,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                       
                       ! derivative with respect to the spin down density
                       IF ( densd .LT. 1.0d0 ) THEN
                               hh = densd*0.00010d0
                       ELSE
                               hh = 0.00010d0
                       ENDIF
                       CALL vxcalt(CORRLEVEL,densu,densd+hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                       CALL vxcalt(CORRLEVEL,densu,densd-hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                       CALL vxcalt(CORRLEVEL,densu,densd+2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                       CALL vxcalt(CORRLEVEL,densu,densd-2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                       dV(1,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                       dV(2,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )

                CASE('PBE')
                        
                        ! derivative with respect to the spin up density
                        IF ( densu .LT. 1.0d0 ) THEN
                               hh = densu*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu+hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu-hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu+2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu-2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the spin down density
                        IF ( densd .LT. 1.0d0 ) THEN
                               hh = densd*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd+hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd-hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd+2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd-2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the length of the total gradient
                        IF ( lgrh .LT. 1.0d0 ) THEN
                                hh = lgrh*0.00010d0
                        ELSE
                                hh = 0.00010d0
                        ENDIF
                        !print*,'(4)',lgrh
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh+hh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh-hh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh+2.0d0*hh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh-2.0d0*hh,ldensu,ldensd,Vm2)
                        dV(1,1,5) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(1,2,5) = dV(1,1,5)
                        dV(2,1,5) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        dV(2,2,5) = dV(2,1,5)

                CASE('B3LYP')
                        !------------------------------------------------------------
                        ! remember that 20% of the non-local exact exchange potential
                        ! is still missing from the B3LYP exchange-correlation.
                        !------------------------------------------------------------
                        
                        ! derivative with respect to the spin up density
                        IF ( densu .LT. 1.0d0 ) THEN
                               hh = densu*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu+hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu-hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu+2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu-2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the spin down density
                        IF ( densd .LT. 1.0d0 ) THEN
                               hh = densd*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd+hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd-hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd+2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd-2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the x-component of the spin up density gradient
                        IF ( abs(gdensu(1)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensu(1))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+hh*xv,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-hh*xv,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+2.0d0*hh*xv,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-2.0d0*hh*xv,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,2) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,2) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the x-component of the spin down density gradient
                        IF ( abs(gdensd(1)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensd(1))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+hh*xv,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-hh*xv,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+2.0d0*hh*xv,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-2.0d0*hh*xv,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,2) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,2) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the y-component of the spin up density gradient
                        IF ( abs(gdensu(2)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensu(2))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+hh*yv,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-hh*yv,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+2.0d0*hh*yv,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-2.0d0*hh*yv,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,3) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,3) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the y-component of the spin down density gradient
                        IF ( abs(gdensd(2)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensd(2))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+hh*yv,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-hh*yv,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+2.0d0*hh*yv,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-2.0d0*hh*yv,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,3) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,3) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the z-component of the spin up density gradient
                        IF ( abs(gdensu(3)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensu(3))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+hh*zv,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-hh*zv,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu+2.0d0*hh*zv,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu-2.0d0*hh*zv,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,4) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,4) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the the z-component of the spin down density gradient
                        IF ( abs(gdensd(3)) .LT. 1.0d0 ) THEN
                               hh = abs(gdensd(2))*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+hh*zv,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-hh*zv,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd+2.0d0*hh*zv,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd-2.0d0*hh*zv,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,4) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,4) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the laplacian of the spin up density gradient
                        IF ( abs(ldensu) .LT. 1.0d0 ) THEN
                               hh = abs(ldensu)*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu+hh,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu-hh,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu+2.0d0*hh,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu-2.0d0*hh,ldensd,Vm2)
                        dV(1,1,6) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,6) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the laplacian of the spin down density gradient
                        IF ( abs(ldensd) .LT. 1.0d0 ) THEN
                               hh = abs(ldensd)*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd+hh,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd-hh,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd+2.0d0*hh,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd-2.0d0*hh,Vm2)
                        dV(1,2,6) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,6) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                CASE DEFAULT
                        ! PBE = default 

                        ! derivative with respect to the spin up density
                        IF ( densu .LT. 1.0d0 ) THEN
                               hh = densu*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu+hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu-hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu+2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu-2.0d0*hh,densd,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,1,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the spin down density
                        IF ( densd .LT. 1.0d0 ) THEN
                               hh = densd*0.00010d0
                        ELSE
                               hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd+hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd-hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd+2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd-2.0d0*hh,gdensu,gdensd,lgrh,ldensu,ldensd,Vm2)
                        dV(1,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(2,2,1) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        
                        ! derivative with respect to the length of the total gradient
                        IF ( lgrh .LT. 1.0d0 ) THEN
                                hh = lgrh*0.00010d0
                        ELSE
                                hh = 0.00010d0
                        ENDIF
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh+hh,ldensu,ldensd,Vp1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh-hh,ldensu,ldensd,Vm1)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh+2.0d0*hh,ldensu,ldensd,Vp2)
                        CALL vxcalt(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh-2.0d0*hh,ldensu,ldensd,Vm2)
                        dV(1,1,5) = (1.0d0/(12.0d0*hh))*(Vm2(1) - 8.0d0*Vm1(1) + 8.0d0*Vp1(1) - Vp2(1) )
                        dV(1,2,5) = dV(1,1,5)
                        dV(2,1,5) = (1.0d0/(12.0d0*hh))*(Vm2(2) - 8.0d0*Vm1(2) + 8.0d0*Vp1(2) - Vp2(2) )
                        dV(2,2,5) = dV(2,1,5)

        END SELECT

END SUBROUTINE dvxc

