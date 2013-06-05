SUBROUTINE nucgradvxc(CORRLEVEL,NATOMS,BAS,Pup,Pdown,gradS,r,gradV)
        ! this subroutine calculates the nuclear gradient of the exchange correlation 
        ! potential (for both spin chanells).
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NATOMS
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION, INTENT(OUT) :: gradV(NATOMS,2,3)
        DOUBLE PRECISION, EXTERNAL :: rho
        DOUBLE PRECISION :: densu,densd,gdensu(3),gdensd(3),gdens(3),ldensu,ldensd
        DOUBLE PRECISION :: P(BAS%NBAS,BAS%NBAS),lgrh,dV(2,2,6)
        DOUBLE PRECISION :: grhou(NATOMS,3),ggrhou(NATOMS,3,3),glaplrhu(NATOMS,3) ! spin-up chanell
        DOUBLE PRECISION :: grhod(NATOMS,3),ggrhod(NATOMS,3,3),glaplrhd(NATOMS,3) ! spin-down chanell
        INTEGER :: I,J,n

        ! Calculating the spin up and down charge densities
        densu = rho(BAS,Pup,r)
        densd = rho(BAS,Pdown,r)
        
        ! Calculating the gradient of the total charge density
        IF ( CORRLEVEL .EQ. 'PBE' ) THEN
                P = Pup + Pdown
                CALL gradrho(BAS,P,r,gdens)
                lgrh = sqrt(DOT_PRODUCT(gdens,gdens))
                !print*,'(5)',lgrh
        ELSE
                lgrh = 1.0d0
        ENDIF
        
        ! Calculatin the gradients of the spin up and down charge densities 
        IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                CALL gradrho(BAS,Pup,r,gdensu)
                CALL gradrho(BAS,Pdown,r,gdensd)
        ELSE
                gdensu = 0.0d0
                gdensd = 0.0d0
        ENDIF
        
        ! Calculating the laplacians of the spin up and down charge densities
        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                ldensu = laprho(BAS,Pup,r)
                ldensd = laprho(BAS,Pdown,r)
        ELSE
                ldensu = 0.0d0
                ldensd = 0.0d0
        ENDIF
        
        
        ! Calculating the various derivatives of the exchange correlation potential
        CALL dvxc(CORRLEVEL,densu,densd,gdensu,gdensd,lgrh,ldensu,ldensd,dV)
        
        ! Calculating the nuclear gradient of the spin up and down charge densities
        CALL nucgradrho(NATOMS,BAS,Pup,gradS,r,grhou)
        CALL nucgradrho(NATOMS,BAS,Pdown,gradS,r,grhod)

        ! Calculating the nuclear gradient of the gradient of the spin up and down charge densities
        IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                CALL nucgradgradrho(NATOMS,BAS,Pup,gradS,r,ggrhou)
                CALL nucgradgradrho(NATOMS,BAS,Pdown,gradS,r,ggrhod)
        ELSE
                ggrhou = 0.0d0
                ggrhod = 0.0d0
        ENDIF

        ! Calculating the nuclear gradients of the laplacians for the spin up and down charge densities
        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                CALL nucgradlaprho(NATOMS,BAS,Pup,gradS,r,glaplrhu)
                CALL nucgradlaprho(NATOMS,BAS,Pdown,gradS,r,glaplrhd)
        ELSE
                glaplrhu = 0.0d0
                glaplrhd = 0.0d0
        ENDIF

        ! The nuclear gradient of the spin-up and down chanell of the potential is calculated here:
        ! See Eqn (5) p. 123 in the black note book and description of output from dvxc.f90
        ! n=1 <=> Spin up and n=2 <=> spin down
        gradV(:,:,:) = 0.0d0

        DO J=1,NATOMS
                DO n=1,2
                        ! Partial derivatives involving the charge density
                        gradV(J,n,:) = dV(n,1,1)*grhou(J,:) + dV(n,2,1)*grhod(J,:)
                        ! Partial derivatives involving the x-component of the gradient of the charge density
                        gradV(J,n,:) = gradV(J,n,:) + dV(n,1,2)*ggrhou(J,1,:) + dV(n,2,2)*ggrhod(J,1,:)
                        ! Partial derivatives involving the y-component of the gradient of the charge density
                        gradV(J,n,:) = gradV(J,n,:) + dV(n,1,3)*ggrhou(J,2,:) + dV(n,2,3)*ggrhod(J,2,:)
                        ! Partial derivatives involving the z-component of the gradient of the charge density
                        gradV(J,n,:) = gradV(J,n,:) + dV(n,1,4)*ggrhou(J,3,:) + dV(n,2,4)*ggrhod(J,3,:)
                        ! Partial derivatives involving the length of the total density
                        gradV(J,n,:) = gradV(J,n,:) + dV(n,1,5)*(1.0d0/lgrh)*( ggrhou(J,1,:)+ggrhod(J,1,:)  +  ggrhou(J,2,:)+ggrhod(J,2,:)  +  ggrhou(J,3,:)+ggrhod(J,3,:) )
                        ! Partial derivatives involving the laplacian of the gradiant of the charge density
                        gradV(J,n,:) = gradV(J,n,:) + dV(n,1,6)*glaplrhu(J,:) + dV(n,2,6)*glaplrhd(J,:)
                ENDDO
        ENDDO
END SUBROUTINE nucgradvxc

