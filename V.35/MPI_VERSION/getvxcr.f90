SUBROUTINE getvxcr(CORRLEVEL,NATOMS,BAS,Pup,Pdown,r,Vxcr)
        !===============================================================!
        ! This subroutine calculate the r-dependent matrix elements     !
        ! Phi_{i}(r)*Vxc(r)*Phi_{j}(r) and the nuclear gradient         !
        ! of Phi_{i}(r)*Vxc(r)*Phi_{j}(r). After integration over the   !
        ! entire 3-DIM space with some quadrature these will become     !
        ! the exchange correlation potential matrix elements to be      !
        ! plugged into the Kohn-sham hamiltonian, and the corresponding !
        ! force expression for the exchange correlation contribution    !
        ! to the interatomic forces.                                    !
        !===============================================================!
        USE datatypemodule
        USE exchcorrmodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NATOMS
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),r(3)
        DOUBLE PRECISION, INTENT(OUT) :: Vxcr(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: V(2),gradV(NATOMS,2,3),gradI(3),gradJ(3),FI,FJ,gFI(3),gFJ(3),gdens(3),lengthgd,gg,ggb3lyp(3),densu,densd
        DOUBLE PRECISION :: lengthgdu,lengthgdd,gdensu(3),gdensd(3),ALLFI(BAS%NBAS),gALLFI(BAS%NBAS,3)
        DOUBLE PRECISION, EXTERNAL :: basfunkval,rho
        INTEGER :: I,J,n,m

        ! Calculating the exchange correlation potential.
        ! In the case of PBE or B3LYP, this corresponds to the 
        ! tilde part of the potentail described at page 158, just 
        ! after Eqn (8.14) in Martin's "Electronic Structure"
        
        !CALL vxc(CORRLEVEL,BAS,Pup,Pdown,r,V)
        densu = rho(BAS,Pup,r)
        densd = rho(BAS,Pdown,r)
        CALL vxc2(CORRLEVEL,BAS,Pup,Pdown,densu,densd,r,V)
        CALL gradrho2(BAS,Pup+Pdown,r,gdens,ALLFI,gALLFI)
        ! Calculating the derivatives of the Exchange-Correlation energies 
        ! with respect to the length of the gradient of the respective
        ! densities. Only needed for PBE and B3LYP

        IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
               
                !densu = rho(BAS,Pup,r)
                !densd = rho(BAS,Pdown,r)
                !CALL gradrho(BAS,Pup+Pdown,r,gdens)
                !CALL gradrho2(BAS,Pup+Pdown,r,gdens,ALLFI,gALLFI)
                lengthgd = sqrt(DOT_PRODUCT(gdens,gdens))
                
                IF ( CORRLEVEL .EQ. 'PBE' ) THEN
                        gg = gVPBE(densu,densd,lengthgd)
                ENDIF

                IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL gradrho(BAS,Pup,r,gdensu)
                        CALL gradrho(BAS,Pdown,r,gdensd)
                        
                        lengthgdu = sqrt(DOT_PRODUCT(gdensu,gdensu))
                        lengthgdd = sqrt(DOT_PRODUCT(gdensd,gdensd))
                        
                        CALL gVB3LYP(densu,densd,gdensu,gdensd,ggb3lyp)
                ENDIF

        ENDIF
        
        ! Calculating the the spatial dependent matrix elements
        ! Phi_{i}(r)*Vxc(r)*Phi_{j}(r) and the corresponding nuclear gradients
        ! See page 123 in the black note book.

        Vxcr = 0.0d0

        DO I=1,BAS%NBAS
                DO J=I,BAS%NBAS
                        !-------------------------------------------------
                        ! The matrix elements Phi_{i}(r)*Vxc(r)*Phi_{j}(r)
                        !-------------------------------------------------
                        !FI = basfunkval(BAS%PSI(I),r)
                        !FJ = basfunkval(BAS%PSI(J),r)
                        FI = ALLFI(I)
                        FJ = ALLFI(J)
                        Vxcr(:,I,J) = FI*V*FJ
                        ! Adding the operator part concerning the gradient
                        ! correction. See for instance Eqn (8.14), p. 158 in R. Martin's
                        ! "Electronic structure" (Basic theory and practical Methods)
                        IF ( CORRLEVEL .EQ. 'PBE' ) THEN
                                !CALL gradbasfunkval(BAS%PSI(I),r,gFI)
                                !CALL gradbasfunkval(BAS%PSI(J),r,gFJ)
                                gFI = gALLFI(I,:)
                                gFJ = gALLFI(J,:)
                                ! Adding the operator part of the gradient potential (Martin Eqn. (8.14) ):
                                Vxcr(:,I,J) = Vxcr(:,I,J) + (gg/lengthgd)*( FI*DOT_PRODUCT(gdens,gFJ) + DOT_PRODUCT(gFI,gdens)*FJ )
                        ENDIF

                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                !CALL gradbasfunkval(BAS%PSI(I),r,gFI)
                                !CALL gradbasfunkval(BAS%PSI(J),r,gFJ)
                                gFI = gALLFI(I,:)
                                gFJ = gALLFI(J,:)
                                
                                ! Adding the operator parts of the gradient potential (Martin Eqn. (8.14) ):
                                Vxcr(1,I,J) = Vxcr(1,I,J) + (ggb3lyp(1)/lengthgdu )*( FI*DOT_PRODUCT(gdensu,gFJ) + DOT_PRODUCT(gFI,gdensu)*FJ )
                                Vxcr(2,I,J) = Vxcr(2,I,J) + (ggb3lyp(2)/lengthgdd )*( FI*DOT_PRODUCT(gdensd,gFJ) + DOT_PRODUCT(gFI,gdensd)*FJ )
                                Vxcr(:,I,J) = Vxcr(:,I,J) + (ggb3lyp(3)/lengthgd  )*( FI*DOT_PRODUCT(gdens,gFJ)  + DOT_PRODUCT(gFI,gdens)*FJ  )
                        ENDIF

                        IF ( I .NE. J ) Vxcr(:,J,I) = Vxcr(:,I,J)
                ENDDO
        ENDDO
END SUBROUTINE getvxcr

