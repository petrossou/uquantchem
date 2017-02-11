SUBROUTINE relax(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW)
        ! This routine relaxes the nuclear positions using the conjugate
        ! gradient method as described in the book "Scientific computing",
        ! second edition by M. T. Heath p. 283, using the flawor of Polak and
        ! Ribiere.
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART
        LOGICAL, INTENT(IN) :: APPROXEE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: DR,MIX,Tol,PRYSR(25,25),PRYSW(25,25),Ftol
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),NucE
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,DRMO,Ebis,Eprim,nom
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST
        INTEGER :: I,J,II,JJ,KK,NLITER,JSTART,JEND
        REAL :: RANDOM
        
        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        NLITER = 10
        !ST = .FALSE.

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS )
                
                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. 1.0E-4 ) ST = .FALSE.

         !=======================================================================================================================
          
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R 
                ENDDO
       
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,.FALSE., SCRATCH )
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                ENDIF
      
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,.FALSE.,SCRATCH)
                ENDIF
     
                ! Calculating forces on atoms:
                CALL forces(NATOMS,Ne,NB,NRED,Cup,Cdown,EHFeigenup,EHFeigendown,ATOMS,gradS,gradT,gradV,gradIntsv,force)
       
      !=======================================================================================================================
                
                E0 = ETOT
                
                IF ( I .EQ. 0 ) THEN
                        grad = -force
                ELSE
                        const1 = 0.0d0
                        const2 = 0.0d0
                        DO J=1,NATOMS
                                const1 = const1 +  DOT_PRODUCT(force(J,:)-forceold(J,:),force(J,:))
                                const2 = const2 +  DOT_PRODUCT(forceold(J,:),forceold(J,:))
                        ENDDO
                        ! Uppdating the gradient according to Polack and Ribiere
                        DO J=1,NATOMS
                                IF ( ST ) THEN
                                        ! steepest deecent
                                        grad(J,:) = force(J,:) 
                                ELSE
                                        ! Conjugate gradient
                                        grad(J,:) = force(J,:) + (const1/const2)*gradold(J,:)
                                ENDIF

                        ENDDO
                ENDIF

                forceold = force
                gradold  = grad
                
       !=========================================================================================================================================
       ! Here we search for the energy minimum along the gradient grad.
       !=========================================================================================================================================
                DRMO    = 10.0d0
                Eprim   = 10.0d0
                DRM     = 0.0d0
                KK      = 0
                LSEARCH = .TRUE.

                DO WHILE ( LSEARCH )
                        IF ( DABS(Eprim) .LT. DR**2 .OR. KK .EQ. NLITER ) LSEARCH = .FALSE.
                        
                        JSTART = 1
                        JEND   = 5
                        
                        IF ( KK .EQ. 0 ) JSTART = 2
                        
                        IF ( .not. LSEARCH ) THEN
                                JSTART = 1
                                JEND   = 1
                        ENDIF
                        
                        DO J=JSTART,JEND
                                ! Calculating new positions on atoms:
                                DO JJ=1,NATOMS
                                        g = grad(JJ,:)
                                        IF ( J .EQ. 1 ) ATOMS(JJ)%R = R0(JJ,:) + g*DRM
                                        IF ( J .EQ. 2 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM-DR)
                                        IF ( J .EQ. 3 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM+DR)
                                        IF ( J .EQ. 4 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM-2*DR)
                                        IF ( J .EQ. 5 ) ATOMS(JJ)%R = R0(JJ,:) + g*(DRM+2*DR)

                                ENDDO
                
                                ! Calculating the nuclear-nuclear repulsion energy of the
                                ! updated nuclear configuration
                                NucE = 0.0d0
                                DO II=1,NATOMS
                                        DO JJ=II+1,NATOMS
                                                Rn = ATOMS(II)%R - ATOMS(JJ)%R
                                                NucE = NucE + ATOMS(II)%Z*ATOMS(JJ)%Z/sqrt(DOT_PRODUCT(Rn,Rn))
                                        ENDDO
                                ENDDO

                                ! Updating the basis set:
                                CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
                                ! Here the normalization of the basis functions is performed
                                CALL normalize(BAS)
                                ! Calculating the matrices in the basis of the updated positions:
                                CALL overlap(NATOMS,BAS,S,gradS)
                                CALL kinetic(NATOMS,BAS,T,gradT)
                                CALL potential(BAS,NATOMS,ATOMS,V,gradV)
                
                                H0 = T + V

                                ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                IF ( LSEARCH ) THEN
                                        CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,Tol/100,.FALSE.)
                                ELSE
                                        CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,Tol/100,.TRUE.)
                                ENDIF
                                
                                ! Calculating the energies E0,E1,E2 in order to be used to obtain approximations of E''(R) and E'(R)
                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,.FALSE., SCRATCH )
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                ELSE IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,.FALSE.,SCRATCH)
                                ENDIF
                      
                                IF ( J .EQ. 1 ) E0 = ETOT
                                IF ( J .EQ. 2 ) E1 = ETOT
                                IF ( J .EQ. 3 ) E2 = ETOT
                                IF ( J .EQ. 4 ) E3 = ETOT
                                IF ( J .EQ. 5 ) E4 = ETOT
                                
                                ! Calculating the new position in the line-searc for the minimum:
                                IF ( J .EQ. 5 ) THEN 
                                    ! Here we use the 5-point finite differnce
                                    ! formulas presented in Koonin's book
                                    ! Computaional Physics. page 6, table 1.2
                                    Eprim = (1.0d0/12.0d0)*(E3 - 8*E1 + 8*E2 - E4 )/DR
                                    Ebis  = (1.0d0/12.0d0)*(-E3 + 16*E1 - 30*E0 + 16*E2 - E4 )/(DR**2)
                                    IF ( Eprim .NE. 0.0d0 ) THEN
                                            DRMO = DRM
                                            IF ( Ebis .NE. 0.0d0 ) DRM = DRM - Eprim/dabs(Ebis)
                                    ELSE
                                            DRMO = DRM
                                    ENDIF
                                    ! Safety cutoff:
                                    IF ( DABS(DRM) .GT. 1.0d0 ) DRM = 1.0d0
                                ENDIF
                        
                        ENDDO
                        KK = KK + 1
                ENDDO

       !=======================================================================================================================
       ! End of line search 
       !=======================================================================================================================

                IF ( I .EQ. 0 ) THEN 
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,'(A79)')'                             Relaxing the nuclear positions:                   '
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,*)' '
                        WRITE(*,'(A114)')'  ITER           E [au]                        DE [au]                     <|F|> [au]                      DR [au]'
                ENDIF
                
                DE = (ETOT-EOLD)
                EOLD = E0
                leng = 0.0d0
                DO II=1,NATOMS
                        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
                
                WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,I4,E30.20)')I,E0,DE,leng,DRM,KK,DABS(DRM-DRMO)

                I = I + 1
                SCRATCH = .FALSE.
      ENDDO
        
      IF ( leng .LT. FTol ) THEN
              WRITE(*,*)' '
              WRITE(*,*)'                 The nuclear positions have been relaxed with respect'
              WRITE(*,*)'                 to the interatomic forces. An force convergence of '
              WRITE(*,'(A29,E10.5,A24)')' F  < ',FTol,' [au] has been obtained.'
              WRITE(*,*)' ' 
              WRITE(*,*)'                The relaxed positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ELSE
              WRITE(*,*)' '
              WRITE(*,*)'          Failure to relax the nuclear positions within prescribed force tolerance '
              WRITE(*,*)' ' 
              WRITE(*,*)'                     The last positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ENDIF

      STOP

      END SUBROUTINE relax
