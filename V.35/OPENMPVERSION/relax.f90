SUBROUTINE relax(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,&
           & NLSPOINTS,PORDER,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY, &
           & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,EDIR,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX,LIMPRECALC)
        ! This routine relaxes the nuclear positions using the conjugate
        ! gradient method as described in the book "Scientific computing",
        ! second edition by M. T. Heath p. 283, using the flawor of Polak and
        ! Ribiere.
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        INTEGER, INTENT(IN) :: NB,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART,NLSPOINTS,PORDER,PULAY,EDIR,NBAUX,LIMPRECALC
        INTEGER*8, INTENT(IN) :: NRED
        LOGICAL, INTENT(IN) :: APPROXEE,ADEF,RIAPPROX
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL,EPROFILE
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS),ATOMSAUX(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS,BASAUX
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25),Ftol,LQ(LORDER,3),CGQ(CGORDER,2),EETOL,ETEMP,EFIELDMAX
        DOUBLE PRECISION, INTENT(INOUT) :: DR,VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX),gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,NB,NB,NBAUX)
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD,mu,ENTROPY,OMEGA,TIME
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),NucE
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3),DPTENSORT(3,NB,NB),DPTENSOR(NB,NB)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,DRMOLD
        DOUBLE PRECISION :: Pup(Nb,NB),Pdown(NB,NB),P(NB,NB),NucFieldE
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST,RIAPPROXX,SLUTA
        INTEGER :: I,J,II,JJ,KK,NPOINTS,JSTART,JEND,POLYORDER,NSCF
        DOUBLE PRECISION, ALLOCATABLE :: ENERGY(:),DRF(:)
        REAL :: RANDOM

        IF ( RIAPPROX .AND. NB .LE. LIMPRECALC ) THEN
                RIAPPROXX = .FALSE.
        ELSE
                RIAPPROXX = RIAPPROX
        ENDIF

        OMEGA = 0.0d0
        TIME = 0.0d0

        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        NPOINTS = NLSPOINTS
        POLYORDER = PORDER
        ST = .TRUE.
        E0 = 0.0d0

        ALLOCATE(ENERGY(NPOINTS),DRF(NPOINTS))

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS )
                
                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. 1.d0 .AND. I .GT. 1 )  ST = .FALSE.

         !=======================================================================================================================
          
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R 
                ENDDO
       
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,&
                        & DIISORD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE.,NBAUX,VRI,WRI,RIAPPROXX )
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                        Pup = P/2.0d0
                        Pdown = P/2.0d0
                ENDIF
      
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,&
                        & MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY,NBAUX,VRI,WRI,RIAPPROXX)
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY,NBAUX,VRI,WRI,RIAPPROXX)
                ENDIF

                ! Calculating forces on atoms:
                CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,&
                & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELDMAX,OMEGA,TIME,ADEF,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROXX)
       
                leng = 0.0d0
                DO II=1,NATOMS
                        leng = leng + sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
      !=======================================================================================================================
                EOLD = E0
                E0 = ETOT
                DE = ETOT-EOLD

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
                
                ! Calculating the norm of the gradient:
                normgrad = 0.0d0
                DO J=1,NATOMS
                        normgrad = normgrad + DOT_PRODUCT(grad(J,:),grad(J,:))
                ENDDO
                normgrad = (1.0d0/(1.0d0*NATOMS))*sqrt(normgrad)
                !normgrad = 1.0d0

                forceold = force
                gradold  = grad
                
       !=========================================================================================================================================
       ! Here we search for the energy minimum along the gradient grad.
       !=========================================================================================================================================
                
                LSEARCH = .TRUE.
                SLUTA = .FALSE.
                ENERGY(1) = E0
                DRF(1) = 0.0d0

                DO WHILE ( LSEARCH .AND. leng .GT. FTOL )
                        JSTART = 2
                        JEND   = NPOINTS 
                        
                        IF ( SLUTA ) LSEARCH = .FALSE.
                        
                        IF ( .not. LSEARCH ) THEN
                                JSTART = 1
                                JEND   = 1
                                DRM = DABS(DRM) ! After all the step should be in the direction of the gradient
                                IF (  I .GT. 1  ) DR = 2.0d0*DRM
                                IF (  DABS(DR) .LT. 1.0E-4 ) DR = 1.0E-4
                        ENDIF
                        
                        DO J=JSTART,JEND
                                ! Calculating new positions on atoms:
                                DO JJ=1,NATOMS
                                        g = grad(JJ,:)/normgrad
                                        IF ( LSEARCH ) THEN
                                                DRF(J) = (J-1)*DR/(NPOINTS-1)
                                                ATOMS(JJ)%R = R0(JJ,:) + g*(J-1)*DRF(J)
                                                IF ( RIAPPROX ) ATOMSAUX(JJ)%R = ATOMS(JJ)%R
                                        ELSE
                                                ATOMS(JJ)%R = R0(JJ,:) + g*DRM
                                                IF ( RIAPPROX ) ATOMSAUX(JJ)%R = ATOMS(JJ)%R
                                        ENDIF
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
                                
                                ! Adding the potential energy emanating from external E-field and the atomic nuclea

                                IF ( ADEF ) THEN
                                        NucFieldE = 0.0d0
                                        DO KK=1,NATOMS
                                                IF ( EDIR .EQ. 1 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(2)
                                                ELSE IF ( EDIR .EQ. 2 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(3)
                                                ELSE IF ( EDIR .EQ. 3 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(1)
                                                ENDIF
                                        ENDDO
                                        NucE = NucE + NucFieldE
                                ENDIF
                        
                                ! Updating the basis set:
                                CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
                                ! Here the normalization of the basis functions is performed
                                CALL normalize(BAS)
                                IF ( RIAPPROX) THEN
                                        ! Updating the aux-basis set:
                                        CALL gettotalbasis(NATOMS,ATOMSAUX,NBAUX,BASAUX,.FALSE.)
                                        ! Here the normalization of the aux-basis functions is
                                        ! performed
                                        CALL normalize(BASAUX)
                                ENDIF
                                ! Calculating the matrices in the basis of the updated positions:
                                CALL overlap(NATOMS,BAS,S,gradS)
                                CALL kinetic(NATOMS,BAS,T,gradT)
                                CALL potential(BAS,NATOMS,ATOMS,V,gradV)
                
                                H0 = T + V

                                IF ( ADEF ) THEN
                                        CALL dipoletensor(BAS,DPTENSORT)
                                        IF ( EDIR .EQ. 1 ) DPTENSOR = DPTENSORT(2,:,:)
                                        IF ( EDIR .EQ. 2 ) DPTENSOR = DPTENSORT(3,:,:)
                                        IF ( EDIR .EQ. 3 ) DPTENSOR = DPTENSORT(1,:,:)
                                        H0 = H0 + DPTENSOR*EFIELDMAX
                                ENDIF

                                ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                IF ( LSEARCH ) THEN
                                        IF ( .not. RIAPPROX ) THEN
                                                CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,EETOL,.FALSE.,Pup+Pdown)
                                        ELSE
                                                CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.FALSE.)
                                                CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.FALSE.)
                                                IF ( NB .LE. LIMPRECALC ) CALL eeintsRI(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,&
                                                & PRYSR,PRYSW,APPROXEE,EETOL,.FALSE.,Pup+Pdown,NBAUX,VRI,WRI,gradVRI,gradWRI)

                                        ENDIF

                                ELSE
                                        IF ( .not. RIAPPROX ) THEN
                                                CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,EETOL,.TRUE.,Pup+Pdown)
                                        ELSE
                                                CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.TRUE.)
                                                CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.TRUE.)
                                                IF ( NB .LE. LIMPRECALC ) CALL eeintsRI(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,&
                                                & PRYSR,PRYSW,APPROXEE,EETOL,.TRUE.,Pup+Pdown,NBAUX,VRI,WRI,gradVRI,gradWRI)
                                        ENDIF
                                ENDIF
                                
                                ! Calculating the energies E0,E1,E2 in order to be used to obtain approximations of E''(R) and E'(R)
                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,&
                                        & MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE.,NBAUX,VRI,WRI,RIAPPROXX )
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                        Pup = P/2.0d0
                                        Pdown = P/2.0d0
                                ELSE IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                
                                CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,&
                                & MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY,NBAUX,VRI,WRI,RIAPPROXX)
                                
                                ELSE IF ( ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) .AND. LSEARCH ) THEN
                                
                                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY,NBAUX,VRI,WRI,RIAPPROXX)
                                
                                ENDIF

                                IF ( LSEARCH ) ENERGY(J) = ETOT

                                ! Calculating the new position in the line-searc for the minimum
                                ! by fitting the NPOINTS energies to a POLYORDER:th  order polynomial, 
                                ! and by minimizing this polynomial in the intervall [0,DR]
                                IF ( J .EQ. NPOINTS .AND. LSEARCH ) THEN
                                    DRMOLD = DRM
                                    CALL linesearchmin(POLYORDER,ENERGY,DRF,NPOINTS,DR,DRM)
                                    ! Safety cutoff:
                                    IF ( DRM .LT. 1.0E-10 ) DRM = DRMOLD
                                    IF ( DABS(DRM) .GT. DR ) DRM = DR
                                    SLUTA = .TRUE.
                                ENDIF
                        ENDDO
                        !-------------------------------------------------
                        ! Saving the computed new atomic positions to file
                        !-------------------------------------------------
                        OPEN(16,FILE='ATOMPOSITIONS.dat',ACTION='WRITE')
                        DO JJ=1,NATOMS
                             WRITE(16,'(A4,I4,3(F15.10))')'ATOM',ATOMS(JJ)%Z,ATOMS(JJ)%R 
                        ENDDO
                        CLOSE(16)
                        !---------------------------------------------------------------
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
                
                !leng = 0.0d0
                !DO II=1,NATOMS
                !        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                !ENDDO
                
                !leng = leng/NATOMS
                
                WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)')I,E0,DE,leng,DRM

                I = I + 1
                IF ( ABS(DR) .LT. 0.01 ) SCRATCH = .FALSE.
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

      !STOP

      END SUBROUTINE relax
