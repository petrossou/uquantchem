SUBROUTINE relaxf(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY, &
           & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL)
        ! This routine relaxes the nuclear positions using the conjugate
        ! gradient method as described in the book "Scientific computing",
        ! second edition by M. T. Heath p. 283, using the flawor of Polak
        ! and Ribiere, by searchig for the atomic configuration where all the
        ! forces are zero, by means of a Newton-Raphson search on the forces.      
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART,PULAY
        LOGICAL, INTENT(IN) :: APPROXEE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25),Ftol,LQ(LORDER,3),CGQ(CGORDER,2),EETOL
        DOUBLE PRECISION, INTENT(INOUT) :: DR
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),NucE
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,DRMOLD
        DOUBLE PRECISION :: Pup(Nb,NB),Pdown(NB,NB),P(NB,NB),p1,p2,hstep,dF(NATOMS,3),dF2(NATOMS,3),forcef(2,NATOMS,3)
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST,SLUT
        INTEGER :: I,J,II,JJ,KK,NPOINTS,JSTART,JEND,POLYORDER,NSCF
        DOUBLE PRECISION :: ENERGY(2),DRF(2)
        REAL :: RANDOM
       
        hstep = 0.00001 
        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        ST = .TRUE.
        E0 = 0.0d0

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS )
                
                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. 10.0d0*FTOL )  ST = .FALSE.

         !=======================================================================================================================
          
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R 
                ENDDO
       
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE. )
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                        Pup = P/2.0d0
                        Pdown = P/2.0d0
                ENDIF
      
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.)
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.)
                ENDIF

                ! Calculating forces on atoms:
                CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,&
                & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force)
       
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
                ENERGY(1) = E0
                SLUT = .FALSE.
                J = 0

                DO WHILE ( LSEARCH )
                        
                        JSTART = 1
                        JEND   = 2
                        
                        IF ( J .EQ. 2  ) THEN
                                JSTART = 1
                                JEND   = 1
                        ENDIF
                       
                        IF ( SLUT ) LSEARCH = .FALSE.

                        DO J=JSTART,JEND
                                ! Calculating new positions on atoms:
                                DO JJ=1,NATOMS
                                        g = grad(JJ,:)/normgrad
                                        IF ( LSEARCH ) THEN
                                                IF ( J .EQ. 1 ) ATOMS(JJ)%R = R0(JJ,:) - g*hstep
                                                IF ( J .EQ. 2 ) ATOMS(JJ)%R = R0(JJ,:) + g*hstep
                                        ELSE
                                                ATOMS(JJ)%R = R0(JJ,:) + g*DRM
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
                                CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,EETOL,.TRUE.,Pup+Pdown)
                                
                                ! Calculating the energies E0,E1,E2 in order to be used to obtain approximations of E''(R) and E'(R)
                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE. )
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                        Pup = P/2.0d0
                                        Pdown = P/2.0d0
                                ELSE IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                
                                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.)
                                
                                ELSE IF ( ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) .AND. LSEARCH ) THEN
                                
                                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.)
                                
                                ENDIF
                                
                                ! Calculating forces on atoms:
                                IF ( LSEARCH ) THEN
                                        CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,&
                                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force)
                                        ENERGY(J) = ETOT
                                        forcef(J,:,:) = force
                                ENDIF

                                ! Calculating the new position in the line-searc for the minimum
                                ! using Newton-Raphson along the gradient, g, 
                                IF ( J .EQ. 2 .AND. LSEARCH ) THEN
                                        ! calculating the derivative and the second derivative
                                        ! of the force along the gradient g.
                                        dF = (0.50d0/hstep)*( forcef(2,:,:) - forcef(1,:,:))
                                        dF2 = (1.0d0/(hstep**2))*( forcef(2,:,:) + forcef(1,:,:) -2.0d0*forceold )
                                        ! Calculating the projections along g of the above derivatives
                                        p1 = 0.0d0
                                        p2 = 0.0d0
                                        DRM = 2.0d0*DR
                                        DO JJ=1,NATOMS
                                                g = grad(JJ,:)/normgrad
                                                p1 = p1 + DOT_PRODUCT(dF(JJ,:),g)
                                                p2 = p2 + DOT_PRODUCT(dF2(JJ,:),g)
                                        ENDDO
                                        IF ( ( p1 .LT. 0.0d0 .AND. p2 .LT. 0.0d0) .OR. ( p1 .GT. 0.0d0 .AND. p2 .LT.  0.0d0 ) ) THEN
                                                IF ( p1 .LT. 0.0d0 ) DRM = - p1/p2 + sqrt( (p1/p2)**2 - 2*normgrad/p2 )
                                                IF ( p1 .GT. 0.0d0 ) DRM = - p1/p2 - sqrt( (p1/p2)**2 - 2*normgrad/p2 )
                                        ELSE
                                                IF ( p1 .NE. 0.0d0 ) THEN
                                                        DRM = -normgrad/p1
                                                ENDIF
                                                IF ( p2 .LT. 0.0d0 .AND. p1 .EQ.  0.0d0 ) THEN
                                                        DRM = sqrt(-2.0d0*normgrad/p2)
                                                ENDIF
                                                IF ( DABS(DRM) .GT. DR ) DRM = DR*DRM/DABS(DRM)
                                        ENDIF
                                        SLUT = .TRUE.
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
                
                leng = 0.0d0
                DO II=1,NATOMS
                        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
                
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

      STOP

      END SUBROUTINE relaxf
