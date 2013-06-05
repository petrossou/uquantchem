SUBROUTINE relax(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
& IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL)
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Qstart,Qend,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend)
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART,Istart,Iend,numprocessors,id,NLSPOINTS,PORDER,PULAY
        INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
        LOGICAL, INTENT(IN) :: APPROXEE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), IntsvR(Istart:Iend)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25),FTol,LQ(LORDER,3),CGQ(CGORDER,2),EETOL
        DOUBLE PRECISION, INTENT(INOUT)  :: DR
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION :: Pup(NB,NB),Pdown(NB,NB),P(NB,NB)
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsvR(NATOMS,3,Istart:Iend),NucE
        DOUBLE PRECISION :: leng, f(3),g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),DR0,DR1,DR2,R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3),DRMOLD
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST
        INTEGER :: I,J,II,JJ,KK,NLITER,JSTART,JEND,NPOINTS,POLYORDER,NSCF,ierr
        DOUBLE PRECISION, ALLOCATABLE :: ENERGY(:),DRF(:)

        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        NLITER = 10
        ST = .TRUE.
        NPOINTS = NLSPOINTS
        POLYORDER = PORDER

        ALLOCATE(ENERGY(NPOINTS),DRF(NPOINTS))

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS ) 

                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. 10.0d0*FTOL ) ST = .FALSE.
     
        !=======================================================================================================================
                
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R
                ENDDO
  
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                        Pup = P/2.0d0
                        Pdown = P/2.0d0
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX, &
                                  & DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                ENDIF
                IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR.  CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,&
                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                ENDIF

                ! Calculating forces on atoms:
                CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY,&
                & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL)
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
                                const1 = const1 + DOT_PRODUCT(force(J,:)-forceold(J,:),force(J,:))
                                const2 = const2 + DOT_PRODUCT(forceold(J,:),forceold(J,:))
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

                forceold = force
                gradold  = grad

        !=========================================================================================================================================
        ! Here we search for the energy minimum along the gradient grad.
        !=========================================================================================================================================
                
                
                LSEARCH = .TRUE.
                ENERGY(1) = E0
                DRF(1) = 0.0d0

                DO WHILE ( LSEARCH )
                        JSTART = 2
                        JEND   = NPOINTS

                        IF ( J .EQ. NPOINTS+1 ) LSEARCH = .FALSE.

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
                                IF ( LSEARCH ) THEN
                                        CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,.FALSE.,id,Pup+Pdown)
                                ELSE
                                        CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,.TRUE.,id,Pup+Pdown)
                                ENDIF

                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1, &
                                                & numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                        Pup = P/2.0d0
                                        Pdown = P/2.0d0
                                ENDIF
                
                                IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1, & 
                                        & numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                                ENDIF
                                
                                IF ( ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) .AND. LSEARCH ) THEN
                                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,&
                                        & Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                                ENDIF

                                IF ( LSEARCH ) ENERGY(J) = ETOT
                       
                                ! Calculating the new position in the line-searc for the minimum
                                ! by fitting the NPOINTS energies to a  POLYORDER:th  order polynomial,
                                ! and by minimizing this polynomial in the  intervall [0,DR]
                                IF ( J .EQ. NPOINTS .AND. LSEARCH ) THEN
                                        DRMOLD = DRM
                                        CALL linesearchmin(POLYORDER,ENERGY,DRF,NPOINTS,DR,DRM)
                                        ! Safety cutoff:
                                        IF ( DRM .LT. 1.0E-10 ) DRM = DRMOLD
                                        IF ( DABS(DRM) .GT. DR ) DRM = DR
                                ENDIF

                          ENDDO

                          IF ( id .EQ. 0 ) THEN
                                  !-------------------------------------------------
                                  ! Saving the computed new atomic positions to file
                                  !-------------------------------------------------
                                  OPEN(16,FILE='ATOMPOSITIONS.dat',ACTION='WRITE')
                                  DO JJ=1,NATOMS
                                        WRITE(16,'(A4,I4,3(F15.10))')'ATOM',ATOMS(JJ)%Z,ATOMS(JJ)%R
                                  ENDDO
                                  CLOSE(16)
                                  !---------------------------------------------------------------
                          ENDIF
                          !================================================================================
                          ! Here we make sure the the mean length of the force
                          ! is consistent over all threads to avoid deadlock
                          !================================================================================
                          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                          CALL MPI_BCAST(leng,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                          !================================================================================
                  ENDDO
          !=========================================================================================================================================

               IF ( I .EQ. 0 .AND. id .EQ. 0 ) THEN 
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
                
                IF ( id .EQ. 0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)')I,E0,DE,leng,DRM

                I = I + 1
                SCRATCH = .FALSE.

      ENDDO
        
      IF ( leng .LT. FTol .AND. id .EQ. 0 ) THEN
              WRITE(*,*)' '
              WRITE(*,'(A69)')'                 The nuclear positions have been relaxed with respect'
              WRITE(*,'(A68)')'                 to the interatomic forces. An force convergence of '
              WRITE(*,'(A27,E12.5,A24)')' F  < ',FTol,' [au] has been obtained.'
              WRITE(*,*)' ' 
              WRITE(*,'(A70)')'                The relaxed positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ELSE
              IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'          Failure to relax the nuclear positions within prescribed force tolerance '
                        WRITE(*,*)' ' 
                        WRITE(*,'(A71)')'                     The last positions of the nuclea are the following:'
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
                        DO J=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
                        ENDDO
                ENDIF
      ENDIF

      STOP

      END SUBROUTINE relax
