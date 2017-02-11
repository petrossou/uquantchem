SUBROUTINE relax(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW)
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART
        LOGICAL, INTENT(IN) :: APPROXEE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: DR,MIX,Tol,PRYSR(25,25),PRYSW(25,25)
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),NucE
        DOUBLE PRECISION :: leng, f(3), Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB)
        LOGICAL :: CFORCE,SCRATCH
        INTEGER :: I,J,II,JJ

        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        EOLD = 0.0d0
        SCRATCH = .TRUE.

        DO WHILE ( DE .GT. Tol .AND. I .LT. NSTEPS ) 
                
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,.FALSE.,SCRATCH)
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,.FALSE.,SCRATCH)
                ENDIF
                
                ! Calculating forces on atoms:
                CALL forces(NATOMS,Ne,NB,NRED,Cup,Cdown,EHFeigenup,EHFeigendown,ATOMS,gradS,gradT,gradV,gradIntsv,force)

                ! Calculating new positions on atoms:
                DO J=1,NATOMS
                        f = force(J,:)
                        leng = sqrt(DOT_PRODUCT(f,f))
                        IF ( leng .GT. 1.0d0 ) f = f/leng
                        
                        ! Updating the atomic positions:
                        ATOMS(J)%R = ATOMS(J)%R + f*DR
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
                CALL eeints(NATOMS,NATOMS,BAS,Intsv,gradIntsv,NRED,NRED,PRYSR,PRYSW,APPROXEE,Tol,CFORCE)

                IF ( I .EQ. 0 ) THEN 
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,'(A79)')'                             Relaxing the nuclear positions:                   '
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,*)' '
                        WRITE(*,'(A85)')'  ITER           E [au]                        DE [au]                     <|F|> [au]'
                ENDIF
                
                DE = DABS(ETOT-EOLD)
                EOLD = ETOT
                leng = 0.0d0
                DO II=1,NATOMS
                        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
                
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,ETOT,DE,leng

                I = I + 1
                SCRATCH = .FALSE.

      ENDDO
        
      IF ( DE .LT. Tol ) THEN
              WRITE(*,*)' '
              WRITE(*,'(A69)')'                 The nuclear positions have been relaxed with respect'
              WRITE(*,'(A68)')'                 to the interatomic forces. An energy convergence of'
              WRITE(*,'(A27,E12.5,A24)')' DE < ',Tol,' [au] has been obtained.'
              WRITE(*,*)' ' 
              WRITE(*,'(A70)')'                The relaxed positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ELSE
              WRITE(*,*)' '
              WRITE(*,'(A84)')'          Failure to relax the nuclear positions within prescribed energy tolerance'
              WRITE(*,*)' ' 
              WRITE(*,'(A71)')'                     The last positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ENDIF

      STOP

      END SUBROUTINE relax
