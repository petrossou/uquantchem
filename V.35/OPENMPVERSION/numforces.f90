SUBROUTINE numforce(S,H0,Intsv,NB,NRED,Ne,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,DRF,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
           & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,NEPERIOD,EDIR,EPROFILE,EFIELDMAX,OMEGA,TIME,ADEF,Pup,Pdown,PRINTOUT,NBAUX,ATOMSAUX,BASAUX,VRI,WRI,RIAPPROX,force)
        ! This routine calculatesd the interatomic forces by finite diffrences
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        INTEGER, INTENT(IN) :: NB,Ne,NATOMS,DIISORD,DIISSTART,NEPERIOD,EDIR,NBAUX
        INTEGER*8, INTENT(IN) :: NRED
        LOGICAL, INTENT(IN) :: APPROXEE,ADEF,PRINTOUT,RIAPPROX
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL,EPROFILE
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), Intsv(NRED),force(NATOMS,3)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS),ATOMSAUX(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS,BASAUX
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25),Ftol,LQ(LORDER,3),CGQ(CGORDER,2),EETOL,ETEMP,EFIELDMAX,OMEGA,TIME
        DOUBLE PRECISION, INTENT(IN) :: DRF
        DOUBLE PRECISION, INTENT(INOUT) :: VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
        DOUBLE PRECISION :: gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,NB,NB,NBAUX)
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD,mu,ENTROPY
        DOUBLE PRECISION :: NucE,DPTENSOR(NB,NB),DPTENSORT(3,NB,NB),EFIELD
        DOUBLE PRECISION :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(1,3,1)
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,DRMOLD
        DOUBLE PRECISION, INTENT(INOUT)  :: Pup(Nb,NB),Pdown(NB,NB)
        DOUBLE PRECISION :: Pup0(Nb,NB),Pdown0(NB,NB)
        DOUBLE PRECISION :: ENERGY(NATOMS,3,2)
        LOGICAL :: SCRATCH,LSEARCH
        INTEGER :: I,J,II,JJ,KK,NSCF,DIISORDD
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

        SCRATCH = .FALSE.
        DIISORDD = 0

         !=======================================================================================================================
          
         DO J=1,NATOMS
              R0(J,:) = ATOMS(J)%R 
         ENDDO
        force = 0.0d0

        ! Calculating the forces on the atoms by numerical differentiation:
        DO I=1,NATOMS
                DO J=1,3
                        DO KK=0,1
                                DO JJ=1,NATOMS
                                        ATOMS(JJ)%R = R0(JJ,:)
                                        IF ( RIAPPROX ) ATOMSAUX(JJ)%R = R0(JJ,:)
                                ENDDO
                                Pup0 = Pup
                                Pdown0 = Pdown
                                ATOMS(I)%R(J) = R0(I,J) + (2*KK - 1)*DRF
                                IF ( RIAPPROX ) ATOMSAUX(I)%R(J) = R0(I,J) + (2*KK - 1)*DRF
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
                                ! In the case we have an electric field present:
                                IF ( ADEF ) THEN
                                        IF ( EPROFILE .EQ. 'AC' ) THEN
                                                CALL dipoletensorinhom(BAS,EDIR,OMEGA,TIME,DPTENSOR)
                                        ELSE
                                                CALL dipoletensor(BAS,DPTENSORT)
                                                IF ( EDIR .EQ. 1 ) DPTENSOR = DPTENSORT(2,:,:)
                                                IF ( EDIR .EQ. 2 ) DPTENSOR = DPTENSORT(3,:,:)
                                                IF ( EDIR .EQ. 3 ) DPTENSOR = DPTENSORT(1,:,:)
                                        ENDIF

                                        IF ( EPROFILE .EQ. 'AC' .OR. EPROFILE .EQ. 'HO' ) THEN
                                                IF ( TIME .GE. 0.0d0 .AND. TIME .LE.  2.0d0*pi/OMEGA ) THEN
                                                        EFIELD = EFIELDMAX*(OMEGA*TIME/(2.0d0*pi))
                                                ENDIF

                                                IF ( TIME .GT. 2.0d0*pi/OMEGA .AND. TIME .LE.  (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                                        EFIELD = EFIELDMAX
                                                ENDIF

                                                IF ( TIME .GT. (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) .AND.  TIME .LT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                                        EFIELD = EFIELDMAX*( (NEPERIOD + 2.0d0) - OMEGA*TIME/(2.0d0*pi) )
                                                ENDIF

                                                IF ( TIME .GT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) EFIELD = 0.0d0

                                                IF ( EPROFILE .EQ. 'AC' ) THEN
                                                        DPTENSOR = DPTENSOR*EFIELD
                                                ELSE IF ( EPROFILE .EQ. 'HO' ) THEN
                                                        EFIELD = EFIELD*sin(OMEGA*TIME)
                                                        DPTENSOR = DPTENSOR*EFIELD
                                                ENDIF
                                        ELSE IF ( EPROFILE .EQ. 'ST' ) THEN
                                                EFIELD = EFIELDMAX
                                                DPTENSOR = DPTENSOR*EFIELD
                                        ENDIF

                                        H0 = H0 - DPTENSOR
                                ENDIF      

                                ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                IF ( .not. RIAPPROX ) THEN
                                        CALL eeints(NATOMS,1,BAS,Intsv,gradIntsv,NRED,1,PRYSR,PRYSW,APPROXEE,EETOL,.FALSE.,Pup0+Pdown0)
                                ELSE
                                        CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.FALSE.)
                                        CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.FALSE.)
                                ENDIF
                                ! Calculating the energies E0,E1,E2 in order to be used to obtain approximations of E''(R) and E'(R)
                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,Pup0+Pdown0,MIX,&
                                        & DIISORDD,DIISSTART,NSCF,-1,.FALSE., SCRATCH,.FALSE.,NBAUX,VRI,WRI,RIAPPROX )
                                ELSE IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                                        CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,&
                                        & ETOT,Cup,Cdown,Pup0,Pdown0,MIX,DIISORDD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY,NBAUX,VRI,WRI,RIAPPROX)
                                ELSE IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                        & ETOT,Cup,Cdown,Pup0,Pdown0,MIX,DIISORDD,DIISSTART,NSCF,-1,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY,NBAUX,VRI,WRI,RIAPPROX)
                                ENDIF

                                !----------------------------------------------------------------------
                                ! Storing the energy to be used in the finite diffrence ! calculation
                                !----------------------------------------------------------------------
                                ENERGY(I,J,KK+1) = ETOT
                        ENDDO
                        !----------------------------------------------------------------
                        ! Calculating the forces by finite differences:
                        !----------------------------------------------------------------
                        force(I,J) = -(1.0/(2.0d0*DRF))*( ENERGY(I,J,2) - ENERGY(I,J,1) )
                        IF (PRINTOUT ) THEN 
                                IF ( I .EQ. 1 .AND. J .EQ. 1 ) THEN
                                        WRITE(*,'(A91)')'   ======================================================================================='
                                        WRITE(*,'(A91)')'                       Calculating forces by finite differences                           '
                                        WRITE(*,'(A91)')'   ======================================================================================='
                                        WRITE(*,'(A91)')' '
                                        WRITE(*,'(A91)')'   ATOM COORD        E(R-DRF) [au]                 E(R+DRF) [au]                    F [au]'
                                ENDIF
                                WRITE(*,'(I6,I6,E30.20,E30.20,E30.20,E30.20)')I,J,ENERGY(I,J,1),ENERGY(I,J,2),force(I,J)
                        ENDIF
                ENDDO
        ENDDO
        IF (PRINTOUT ) WRITE(*,'(A91)')' '
        DO I=1,NATOMS
                !----------------------------------------------------------------
                ! Resetting the atomic positions to their original configuration
                !----------------------------------------------------------------
                ATOMS(I)%R = R0(I,:)
        ENDDO
      END SUBROUTINE numforce
