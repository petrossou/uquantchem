SUBROUTINE moleculardynamics(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,DT,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                             & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id)
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        INTEGER, INTENT(IN) :: NB,NRED,Ne,NATOMS,NTIMESTEPS,DIISORD,DIISSTART,SAMPLERATE,Istart,Iend,numprocessors,id
        INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
        LOGICAL, INTENT(IN) :: APPROXEE,WRITEONFLY,MOVIE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB), IntsvR(Istart:Iend)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25)
        DOUBLE PRECISION, INTENT(IN) :: DT,TEMPERATURE
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsvR(NATOMS,3,Istart:Iend),NucE
        DOUBLE PRECISION :: leng, f(3), g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3)
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,DRMOLD,VCM(3),MTOT
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST
        INTEGER :: I,J,II,JJ,KK,NPOINTS,JSTART,JEND,MM,clock,ierr
        DOUBLE PRECISION :: ENERGY(NTIMESTEPS),VEL(NTIMESTEPS,NATOMS,3),R(NTIMESTEPS,NATOMS,3),EKIN(NTIMESTEPS),EPOT(NTIMESTEPS),TEMP(NTIMESTEPS)
        REAL :: RAND
        DOUBLE PRECISION, PARAMETER :: KB = 0.000003166811524300030d0 !  Boltzmanns constant in au/K
        
        
        IF ( id .EQ. 0 ) OPEN(100,FILE='MOLDYNENERGY.dat',ACTION='WRITE')
         
        IF ( MOVIE .AND. id .EQ. 0 ) OPEN(200,FILE='MOLDYNMOVIE.xsf',ACTION='WRITE')
        
        I = 1
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        ST = .TRUE.
        E0 = 0.0d0
        KK = 0
        VCM = 0.0d0
        MTOT = 0.0d0

        ! Initialization of positions and velocities:
        DO J=1,NATOMS
                R(1,J,:) = ATOMS(J)%R
                CALL RANDOM_NUMBER(RAND)
                f(1) = RAND
                CALL RANDOM_NUMBER(RAND)
                f(2) = RAND
                CALL RANDOM_NUMBER(RAND)
                f(3) = RAND
                IF ( MOD(J,2) .NE. 0 ) THEN
                        VEL(1,J,:) = (f/sqrt(DOT_PRODUCT(f,f)))*sqrt(2.0d0*3.0d0*KB*TEMPERATURE/ATOMS(J)%M)
                ELSE
                        VEL(1,J,:) = -(VEL(1,J-1,:)/sqrt(DOT_PRODUCT(VEL(1,J-1,:),VEL(1,J-1,:))))*sqrt(2.0d0*3.0d0*KB*TEMPERATURE/ATOMS(J)%M)
                ENDIF
                ! Calculating center of mass velocity (to be subtracted later)
                VCM = VCM + VEL(1,J,:)*ATOMS(J)%M
                MTOT = MTOT + ATOMS(J)%M
        ENDDO
        VCM = VCM/MTOT

        ! subtracting the center of mass velocity from all the atoms:
        DO J=1,NATOMS
                VEL(1,J,:)  = VEL(1,J,:) - VCM
        ENDDO

        
        DO WHILE ( I .LT. NTIMESTEPS )
                
        
         !=======================================================================================================================
                
                IF ( I .EQ. 1 ) THEN
                        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                                CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                                EHFeigenup = EHFeigen
                                EHFeigendown = EHFeigen
                                Cup = EIGENVECT
                                Cdown = EIGENVECT
                        ENDIF
      
                        IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                               CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                        ENDIF
     
                        ! Calculating forces on atoms:
                        CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,EHFeigenup,EHFeigendown,ATOMS,gradS,gradT,gradV,gradIntsvR,force,numprocessors,id)
                 ENDIF
                 
                 ! Saveing the force of the current time-step:
                 forceold = force

                 !=======================================================================================================================
                 ! The velocity Verlet algorithm as presented in J.M Thijssen's
                 ! "Computational Physics" P.181. Here we update the positions:
                 !=======================================================================================================================
                  DO J=1,NATOMS
                          R(I+1,J,:) = R(I,J,:) + DT*VEL(I,J,:) + (DT**2)*force(J,:)/(2.0d0*ATOMS(J)%M)
                  ENDDO

                 ! Uppdating the atomic positions to the positions of the
                 ! current time-step:
                 DO J=1,NATOMS
                        ATOMS(J)%R = R(I+1,J,:)
                 ENDDO

                 !-----------------------------------------------------------
                 ! Calculating the nuclear-nuclear repulsion energy of the
                 ! updated nuclear configuration
                 !-----------------------------------------------------------
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
                 CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,.TRUE.,id)

                 !=============================================
                 ! Calculating the energies and forces for t+Dt
                 !=============================================
                 IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                         CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                         EHFeigenup = EHFeigen
                         EHFeigendown = EHFeigen
                         Cup = EIGENVECT
                         Cdown = EIGENVECT
                 ENDIF
      
                 IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,MIX,DIISORD,DIISSTART,numprocessors,id,.FALSE.,SCRATCH)
                 ENDIF

                 ! Calculating forces on atoms:
                 CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,EHFeigenup,EHFeigendown,ATOMS,gradS,gradT,gradV,gradIntsvR,force,numprocessors,id)
                 
                 
                 !=======================================================================================================================
                 ! The velocity Verlet algorithm as presented in J.M Thijssen's
                 ! "Computational Physics" P.181. Here we update the velocities:
                 !=======================================================================================================================
                  DO J=1,NATOMS
                          VEL(I+1,J,:) = VEL(I,J,:) + (DT/2.0d0)*(force(J,:) + forceold(J,:) )/ATOMS(J)%M
                  ENDDO
                 
                 !===============================================
                 ! Calculating the diferent energy contributions:
                 !===============================================

                 ! The kinetic energy:
                 EKIN(I+1) = 0.0d0
                 DO J=1,NATOMS
                        EKIN(I+1) = EKIN(I+1) + (ATOMS(J)%M/2.0d0)*DOT_PRODUCT(VEL(I+1,J,:),VEL(I+1,J,:))
                 ENDDO

                 ! The potential energy
                 EPOT(I+1) = ETOT
                 
                 ! The total energy
                 ENERGY(I+1) = EPOT(I+1) + EKIN(I+1)

                 !=================================================
                 ! Calculating the temperature:
                 !=================================================

                 TEMP(I+1) = EKIN(I+1)/(3.0d0*NATOMS*KB)
                
                 IF ( I .EQ. 1 .AND. id .EQ. 0 ) THEN
                         WRITE(100,'(A124)')'        ITER           E [au]                        EK [au]                       EP [au]                        T [K]     '
                         IF (WRITEONFLY) WRITE(*,'(A124)')'        ITER           E [au]                        EK [au]                       EP [au]                        T [K]     '
                         IF ( MOVIE ) WRITE(200,'(A9,I10)')'ANIMSTEPS',INT(NTIMESTEPS*1.0/(SAMPLERATE*1.0))
                 ENDIF

                 ! increasing the time index one step
                 I = I + 1
                 
                 !=============================================================
                 ! Enabaling the use of the density matrix of the previous step 
                 ! to be used as a starting guess for the next time step
                 !=============================================================
                 SCRATCH = .FALSE.

                 !==============================
                 ! Saving the energies to file:
                 !==============================
                 IF (WRITEONFLY .AND. id .EQ. 0 )  THEN
                         WRITE(100,'(I10,E30.20,E30.20,E30.20,E30.20)')I,ENERGY(I),EKIN(I),EPOT(I),TEMP(I)
                         WRITE(*,'(I10,E30.20,E30.20,E30.20,E30.20)')I,ENERGY(I),EKIN(I),EPOT(I),TEMP(I)
                 ENDIF
                 !==================================================
                 ! saving a xsf file used by xcrysden to make movies
                 !==================================================
                 IF ( MOVIE .AND. MOD(I,SAMPLERATE) .EQ. 0 .AND. id .EQ. 0 ) THEN
                        KK = KK + 1
                        IF (WRITEONFLY) THEN
                                WRITE(200,'(A5,I10)')'ATOMS',KK
                                DO J=1,NATOMS
                                        f = R(I,J,:)*0.52917720859 ! positions saved in angstroem
                                        g = VEL(I,J,:)*0.52917720859 ! velocities in angstroem/au
                                        WRITE(200,'(I4,6(E30.20))')ATOMS(J)%Z,f,g
                                ENDDO
                       ENDIF
                 ENDIF
                 !===============================
                 !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      ENDDO
      

      !===================================================
      ! If specified not to save on-the-fly, we save after
      ! the calculation has been finished 
      !===================================================
      KK = 0
      IF ( .not. WRITEONFLY .AND. id .EQ. 0 ) THEN
                DO J=2,NTIMESTEPS
                        WRITE(100,'(I10,E30.20,E30.20,E30.20,E30.20)')J,ENERGY(J),EKIN(J),EPOT(J),TEMP(J)
                        IF ( MOVIE .AND. MOD(J,SAMPLERATE) .EQ. 0 ) THEN
                                KK = KK+1
                                WRITE(200,'(A5,I10)')'ATOMS',KK
                                DO JJ=1,NATOMS
                                       f = R(J,JJ,:)*0.52917720859 ! positions saved in angstroem
                                       g = VEL(J,JJ,:)*0.52917720859 ! velocities in angstroem/au
                                       WRITE(200,'(I4,6(E30.20))')ATOMS(JJ)%Z,f,g
                                ENDDO
                        ENDIF
               ENDDO
      ENDIF
      IF ( id .EQ. 0 ) THEN
        close(100)
        IF ( MOVIE ) close(200)
      ENDIF
      STOP
      
      END SUBROUTINE moleculardynamics
