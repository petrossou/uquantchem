SUBROUTINE dqmc(N,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,beta,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,EHF,nucE,ESIGMA,EDQMC,RESAMPVAL,NRECALC,NOREDIST,REDISTRIBUTIONFREQ,numprocessors,id)
        ! This subroutine calculates the ground state energy using Diffusion
        ! Quantum Monte Carlo with importance sampling. The trial function is
        ! that of a slaterdeterminant times a Jastrow factor ( see DQMC motes
        ! p.7 Eqn(16).
        ! N = Total number of electrons, NATOMS= Total number of atoms in
        ! molecule, Cup,Cdown = Hartree Fock basisset expansion coefficients,
        ! BAS = basis set used to solve HF-equations, ATOMS = (ATOMIC POSITION
        ! AND SPECIES ), b,c = Jastrow parameters, SAMPLERATE = elvery
        ! SAMPLERATE:th walker configuration is sampled to avoid sampling of
        ! correlated walker distributions. NREPLICAS= number of initial walkers,
        ! TIMESTEP = time step of the diffusion, TEND end time of diffusion
        ! walk, TSTART= we start calculating the energy mean value of the
        ! distribution here, EDQMC= total energy, ESIGMA = mean square deviation
        ! of calculated total energy. beta = energy updating factor
        ! EHF = Hartee Fock energy,ESIGMA = standard deviation of total energy
        ! EDQMC = total energy.
        USE datatypemodule
        USE random
        IMPLICIT NONE
        INCLUDE "mpif.h"
        INTEGER, INTENT(IN)  :: N,NATOMS,NRECALC,numprocessors,REDISTRIBUTIONFREQ
        TYPE(BASIS), INTENT(IN) :: BAS
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: TIMESTEP,TEND,TSTART,beta,EHF,nucE,BJASTROW,CJASTROW,RESAMPVAL
        INTEGER, INTENT(IN) :: SAMPLERATE,NREPLICAS,id
        LOGICAL, INTENT(IN) :: NOREDIST
        DOUBLE PRECISION, INTENT(OUT) :: EDQMC, ESIGMA
        DOUBLE PRECISION, EXTERNAL :: trialfnk,EL,laplacetrialfnk,laplacehforbitalval,hforbitalval
        EXTERNAL :: guideforce
        DOUBLE PRECISION :: SIGMA,VECT1(3),VECT2(3)
        DOUBLE PRECISION :: EGUESS,DT,V
        DOUBLE PRECISION :: P1,P2,b,c,TEST
        INTEGER :: I,J,K,NTIMESTEPS,NUMBEROFNEWWALKERS,COUN,N3,NREP2,NSTART,NEND,M,clock,ierr,NWALKTOT,NRT,scount,scountx,NREDIST,II
        INTEGER, ALLOCATABLE :: NWALK(:,:),seed(:),rcounts(:),displs(:),rcountsx(:),displsx(:),COUNC(:),offset(:)
        DOUBLE PRECISION, ALLOCATABLE :: ET(:),ER(:),E(:),force(:),xt(:),xtrial(:),xtemp(:,:),xold(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: xnew(:,:),x1(:),x2(:),xstart(:,:),EC(:),xsend(:,:),xresiev(:,:),ERMEAN(:)
        REAL :: RAND,ratio,testa,ita,s
        LOGICAL :: REDISTRIBUTED,DISTR

        REDISTRIBUTED = .FALSE.
        NREDIST = 0
        ! Allocating some variables:
        N3 = 3*N
        NREP2 = 2*NREPLICAS
        NWALKTOT = NREPLICAS
        ALLOCATE(xt(N3),xtrial(N3),xstart(N3,NREP2),offset(numprocessors))
        ALLOCATE(xtemp(N3,NREP2))
        ALLOCATE(force(N3),x1(N3),x2(N3))
        b = BJASTROW    ! When b = 0, the below Metropolis generated distribution should give the
                        ! Hartree Fock energy, when monte carlo integration is performed, see ER(1)
        c = CJASTROW
       
        NTIMESTEPS = INT(TEND/TIMESTEP)
        
        ALLOCATE(NWALK(numprocessors,NTIMESTEPS),ET(NTIMESTEPS),ER(NTIMESTEPS),E(NTIMESTEPS),ERMEAN(NTIMESTEPS))
        DT = TIMESTEP

        SIGMA = sqrt(DT)
        
        ! USED FOR THE MPI_GATHERV communication
        ALLOCATE(rcounts(numprocessors),displs(numprocessors),COUNC(numprocessors),EC(numprocessors))
        ALLOCATE(rcountsx(numprocessors),displsx(numprocessors))
        
        DO K=1,numprocessors
             rcounts(K) = 1
             displs(K) = K-1
         ENDDO
        
        !----------------------------------------------------------------
        ! Calculating the initial distribution corresponding to the trial
        ! function, using the Metroplplis algorithm:
        !---------------------------------------------------------------
        
        K = 0
        DO J=1,NATOMS
                DO I=1,ATOMS(J)%Z
                        K = K+1
                        xt(3*(K-1) + 1 ) = ATOMS(J)%R(1) + random_normal() 
                        xt(3*(K-1) + 2 ) = ATOMS(J)%R(2) + random_normal() 
                        xt(3*(K-1) + 3 ) = ATOMS(J)%R(3) + random_normal() 
                ENDDO
        ENDDO
        
        K = 0
        ER(:) = 0.0d0
        ERMEAN(:) = 0.0d0
        DO I=1,NREPLICAS*SAMPLERATE
                !-----------
                ! Trial move
                !-----------
                DO J=1,N
                        xtrial(3*(J-1) + 1 ) = xt(3*(J-1) + 1 ) + SIGMA*random_normal()
                        xtrial(3*(J-1) + 2 ) = xt(3*(J-1) + 2 ) + SIGMA*random_normal()
                        xtrial(3*(J-1) + 3 ) = xt(3*(J-1) + 3 ) + SIGMA*random_normal()
                ENDDO
                !----------------
                ! Metropolis step:
                !----------------
                
                ratio = REAL((trialfnk(N,N3,Cup,Cdown,BAS,b,c,xtrial)/trialfnk(N,N3,Cup,Cdown,BAS,b,c,xt))**2)
                
                CALL RANDOM_NUMBER(testa)
                !------------------------
                ! ACCEPT OR REJECT MOVE ?
                !------------------------
                IF ( ratio .GE. 1.0 .OR. testa .LT. ratio ) xt = xtrial

                IF ( MOD(I,SAMPLERATE) .EQ. 0 ) THEN
                        K = K+1
                        ER(1) = ER(1) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xt)
                        xstart(:,K) = xt
                ENDIF
                
        ENDDO
        
        ! The starting guess:
        ER(1) = ER(1)/K
        ERMEAN(1) = ER(1)
        II = 1
        E(1) = ER(1)
        ET(1) = ER(1)

        !--------------------------------------------------
        ! Dividing up the walkers on the different threads:
        !--------------------------------------------------
       
        ! The number of replicas/processor:
        
        NRT = INT(NREPLICAS/numprocessors)
        
        DO I=1,numprocessors
           NWALK(I,1) = NRT
        ENDDO

        DO I=1,MOD(NREPLICAS,numprocessors)
            NWALK(I,1) = NWALK(I,1) + 1
        ENDDO
        
        offset(1) = 0
        DO I=2,numprocessors
          offset(I) = NWALK(I-1,1)+offset(I-1)
        ENDDO

        ALLOCATE(xold(N3,5*NRT),xnew(N3,5*NRT))
        
        xold(:,:) = 0.0d0

        DO I=1,NWALK(id+1,1)
               xold(:,I) = xstart(:,offset(id+1)+I)
        ENDDO
        
        IF ( id .EQ. 0 ) OPEN(12,FILE='ENERGYDQMC.dat',ACTION='WRITE')
        !---------------------------------------------------
        ! Here the Quantum Diffusion Monte Carlo loop starts
        !---------------------------------------------------
        IF ( id .EQ. 0 ) THEN
         print*,'   =========================================================='
         print*,'          Entering the Diffusion Monte Carlo loop            '
         print*,'   =========================================================='
         print*,' '
        ENDIF
        DO I=2,NTIMESTEPS
              IF ( id .EQ. 0 ) WRITE(12,'(I20,F30.20,F30.20,F30.20,I20)')I,ET(I-1)+nucE,ER(I-1)+nucE,E(I-1)+nucE,NWALKTOT
             10 CONTINUE
                COUN = 0
                E(I) = 0.0d0
                !$OMP DO
                DO J=1,NWALK(id+1,I-1)
                        !------------------------------------
                        ! Seeding the random number generator
                        !------------------------------------
                        IF ( I .EQ. 2 .AND. J .EQ. 1 ) THEN
                                call random_seed(size=M)
                                ALLOCATE(seed(M))
                                CALL SYSTEM_CLOCK(COUNT=clock)
                                seed = clock + (37+id)* (/ (i - 1, i = 1, M) /)
                                CALL RANDOM_SEED(PUT = seed )
                                !DEALLOCATE(seed)
                        ENDIF
                       
                        !------------------------------------------
                        ! End of seeding.
                        !------------------------------------------

                        !------------------------------
                        ! The diffusion step/trial move
                        !------------------------------
                        
                        CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,xold(:,J),force)
                        
                        DO K=1,N
                                xtemp(3*(K-1) + 1,J ) = xold(3*(K-1) + 1,J ) + force(3*(K-1) + 1)*DT/2.0d0 + SIGMA*random_normal()
                                xtemp(3*(K-1) + 2,J ) = xold(3*(K-1) + 2,J ) + force(3*(K-1) + 2)*DT/2.0d0 + SIGMA*random_normal()
                                xtemp(3*(K-1) + 3,J ) = xold(3*(K-1) + 3,J ) + force(3*(K-1) + 3)*DT/2.0d0 + SIGMA*random_normal()
                        ENDDO
                        !----------------
                        ! Metropolis step
                        !----------------

                        x1 = xtemp(:,J) - xold(:,J) - force*DT/2.0d0
                        P1 = exp(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(x1,x1))*trialfnk(N,N3,Cup,Cdown,BAS,b,c,xold(:,J))**2
                        
                        CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,xtemp(:,J),force)
                        

                        x2 = xold(:,J) - xtemp(:,J) - force*DT/2.0d0
                        P2 = exp(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(x2,x2))*trialfnk(N,N3,Cup,Cdown,BAS,b,c,xtemp(:,J))**2
                        
                        CALL RANDOM_NUMBER(ita)

                        IF ( ita .GT. REAL(P2/P1) ) THEN
                                xtemp(:,J) = xold(:,J)
                                COUN = COUN + 1
                                IF ( COUN > 5*NRT ) THEN
                                        WRITE(*,*)'======================================================================'
                                        WRITE(*,*)'            ABORTING SINCE NUMBER OF WALKERS > 5*NREPLICAS.'
                                        WRITE(*,*)' Try INCREASING NREPLICAS, CHANGING THE BETA PARAMETER OR THE TIMESTEP'
                                        WRITE(*,*)'           OR SETTING THE PARAMETER REDISTRIBUTIONFREQ > 0            '
                                        WRITE(*,*)'======================================================================'
                                        STOP
                                ENDIF
                                xnew(:,COUN) = xtemp(:,J)
                                E(I) = E(I) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(:,J))
                        ENDIF
                        
                        !---------------------------------------------------------------------------------
                        ! Birth/death step if there is a new step accepted by  the Metropolis algorithm
                        ! we check if there is going to a birth or death of particles at the new position.
                        !---------------------------------------------------------------------------------

                        IF ( ita .LT. REAL(P2/P1) .OR. REAL(P2/P1) .GE. 1.0d0 ) THEN
                                V = 0.50d0*( EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(:,J)) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xold(:,J)) )
                                CALL RANDOM_NUMBER(RAND)
                                s = REAL(exp(-DT*(V-ER(I-1)))) + RAND
                                IF ( s .GE. 1.0 ) THEN
                                        NUMBEROFNEWWALKERS = floor(s)
                                        DO K=1,NUMBEROFNEWWALKERS
                                                COUN = COUN + 1
                                                IF ( COUN > 5*NRT ) THEN
                                                        WRITE(*,*)'======================================================================'
                                                        WRITE(*,*)'            ABORTING SINCE NUMBER OF WALKERS > 5*NREPLICAS.'
                                                        WRITE(*,*)' Try INCREASING NREPLICAS, CHANGING THE BETA PARAMETER OR THE TIMESTEP'
                                                        WRITE(*,*)'           OR SETTING THE PARAMETER REDISTRIBUTIONFREQ > 0            '
                                                        WRITE(*,*)'======================================================================'
                                                        STOP
                                                ENDIF
                                                xnew(:,COUN) = xtemp(:,J)
                                                E(I) = E(I) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(:,J))
                                        ENDDO
                                ENDIF
                        ENDIF
                ENDDO 
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                
                ! Collecting the number of walkers COUN at every thread into the array COUNC, so that 
                ! COUNC(id+1) = The number of walkers/replicas at thread id. COUNC is located at thread id = 0.
                call MPI_GATHERV(COUN,1,MPI_INTEGER,COUNC,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

                ! Collecting the Energy E(I) at every thread into the array EC, so that 
                ! EC(id+1) = The energy E(I) calulated at thread id. EC is located at thread id = 0.
                call MPI_GATHERV(E(I),1,MPI_DOUBLE_PRECISION,EC,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                
                ! Sending COUNC and EC to all the other threads
                CALL MPI_BCAST(COUNC,numprocessors, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(EC,numprocessors, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                
                NWALK(id+1,I) = COUN
                NWALKTOT = SUM(COUNC)
                E(I) = SUM(EC)/NWALKTOT
                ET(I) = (ET(I-1)*(I-1) + E(I))/I

                ! If a too pathological energy E(I) has been accepted 
                ! discard the new configuration of walkers and go 
                ! back and recalculate the population of walkers at 
                ! the same time step.
                
                IF ( abs((E(I-1)-E(I))/E(I-1)) .GT. RESAMPVAL ) GOTO 10 

                IF ( MOD(I,NRECALC) .EQ. 0 ) THEN
                        ER(I) = SUM(ERMEAN)/II + beta*log(1.0d0*NREPLICAS/(1.0d0*NWALKTOT))
                        II = II + 1
                        ERMEAN(II) = ER(I)
                ELSE
                        ER(I) = ER(I-1)
                ENDIF
                
                xold(:,:) = xnew(:,:)

                !------------------------------------------------------------------------------------
                !                              LOAD BALANCING STEP 
                !------------------------------------------------------------------------------------
                ! If all the walkers at one node are dead, or if one node runs the danger of 
                ! getting over-poulated by walkers/replicas we here redistribute the walkers eavenly
                !------------------------------------------------------------------------------------
                IF ( REDISTRIBUTIONFREQ .NE. 0 ) THEN
                       IF( MOD(I,REDISTRIBUTIONFREQ) .EQ. 0 ) THEN 
                              DISTR = .TRUE.
                       ELSE
                              DISTR = .FALSE.
                       ENDIF
                ELSE
                       DISTR = .FALSE.
                ENDIF

                IF ( ( MINVAL(COUNC,1) .EQ. 0 .OR. MAXVAL(COUNC,1) .GE. 2*NRT .OR. DISTR ) .AND. .not. NOREDIST ) THEN
                   
                   displsx(1) = 0
                   DO K=1,numprocessors
                      rcountsx(K) = N3*COUNC(K)
                      IF ( I .GE. 2 ) displsx(K) = displsx(K-1) + N3*COUNC(K-1)
                   ENDDO

                   IF ( id .EQ. 0 ) THEN
                      REDISTRIBUTED = .TRUE.
                      NREDIST = NREDIST + 1
                   ENDIF

                   ! Collecting all the walker-coordinates in xnew at the different threads and put them into xstart located at thread 0 (id=0)
                   CALL MPI_GATHERV(xnew(:,1:COUNC(id+1)),rcountsx(id+1),MPI_DOUBLE_PRECISION,xstart,rcountsx,displsx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                   
                   ! Sending the array of walkers, xstart, located at thread 0 (id = 0 ) to all other threads.
                   CALL MPI_BCAST(xstart,N3*NREP2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

                   ! The new number of replicas/processor:
                   DO K=1,numprocessors
                      NWALK(K,I) = INT(NWALKTOT/numprocessors)
                   ENDDO

                   DO K=1,MOD(NWALKTOT,numprocessors)
                      NWALK(K,I) = NWALK(K,I) + 1
                   ENDDO
        
                   offset(1) = 0
                   DO K=2,numprocessors
                      offset(K) = NWALK(K-1,I)+offset(K-1)
                   ENDDO
                   
                   ! Re-distributing the walkers/replicas on the different threads:
                   xold(:,:) = 0.0d0
                   DO K=1,NWALK(id+1,I)
                       xold(:,K) = xstart(:,offset(id+1)+K)
                   ENDDO
               ELSE
                  !----------------------------------------------------------------
                  !                     END OF LOAD BALANCING
                  !----------------------------------------------------------------
                  xold(:,:) = xnew(:,:)
              ENDIF

        ENDDO
        !-----------------------------
        ! Calculating the total energy
        !-----------------------------
        
        NSTART = INT(TSTART/TIMESTEP)
        NEND = INT(TEND/TIMESTEP)
        
        IF (NSTART .EQ. 0 ) NSTART = 1
        
        EDQMC = 0.0d0
        II = 0
        
        DO I=NSTART,NEND
             EDQMC = EDQMC + E(I) 
             II = II + 1
             ET(I) = EDQMC/II
        ENDDO

        EDQMC = EDQMC/II
        
        !-----------------------------------
        ! Calculating the standard deviation
        ! of the total energy
        !-----------------------------------
        
        ESIGMA = 0.0d0
        DO I=NSTART,NEND
                ESIGMA =  ESIGMA + (ET(I)-EDQMC)**2
        ENDDO

        ESIGMA = sqrt(ESIGMA/(1.0d0*(NEND-NSTART+1)))
        IF ( id .EQ. 0 ) THEN
         WRITE(*,*)'======================================================================='
         WRITE(*,*)'        TOTAL ENERGY OF DIFFUSION QUANTUM MONTECARLO CALCULATION  '
         WRITE(*,*)'======================================================================='
         WRITE(*,*)' '
         WRITE(*,'(A9,F30.20,A5,F23.20,A5)')'E=',EDQMC+nucE,' +/- ',ESIGMA,' au'
         WRITE(*,*)' '
         WRITE(*,*)'======================================================================='
         IF ( REDISTRIBUTED ) THEN
            WRITE(*,*)' '
            IF ( NREDIST .LT. 10000 ) THEN 
               WRITE(*,'(A40,I6,A6)')'The walkers/replicas were redistributed ',NREDIST,' TIMES'
               WRITE(*,*)' '
            ELSE
               WRITE(*,*)'The walkers/replicas were redistributed > 10000 TIMES'
               WRITE(*,*)' '
               WRITE(*,*)'======================================================================='
            ENDIF
         ENDIF
         OPEN(13,FILE='DQMCRESULTS.dat',ACTION='WRITE')
         DO I=NSTART,NEND
            WRITE(13,'(I20,F30.20,F30.20,F30.20)')I,ET(I)+nucE,ER(I)+nucE,E(I)+nucE
         ENDDO
         CLOSE(12)
         CLOSE(13)
        ENDIF
END SUBROUTINE dqmc
