SUBROUTINE dqmc(N,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,beta,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,EHF,nucE,ESIGMA,EDQMC)
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
        INTEGER, INTENT(IN)  :: N,NATOMS
        TYPE(BASIS), INTENT(IN) :: BAS
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: TIMESTEP,TEND,TSTART,beta,EHF,nucE,BJASTROW,CJASTROW
        INTEGER, INTENT(IN) :: SAMPLERATE,NREPLICAS
        DOUBLE PRECISION, INTENT(OUT) :: EDQMC, ESIGMA
        DOUBLE PRECISION, EXTERNAL :: trialfnk,EL,laplacetrialfnk,laplacehforbitalval,hforbitalval
        EXTERNAL :: guideforce
        DOUBLE PRECISION :: SIGMA,VECT1(3),VECT2(3)
        DOUBLE PRECISION :: EGUESS,DT,V
        DOUBLE PRECISION :: P1,P2,b,c,TEST
        INTEGER :: I,J,K,NTIMESTEPS,COUN,NUMBEROFNEWWALKERS,N3,NREP2,NSTART,NEND,M,clock
        INTEGER, ALLOCATABLE :: NWALK(:),seed(:)
        DOUBLE PRECISION, ALLOCATABLE :: ET(:),ER(:),E(:),force(:),xt(:),xtrial(:),xtemp(:,:),xold(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: xnew(:,:),x1(:),x2(:)
        REAL :: RAND,ratio,testa,ita,s
       

        ! Allocating some variables:
        N3 = 3*N
        NREP2 = 2*NREPLICAS
        ALLOCATE(xt(N3),xtrial(N3),xold(NREP2,N3))
        ALLOCATE(xnew(NREP2,N3),xtemp(NREP2,N3))
        
        b = BJASTROW    ! When b = 0, the below Metropolis generated distribution should give the
                        ! Hartree Fock energy, when monte carlo integration is performed, see ER(1)
        c = CJASTROW
       
        NTIMESTEPS = INT(TEND/TIMESTEP)
        
        ALLOCATE(NWALK(NTIMESTEPS),ET(NTIMESTEPS),ER(NTIMESTEPS),E(NTIMESTEPS))
        DT = TIMESTEP

        SIGMA = sqrt(DT)
        
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
        
        xold(:,:) = 0.0d0
        K = 0
        ER(1) = 0
        
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
                        xold(K,:) = xt
                ENDIF
                
        ENDDO
        
        ! The starting guess:
        ER(1) = ER(1)/K
        E(1) = ER(1)
        ET(1) = ER(1)
        NWALK(1) = NREPLICAS
        
        OPEN(12,FILE='ENERGYDQMC.dat',ACTION='WRITE')
        !---------------------------------------------------
        ! Here the Quantum Diffusion Monte Carlo loop starts
        !---------------------------------------------------
        print*,'   =========================================================='
        print*,'          Entering the Diffusion Monte Carlo loop            '
        print*,'   =========================================================='
        print*,' '
        
        DO I=2,NTIMESTEPS
                COUN = 0
                E(I) = 0.0d0
                WRITE(12,'(I20,F30.20,F30.20,I20)')I,ET(I-1)+nucE,E(I-1)+nucE,NWALK(I-1)
                !$OMP PARALLEL PRIVATE(x1,x2,P1,P2,force,ita,NUMBEROFNEWWALKERS,s,V,K,RAND,J,seed,clock) &
                !$OMP & SHARED(COUN,xnew,I,E,ER)
                ALLOCATE(force(N3),x1(N3),x2(N3))
                !$OMP DO
                DO J=1,NWALK(I-1)
                        !------------------------------------
                        ! Seeding the random number generator
                        !------------------------------------
                        IF ( I .EQ. 2 ) THEN
                                call random_seed(size=M)
                                ALLOCATE(seed(M))
                                CALL SYSTEM_CLOCK(COUNT=clock)
                                seed = clock + 37 * (/ (i - 1, i = 1, M) /)
                                CALL RANDOM_SEED(PUT = seed )
                                DEALLOCATE(seed)
                        ENDIF
                        !------------------------------------------
                        ! End of seeding.
                        !------------------------------------------

                        !------------------------------
                        ! The diffusion step/trial move
                        !------------------------------
                        
                        CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,xold(J,:),force)
                        
                        DO K=1,N
                                xtemp(J,3*(K-1) + 1 ) = xold(J,3*(K-1) + 1 ) + force(3*(K-1) + 1)*DT/2.0d0 + SIGMA*random_normal()
                                xtemp(J,3*(K-1) + 2 ) = xold(J,3*(K-1) + 2 ) + force(3*(K-1) + 2)*DT/2.0d0 + SIGMA*random_normal()
                                xtemp(J,3*(K-1) + 3 ) = xold(J,3*(K-1) + 3 ) + force(3*(K-1) + 3)*DT/2.0d0 + SIGMA*random_normal()
                        ENDDO
                        !----------------
                        ! Metropolis step
                        !----------------

                        x1 = xtemp(J,:) - xold(J,:) - force*DT/2.0d0
                        P1 = exp(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(x1,x1))*trialfnk(N,N3,Cup,Cdown,BAS,b,c,xold(J,:))**2
                        
                        CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,xtemp(J,:),force)
                        

                        x2 = xold(J,:) - xtemp(J,:) - force*DT/2.0d0
                        P2 = exp(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(x2,x2))*trialfnk(N,N3,Cup,Cdown,BAS,b,c,xtemp(J,:))**2
                        
                        CALL RANDOM_NUMBER(ita)
                        
                        IF ( ita .GT. REAL(P2/P1) ) THEN
                                xtemp(J,:) = xold(J,:)
                                COUN = COUN + 1
                                IF ( COUN > 2*NREPLICAS ) THEN
                                        WRITE(*,*)'======================================================================'
                                        WRITE(*,*)'            ABORTING SINCE NUMBER OF WALKERS > 2*NREPLICAS.'
                                        WRITE(*,*)' Try INCREASING NREPLICAS, CHANGING THE BETA PARAMETER OR THE TIMESTEP'
                                        WRITE(*,*)'======================================================================'
                                        STOP
                                ENDIF
                                xnew(COUN,:) = xtemp(J,:)
                                E(I) = E(I) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(J,:))
                        ENDIF
                        
                        !---------------------------------------------------------------------------------
                        ! Birth/death step if there is a new step accepted by  the Metropolis algorithm
                        ! we check if there is going to a birth or death of particles at the new position.
                        !---------------------------------------------------------------------------------

                        IF ( ita .LT. REAL(P2/P1) .OR. REAL(P2/P1) .GE. 1.0d0 ) THEN
                                V = 0.50d0*( EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(J,:)) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xold(J,:)) )
                                CALL RANDOM_NUMBER(RAND)
                                s = REAL(exp(-DT*(V-ER(I-1)))) + RAND
                                IF ( s .GE. 1.0 ) THEN
                                        NUMBEROFNEWWALKERS = floor(s)
                                        DO K=1,NUMBEROFNEWWALKERS
                                                COUN = COUN + 1
                                                IF ( COUN > 2*NREPLICAS ) THEN
                                                        WRITE(*,*)'======================================================================'
                                                        WRITE(*,*)'            ABORTING SINCE NUMBER OF WALKERS > 2*NREPLICAS.'
                                                        WRITE(*,*)' Try INCREASING NREPLICAS, CHANGING THE BETA PARAMETER OR THE TIMESTEP'
                                                        WRITE(*,*)'======================================================================'
                                                        STOP
                                                ENDIF
                                                xnew(COUN,:) = xtemp(J,:)
                                                E(I) = E(I) + EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xtemp(J,:))
                                        ENDDO
                                ENDIF
                        ENDIF
                ENDDO 
                !$OMP END DO
                !$OMP BARRIER
                DEALLOCATE(force,x1,x2)
                !$OMP END PARALLEL
                NWALK(I) = COUN
                E(I) = E(I)/COUN
                ET(I) = (ET(I-1)*(I-1) + E(I))/I
                ER(I) = ER(I-1) + beta*log(1.0d0*NREPLICAS/(1.0d0*NWALK(I)))
                xold = xnew
        ENDDO
        !-----------------------------
        ! Calculating the total energy
        !-----------------------------
        
        EDQMC = ET(NTIMESTEPS)
        
        NSTART = INT(TSTART/TIMESTEP)
        NEND = INT(TEND/TIMESTEP)
        
        IF (NSTART .EQ. 0 ) NSTART = 1
        
        !-----------------------------------
        ! Calculating the standard deviation
        ! of the total energy
        !-----------------------------------
        
        ESIGMA = 0.0d0
        DO I=NSTART,NEND
                ESIGMA =  ESIGMA + (ET(I)-EDQMC)**2
        ENDDO

        ESIGMA = sqrt(ESIGMA/(1.0d0*(NEND-NSTART+1)))
        
        WRITE(*,*)'======================================================================='
        WRITE(*,*)'        TOTAL ENERGY OF DIFFUSION QUANTUM MONTECARLO CALCULATION  '
        WRITE(*,*)'======================================================================='
        WRITE(*,*)' '
        WRITE(*,'(A9,F30.20,A5,F23.20,A5)')'E=',EDQMC+nucE,' +/- ',ESIGMA,' au'
        WRITE(*,*)' '
        WRITE(*,*)'======================================================================='
        CLOSE(12)
END SUBROUTINE dqmc
