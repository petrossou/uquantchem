SUBROUTINE dqmc(N,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,beta,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,EHF,nucE,ESIGMA,EDQMC,NPERSIST, & 
&               NRECALC,NOREDIST,REDISTRIBUTIONFREQ,CUTTOFFFACTOR,NVMC,VMCCALC,WRITEDENS,MESH,LIMITS,numprocessors,id)
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
        INTEGER, INTENT(IN)  :: N,NATOMS,NRECALC,numprocessors,REDISTRIBUTIONFREQ,NVMC,MESH(3)
        TYPE(BASIS), INTENT(IN) :: BAS
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS),LIMITS(3)
        DOUBLE PRECISION, INTENT(IN) :: TIMESTEP,TEND,TSTART,beta,EHF,nucE,BJASTROW,CJASTROW,CUTTOFFFACTOR
        INTEGER, INTENT(IN) :: SAMPLERATE,NPERSIST,id
        INTEGER, INTENT(INOUT) :: NREPLICAS
        LOGICAL, INTENT(IN) :: NOREDIST,VMCCALC,WRITEDENS
        DOUBLE PRECISION, INTENT(OUT) :: EDQMC, ESIGMA
        DOUBLE PRECISION, EXTERNAL :: trialfnk,EL,laplacetrialfnk,laplacehforbitalval,hforbitalval
        EXTERNAL :: guideforce
        DOUBLE PRECISION :: SIGMA,VECT1(3),VECT2(3),ESIGMAM,CORRECTION,fmean,ffmean,fpfmean,CCORR
        DOUBLE PRECISION :: EGUESS,DT,V,TAUEFF,ENT1,ENT2,EVMC,DELTAVOL,ETDUMMY
        DOUBLE PRECISION :: P1,P2,b,c,TEST,ET,ERMEAN,PSIT1,PSIT2,TOTALVIKT,EVMCCOLLECT(numprocessors)
        INTEGER :: I,J,K,NTIMESTEPS,NUMBEROFNEWWALKERS,COUN,N3,NREP2,NSTART,NEND,M,clock,ierr,NRT,scount,scountx,NREDIST,II,NACC,NTOOLD,NMETROPOLIS
        INTEGER :: NWALKTOT,NWALKTOTNEW,NMERGED,NWALKTOTKILLED
        INTEGER :: xcoord,ycoord,zcoord,JJ,status(MPI_STATUS_SIZE),NBEGIN,NWALKTOTSAVED
        INTEGER, ALLOCATABLE :: NWALK(:),NWALKNEW(:),seed(:),rcounts(:),displs(:),rcountsx(:),displsx(:),COUNC(:),offset(:),CNACC(:),NREJECT(:),NREJECTNEW(:),NREJECTKILLEDC(:)
        INTEGER, ALLOCATABLE :: displsreject(:),rcountsreject(:),NREJECTUPPDATE(:),NKILLWC(:),rcountskillx(:),rcountskill(:),displskillx(:),displskill(:),NREJECTKILLED(:)
        INTEGER, ALLOCATABLE :: NREJECTMERGED(:),NREJECTNEWC(:),NREJECTKILLEDOLD(:),EKILLEDOLD(:),NREJECTSAVED(:)
        DOUBLE PRECISION, ALLOCATABLE :: ER(:),E(:),force(:),xt(:),xtrial(:),xtemp(:,:),xold(:,:),xkilledc(:,:),EKILLEDC(:),WKILLEDC(:),OLDWEIGHT(:),xkilledold(:,:)
        DOUBLE PRECISION, ALLOCATABLE :: xoldSAVED(:,:),xkilledSAVED(:,:),WEIGHTSAVED(:),WKILLEDSAVED(:),EKILLEDSAVED(:)
        INTEGER, ALLOCATABLE :: NWALKSAVED(:),NREJECTKILLEDSAVED(:)
        INTEGER :: NKILLWSAVED
        DOUBLE PRECISION, ALLOCATABLE :: xnew(:,:),x1(:),x2(:),x3(:),x4(:),xstart(:,:),EC(:),xsend(:,:),xresiev(:,:),EG(:),EM(:),xmerged(:,:),WEIGHTMERGED(:),ECORR(:)
        REAL :: RAND,ratiov,testa,ita,s
        LOGICAL :: REDISTRIBUTED,DISTR,DEBUGGRANDOM,ALLDEAD,LOCKED
        LOGICAL :: first,NEXTSTEP,SINNGULAR1,SINNGULAR2,RESTART,finns
        DOUBLE PRECISION :: rvect(3),vvect(3),vlength,zvect(3),vect(3),distance,rnuc(3),z,zetaa,aa,vmeanvect(3),vmeanz,vmeanrho,rhovect(3),ratio
        DOUBLE PRECISION :: zbis, rhobis,dr(3),rprim(3),RANDDIR(3),ppobability,SR,SRP,rvectp(3),G1,G2,WTOT
        DOUBLE PRECISION, ALLOCATABLE :: VBIG(:), VBIGMEAN(:),VBIGP(:), VBIGMEANP(:),WEIGHT(:),NEWWEIGHT(:),xkilled(:,:),WTOTC(:),EKILLED(:),WKILLED(:),WKILLEDOLD(:)
        DOUBLE PRECISION, ALLOCATABLE :: NEWWEIGHTC(:),xstartc(:,:),density(:,:,:),densityc(:,:,:)
        INTEGER :: KK,atomnr,NKILLW,NKILLOLD,NEWARRAYSIZE,NSTUCK,NPOPMIN,ISAVED
        REAL :: qtilde,qtildep
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        
        NSTUCK = 0
        REDISTRIBUTED = .FALSE.
        DEBUGGRANDOM = .FALSE.
        NREDIST = 0
        first = .TRUE.
        ALLDEAD = .FALSE.
        LOCKED = .FALSE.

        ! Allocating some variables:
        N3 = 3*N
        NREP2 = 2*NREPLICAS
        NPOPMIN = 2*NREPLICAS
        NTOOLD = 0

        ALLOCATE(xt(N3),xtrial(N3),xstart(N3,NREP2),xstartc(N3,NREP2),offset(numprocessors),NKILLWC(numprocessors),NREJECTMERGED(NREP2),OLDWEIGHT(NREP2))
        ALLOCATE(xtemp(N3,NREP2),NREJECT(NREP2),NREJECTNEW(NREP2),NREJECTNEWC(NREP2),NREJECTUPPDATE(NREP2),xkilled(N3,NREP2),xkilledSAVED(N3,NREP2),xkilledc(N3,NREP2),xmerged(N3,NREP2))
        ALLOCATE(force(N3),x1(N3),x2(N3),x3(N3),x4(N3),WEIGHT(NREP2),WEIGHTSAVED(NREP2),NEWWEIGHT(NREP2),NEWWEIGHTC(NREP2),EKILLEDC(NREP2),WKILLEDC(NREP2))
        ALLOCATE(VBIG(N3),VBIGMEAN(N3),VBIGP(N3),VBIGMEANP(N3),EKILLED(NREP2),WKILLED(NREP2),WKILLEDSAVED(NREP2),WEIGHTMERGED(NREP2),NREJECTKILLED(NREP2),NREJECTKILLEDC(NREP2))
        ALLOCATE(xkilledold(N3,NREP2),WKILLEDOLD(NREP2),NREJECTKILLEDOLD(NREP2),EKILLEDOLD(NREP2),NREJECTKILLEDSAVED(NREP2),EKILLEDSAVED(NREP2),NREJECTSAVED(NREP2))
        
        IF ( WRITEDENS ) THEN
                ALLOCATE(density(MESH(1),MESH(2),MESH(3)))
                density = 0.0d0
                IF ( id .EQ. 0 ) ALLOCATE(densityc(MESH(1),MESH(2),MESH(3)))
        ENDIF
        
        b = BJASTROW    ! When b = 0, the below Metropolis generated distribution should give the
                        ! Hartree Fock energy, when monte carlo integration is performed, see ER(1)
        c = CJASTROW
        
        NTIMESTEPS = INT(TEND/TIMESTEP)
        NREJECT(:) = 0
        NREJECTNEW(:) = 0
        NREJECTUPPDATE(:) = 0
        ALLOCATE(NWALK(numprocessors),NWALKNEW(numprocessors),ER(NTIMESTEPS),E(NTIMESTEPS),NWALKSAVED(numprocessors))
        DT = TIMESTEP
        NWALKTOT = NREPLICAS
        SIGMA = sqrt(DT)
        TAUEFF = DT
        ! USED FOR THE MPI_GATHERV communication
        ALLOCATE(rcounts(numprocessors),displs(numprocessors),COUNC(numprocessors),EC(numprocessors),CNACC(numprocessors),WTOTC(numprocessors))
        ALLOCATE(rcountsx(numprocessors),displsx(numprocessors),rcountsreject(numprocessors),displsreject(numprocessors),rcountskillx(numprocessors))
        ALLOCATE(rcountskill(numprocessors),displskillx(numprocessors),displskill(numprocessors))
        
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
                        xt(3*(K-1) + 1 ) = ATOMS(J)%R(1) + random_normall() 
                        xt(3*(K-1) + 2 ) = ATOMS(J)%R(2) + random_normall() 
                        xt(3*(K-1) + 3 ) = ATOMS(J)%R(3) + random_normall() 
                ENDDO
        ENDDO
        
        K = 0
        ER(:) = 0.0d0
        ERMEAN = 0.0d0

        !-------------------------------------------------
        ! Here we check if a we are continuing a dqmc
        ! calculation or performing a new one from scratch
        !--------------------------------------------------
        RESTART = .FALSE.
        NBEGIN = 1
        NWALKTOTNEW = NREPLICAS

        inquire(file='DQMCRESTART.dat',exist=finns)
        IF ( finns .AND. .not. VMCCALC ) RESTART = .TRUE.

IF ( .not. RESTART ) THEN  !-start of restart if-conditional
       call random_seed(size=M)
       ALLOCATE(seed(M))
       CALL SYSTEM_CLOCK(COUNT=clock)
       seed = clock + (37+id)* (/ (i - 1, i = 1, M) /)
       CALL RANDOM_SEED(PUT = seed )
                                
        !================================
        ! HERE THE VMC CALCULATION STARTS
        !================================
        EVMC = 0.0d0
        NMETROPOLIS = NREPLICAS*SAMPLERATE
        IF ( VMCCALC ) THEN
                NMETROPOLIS = NVMC*SAMPLERATE
        ENDIF
        IF ( .not. VMCCALC .AND. NVMC .GT. NREPLICAS ) THEN
                NMETROPOLIS = NVMC*SAMPLERATE
        ENDIF
        ALLOCATE(ECORR(NMETROPOLIS))

        DO I=1,NMETROPOLIS
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

                ratiov = REAL((trialfnk(N,N3,Cup,Cdown,BAS,b,c,xtrial)/trialfnk(N,N3,Cup,Cdown,BAS,b,c,xt))**2)

                CALL RANDOM_NUMBER(testa)
                !------------------------
                ! ACCEPT OR REJECT MOVE ?
                !------------------------
                IF ( ratiov .GE. 1.0 .OR. testa .LT. ratiov ) xt = xtrial

                IF ( MOD(I,SAMPLERATE) .EQ. 0  ) THEN
                        K = K+1
                        ECORR(K) = EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,xt)
                        EVMC = EVMC + ECORR(K)
                        !---------------------------------------------
                        ! Collecting data for the density calculation:
                        !---------------------------------------------
                        IF ( WRITEDENS .AND. VMCCALC ) THEN
                                DO J=1,N
                                        xcoord =  NINT((xt(3*(J-1) + 1)+LIMITS(1))*(MESH(1)-1)/(2.0d0*LIMITS(1)))
                                        ycoord =  NINT((xt(3*(J-1) + 2)+LIMITS(2))*(MESH(2)-1)/(2.0d0*LIMITS(2)))
                                        zcoord =  NINT((xt(3*(J-1) + 3)+LIMITS(3))*(MESH(3)-1)/(2.0d0*LIMITS(3)))
                                ENDDO
                                IF ( xcoord .GE. 1 .AND. xcoord .LE. MESH(1) .AND. ycoord .GE. 1 .AND. ycoord .LE. MESH(2) .AND. zcoord .GE. 1 .AND. zcoord .LE. MESH(3) ) THEN
                                        density(xcoord,ycoord,zcoord) = density(xcoord,ycoord,zcoord) + 1.0d0
                                ENDIF
                        ENDIF
                        !----------------------------------------------
                        ! End of collecting chargedensity information
                        !----------------------------------------------
                        IF ( NMETROPOLIS - I .LE. (NREPLICAS+1)*SAMPLERATE ) THEN
                                KK = KK + 1
                                IF ( KK .LE. NREPLICAS ) THEN
                                        xstart(:,KK) = xt
                                ENDIF
                        ENDIF
                ENDIF

        ENDDO
        DEALLOCATE(seed)
        
        !================================
        ! HERE THE VMC CALCULATION ENDS
        !================================
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_GATHER(EVMC, 1, MPI_DOUBLE_PRECISION, EVMCCOLLECT(1+id), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(EVMCCOLLECT,numprocessors, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        ! the thread = 0 configuration has to serve as the starting  configuration
        CALL MPI_BCAST(xstart,N3*NREP2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        EVMC = SUM(EVMCCOLLECT)/(K*numprocessors)
        
        ! Accumumulating the charge-density:
        IF (WRITEDENS .AND. VMCCALC ) THEN
                IF ( id .NE. 0 ) THEN
                        call MPI_SEND( density, MESH(1)*MESH(2)*MESH(3), MPI_DOUBLE_PRECISION, 0, id, MPI_COMM_WORLD, ierr )
                ELSE
                        DO I=1,numprocessors-1
                                CALL MPI_RECV(densityc, MESH(1)*MESH(2)*MESH(3), MPI_DOUBLE_PRECISION, I, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                                density = density + densityc
                        ENDDO
                ENDIF
        ENDIF
        IF ( VMCCALC .AND. id .EQ. 0 ) THEN
                ! CALCULATING THE CORRELATION FUNCTION BETWEEN CONSECUTIVE 
                ! METROPOLIS CONFIGURATIONS. SEE KOONIN PAGE 214, EQN (8.17)
                fmean = 0.0d0
                ffmean = 0.0d0
                fpfmean = 0.0d0
                DO I=1,K-1
                        fpfmean = fpfmean + ECORR(I)*ECORR(I+1)
                ENDDO
                fpfmean = fpfmean/(K-1)
                fmean = SUM(ECORR)/K
                ffmean = DOT_PRODUCT(ECORR,ECORR)/K
                CCORR = ( fpfmean - fmean**2)/(ffmean-fmean**2)
                WRITE(*,*)'======================================================================='
                WRITE(*,*)'        TOTAL ENERGY OF VARIATIONAL QUANTUM MONTECARLO CALCULATION  '
                WRITE(*,*)'======================================================================='
                WRITE(*,*)' '
                WRITE(*,'(A14,F30.20,A5)')' E =  ',EVMC+nucE,'au'
                WRITE(*,*)' '
                WRITE(*,*)'======================================================================='
                WRITE(*,*)'      CORRELATION BETWEEN CONSECUTIVE METROPOLIS CONFIGURATIONS:       '
                WRITE(*,*)'======================================================================='
                WRITE(*,*)' '
                WRITE(*,'(A14,F30.20)')' C(1) =  ',CCORR
                WRITE(*,*)' '
                WRITE(*,*)'======================================================================='
                IF ( WRITEDENS ) THEN
                        ! RENORMALIZING THE CHARGEDENSITY
                        DELTAVOL = 8.0d0*(LIMITS(1)/(1.0d0*(MESH(1)-1)))*(LIMITS(2)/(1.0d0*(MESH(2)-1)))*(LIMITS(3)/(1.0d0*(MESH(3)-1)))
                        density = density*N/( DELTAVOL*SUM(reshape(density,(/MESH(1)*MESH(2)*MESH(3)/))))
                        WRITE(*,*)'================================================================'
                        WRITE(*,*)'Charge density has been calculated and is being  saved to file. '
                        WRITE(*,*)'================================================================'
                        2 FORMAT(I6,I6,I6,F30.20)
                        OPEN(9,FILE='CHARGEDENS.dat',ACTION='WRITE')
                        DO I=1,MESH(1)
                                DO J=1,MESH(2)
                                        DO K=1,MESH(3)
                                                WRITE(9,FMT=2)I,J,K,density(I,J,K)
                                        ENDDO
                                ENDDO
                        ENDDO
                ENDIF

                STOP
        ELSE
                IF ( id .EQ. 0 ) WRITE(*,'(A14,F30.20,A5)')'E(VMC) =',EVMC+nucE,'au'
        ENDIF
        !-------------------------------------------------------------
        ! The continuation of a previously saved dqmc run is started
        ! from this point in the code.
        !-------------------------------------------------------------
ELSE   !-middle part of restart if-conditional
        OPEN(14,FILE='DQMCRESTART.dat',STATUS='OLD',ACTION='READ')
        READ(14,*)NBEGIN,NWALKTOTNEW
        READ(14,*)E(NBEGIN),ER(NBEGIN),ET,EVMC
        DO K=1,NWALKTOTNEW
                READ(14,*)xstart(:,K)
        ENDDO
        DO K=1,NWALKTOTNEW
                READ(14,'(I20,F30.20)')NREJECTNEWC(K),NEWWEIGHTC(K)
        ENDDO
        close(14)


ENDIF !-final part of restart if-conditional
      ! The starting guess:
      
      NBEGIN = NBEGIN + 1
      NREPLICAS = NWALKTOTNEW
      NWALKTOT = NWALKTOTNEW
      
      ER(1) = EVMC
        
        IF ( ER(1) .GT.  EHF ) ER(1) = EHF
           
        ERMEAN = ER(1)
        II = 1
        E(1) = ER(1)
        IF ( .not. RESTART ) ET = ER(1)

        !--------------------------------------------------
        ! Dividing up the walkers on the different threads:
        !--------------------------------------------------
        
        ! The number of replicas/processor:
        
        NRT = INT(NREPLICAS/numprocessors)
        
        DO I=1,numprocessors
           NWALK(I) = NRT
        ENDDO

        DO I=1,MOD(NREPLICAS,numprocessors)
            NWALK(I) = NWALK(I) + 1
        ENDDO
        
        offset(1) = 0
        DO I=2,numprocessors
          offset(I) = NWALK(I-1)+offset(I-1)
        ENDDO

        ALLOCATE(xold(N3,10*NRT),xnew(N3,10*NRT),xoldSAVED(N3,10*NRT))
        
        xold(:,:) = 0.0d0

        DO I=1,NWALK(id+1)
               xold(:,I) = xstart(:,offset(id+1)+I)
               IF ( RESTART ) THEN
                       WEIGHT(I)  = NEWWEIGHTC(offset(id+1)+I)
                       NREJECT(I) = NREJECTNEWC(offset(id+1)+I)
               ELSE
                       WEIGHT(I) = 1.0d0
                       NREJECT(I) = 0
               ENDIF
        ENDDO
        
        NWALKTOTKILLED = 0

        IF ( id .EQ. 0 ) THEN
                IF ( .not. RESTART ) THEN 
                        OPEN(12,FILE='ENERGYDQMC.dat',ACTION='WRITE')
                ELSE
                        ! Continue to write to the file
                        OPEN(12,FILE='ENERGYDQMC.dat',ACCESS="SEQUENTIAL")
                        DO K=2,NBEGIN-1
                                READ(12,'(I20,F30.20,F30.20,F30.20,I20,I10,I10)')I,ETDUMMY,ER(K-1),E(K-1),NWALKTOT,NWALKTOTKILLED,NMERGED
                                ER(K-1) = ER(K-1) - nucE
                                E(K-1)  =  E(K-1) - nucE
                        ENDDO
                ENDIF
        ENDIF
        !---------------------------------------------------
        ! Here the Quantum Diffusion Monte Carlo loop starts
        !---------------------------------------------------
        IF ( id .EQ. 0 ) THEN
         print*,'   =========================================================='
         print*,'          Entering the Diffusion Monte Carlo loop            '
         print*,'   =========================================================='
         print*,' '
        ENDIF
        
        !DEBUGGRANDOM = .TRUE.
        NKILLW = 0
        NMERGED = 0
        DO I=NBEGIN,NTIMESTEPS
                IF ( id .EQ. 0 ) WRITE(12,'(I20,F30.20,F30.20,F30.20,I20,I10,I10)')I,ET+nucE,ER(I-1)+nucE,E(I-1)+nucE,NWALKTOT,NWALKTOTKILLED,NMERGED
             20 CONTINUE    
                COUN = 0
                COUNC = 0
                NKILLWC = 0
                WTOT = 0.0d0
                NACC = 0 
                E(I) = 0.0d0
                DO J=1,NWALK(id+1)
                        !------------------------------------
                        ! Seeding the random number generator
                        !------------------------------------
                        IF ( .not. RESTART ) THEN
                                IF ( I .EQ. 2 .AND. J .EQ. 1 ) THEN
                                        call random_seed(size=M)
                                        ALLOCATE(seed(M))
                                        CALL SYSTEM_CLOCK(COUNT=clock)
                                        seed = clock + (37+id)* (/ (i - 1, i = 1, M) /)
                                        CALL RANDOM_SEED(PUT = seed )
                                        !DEALLOCATE(seed)
                                ENDIF
                        ELSE
                                IF ( I .EQ. NBEGIN .AND. J .EQ. 1 ) THEN
                                        call random_seed(size=M)
                                        ALLOCATE(seed(M))
                                        CALL SYSTEM_CLOCK(COUNT=clock)
                                        seed = clock + (37+id)* (/ (i - 1, i = 1, M) /)
                                        CALL RANDOM_SEED(PUT = seed )
                                        !DEALLOCATE(seed)
                                ENDIF
                        ENDIF
                       
                        !------------------------------------------
                        ! End of seeding.
                        !------------------------------------------
                        !--------------------------------------------------------------
                        ! All the numbers in parenthesis to the right are in  reference 
                        ! to the numbering of the steps in Umrigars algorithm, 
                        ! of the so called improved algorithm on pages 2886-2887 in
                        ! J. Chem. Phys. 99, 2865-2890, (1993)
                        !--------------------------------------------------------------
                        !===========================================
                        ! The the forward diffusion step/trial move
                        !===========================================
                        x3 = xold(:,J)
                        PSIT1 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,x3)
                        
                        IF ( PSIT1 .EQ. 0.0d0 ) SINNGULAR1 = .TRUE.
                        IF ( PSIT1 .NE. 0.0d0 ) SINNGULAR1 = .FALSE.

                        P1 = PSIT1**2 ! Initializing the forward move greens function G(R',R) (8)
                        CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,x3,force) ! Calculating the forward force/drift-velocity (10)
                        VBIG = force
                        !--------------------
                        ! loop over electrons
                        !--------------------
                        DO K=1,N
                                rvect(1) = x3(3*(K-1) + 1)
                                rvect(2) = x3(3*(K-1) + 2)
                                rvect(3) = x3(3*(K-1) + 3)

                                vvect(1) = force(3*(K-1) + 1)
                                vvect(2) = force(3*(K-1) + 2)
                                vvect(3) = force(3*(K-1) + 3)

                                vlength  = sqrt(DOT_PRODUCT(vvect,vvect))

                                zvect = 0.0d0

                                CALL findclosestatom(NATOMS,ATOMS,rvect,z,zvect,atomnr,rnuc)     ! (11)

                                CALL preparedifusion(atomnr,rnuc,DT,vlength,vvect,zvect,z,zetaa,vmeanvect,aa,vmeanz,vmeanrho,rhovect,zbis,rhobis,dr,qtilde) ! (12-19)

                                CALL RANDOM_NUMBER(RAND)
                                IF ( RAND .LT. 1.0d0 - qtilde ) THEN
                                        ! Calculating one of the two possible forward trial positions r' (20-22)
                                        rprim(1) = dr(1) + SIGMA*random_normall()
                                        rprim(2) = dr(2) + SIGMA*random_normall()
                                        rprim(3) = dr(3) + SIGMA*random_normall()
                                ELSE
                                        ! Calculating the other of the two possible forward trial positions r' (20-22)
                                        RANDDIR(1) = random_normall()
                                        RANDDIR(2) = random_normall()
                                        RANDDIR(3) = random_normall()
                                        rprim = rnuc + ( RANDDIR/sqrt(DOT_PRODUCT(RANDDIR,RANDDIR)))*random_chisq(6,first)/(4.0d0*zetaa)
                                        IF ( first ) first = .FALSE.
                                ENDIF
                                ! debugging:
                                !IF ( zetaa .NE. zetaa .OR. isnan(DOT_PRODUCT(rprim-dr,rprim-dr)) .OR. isnan(DOT_PRODUCT(rprim-rnuc,rprim-rnuc))  ) THEN
                                !   WRITE(*,'(A20,4(F15.10))')'ERROR AT POSITION -8',zetaa,DOT_PRODUCT(rprim,rprim),DOT_PRODUCT(rnuc,rnuc),DOT_PRODUCT(dr,dr)
                                !   
                                !ENDIF

                                G1 = ((2*PI*DT)**(-1.50d0))*EXP(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(rprim-dr,rprim-dr)) ! (27)
                                G2 = (zetaa**3/PI)*EXP(-2.0d0*zetaa*sqrt(DOT_PRODUCT(rprim-rnuc,rprim-rnuc))) ! (27)

                                P1 = P1*( G1*(1-qtilde) + G2*qtilde )

                                ! This is where the new R' position is saved
                                xtemp(3*(K-1) + 1,J) = rprim(1)
                                xtemp(3*(K-1) + 2,J) = rprim(2)
                                xtemp(3*(K-1) + 3,J) = rprim(3)

                                ! This is where the mean drift-velocity <V> is saved:
                                VBIGMEAN(3*(K-1) + 1) = vmeanvect(1)
                                VBIGMEAN(3*(K-1) + 2) = vmeanvect(2)
                                VBIGMEAN(3*(K-1) + 3) = vmeanvect(3)
                        ENDDO
               
                        !=============================
                        ! Start of the reversed move:
                        !=============================
               
                        x4 = xtemp(:,J)
                        PSIT2 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,x4)

                        IF ( PSIT2 .EQ. 0.0d0 ) SINNGULAR2 = .TRUE.
                        IF ( PSIT2 .NE. 0.0d0 ) SINNGULAR2 = .FALSE.
                        
                        ! This is basically where the fixed node approximation comes in (29)
                        IF ( PSIT1*PSIT2 .LE. 0.0d0 ) THEN
                                ppobability = 0.0
                                NREJECT(J) = NREJECT(J) + 1
                                xtemp(:,J) = xold(:,J)
                        ELSE
                                P2 = PSIT2**2
                                CALL guideforce(N,N3,Cup,Cdown,BAS,b,c,x4,force) ! Calculating the reverse force/drift-velocity (10)
                                VBIGP = force

                               !--------------------
                               ! loop over electrons
                               !--------------------
                               DO K=1,N
                                    rvect(1) = x3(3*(K-1) + 1)
                                    rvect(2) = x3(3*(K-1) + 2)
                                    rvect(3) = x3(3*(K-1) + 3)

                                    rvectp(1) = x4(3*(K-1) + 1)
                                    rvectp(2) = x4(3*(K-1) + 2)
                                    rvectp(3) = x4(3*(K-1) + 3)

                                    vvect(1) = force(3*(K-1) + 1)
                                    vvect(2) = force(3*(K-1) + 2)
                                    vvect(3) = force(3*(K-1) + 3)

                                    vlength  = sqrt(DOT_PRODUCT(vvect,vvect))

                                    zvect = 0.0d0

                                    CALL findclosestatom(NATOMS,ATOMS,rvectp,z,zvect,atomnr,rnuc)   ! (11)

                                    CALL preparedifusion(atomnr,rnuc,DT,vlength,vvect,zvect,z,zetaa,vmeanvect,aa,vmeanz,vmeanrho,rhovect,zbis,rhobis,dr,qtildep)  ! (12-19)
                                    
                                    ! debugging
                                    !IF ( zetaa .NE. zetaa .OR. isnan(DOT_PRODUCT(rprim-dr,rprim-dr)) .OR. isnan(DOT_PRODUCT(rprim-rnuc,rprim-rnuc))  ) THEN
                                    !   WRITE(*,'(A20,4(F15.10))')'ERROR AT POSITION -7',zetaa,DOT_PRODUCT(rprim,rprim),DOT_PRODUCT(rnuc,rnuc),DOT_PRODUCT(dr,dr)
                                    !ENDIF

                                    G1 = ((2*PI*DT)**(-1.50d0))*EXP(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(rvect-dr,rvect-dr))
                                    G2 = (zetaa**3/PI)*EXP(-2.0d0*zetaa*sqrt(DOT_PRODUCT(rvect-rnuc,rvect-rnuc)))

                                    P2 = P2*( G1*(1-qtildep) + G2*qtildep )

                                    ! This is where the mean drift-velocity <V> is saved:
                                    VBIGMEANP(3*(K-1) + 1) = vmeanvect(1)
                                    VBIGMEANP(3*(K-1) + 2) = vmeanvect(2)
                                    VBIGMEANP(3*(K-1) + 3) = vmeanvect(3)
                               ENDDO

                               ratio = (P2/P1)

                               ! To avoid persistent configuration catastrophes we here
                               ! employ the algorithm of Umrigar et al, Chem. Phys. 99, 2865 (1993)
                               ! For details see page 2872 in the above reference.
                               IF ( NREJECT(J) .GT. NPERSIST ) THEN
                                    ratio = ratio*1.10d0**(NREJECT(J)-NPERSIST)
                                    NTOOLD = NTOOLD+1
                               ENDIF

                               IF ( ratio .GT. 1.0 ) ratio = 1.0

                               CALL RANDOM_NUMBER(ita)

                               IF ( ita .GT. ratio ) THEN
                                    xtemp(:,J) = xold(:,J)
                                    NREJECT(J) = NREJECT(J) + 1
                                    ppobability = ratio
                               ENDIF

                               IF ( ita .LE. ratio ) THEN
                                    NACC = NACC + 1
                                    NREJECT(J) = 0
                                    ppobability = ratio
                               ENDIF

                                !---------------------------------------------
                                ! Collecting data for the density calculation:
                                !---------------------------------------------
                                IF ( WRITEDENS ) THEN
                                        DO JJ=1,N
                                                xcoord =  NINT((xtemp(3*(JJ-1) + 1,J)+LIMITS(1))*(MESH(1)-1)/(2.0d0*LIMITS(1)))
                                                ycoord =  NINT((xtemp(3*(JJ-1) + 2,J)+LIMITS(2))*(MESH(2)-1)/(2.0d0*LIMITS(2)))
                                                zcoord =  NINT((xtemp(3*(JJ-1) + 3,J)+LIMITS(3))*(MESH(3)-1)/(2.0d0*LIMITS(3)))
                                        ENDDO
                                        IF ( xcoord .GE. 1 .AND. xcoord .LE. MESH(1) .AND. ycoord .GE. 1 .AND. ycoord .LE. MESH(2) .AND. zcoord .GE. 1 .AND. zcoord .LE. MESH(3) ) THEN
                                                density(xcoord,ycoord,zcoord) = density(xcoord,ycoord,zcoord) + ppobability
                                        ENDIF
                                ENDIF
                        !----------------------------------------------
                        ! End of collecting chargedensity information
                        !----------------------------------------------
                        ENDIF

                        !----------------------------------------------------------------------------------------
                        ! Birth/death step if there is a new step accepted by  the Metropolis algorithm
                        ! we check if there is going to a birth or death (merge) of particles at the new position.
                        !-----------------------------------------------------------------------------------------
                        x3 = xold(:,J)
                        x4 = xtemp(:,J)
                        ENT1 = EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,x3)
                        ENT2 = EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,x4)
                        
                        IF ( DOT_PRODUCT(VBIGP,VBIGP) .NE. 0.0d0 ) THEN
                                 SRP = ER(I-1) - ET + ( ET - ENT2 )*sqrt(DOT_PRODUCT(VBIGMEANP,VBIGMEANP)/DOT_PRODUCT(VBIGP,VBIGP))      ! (36)
                        ELSE
                                 SRP = ER(I-1) - ENT2 
                        ENDIF
                        
                        IF ( DOT_PRODUCT(VBIG,VBIG) .NE. 0.0d0 ) THEN
                                 SR =  ER(I-1) - ET + ( ET - ENT1 )*sqrt(DOT_PRODUCT(VBIGMEAN,VBIGMEAN)/DOT_PRODUCT(VBIG,VBIG))          ! (36)
                        ELSE
                                 SR = ER(I-1) - ENT1
                        ENDIF

                        V = 0.50d0*(SRP + SR)*ppobability + (1-ppobability)*SR

                        WEIGHT(J) = WEIGHT(J)*exp(TAUEFF*V)   ! (38)

                        NUMBEROFNEWWALKERS = NINT(WEIGHT(J))
                        
                        IF ( WEIGHT(J) .GT. 1.0d0*NREPLICAS ) THEN
                           SINNGULAR2 = .TRUE.
                           SINNGULAR1 = .TRUE.
                           WEIGHT(J) = 0.0d0
                        ENDIF

                        IF ( SINNGULAR2  .OR. SINNGULAR1 ) NUMBEROFNEWWALKERS = 0

                        ! debugging:
                        !IF ( qtilde .NE. qtilde .OR. qtildep .NE. qtildep  ) WRITE(*,*)'ERROR AT POSITION -6',qtilde,qtildep
                        
                        ! debugging:
                        !IF ( P1 .NE. P1 .OR. P2 .NE. P2 .OR. P1 .EQ. 0.0d0 ) WRITE(*,*)'ERROR AT POSITION -5',P1,P2
                        
                        ! debugging:
                        !IF ( ppobability .NE. ppobability .OR. ratio .NE. ratio  ) WRITE(*,*)'ERROR AT POSITION -4',ppobability,ratio
                        
                        ! debugging:
                        !IF ( DOT_PRODUCT(VBIGMEAN,VBIGMEAN) .NE. DOT_PRODUCT(VBIGMEAN,VBIGMEAN) ) WRITE(*,*)'ERROR AT POSITION -3',DOT_PRODUCT(VBIGMEAN,VBIGMEAN)
                        
                        ! debugging:
                        !IF ( DOT_PRODUCT(VBIG,VBIG) .NE. DOT_PRODUCT(VBIG,VBIG) ) WRITE(*,*)'ERROR AT POSITION -2',DOT_PRODUCT(VBIG,VBIG)
                        
                        ! debugging:
                        !IF ( WEIGHT(J) .NE. WEIGHT(J) .OR. WEIGHT(J) .GT. 1000000 ) WRITE(*,*)'ERROR AT POSITION -1',WEIGHT(J),V,NUMBEROFNEWWALKERS
                        
                        ! Here new walkers are born ( See, Chem. Phys. 99, 2865 (1993), page 2867 column 2 ) 
                        IF ( NUMBEROFNEWWALKERS .GE. 1  ) THEN
                                IF ( NUMBEROFNEWWALKERS .GT. 2 ) NUMBEROFNEWWALKERS = 2
                                DO K=1,NUMBEROFNEWWALKERS
                                        IF ( COUN .LT. 2*NREPLICAS ) THEN
                                                COUN = COUN + 1
                                                xnew(:,COUN) = x4
                                                NEWWEIGHT(COUN) = WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                                WTOT = WTOT + WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                                NREJECTUPPDATE(COUN) = NREJECT(J)
                                                E(I) = E(I) + ( ENT2*ppobability + ENT1*(1-ppobability) )*WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                        ENDIF
                                ENDDO
                        ELSE
                                ! Here the non singular rejected walkers are saved in order to be merged in pairs after the loop over walkers
                                ! have been excited (  See, Chem. Phys. 99, 2865 (1993), page 2867, column 2 )
                                IF ( .not. SINNGULAR2  .AND. .not. SINNGULAR1 ) THEN
                                   NKILLW = NKILLW + 1
                                   EKILLED(NKILLW) = ( ENT2*ppobability + ENT1*(1-ppobability) )
                                   WKILLED(NKILLW) = WEIGHT(J)
                                   NREJECTKILLED(NKILLW) = NREJECT(J)
                                   xkilled(:,NKILLW) = x4
                                   ! debugging:
                                   !IF ( WEIGHT(J) .NE. WEIGHT(J) .OR. WEIGHT(J) .GT. 1000000 ) WRITE(*,*)'ERROR AT POSITION 0',WEIGHT(J)
                                ENDIF
                                !----------------------------------------------
                        ENDIF
                ENDDO ! the loop over walkers stops here
                ! debugging:
                !IF ( WTOT .NE. WTOT .OR. WTOT .GT. 1000000 ) WRITE(*,*)'ERROR AT POSITION 1',WTOT

                ! Collecting the number of killed walkers at every node, NKILLW, and sending them to 
                ! the array, NKILLWC, located at thread 0.
                call MPI_GATHERV(NKILLW,1,MPI_INTEGER,NKILLWC,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                ! Collecting the number of walkers COUN at every thread into the array COUNC, so that 
                ! COUNC(id+1) = The number of walkers/replicas at thread id. COUNC is located at thread id = 0.
                call MPI_GATHERV(COUN,1,MPI_INTEGER,COUNC,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                ! Collecting the weights of the  walkers WTOT at every thread into the  array WTOTC, so that 
                ! WTOTC(id+1) = The number of walkers/replicas at thread id. WTOTC is located at thread id = 0.
                call MPI_GATHERV(WTOT,1,MPI_DOUBLE_PRECISION,WTOTC,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                ! Collecting the number of accepted moves NACC at every thread into the array CNACC, so that 
                ! CNACC(id+1) = The number of walkers/replicas at thread id. CNACC is located at thread id = 0.
                call MPI_GATHERV(NACC,1,MPI_INTEGER,CNACC,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                ! Collecting the Energy E(I) at every thread into the array EC, so that 
                ! EC(id+1) = The energy E(I) calulated at thread id. EC is located at thread id = 0.
                call MPI_GATHERV(E(I),1,MPI_DOUBLE_PRECISION,EC,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                ! Sending COUNC, WTOTC, CNACC and EC to all the other threads
                CALL MPI_BCAST(COUNC,numprocessors, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(NKILLWC,numprocessors, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(WTOTC,numprocessors, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(CNACC,numprocessors, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(EC,numprocessors, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                ! here we prepair the gathering of the data in EKILLED,WKILLED and xkilled:
                displskillx(1) = 0
                DO K=1,numprocessors
                    rcountskillx(K) = N3*NKILLWC(K)
                    IF ( K .GE. 2 ) displskillx(K) = displskillx(K-1) + N3*NKILLWC(K-1)
                ENDDO

                displskill(1) = 0
                DO K=1,numprocessors
                      rcountskill(K) = NKILLWC(K)
                      IF ( K .GE. 2 ) displskill(K) = displskill(K-1) + NKILLWC(K-1)
                ENDDO
                ! Collecting all the walker-coordinates in xkilled at the different threads and put them into xkilledc located at thread 0 (id=0)
                CALL MPI_GATHERV(xkilled(:,1:NKILLWC(id+1)),rcountskillx(id+1),MPI_DOUBLE_PRECISION,xkilledc,rcountskillx,displskillx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                ! Collecting all the energies of the walkers that have been killed at the different threads, EKILLED, and put them in EKILLEDC located at the thread 0 (id = 0 )
                CALL MPI_GATHERV(EKILLED(1:NKILLWC(id+1)),rcountskill(id+1),MPI_DOUBLE_PRECISION,EKILLEDC,rcountskill,displskill,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                ! Collecting all the weights of the walkers that have been killed at the different threads, EKILLED, and put them in EKILLEDC located at the thread 0 (id = 0 )
                CALL MPI_GATHERV(WKILLED(1:NKILLWC(id+1)),rcountskill(id+1),MPI_DOUBLE_PRECISION,WKILLEDC,rcountskill,displskill,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                ! Collecting all the number of times a killed walker has been denied a move, NREJECTKILLED, and put them in NREJECTKILLEDC located at the thread id=0
                CALL MPI_GATHERV(NREJECTKILLED(1:NKILLWC(id+1)),rcountskill(id+1),MPI_INTEGER,NREJECTKILLEDC,rcountskill,displskill,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
                
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
               
                ! Sending xkilledc,EKILLEDC and WKILLEDC to all the other threads
                CALL MPI_BCAST(xkilledc,N3*NREP2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                
                CALL MPI_BCAST(EKILLEDC,NREP2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(WKILLEDC,NREP2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(NREJECTKILLEDC,NREP2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                NWALKNEW(id+1) = COUN
                NWALKTOTNEW = SUM(COUNC)
                TOTALVIKT = SUM(WTOTC)
                NWALKTOTKILLED = SUM(NKILLWC)
                
                ! debugging:
                !IF ( id .EQ. 0 .AND. TOTALVIKT .NE. TOTALVIKT .OR. TOTALVIKT .GE. 1000000)WRITE(*,*)'ERROR AT POSITION 2'
                
                K = 1
                NMERGED = 0
                !-----------------------------------------------
                ! Merging pairs of walkers that have been killed
                !-----------------------------------------------
                DO WHILE ( K .LE. NWALKTOTKILLED )
                   IF ( K+1 .LE. NWALKTOTKILLED ) THEN
                        NMERGED = NMERGED +1
                        WEIGHTMERGED(NMERGED) = WKILLEDC(K) + WKILLEDC(K+1)
                        TOTALVIKT = TOTALVIKT + WEIGHTMERGED(NMERGED)
                        IF ( id .EQ. 0 ) CALL RANDOM_NUMBER(RAND)
                        CALL MPI_BCAST(RAND,1, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        IF ( RAND .LE. WKILLEDC(K+1)/WEIGHTMERGED(NMERGED)  ) THEN
                             xmerged(:,NMERGED) = xkilledc(:,K+1)
                             NREJECTMERGED(NMERGED) = NREJECTKILLEDC(K+1)
                             E(I) = SUM(EC) + EKILLEDC(K+1)*WEIGHTMERGED(NMERGED)
                        ELSE
                             xmerged(:,NMERGED) = xkilledc(:,K)
                             NREJECTMERGED(NMERGED) = NREJECTKILLEDC(K)
                             E(I) = SUM(EC) + EKILLEDC(K)*WEIGHTMERGED(NMERGED)
                        ENDIF
                    ENDIF
                    K= K+2
                ENDDO
                
                ! debugging:
                IF ( id .EQ. 0 .AND. ( TOTALVIKT .NE. TOTALVIKT .OR. TOTALVIKT .GT. 1000000 ) )WRITE(*,*)'ERROR AT POSITION 3',TOTALVIKT

                NWALKTOTNEW = NWALKTOTNEW + NMERGED
                
                IF ( NMERGED .NE. 0 ) THEN
                   E(I) = E(I)/TOTALVIKT
                ELSE
                   E(I) = SUM(EC)/TOTALVIKT
                ENDIF
                
                NKILLW = 0
                xkilled = 0.0d0
                WKILLED = 0.0d0
                EKILLED = 0.0d0
                NREJECTKILLED = 0
                
                IF ( MOD(NWALKTOTKILLED,2) .NE. 0 ) THEN
                    ! We put the walkers weights and number of killed walkers that haven't been merged in thread id=0
                    IF ( id .EQ. 0 ) THEN
                         NKILLW = 1
                         xkilled(:,NKILLW) = xkilledc(:,NWALKTOTKILLED)
                         WKILLED(NKILLW) = WKILLEDC(NWALKTOTKILLED)
                         EKILLED(NKILLW) = EKILLEDC(NWALKTOTKILLED)
                         NREJECTKILLED(NKILLW) = NREJECTKILLEDC(NWALKTOTKILLED)
                    ENDIF
                ENDIF
                !------------------------------------------------------------
                ! Distributing the merged walkers o the computational threads
                !------------------------------------------------------------
                IF ( NMERGED .NE. 0 .AND. (ET-E(I))/abs(ET) .LT. CUTTOFFFACTOR  ) THEN
                    DO K=1,numprocessors
                        NWALK(K) = INT(NMERGED/numprocessors)
                    ENDDO

                    DO K=1,MOD(NMERGED,numprocessors)
                        NWALK(K) = NWALK(K) + 1
                    ENDDO
        
                    offset(1) = 0
                    DO K=2,numprocessors
                        offset(K) = NWALK(K-1)+offset(K-1)
                    ENDDO

                    ! Distributing the merged walkers/replicas
                    DO K=1,NWALK(id+1)
                       xnew(:,NWALKNEW(id+1)+K) = xmerged(:,offset(id+1)+K)
                       NEWWEIGHT(NWALKNEW(id+1)+K) = WEIGHTMERGED(offset(id+1)+K)
                       NREJECTUPPDATE(NWALKNEW(id+1)+K) =  NREJECTMERGED(offset(id+1)+K)
                    ENDDO
                ENDIF
                
                !----------------------------------------
                ! End of redistribution of merged walkers
                !----------------------------------------
                
                IF ( NWALKTOTNEW .EQ. 0 ) THEN
                        IF ( id .EQ. 0 ) THEN
                           WRITE(*,'(A84)')'===================================================================================='
                           WRITE(*,'(A84)')' RESTARTING FROM MOST RECENTLY SAVED POPULATION SINCE ALL WALKERS HAVE BEEN KILLED! '
                           WRITE(*,'(A84)')'===================================================================================='
                        ENDIF
                        ALLDEAD = .TRUE.
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        GOTO 30
                ENDIF

                TAUEFF = DT*(SUM(CNACC)*1.0d0/NWALKTOT*1.0d0)
                
                IF ( TAUEFF .EQ. 0.0d0 ) TAUEFF = DT

                IF ( MOD(I,NRECALC) .EQ. 0 ) THEN
                      ER(I) = ET - (1.0d0/(1.0d0*NRECALC*TAUEFF))*log(1.0d0*TOTALVIKT/(1.0d0*NREPLICAS))
                ELSE
                      IF ( TOTALVIKT .GT. 0.0d0 ) THEN
                          CORRECTION = beta*(1.0d0/(1.0d0*NRECALC*TAUEFF))*log(1.0d0*TOTALVIKT/(1.0d0*NREPLICAS))
                          ER(I) = SUM(ER)/(I-1) - CORRECTION
                      ELSE
                          ER(I)  = ER(I-1)
                      ENDIF 
                ENDIF
             30 CONTINUE
                ! Here if the configuration of walkers is pathological 
                ! we redo the time-step.
                IF ( abs(ET-E(I))/abs(ET) .LT. CUTTOFFFACTOR .AND. .not. ALLDEAD ) THEN
                        ET = (ET*(I-1) + E(I))/I
                        xold(:,:) = xnew(:,:)
                        NREJECT(:) = NREJECTUPPDATE(:)
                        NWALKTOT = NWALKTOTNEW
                        IF ( NMERGED .NE. 0 ) THEN
                               DO K=1,numprocessors
                                      COUNC(K) = COUNC(K) + NWALK(K)
                               ENDDO
                        ENDIF
                        NWALK(id+1) = COUNC(id+1)
                        WEIGHT = NEWWEIGHT
                        !----------------------------------------------
                        ! Saved if we want to recalculate the time step
                        !----------------------------------------------
                        ! saving the info about the killed walkers:
                        !----------------------------------------------
                        NKILLOLD = NKILLW
                        xkilledold = xkilled
                        WKILLEDOLD = WKILLED
                        EKILLEDOLD = EKILLED
                        NREJECTKILLEDOLD(:) = NREJECTUPPDATE(:)
                        !-----------------------------------------------------------------------------
                        ! If the population get stuck in a pathological configuration ( NSTUCK > 100 )
                        ! the wakers has to be restarted from a good configuration
                        !-----------------------------------------------------------------------------
                        IF ( NWALKTOT .GE. NREPLICAS .AND. NWALKTOT .LE. NPOPMIN .AND. NSTUCK .EQ. 0 .AND. .not. LOCKED ) THEN
                               ISAVED = I
                               xoldSAVED = xold
                               NREJECTSAVED = NREJECT
                               WEIGHTSAVED = WEIGHT
                               NWALKSAVED(id+1) = NWALK(id+1)
                               NWALKTOTSAVED = NWALKTOT
                               ! restoring the info about the killed walkers:
                               NKILLWSAVED = NKILLOLD
                               xkilledSAVED = xkilledold
                               WKILLEDSAVED = WKILLEDOLD
                               EKILLEDSAVED = EKILLEDOLD
                               NREJECTKILLEDSAVED = NREJECTKILLEDOLD
                               NPOPMIN = NWALKTOT
                               !IF ( NPOPMIN .EQ. NREPLICAS ) LOCKED = .TRUE.
                        ENDIF
                        NSTUCK = 0
                        NEXTSTEP = .TRUE.
                ELSE 
                        ! Since we are re-calculating the time step some variables
                        ! need to be restored:
                        NSTUCK = NSTUCK+1
                        
                        !-For debugging purposes:
                        !IF ( id .EQ. 0 ) THEN
                        !    WRITE(*,'(A9,I5,E30.20,F30.20,F30.20,I5)')'Current: ',NWALKTOTNEW,TOTALVIKT,E(I),ER(I),NSTUCK
                        !    WRITE(*,'(A9,I5,F30.20,F30.20,I5)')'  Saved: ',NWALKTOTSAVED,E(ISAVED),ER(ISAVED),NSTUCK
                        !ENDIF
                            
                        IF ( NSTUCK .LT. 100 .AND. .not. ALLDEAD ) THEN
                           ! restoring the info about the killed walkers:
                           NKILLW = NKILLOLD
                           xkilled = xkilledold
                           WKILLED = WKILLEDOLD
                           EKILLED = EKILLEDOLD
                           NREJECTKILLED = NREJECTKILLEDOLD
                        ELSE
                           xold = xoldSAVED
                           NREJECT = NREJECTSAVED
                           WEIGHT = WEIGHTSAVED
                           NWALK(id+1) = NWALKSAVED(id+1)
                           NWALKTOT = NWALKTOTSAVED
                           ! restoring the info about the killed walkers:
                           NKILLW = NKILLWSAVED
                           xkilled = xkilledSAVED
                           WKILLED = WKILLEDSAVED
                           EKILLED = EKILLEDSAVED
                           NREJECTKILLED = NREJECTKILLEDSAVED
                           ER(I-1) = ER(ISAVED)
                           E(I-1) = E(ISAVED)
                           ALLDEAD = .FALSE.
                           NSTUCK = 0
                        ENDIF
                        NEXTSTEP = .FALSE.
                        GOTO 10
                ENDIF
                
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
                !IF ( SUM(WEIGHT) .LE. 1.0d0 ) THEN 
                !   print*,I,'LIGHT AS A FEATHER!'
                !   STOP
                !ENDIF
                IF ( ( MINVAL(COUNC,1) .EQ. 1 .OR. MAXVAL(COUNC,1) .GE. 2*NRT .OR. DISTR ) .AND. .not. NOREDIST ) THEN
                   displsx(1) = 0
                   DO K=1,numprocessors
                      rcountsx(K) = N3*COUNC(K)
                      IF ( K .GE. 2 ) displsx(K) = displsx(K-1) + N3*COUNC(K-1)
                   ENDDO

                   displsreject(1) = 0
                   DO K=1,numprocessors
                      rcountsreject(K) = COUNC(K)
                      IF ( K .GE. 2 ) displsreject(K) = displsreject(K-1) + COUNC(K-1)
                   ENDDO
                   
                   IF ( id .EQ. 0 ) THEN
                      REDISTRIBUTED = .TRUE.
                      NREDIST = NREDIST + 1
                   ENDIF
                   
                   ! If the total number of walkers grow beyond wahat was initially allocated
                   ! then we have to reallcoate the size of the arrays
                   IF ( SUM(COUNC) .GT. NREP2 ) THEN 
                       NEWARRAYSIZE = SUM(COUNC)
                       DEALLOCATE(xstartc,NEWWEIGHTC,NREJECTNEWC)
                       ALLOCATE(xstartc(N3,NEWARRAYSIZE),NEWWEIGHTC(NEWARRAYSIZE),NREJECTNEWC(NEWARRAYSIZE))
                   ELSE
                       NEWARRAYSIZE = NREP2
                   ENDIF

                   !IF ( COUNC(id+1) .NE. 0 ) THEN
                     ! Collecting all the walker-coordinates in xnew at the different threads and put them into xstart located at thread 0 (id=0)
                     CALL MPI_GATHERV(xnew(:,1:COUNC(id+1)),rcountsx(id+1),MPI_DOUBLE_PRECISION,xstartc,rcountsx,displsx,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                     ! Collecting all the walker-weights in NEWWEIGHT at the different threads and put them into WEIGHT located at thread 0 (id=0)
                     CALL MPI_GATHERV(WEIGHT(1:COUNC(id+1)),rcountsreject(id+1),MPI_DOUBLE_PRECISION,NEWWEIGHTC,rcountsreject,displsreject,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
                     ! Collecting all the walker-rejaction numbers in NREJECT at the different threads and put them into NREJECTNEW located at thread 0 (id=0)
                     CALL MPI_GATHERV(NREJECT(1:COUNC(id+1)),rcountsreject(id+1),MPI_INTEGER,NREJECTNEWC,rcountsreject,displsreject,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
                   !ENDIF
                   
                   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr) 

                   ! Sending the array of walkers, xstart, located at thread 0 (id = 0 ) to all other threads.
                   CALL MPI_BCAST(xstartc,N3*NEWARRAYSIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)  
                   ! Sending the array of walker-weights, WEIGHT, located at thread 0 (id = 0 ) to all other threads.
                   CALL MPI_BCAST(NEWWEIGHTC,NEWARRAYSIZE, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   ! Sending the array of walker-rejection numbers, NREJECTNEW, located at thread 0 (id = 0 ) to all other threads.
                   CALL MPI_BCAST(NREJECTNEWC,NEWARRAYSIZE, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
                   
                   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                   ! The new number of replicas/processor:
                   DO K=1,numprocessors
                      NWALK(K) = INT(NWALKTOTNEW/numprocessors)
                   ENDDO

                   DO K=1,MOD(NWALKTOTNEW,numprocessors)
                      NWALK(K) = NWALK(K) + 1
                   ENDDO
        
                   offset(1) = 0
                   DO K=2,numprocessors
                      offset(K) = NWALK(K-1)+offset(K-1)
                   ENDDO
                   
                   ! Re-distributing the walkers/replicas, the rejection numbers and the walker-weights on the different threads:
                   xold(:,:) = 0.0d0
                   DO K=1,NWALK(id+1)
                       xold(:,K) = xstartc(:,offset(id+1)+K)
                       NREJECT(K) = NREJECTNEWC(offset(id+1)+K)
                       WEIGHT(K) =  NEWWEIGHTC(offset(id+1)+K)
                   ENDDO
                   IF ( id .EQ. 0 .AND. DISTR .AND. NSTUCK .EQ. 0 ) THEN 
                                4 FORMAT(3000(F30.20))
                                OPEN(15,FILE='DQMCRESTART.dat',ACTION='WRITE')
                                WRITE(15,*)I,NWALKTOTNEW
                                WRITE(15,'(F30.20,F30.20,F30.20,F30.20)')E(I),ER(I),ET,EVMC
                                DO K=1,NWALKTOTNEW
                                        WRITE(15,FMT=4)xstartc(:,K)
                                ENDDO
                                DO K=1,NWALKTOTNEW
                                        WRITE(15,'(I20,F30.20)')NREJECTNEWC(K),NEWWEIGHTC(K)
                                ENDDO
                                close(15)
                   ENDIF

               ELSE
                  !----------------------------------------------------------------
                  !                     END OF LOAD BALANCING
                  !----------------------------------------------------------------
                  xold(:,:) = xnew(:,:)
              ENDIF
          10 CONTINUE
              CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              IF ( .not. NEXTSTEP ) GOTO 20
        ENDDO
        IF ( id .EQ. 0 ) THEN
           !-----------------------------
           ! Calculating the total energy
           !-----------------------------
        
           NSTART = INT(TSTART/TIMESTEP)
           NEND = INT(TEND/TIMESTEP)
        
           IF (NSTART .EQ. 0 ) NSTART = 1
        
           EDQMC = 0.0d0
           II = 0

           ALLOCATE(EG(NTIMESTEPS),EM(NTIMESTEPS))
     
           EG(:) = 0.0d0
           EM(:) = 0.0d0

           DO I=NSTART,NEND
               EDQMC = EDQMC + E(I) 
               II = II + 1
               EM(II) = EDQMC/II
               IF ( II .GT. 1 ) THEN 
                       EG(II) = ( EG(II-1)*(II-1) +EM(II) )/II
               ELSE
                       EG(II) = EM(II)
               ENDIF
           ENDDO

           EDQMC = EDQMC/II
        
           !-----------------------------------
           ! Calculating the standard deviation
           ! of the total energy
           !-----------------------------------
        
           ESIGMA = 0.0d0
           ESIGMAM = 0.0d0
           DO I=1,II
                 ESIGMA  =  ESIGMA + (EM(I)-EDQMC)**2
                 ESIGMAM =  ESIGMAM + (EM(I)-EG(II))**2
           ENDDO

           ESIGMA = sqrt(ESIGMA/(1.0d0*(NEND-NSTART+1)))
           ESIGMAM = sqrt(ESIGMAM/(1.0d0*(NEND-NSTART+1)))

           WRITE(*,*)'======================================================================='
           WRITE(*,*)'        TOTAL ENERGY OF DIFFUSION QUANTUM MONTECARLO CALCULATION  '
           WRITE(*,*)'======================================================================='
           WRITE(*,*)' '
           WRITE(*,'(A14,F30.20,A5,F23.20,A5)')' <E(LOCAL)>=  ',EDQMC+nucE,' +/- ',ESIGMA,' au'
           WRITE(*,*)' '
           WRITE(*,'(A14,F30.20,A5,F23.20,A5)')'<<E(LOCAL)>>= ',EG(II)+nucE,' +/- ',ESIGMAM,' au'
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
           DO I=1,II
             WRITE(13,'(I20,F30.20,F30.20,F30.20,F30.20)')NSTART+I-1,EM(I)+nucE,EG(I)+nucE,ER(I+NSTART-1)+nucE,E(I+NSTART-1)+nucE
           ENDDO
           CLOSE(12)
           CLOSE(13)
       ENDIF
       IF ( NTOOLD .NE. 0 ) THEN
          WRITE(*,*)' '
          WRITE(*,'(A3,I6,A46,I4,A11)')'At ',NTOOLD,' TIMES there were walkers stuck for more than ' ,NPERSIST, 'generations'
       ENDIF
       ! Accumumulating the charge-density:
       IF (WRITEDENS ) THEN
                IF ( id .NE. 0 ) THEN
                        call MPI_SEND( density, MESH(1)*MESH(2)*MESH(3), MPI_DOUBLE_PRECISION, 0, id, MPI_COMM_WORLD, ierr )
                ELSE
                        DO I=1,numprocessors-1
                                CALL MPI_RECV(densityc, MESH(1)*MESH(2)*MESH(3), MPI_DOUBLE_PRECISION, I, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
                                density = density + densityc
                        ENDDO
                ENDIF
                IF ( id .EQ. 0 ) THEN
                        ! RENORMALIZING THE CHARGEDENSITY
                        DELTAVOL = 8.0d0*(LIMITS(1)/(1.0d0*(MESH(1)-1)))*(LIMITS(2)/(1.0d0*(MESH(2)-1)))*(LIMITS(3)/(1.0d0*(MESH(3)-1)))
                        density = density*N/( DELTAVOL*SUM(reshape(density,(/MESH(1)*MESH(2)*MESH(3)/))))
                        WRITE(*,*)'================================================================'
                        WRITE(*,*)'Charge density has been calculated and is being  saved to file. '
                        WRITE(*,*)'================================================================'
                        3 FORMAT(I6,I6,I6,F30.20)
                        OPEN(9,FILE='CHARGEDENS.dat',ACTION='WRITE')
                        DO I=1,MESH(1)
                                DO J=1,MESH(2)
                                        DO K=1,MESH(3)
                                                WRITE(9,FMT=3)I,J,K,density(I,J,K)
                                        ENDDO
                                ENDDO
                        ENDDO
                ENDIF
       ENDIF
CONTAINS
  FUNCTION random_normall() RESULT(fn_val)

  ! Adapted from the following Fortran 77 code
  !      ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM.
  !      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
  !      VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435.

  !  The function random_normal() returns a normally distributed pseudo-random
  !  number with zero mean and unit variance.

  !  The algorithm uses the ratio of uniforms method of A.J. Kinderman
  !  and J.F. Monahan augmented with quadratic bounding curves.

  REAL :: fn_val
  !     Local variables
  REAL :: sss = 0.449871, ttt = -0.386595, aaa = 0.19600, bbb = 0.25472, rrr1 = 0.27597, rrr2 = 0.27846, uuu, vvv, xxx, yyy, qqq

  !     Generate P = (u,v) uniform in rectangle enclosing acceptance region

  DO
    CALL RANDOM_NUMBER(uuu)
    CALL RANDOM_NUMBER(vvv)
  IF ( DEBUGGRANDOM )  print*,id,uuu
   vvv = 1.7156 * (vvv - 0.5 )

  !     Evaluate the quadratic form
    xxx = uuu - sss
    yyy = ABS(vvv) - ttt
    qqq = xxx**2 + yyy*(aaa*yyy - bbb*xxx)

  !     Accept P if inside inner ellipse
    IF (qqq < rrr1) EXIT
  !     Reject P if outside outer ellipse
    IF (qqq > rrr2) CYCLE
  !     Reject P if outside acceptance region
    IF (vvv**2 < -4.0*LOG(uuu)*uuu**2) EXIT
  END DO

  !     Return ratio of P's coordinates as the normal deviate
  fn_val = vvv/uuu
  RETURN
  END FUNCTION random_normall

END SUBROUTINE dqmc
