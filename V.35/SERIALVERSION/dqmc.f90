SUBROUTINE dqmc(N,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,beta,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,EHF,nucE,ESIGMA,EDQMC,NPERSIST,NRECALC,CUTTOFFFACTOR,NVMC,VMCCALC, &
                & WRITEDENS,MESH,LIMITS)
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
        INTEGER, INTENT(IN)  :: N,NATOMS,NRECALC,NVMC
        TYPE(BASIS), INTENT(IN) :: BAS
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS),LIMITS(3)
        DOUBLE PRECISION, INTENT(IN) :: TIMESTEP,TEND,TSTART,beta,EHF,nucE,BJASTROW,CJASTROW,CUTTOFFFACTOR
        INTEGER, INTENT(IN) :: SAMPLERATE,NREPLICAS,NPERSIST,MESH(3)
        LOGICAL, INTENT(IN) :: VMCCALC,WRITEDENS
        DOUBLE PRECISION, INTENT(OUT) :: EDQMC, ESIGMA
        DOUBLE PRECISION, EXTERNAL :: trialfnk,EL,laplacetrialfnk,laplacehforbitalval,hforbitalval
        EXTERNAL :: guideforce,findclosestatom,preparedifusion
        DOUBLE PRECISION :: SIGMA,VECT1(3),VECT2(3),ESIGMAM,fmean,ffmean,fpfmean,CCORR,CORRECTION
        DOUBLE PRECISION :: EGUESS,DT,V,ET,ENT1,ENT2,PSIT1,PSIT2,EVMC,PSI1,PSI2
        DOUBLE PRECISION :: P1,P2,b,c,TEST,TAUEFF,ERMEAN,DELTAVOL
        INTEGER :: I,J,K,NTIMESTEPS,COUN,NUMBEROFNEWWALKERS,N3,NREP2,NSTART,NEND,M,clock,II,NACC,NTOOLD,NWALK,NWALKNEW,NWALKSAVED,NMETROPOLIS
        INTEGER, ALLOCATABLE :: seed(:),NREJECT(:),NREJECTUPPDATE(:),NREJECTSAVED(:)
        DOUBLE PRECISION, ALLOCATABLE :: ER(:),E(:),force(:),xt(:),xtrial(:),xtemp(:,:),xold(:,:),xkilledOLD(:),xkilledSAVED(:)
        DOUBLE PRECISION, ALLOCATABLE :: xnew(:,:),xsaved(:,:),x1(:),x2(:),EG(:),EM(:),x3(:),x4(:),ECORR(:),density(:,:,:)
        REAL :: RAND,ratiov,testa,ita,s
        INTEGER :: xcoord,ycoord,zcoord,JJ,ISAVED,NPOPMIN,NSTUCK
        LOGICAL :: first,ALLDEAD,LOCKED,NONPATHOLOGICAL
        DOUBLE PRECISION :: rvect(3),vvect(3),vlength,zvect(3),vect(3),distance,rnuc(3),z,zetaa,aa,vmeanvect(3),vmeanz,vmeanrho,rhovect(3)
        DOUBLE PRECISION :: zbis, rhobis,dr(3),rprim(3),RANDDIR(3),ppobability,SR,SRP,rvectp(3),G1,G2,WTOT,EKILLED,WKILLED,ratio
        DOUBLE PRECISION, ALLOCATABLE :: VBIG(:), VBIGMEAN(:),VBIGP(:), VBIGMEANP(:),WEIGHT(:),NEWWEIGHT(:),WEIGHTSAVED(:),xkilled(:)
        INTEGER :: KK,atomnr,NKILLW,NKILLWOLD,NKILLWSAVED
        DOUBLE PRECISION :: EKILLEDOLD,WKILLEDOLD,EKILLEDSAVED,WKILLEDSAVED
        REAL :: qtilde,qtildep
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        ! Allocating some variables:
        N3 = 3*N
        NREP2 = 2*NREPLICAS

        NSTUCK = 0
        NPOPMIN = 2*NREPLICAS
        ALLDEAD = .FALSE.
        LOCKED = .FALSE.

        ALLOCATE(force(N3),xt(N3),xtrial(N3),xtemp(NREP2,N3),xold(NREP2,N3),xkilled(N3),xkilledOLD(N3),xkilledSAVED(N3))
        ALLOCATE(xnew(NREP2,N3),xsaved(NREP2,N3),x1(N3),x2(N3),x3(N3),x4(N3),NREJECT(NREP2),NREJECTUPPDATE(NREP2),NREJECTSAVED(NREP2),WEIGHT(NREP2),NEWWEIGHT(NREP2))
        ALLOCATE(VBIG(N3),VBIGMEAN(N3),VBIGP(N3),VBIGMEANP(N3),WEIGHTSAVED(NREP2))

        IF ( WRITEDENS ) THEN
                ALLOCATE(density(MESH(1),MESH(2),MESH(3)))
                density = 0.0d0
        ENDIF
        first = .TRUE.

        b = BJASTROW    ! When b = 0, the below Metropolis generated distribution should give the
                        ! Hartree Fock energy, when monte carlo integration is performed, see ER(1)
        c = CJASTROW
       
        !------------------------------------
        ! Seeding the random number generator
        !------------------------------------
        call random_seed(size=M)
        ALLOCATE(seed(M))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * (/ (i - 1, i = 1, M) /)
        CALL RANDOM_SEED(PUT = seed )
        !------------------------------------------
        ! End of seeding.
        !------------------------------------------
        
        NTIMESTEPS = INT(TEND/TIMESTEP)
        NREJECT(:) = 0
        NREJECTUPPDATE(:) = 0 
        NTOOLD = 0
        ALLOCATE(ER(NTIMESTEPS),E(NTIMESTEPS),EG(NTIMESTEPS),EM(NTIMESTEPS))
        DT = TIMESTEP

        SIGMA = sqrt(DT)
        TAUEFF = DT

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
        KK = 0
        ER(:) = 0.0d0
        EG(:) = 0.0d0
        ERMEAN = 0.0d0
        WEIGHT = 0.0d0

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
                PSI1 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,xtrial)
                PSI2 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,xt)
                ratiov = REAL((PSI1/PSI2)**2)
               
                CALL RANDOM_NUMBER(testa)
                !------------------------
                ! ACCEPT OR REJECT MOVE ?
                !------------------------
                IF ( ratiov .GE. 1.0 .OR. testa .LT. ratiov ) THEN
                        xt = xtrial
                        PSI1 = PSI2
                ENDIF

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
                        !--------------------------
                        ! End of density collection
                        !--------------------------
                        IF ( NMETROPOLIS - I .LE. (NREPLICAS+1)*SAMPLERATE ) THEN
                                KK = KK + 1
                                IF ( KK .LE. NREPLICAS ) THEN
                                        xold(KK,:) = xt
                                        WEIGHT(KK) = 1.0d0
                                ENDIF
                        ENDIF
                ENDIF
                
        ENDDO
        
        !================================
        ! HERE THE VMC CALCULATION ENDS
        !================================
        EVMC = EVMC/K
        IF ( VMCCALC ) THEN
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
                        density = density*N/( DELTAVOL*SUM(reshape(density,(/MESH(1)*MESH(2)*MESH(3)/))) )
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
                        CLOSE(9)
                ENDIF
                STOP
        ELSE
                WRITE(*,'(A14,F30.20,A5)')'E(VMC) =',EVMC+nucE,'au'
        ENDIF

        ! The starting guess:
        ER(1) = EVMC
        IF ( ER(1) .GT. EHF ) ER(1)  = EHF
        ERMEAN = ER(1)
        EG(1) = ER(1)
        II = 1

        E(1) = ER(1)
        ET = ER(1)
        NWALK = NREPLICAS
        
        OPEN(12,FILE='ENERGYDQMC.dat',ACTION='WRITE')
        !---------------------------------------------------
        ! Here the Quantum Diffusion Monte Carlo loop starts
        !---------------------------------------------------
        print*,'   =========================================================='
        print*,'          Entering the Diffusion Monte Carlo loop            '
        print*,'   =========================================================='
        print*,' '
        !WRITE(*,'(A9,A24,A24,A27)')'N','<E> [au]','E [au]','NREPLICAS'
        NKILLW = 0
        DO I=2,NTIMESTEPS
                WRITE(12,'(I20,F30.20,F30.20,F30.20,I20)')I,ET+nucE,ER(I-1)+nucE,E(I-1)+nucE,NWALK
            10  CONTINUE   
                COUN = 0
                WTOT = 0.0d0
                E(I) = 0.0d0
                NACC = 0
                DO J=1,NWALK
                        !--------------------------------------------------------------
                        ! All the numbers in parenthesis to the right are in  reference 
                        ! to the numbering of the steps in Umrigars algorithm, 
                        ! of the so called improved algorithm on pages 2886-2887 in
                        ! J. Chem. Phys. 99, 2865-2890, (1993)
                        !--------------------------------------------------------------
                        !===========================================
                        ! The the forward diffusion step/trial move
                        !===========================================
                        x3 = xold(J,:)
                        PSI1 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,x3)
                        PSIT1 = PSI1
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
                                        rprim(1) = dr(1) + SIGMA*random_normal()
                                        rprim(2) = dr(2) + SIGMA*random_normal()
                                        rprim(3) = dr(3) + SIGMA*random_normal()
                                ELSE
                                        ! Calculating the other of the two possible forward trial positions r' (20-22)
                                        RANDDIR(1) = random_normal()
                                        RANDDIR(2) = random_normal()
                                        RANDDIR(3) = random_normal()
                                        rprim = rnuc + ( RANDDIR/sqrt(DOT_PRODUCT(RANDDIR,RANDDIR)))*random_chisq(6,first)/(4.0d0*zetaa)
                                        IF ( first ) first = .FALSE.
                                ENDIF
                                G1 = ((2*PI*DT)**(-1.50d0))*EXP(-(1.0d0/(2.0d0*DT))*DOT_PRODUCT(rprim-dr,rprim-dr)) ! (27)
                                G2 = (zetaa**3/PI)*EXP(-2.0d0*zetaa*sqrt(DOT_PRODUCT(rprim-rnuc,rprim-rnuc))) ! (27)

                                P1 = P1*( G1*(1-qtilde) + G2*qtilde )

                                ! This is where the new R' position is saved
                                xtemp(J,3*(K-1) + 1) = rprim(1)
                                xtemp(J,3*(K-1) + 2) = rprim(2)
                                xtemp(J,3*(K-1) + 3) = rprim(3)
                        
                                ! This is where the mean drift-velocity <V> is saved:
                                VBIGMEAN(3*(K-1) + 1) = vmeanvect(1)
                                VBIGMEAN(3*(K-1) + 2) = vmeanvect(2)
                                VBIGMEAN(3*(K-1) + 3) = vmeanvect(3)
                        ENDDO
                        
                        !=============================
                        ! Start of the reversed move:
                        !=============================
                        x4 = xtemp(J,:)
                        PSI2 = trialfnk(N,N3,Cup,Cdown,BAS,b,c,x4)
                        PSIT2 = PSI2
                        ! This is basically where the fixed node approximation comes in (29)
                        IF ( PSIT1*PSIT2 .LE. 0.0d0 ) THEN
                                ppobability = 0.0
                                NREJECT(J) = NREJECT(J) + 1
                                xtemp(J,:) = xold(J,:)   
                                GOTO 20
                        ENDIF 
                        
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
                                xtemp(J,:) = xold(J,:)
                                NREJECT(J) = NREJECT(J) + 1
                                ppobability = ratio
                        ENDIF
                        
                        IF ( ita .LE. ratio ) THEN 
                                NACC = NACC + 1
                                NREJECT(J) = 0
                                ppobability = ratio
                                PSI1 = PSI2
                        ENDIF
                        !---------------------------------------------
                        ! Collecting data for the density calculation:
                        !---------------------------------------------
                        IF ( WRITEDENS ) THEN
                                DO JJ=1,N
                                        xcoord =  NINT((xtemp(J,3*(JJ-1) + 1)+LIMITS(1))*(MESH(1)-1)/(2.0d0*LIMITS(1)))
                                        ycoord =  NINT((xtemp(J,3*(JJ-1) + 2)+LIMITS(2))*(MESH(2)-1)/(2.0d0*LIMITS(2)))
                                        zcoord =  NINT((xtemp(J,3*(JJ-1) + 3)+LIMITS(3))*(MESH(3)-1)/(2.0d0*LIMITS(3)))
                                ENDDO
                                IF ( xcoord .GE. 1 .AND. xcoord .LE. MESH(1) .AND. ycoord .GE. 1 .AND. ycoord .LE. MESH(2) .AND. zcoord .GE. 1 .AND. zcoord .LE. MESH(3) ) THEN
                                        density(xcoord,ycoord,zcoord) = density(xcoord,ycoord,zcoord) + ppobability
                                ENDIF
                        ENDIF
                        !--------------------------
                        ! End of density collection
                        !--------------------------
                    20  CONTINUE
                        
                        !----------------------------------------------------------------------------------------
                        ! Birth/death step if there is a new step accepted by  the Metropolis algorithm
                        ! we check if there is going to a birth or death (merge) of particles at the new position.
                        !-----------------------------------------------------------------------------------------
                        x3 = xold(J,:)
                        x4 = xtemp(J,:)
                        ENT1 = EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,x3) 
                        ENT2 = EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,x4)
                        
                        SRP = ER(I-1) - ET + ( ET - ENT2 )*sqrt(DOT_PRODUCT(VBIGMEANP,VBIGMEANP)/DOT_PRODUCT(VBIGP,VBIGP))      ! (36)
                        SR =  ER(I-1) - ET + ( ET - ENT1 )*sqrt(DOT_PRODUCT(VBIGMEAN,VBIGMEAN)/DOT_PRODUCT(VBIG,VBIG))          ! (36)
                        
                        V = 0.50d0*(SRP + SR)*ppobability + (1-ppobability)*SR

                        WEIGHT(J) = WEIGHT(J)*exp(TAUEFF*V)   ! (38)

                        NUMBEROFNEWWALKERS = NINT(WEIGHT(J))
                        NONPATHOLOGICAL = .TRUE.

                        IF ( WEIGHT(J) .GT. 1.0d0*NREPLICAS ) THEN
                                NUMBEROFNEWWALKERS = 0
                                NONPATHOLOGICAL = .FALSE.
                                WEIGHT(J)  = 0.0d0
                        ENDIF
                        
                        ! Here new walkers are born ( See, Chem. Phys. 99, 2865 (1993), page 2867 column 2 )
                        IF ( NUMBEROFNEWWALKERS .GE. 1  ) THEN
                                IF ( NUMBEROFNEWWALKERS .GT. 2 ) NUMBEROFNEWWALKERS = 2
                                DO K=1,NUMBEROFNEWWALKERS
                                        IF ( COUN < 2*NREPLICAS ) THEN
                                                COUN = COUN + 1
                                                xnew(COUN,:) = x4
                                                NEWWEIGHT(COUN) = WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                                WTOT = WTOT + WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                                NREJECTUPPDATE(COUN) = NREJECT(J)
                                                 E(I) = E(I) + ( ENT2*ppobability + ENT1*(1-ppobability) )*WEIGHT(J)/(1.0d0*NUMBEROFNEWWALKERS)
                                        ENDIF
                                ENDDO
                        ELSE
                                IF ( NONPATHOLOGICAL ) THEN
                                        ! Here 2 rejected walkers are merged into one (  See, Chem. Phys. 99, 2865 (1993), page 2867, column 2 )
                                        NKILLW = NKILLW + 1
                                        IF ( NKILLW .EQ. 1 ) THEN
                                                EKILLED = ( ENT2*ppobability + ENT1*(1-ppobability) )
                                                WKILLED = WEIGHT(J)
                                                xkilled = x4
                                        ELSE
                                                WKILLED = WKILLED + WEIGHT(J)
                                                CALL RANDOM_NUMBER(RAND)
                                                IF ( RAND .LE. WEIGHT(J)/WKILLED ) THEN
                                                        xkilled = x4
                                                        EKILLED = ( ENT2*ppobability + ENT1*(1-ppobability) )
                                                ENDIF
                                                IF ( COUN < 2*NREPLICAS   ) THEN
                                                        COUN = COUN + 1
                                                        WTOT = WTOT + WKILLED
                                                        xnew(COUN,:)  = xkilled
                                                        NEWWEIGHT(COUN) = WKILLED
                                                        E(I) = E(I) + EKILLED*WKILLED
                                                        NKILLW = 0
                                                ENDIF
                                                NKILLW = 0
                                        ENDIF
                                ENDIF
                        ENDIF
                ENDDO
                IF ( COUN .EQ. 0 ) THEN
                        WRITE(*,'(A84)')'===================================================================================='
                        WRITE(*,'(A84)')' RESTARTING FROM MOST RECENTLY SAVED POPULATION SINCE ALL WALKERS HAVE BEEN KILLED! '
                        WRITE(*,'(A84)')'===================================================================================='
                        ALLDEAD = .TRUE.
                        GOTO 30
                ENDIF
                NWALKNEW = COUN
                E(I) = E(I)/WTOT

                TAUEFF = DT*(NACC*1.0d0/NWALK*1.0d0)
                
                IF ( TAUEFF .EQ. 0.0d0 ) TAUEFF = DT

                IF ( MOD(I,NRECALC) .EQ. 0  ) THEN
                        ER(I) = ET - (1.0d0/(1.0d0*NRECALC*TAUEFF))*log(1.0d0*WTOT/(1.0d0*NREPLICAS))
                ELSE
                      CORRECTION = beta*(1.0d0/(1.0d0*NRECALC*TAUEFF))*log(1.0d0*WTOT/(1.0d0*NREPLICAS))
                      ER(I) = SUM(ER)/(I-1) - CORRECTION
                ENDIF
             30 CONTINUE
                ! Here if the configuration of walkers is pathological 
                ! we redo the time-step.
                IF ( abs((ET-E(I))/ET) .LT. CUTTOFFFACTOR .AND. .not. ALLDEAD ) THEN
                        ET = (ET*(I-1) + E(I))/I
                        IF ( NWALKNEW .GE. NREPLICAS .AND. NWALKNEW .LE. NPOPMIN .AND. NSTUCK .EQ. 0 .AND. .not. LOCKED ) THEN
                                xsaved = xnew
                                ISAVED = I
                                NREJECTSAVED = NREJECTUPPDATE
                                NWALKSAVED = NWALKNEW
                                WEIGHTSAVED = NEWWEIGHT
                                NPOPMIN = NWALKNEW
                                WKILLEDSAVED = WKILLED
                                EKILLEDSAVED = EKILLED
                                xkilledSAVED = xkilled
                                NKILLWSAVED = NKILLW
                                !IF ( NPOPMIN .EQ. NREPLICAS ) LOCKED = .TRUE.
                        ENDIF
                        xold = xnew
                        NREJECT = NREJECTUPPDATE
                        NWALK = NWALKNEW
                        WEIGHT = NEWWEIGHT
                        ! saving the info about the killed walkers that haven't
                        ! gotten merged back into the population
                        WKILLEDOLD = WKILLED
                        EKILLEDOLD = EKILLED
                        xkilledOLD = xkilled
                        NKILLWOLD = NKILLW
                        NSTUCK = 0
                ELSE
                        ! restoring the info about the killed walkers that
                        ! gotten merged back into the population
                        ER(I) = 0.0d0
                        WKILLED = WKILLEDOLD
                        EKILLED = EKILLEDOLD
                        xkilled = xkilledOLD
                        NKILLW = NKILLWOLD
                        NSTUCK = NSTUCK + 1
                        IF ( NSTUCK .GT. 100 .OR. ALLDEAD ) THEN
                                xold = xsaved
                                NREJECT = NREJECTSAVED
                                NWALK = NWALKSAVED
                                WEIGHT = WEIGHTSAVED
                                ER(I-1) = ER(ISAVED)
                                E(I-1) = E(ISAVED)
                                ! restoring the info about the killed walkers that
                                ! gotten merged back into the population
                                WKILLED = WKILLEDSAVED
                                EKILLED = EKILLEDSAVED
                                xkilled = xkilledSAVED
                                NKILLW = NKILLWSAVED
                                ALLDEAD = .FALSE.
                                NSTUCK = 0
                        ENDIF
                        GOTO 10 
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
        
        OPEN(13,FILE='DQMCRESULTS.dat',ACTION='WRITE')
        DO I=1,II
            WRITE(13,'(I20,F30.20,F30.20,F30.20,F30.20)')NSTART+I-1,EM(I)+nucE,EG(I)+nucE,ER(I+NSTART-1)+nucE,E(I+NSTART-1)+nucE
        ENDDO
        CLOSE(13)
        CLOSE(12)

        IF ( NTOOLD .NE. 0 ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A3,I6,A60)')'At ',NTOOLD,' TIMES there were walkers stuck for more than 50 generations'
        ENDIF
        IF ( WRITEDENS ) THEN
                ! RENORMALIZING THE CHARGEDENSITY
                DELTAVOL = 8.0d0*(LIMITS(1)/(1.0d0*(MESH(1)-1)))*(LIMITS(2)/(1.0d0*(MESH(2)-1)))*(LIMITS(3)/(1.0d0*(MESH(3)-1)))
                density = density*N/( DELTAVOL*SUM(reshape(density,(/MESH(1)*MESH(2)*MESH(3)/))) )
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
                CLOSE(9)
        ENDIF
END SUBROUTINE dqmc
