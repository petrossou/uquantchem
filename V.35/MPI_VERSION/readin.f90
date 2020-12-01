SUBROUTINE readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBER,SPATIALPOSITION,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, & 
& SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,REDISTRIBUTIONFREQ,NOREDIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,CUSPCORR,rc,CORRALCUSP, &
& NVMC,HFORBWRITE,IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol,NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS,&
& EETOL,RELALGO,SOFTSTART,HUCKEL,ZEROSCFTYPE,ETEMP,IORBNR,AORBS,DIISORDEX,DOTDFT,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF,OPTH,MIXEX,DOABSSPECTRUM, &
& DIFFDENS,NSCCORR,MIXTDDFT,SCERR,RIAPPROX,LIMPRECALC,DIAGDG,FIELDDIR,FIELDREAD,DAMPING,SCALARRELC,URHFTORHF)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS,NLINES
      CHARACTER(LEN=20), INTENT(OUT) :: CORRLEVEL,EPROFILE
      INTEGER, INTENT(OUT) :: ATOMICNUMBER(NATOMS),Ne,MESH(3),NEEXC,SAMPLERATE,NREPLICAS,REDISTRIBUTIONFREQ,NRECALC,NPERSIST,DIISORD,DIISSTART,NVMC,IOSA,NSTEPS,NLSPOINTS,PORDER,PULAY,FIXNSCF,DORDER
      INTEGER, INTENT(OUT) :: NLEBEDEV,NCHEBGAUSS,RELALGO,ZEROSCFTYPE,IORBNR(2),AORBS,DIISORDEX,EDIR,NEPERIOD,NSCCORR,LIMPRECALC
      DOUBLE PRECISION,INTENT(OUT) :: SPATIALPOSITION(NATOMS,3),Tol,LIMITS(3),ENEXCM,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,CUTTOFFFACTOR,MIX,rc,DR,FTol,TEMPERATURE,kappa,alpha,EETOL,ETEMP
      DOUBLE PRECISION,INTENT(OUT) :: OMEGA,EFIELDMAX,MIXEX,MIXTDDFT,SCERR,FIELDDIR(3),DAMPING
      LOGICAL, INTENT(OUT) :: WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,NOREDIST,CUSPCORR,CORRALCUSP,HFORBWRITE,CFORCE,RELAXN,WRITEONFLY,MOVIE,MOLDYN,ZEROSCF,XLBOMD
      LOGICAL, INTENT(OUT) :: SOFTSTART,HUCKEL,DOTDFT,ADEF,OPTH,DOABSSPECTRUM,DIFFDENS,RIAPPROX,DIAGDG,FIELDREAD,SCALARRELC,URHFTORHF
      INTEGER :: I,J,TOTALNLINES
      LOGICAL :: finns,LINEHASBEANREAD,ABSSPECTSET
      CHARACTER(LEN=20), ALLOCATABLE :: VARNAME(:)
      CHARACTER(LEN=2) :: DUMMY
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
    
      FIELDREAD = .FALSE.
     
      inquire(file='INPUTFILE',exist=finns)
      IF ( finns ) THEN
        OPEN(10,FILE='INPUTFILE',STATUS='OLD',ACTION='READ')
        TOTALNLINES = NATOMS+NLINES
        ALLOCATE(VARNAME(TOTALNLINES))
        I = 1
        DO I=1,TOTALNLINES
                READ(10,*)VARNAME(I),DUMMY
        ENDDO
        
        REWIND(10)

        ABSSPECTSET = .FALSE.


        ! Here we put some default values:
        URHFTORHF = .FALSE. ! If true the density matrices of an open shell calculation (URHF) Pu and Pd are transformed to be equal
                            ! by the following transformation Pu* = 0.5*(Pu + Pd) , Pd* = Pu*. This is for the moment only
                            ! implemented for URHF calculations to make the convergence of open shell systems more stable. TODO:
                            ! Extend this to DFT calculations.
        SCALARRELC = .FALSE.     ! If true the total energy Hartree-Fock and DFT calculations are done with scalar relatevistic
                                ! correctons, i.e Darwin correction to the nuclear and Hartree potentials including the p**4 mass
                                ! correction term. TODO: Calculate the nuclear forces with this correction as well.
        
        DAMPING = 0.0d0                         ! When doing molecular dynamics calculations one can introduce (fenomenological) dynamical friction,
                                                ! where the friction forces are proportonal to the velocities of the atomic nuclea. The proportionality 
                                                ! constant is given by the parameter DAMPING


        FIELDDIR = (/ 0.0d0, 0.0d0, 1.0d0 /)   ! The direction of the electric field in a TDFT or adiabatic electric field (ADEF = .T.) 
                                               ! calculation. If this derection  is set in the INPUTFILE the  EDIR directive is ignored 
                                               ! and the direction of travel of the electromagnetic field is set to ( cos(TH)*sin(FI),  cos(TH)*cos(FI), - sin(TH) )
                                               ! where TH = arcos(FIELDDIR(3))  and  FI =  arctan(FIELDDIR(1)/FIELDDIR(2))
        DIAGDG = .FALSE.        ! creates a diagonal density matrix and uses it as a startinguess for the SCF calculations as long 
                                ! as the file DENSMATSTARTGUESS.dat is not in  the run-directory or as long as the number of  atoms = 1


        RIAPPROX = .FALSE.      ! If true and if the file BASISFILEAUX is present in the run-directory the Resolution of the Identity approximation 
                                ! will be empoyed.
        LIMPRECALC = 100        ! (ONLY USED WITH RIAPPROX = .TRUE. ) Basis-set size limmit for the pre-calculation of the electron-electron integrals (ij|kl) and their derivatives 
                                ! with respect to nuclear cordinates. If the number of basis, NB, is less or equal to LIMPRECALC then all these integrals
                                ! will be constructed and put into the memory in  the beginning of the calculation. For RIAPPROX = .FALSE, the integrals are
                                ! always pre-calculated

        NSCCORR = 0             ! If NSCCORR > 0, then NSCCORR is the number of iterations that is used at each time-step in a TDDFT/THF calculation
                                ! to solve the equation: P(t+Dt) = exp(0.5*i*Dt*(F[P(t)]+F[P(t+Dt)]))*P(t)*exp(-0.5*i*Dt*(F[P(t)]+F[P(t+Dt)]))*P(t),

        MIXTDDFT = 0.50d0       ! Mixing parameter used in the iterative solution of the above equation
        SCERR = 0.1E-10         ! Convergence criterion for the difference between two consecutive P(t+Dt) obtained from the iterations of the above equation.

        Tol = 1.0E-6
        HUCKEL = .FALSE.
        WRITECICOEF = .FALSE. 
        WRITEDENS   = .FALSE.
        WHOMOLUMO   = .FALSE.
        DIFFDENS  = .FALSE.     ! If true and WRITEDENS = .T. and DOTDFT = .T.  the change in density instead of the density will be
                                ! saved. Decreas of orbital occupation (compared to t=0 case ) will be saved as positive densities
                                ! and increas of orbitals occupation of orbitals that are unuccupied at t= 0 will be saved as  negative contribution to density.
        DOABSSPECTRUM = .FALSE. ! IF DOABSSPECTRUM = .TRUE. the absorbtion spectrum is calculated. DOABSSPECTRUM = .TRUE.
                                ! automatically if (TEND/TIMESTEP) .LT. 10000
        MIXEX = 0.0d0         ! Mixing parameter used in Orthogonality Constrained DFT calculations of Exited States in the case of 
                              ! an MD-calculation. MIXEX = amount of the excited state density matrix that is going to get mixed in
                              ! to the starting-guess of the groundstate density matrix.
        ETEMP = -1.0d0        ! Electronic temperature for calculating temperature smearing of the density matrice and entropy. If ETEMP < 0
                              ! no temperature smearing will be employed.
        SOFTSTART   = .FALSE. ! If true the start up of the XLBOMD and fast-QMMD will be done by doing full scf-convergence the first 50 time-steps using the largest
                              ! disipation possible (i.e lowest ORDER, DORDER = 4 ) and then continiuing with DORDER=4 for additional 50-time steps but with
                              ! ZEROSCF = .TRUE. and finally setting the order, DORDER, to the user specified value and continiuning with XL-BOMD or fast-QMMD (ZEROSCF = .TRUE.)
        DOTDFT   = .FALSE.     ! Weather or not to perform a TDFT calculation
        OMEGA    = 137.0359990740d0*2.0d0*pi/(6000/0.5291772192) ! The default frequency of the electromagnetic field used in the TDFT calculation, corresponding to a 
                               ! Wavelength of 600 nm.
        EDIR     = 1           ! The direction in which the tranverse electromagnetic wave travels. 1 = x (polarized along y), 2 = y (polarized along z ), 3 = z (polarized along x).
        NEPERIOD = 1           ! Number of periods the electromagnetic field in the TDFT calculations retain its full strength.
        EPROFILE = 'HO'        ! Here the profile of the perturbing field is set. DP = dirac puls, AC = Alternating pulse non-homogemeous field in
                               ! direction of travell, HO = Homogeneous field. For EPROFILE =  AC, the wave is an amplitude modulated E-M wave with frequency OMEGA.
        EFIELDMAX = 0.030d0    ! Maxiumum electric field strength used in TDFT calculation. Here the default value of 0.03 a.u. corresponds to an intensity of 3.17*10E+13 W/(cm^(-2))
        IORBNR(1) = 0          ! Index of orbital being excited and treated as a hole. If IORBNR(1) = 0, then no excitation calculation will be done.
                               ! If IORBNR > 0, then a spin-up orbital will be treated as a whole, if IORBNR < 0 a spin-down orbital will be treated as a whole.
        IORBNR(2) = 0         ! Not in use anymore. Rest from version V.31
        AORBS = 0             ! If AORBS > 0, then the spin-up orbital with index |AORBS| willl be saved to the file 'ARBORB.xsf'.
                              ! If AORBS < 0, then the spin-down orbital with index |AORBS| willl be saved to the file 'ARBORB.xsf'.
        ADEF = .FALSE.        ! if true then all the time independent calculations will be performed with a static efield present.
                              ! the parameter EPROFILE the field will have different time evolutions when doing molecular dynamics.
        OPTH = .FALSE.        ! When doing a orthogonality Constrained DFT calculation for excited states, this logical flag decides wheter or not 
                              ! the holes corresponding to the orbital being excited is to be optimized (.TRUE.) or not (.FALSE.)
        CFORCE = .FALSE.
        RELAXN      = .FALSE. ! If true the nuclear positions are relaxed with respect to hartree-fock forces and energy
        NSTEPS      =   500   ! Number of maximum relaxation steps used when relaxing the nuclear positions
        NLSPOINTS   =   10    ! Points used in the polynomial fit along the gradient in the nuclear relaxation
        RELALGO = 1           ! RELALGO = 1, then the relaxation will be done by searching for the configuration that gives FORCES = 0 (relaxf.f90)
                              ! RELALGO = 2, then the relaxation will be done by minimizing the total energy (relax.f90)
        PORDER      =   7     ! PORDER-1 = order of the polynomial used in the least square fit along the direction of the gradient.(nuclear relaxation)
        DR          = 0.5     ! step length used in nuclear relaxation.
        FTol        = 1.0E-4  ! Force tolerance when relaxing molecular structure

        MESH        = (/ 100, 100, 100 /)
        LIMITS      = (/ 5.0d0, 5.0d0, 5.0d0 /)
        LEXCSP      = .FALSE. ! If true any CISD calculation will be limmited to
                              ! to basis sets consisting of slaterdeterminants
                              ! corresponding to the NEEXC highest occupated
                              ! HF-levels and only HF-oneparticle energy levels
                              ! of energy < ENEXCM will be allowed
        MIX = 0.0d0           ! Miximg parameter used in the HF-calculations rho(n)*(1-MIX) + rho(n-1)*MIX -> rho(n+1)
        SPINCONSERVE = .TRUE.
        RESTRICT = .FALSE.
        APPROXEE = .TRUE.     ! If true all the electron-electron tensor elements (i,j|k,l) <= sqrt( (i,j|i,j)*(k,l|kl) ) < EETOL will not be calculated
        EETOL = 1.0E-10       ! Tolerance for the calculation of (i,j|k,l)
        MOLDYN = .FALSE.      ! If true perform molecular dynamics calculation.
        WRITEONFLY = .FALSE.  ! Movie frames and energies are saved on the fly if true, otherwise they are saved after molecular dynamics calculation has finished
        MOVIE = .FALSE.       ! If true, a xfs-file (xcrysden) is saved for animating the molecular dynamics calculation
        TEMPERATURE = 300.0d0 ! Temperature for molecular dynamics calculation
        ZEROSCF = .FALSE.     ! No self consistency is done in molecular dynamics
        ZEROSCFTYPE = 1       ! The type of energy expression used in the zero-scf-MD (fast-QMMD, ZEROSCF = .T.) is selected by this parameter/flag.
                              ! If ZEROSCFTYPE = 1 the standard expression: ENERGY = 2*tr[h*D] + tr[D*G(D)], with  D = theta[mu*I-H(P)] will be used.
                              ! If ZEROSCFTYPE = 2 the expression: ENERGY = 2*tr[h*D] + tr[(2*D-P)*G(P)] will instead be used. Here: G = 2*J - K 
        XLBOMD = .FALSE.      ! If true then A.M. Niklassons method for time-reversible propagation of the density matrix will be used. See
                              ! A.M Niklasson et al in Phys. Rev. Lett. 100,  123004 (2008) and A.M Niklasson et al  in Phys. Rev. B 86, 174308 (2012) 
        PULAY = 4             ! Option for selecting how the Pulay-force is selected. For more details see the subroutine "forces.f90"
        FIXNSCF = -1          ! If FIXNSCF > 0 the number of SCF cykles in a Molecular dynamics run is kept fixed to FIXNSCF
        !===========================================================
        ! Here the XLBOMD default parameters are given:
        !===========================================================
        kappa =  -1.0d0
        alpha =  -1.0d0
        DORDER = 8
        !CN    =  (/ -6.0d0, 14.0d0, -8.0d0, -3.0d0, 4.0d0, -1.0d0 /)
        !===========================================================
        Ne = 0
        SAMPLERATE = 10
        NREPLICAS = 2000
        TIMESTEP = 0.0025
        TEND = 10.0d0
        TSTART = TIMESTEP
        BETA = 0.0d0        ! If BETA .NE. 0.0d0, then the population control energy ER will be continiously recalculated
        BJASTROW = 1.0d0
        CJASTROW = 0.50d0
        NPERSIST = 50
        NRECALC = INT(1.0d0/TIMESTEP)            ! The population control energy ER will be recalculated every RESAMPVAL:th time-step.
        REDISTRIBUTIONFREQ = 0 ! The redistribution frequency for the walkers/replicas in the diffusion quantum monte carlo calculation. 
                               ! if not equal to zero the walkers/replicas are evenly redestributed everey REDISTRIBUTIONFREQ:th time step.
                               ! This parameter only exists in the MPI-version of the program.
        NOREDIST = .FALSE.     ! This parameter only exists in the MPI-version of the program. If True no redistribution of the walkers will 
                               ! will be permitted over the computational grid
        CUTTOFFFACTOR = 1.0d0
        DIISORD = 3     ! The order ( number of previous density matrices ) to
                        ! be used in direct inversion in the iterative subspace
                        ! method DIIS. The integere m in Eqn (3) in
                        ! Chem. Phys. Lett. 73, 393 (1980)
        DIISORDEX = 3   ! Same as above but for the excitation calculation
        DIISSTART = 50  ! The maximum number of iterations allowed without the
                        ! DIIS-mixing sheme starts. After DIISSTART electronic
                        ! iterations the DIIS starts regrdless the size of the 
                        ! density matrix difference vector deltaP ( Eqn. (4) in
                        ! Chem. Phys. Lett. 73, 393 (1980)
        CUSPCORR = .TRUE.       ! If true the basis functions consisting
                                ! of contracted primitive gaussians are 
                                ! corrected at the nuclear cusps. ( see:  J.Chem.Phys. 115, 5362 (2001) )
        rc = 0.10d0             ! The cusp correction comes into play for r <= rc
        CORRALCUSP = .FALSE.    ! If true all basis functions will be cusp
                                ! corrected. If false only basis functions
                                ! constructed from contractions of more than 1
                                ! primitive gaussiam will be cusp corrected.
        NVMC = 1000000          ! Number of metropolis moves used in the
                                ! variational monte carlo runs
        HFORBWRITE = .FALSE.    ! If true a HF-orbital will be saved on disk in
                                ! the file ORBITAL.dat. The index of the orbital 
                                ! being saved is given by IOSA
        IOSA = 0                ! The index of the HF-orbital that is saved if
                                ! HFORBWRITE = .TRUE.
        NLEBEDEV = 3            ! Parameter used to choose the angular mesh in
                                ! the lebedev quadrature ( 1 = 110 points, 2 =
                                ! 170 points, 3 = 194 points, 4 = 230 points, 5
                                ! = 266 points, 6 = 302, 7 = 350 points, 8 = 434
                                ! points, 9 = 590 points, 10 = 770, 11 = 974, 12
                                ! = 1202 points
        NCHEBGAUSS = 100        ! the number of radial points to be used in the
                                ! Gauss-chebyshev quadrature
        J = 0
        DO I=1,TOTALNLINES
                LINEHASBEANREAD = .FALSE.
                IF ( VARNAME(I) .EQ. 'NATOMS' ) THEN
                        READ(10,*)DUMMY
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'TOL') THEN
                        READ(10,*)DUMMY,Tol
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MESH') THEN
                        READ(10,*)DUMMY,MESH
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'LIMITS') THEN
                        READ(10,*)DUMMY,LIMITS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'Ne') THEN
                        READ(10,*)DUMMY,Ne
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CORRLEVEL') THEN 
                        READ(10,*)DUMMY,CORRLEVEL
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'WRITECICOEF') THEN 
                        READ(10,*)DUMMY,WRITECICOEF
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'WRITEDENS') THEN 
                        READ(10,*)DUMMY,WRITEDENS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'WHOMOLUMO') THEN 
                        READ(10,*)DUMMY,WHOMOLUMO
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ATOM' ) THEN
                        J = J+1
                        READ(10,*)DUMMY,ATOMICNUMBER(J),SPATIALPOSITION(J,:)
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'LEXCSP') THEN
                        READ(10,*)DUMMY,LEXCSP
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NEEXC') THEN
                        READ(10,*)DUMMY,NEEXC
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ENEXCM') THEN
                        READ(10,*)DUMMY,ENEXCM
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'SPINCONSERVE') THEN
                        READ(10,*)DUMMY,SPINCONSERVE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'RESTRICT') THEN
                        READ(10,*)DUMMY,RESTRICT
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'APPROXEE') THEN
                        READ(10,*)DUMMY,APPROXEE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'SAMPLERATE' ) THEN
                        READ(10,*)DUMMY,SAMPLERATE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NREPLICAS' ) THEN
                        READ(10,*)DUMMY,NREPLICAS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'TIMESTEP' ) THEN
                        READ(10,*)DUMMY,TIMESTEP
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'TEND' ) THEN
                        READ(10,*)DUMMY,TEND
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'TSTART' ) THEN
                        READ(10,*)DUMMY,TSTART
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'BETA' ) THEN
                        READ(10,*)DUMMY,BETA
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'BJASTROW' ) THEN
                        READ(10,*)DUMMY,BJASTROW
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CJASTROW' ) THEN
                        READ(10,*)DUMMY,CJASTROW
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NPERSIST' ) THEN
                        READ(10,*)DUMMY,NPERSIST
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'REDISTRIBUTIONFREQ' ) THEN
                        READ(10,*)DUMMY,REDISTRIBUTIONFREQ
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NOREDIST' ) THEN
                        READ(10,*)DUMMY,NOREDIST
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NRECALC' ) THEN
                        READ(10,*)DUMMY,NRECALC
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CUTTOFFFACTOR' ) THEN
                        READ(10,*)DUMMY,CUTTOFFFACTOR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MIX' ) THEN
                        READ(10,*)DUMMY,MIX
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DIISORD' ) THEN
                        READ(10,*)DUMMY,DIISORD
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DIISSTART' ) THEN
                        READ(10,*)DUMMY,DIISSTART
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CUSPCORR' ) THEN
                        READ(10,*)DUMMY,CUSPCORR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'rc' ) THEN
                        READ(10,*)DUMMY,rc
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CORRALCUSP' ) THEN
                        READ(10,*)DUMMY,CORRALCUSP
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NVMC' ) THEN
                        READ(10,*)DUMMY,NVMC
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'HFORBWRITE' ) THEN
                        READ(10,*)DUMMY,HFORBWRITE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'IOSA' ) THEN
                        READ(10,*)DUMMY,IOSA
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'CFORCE' ) THEN
                        READ(10,*)DUMMY,CFORCE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'RELAXN' ) THEN
                        READ(10,*)DUMMY,RELAXN
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NSTEPS' ) THEN
                        READ(10,*)DUMMY,NSTEPS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DR' ) THEN
                        READ(10,*)DUMMY,DR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'FTOL' ) THEN
                        READ(10,*)DUMMY,FTol
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NLSPOINTS' ) THEN
                        READ(10,*)DUMMY,NLSPOINTS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'PORDER' ) THEN
                        READ(10,*)DUMMY,PORDER
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'WRITEONFLY' ) THEN
                        READ(10,*)DUMMY,WRITEONFLY
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MOVIE' ) THEN
                        READ(10,*)DUMMY,MOVIE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MOLDYN' ) THEN
                        READ(10,*)DUMMY,MOLDYN
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'TEMPERATURE' ) THEN
                        READ(10,*)DUMMY,TEMPERATURE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ZEROSCF' ) THEN
                        READ(10,*)DUMMY,ZEROSCF
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'XLBOMD' ) THEN
                        READ(10,*)DUMMY,XLBOMD
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'KAPPA' ) THEN
                        READ(10,*)DUMMY,kappa
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ALPHA' ) THEN
                        READ(10,*)DUMMY,alpha
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DORDER' ) THEN
                        READ(10,*)DUMMY,DORDER
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'PULAY' ) THEN
                        READ(10,*)DUMMY,PULAY
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'FIXNSCF' ) THEN
                        READ(10,*)DUMMY,FIXNSCF
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NLEBEDEV' ) THEN
                        READ(10,*)DUMMY,NLEBEDEV
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NCHEBGAUSS') THEN
                        READ(10,*)DUMMY,NCHEBGAUSS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'EETOL') THEN
                        READ(10,*)DUMMY,EETOL
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'RELALGO') THEN
                        READ(10,*)DUMMY,RELALGO
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'HUCKEL') THEN
                        READ(10,*)DUMMY,HUCKEL
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'SOFTSTART') THEN
                        READ(10,*)DUMMY,SOFTSTART
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ZEROSCFTYPE') THEN
                        READ(10,*)DUMMY,ZEROSCFTYPE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ETEMP') THEN
                        READ(10,*)DUMMY,ETEMP
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'IORBNR') THEN
                        READ(10,*)DUMMY,IORBNR(1)
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'AORBS') THEN
                        READ(10,*)DUMMY,AORBS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DIISORDEX') THEN
                        READ(10,*)DUMMY,DIISORDEX
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DOTDFT') THEN
                        READ(10,*)DUMMY,DOTDFT
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'OMEGA') THEN
                        READ(10,*)DUMMY,OMEGA
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'EDIR' .AND. .not. FIELDREAD ) THEN
                        READ(10,*)DUMMY,EDIR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NEPERIOD') THEN
                        READ(10,*)DUMMY,NEPERIOD
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'EPROFILE') THEN
                        READ(10,*)DUMMY,EPROFILE
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'EFIELDMAX') THEN
                        READ(10,*)DUMMY,EFIELDMAX
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'ADEF') THEN
                        READ(10,*)DUMMY,ADEF
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'OPTH') THEN
                        READ(10,*)DUMMY,OPTH
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MIXEX') THEN
                        READ(10,*)DUMMY,MIXEX
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DIFFDENS') THEN
                        READ(10,*)DUMMY,DIFFDENS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'NSCCORR') THEN
                        READ(10,*)DUMMY,NSCCORR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'MIXTDDFT') THEN
                        READ(10,*)DUMMY,MIXTDDFT
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'SCERR') THEN
                        READ(10,*)DUMMY,SCERR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'RIAPPROX') THEN
                        READ(10,*)DUMMY,RIAPPROX
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'LIMPRECALC') THEN
                        READ(10,*)DUMMY,LIMPRECALC
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DOABSSPECTRUM') THEN
                        READ(10,*)DUMMY,DOABSSPECTRUM
                        LINEHASBEANREAD = .TRUE.
                        ABSSPECTSET = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DIAGDG') THEN
                        READ(10,*)DUMMY,DIAGDG
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DAMPING') THEN
                        READ(10,*)DUMMY,DAMPING
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'FIELDDIR') THEN
                        READ(10,*)DUMMY,FIELDDIR
                        LINEHASBEANREAD = .TRUE.
                        FIELDREAD = .TRUE.
                        EDIR = 2
                ENDIF
                IF ( VARNAME(I) .EQ. 'SCALARRELC') THEN
                        READ(10,*)DUMMY,SCALARRELC
                        LINEHASBEANREAD = .TRUE.
                        FIELDREAD = .TRUE.
                        EDIR = 2
                ENDIF
                IF ( VARNAME(I) .EQ. 'URHFTORHF') THEN
                        READ(10,*)DUMMY,URHFTORHF
                        LINEHASBEANREAD = .TRUE.
                        FIELDREAD = .TRUE.
                        EDIR = 2
                ENDIF
                IF ( .not. LINEHASBEANREAD )READ(10,*)DUMMY
        ENDDO
        
        CLOSE(10)
        
        !=====================================================================
        ! If the user has not explicitly stated that he/she wants to calculate
        ! the absorbtion spectrum than this calculation will be done if the 
        ! the fourier-grid is reasonably small, i.e less than 10000
        !======================================================================
        IF ( .not. ABSSPECTSET .AND. INT(TEND/TIMESTEP) .LT. 10000 ) DOABSSPECTRUM = .TRUE.
        
      ELSE
              print*,'ABORTING SINCE FILE NAMED INPUTFILE IS MISSING'
              STOP
      ENDIF
END SUBROUTINE readin
