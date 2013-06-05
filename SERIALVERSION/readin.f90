SUBROUTINE readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBER,SPATIALPOSITION,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, & 
& SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,CUSPCORR,rc,CORRALCUSP,NVMC,HFORBWRITE,IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol, &
& NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS,EETOL,RELALGO,SOFTSTART)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS,NLINES
      CHARACTER(LEN=20), INTENT(OUT) :: CORRLEVEL
      INTEGER, INTENT(OUT) :: ATOMICNUMBER(NATOMS),Ne,MESH(3),NEEXC,SAMPLERATE,NREPLICAS,NRECALC,NPERSIST,DIISORD,DIISSTART,NVMC,IOSA,NSTEPS,NLSPOINTS,PORDER,PULAY,FIXNSCF,DORDER
      INTEGER, INTENT(OUT) :: NLEBEDEV,NCHEBGAUSS,RELALGO
      DOUBLE PRECISION,INTENT(OUT) :: SPATIALPOSITION(NATOMS,3),Tol,LIMITS(3),ENEXCM,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,CUTTOFFFACTOR,MIX,rc,DR,FTol,TEMPERATURE,kappa,alpha
      DOUBLE PRECISION,INTENT(OUT) :: EETOL
      LOGICAL, INTENT(OUT) :: WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,CUSPCORR,CORRALCUSP,HFORBWRITE,CFORCE,RELAXN,WRITEONFLY,MOVIE,MOLDYN,ZEROSCF,XLBOMD,SOFTSTART
      INTEGER :: I,J,TOTALNLINES
      LOGICAL :: finns,LINEHASBEANREAD
      CHARACTER(LEN=20), ALLOCATABLE :: VARNAME(:)
      CHARACTER(LEN=2) :: DUMMY
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

        ! Here we put some default values:
        Tol = 1.0E-6
        WRITECICOEF = .FALSE. 
        WRITEDENS   = .FALSE.
        WHOMOLUMO   = .FALSE.
        SOFTSTART   = .FALSE. ! If true the start up of the XLBOMD and fast-QMMD will be done by doing full scf-convergence the first 50 time-steps using the largest 
                              ! disipation possible (i.e lowest ORDER, DORDER = 4 ) and then continiuing with DORDER=4 for additional 50-time steps but with
                              ! ZEROSCF = .TRUE. and finally setting the order, DORDER, to the user specified value and continiuning with XL-BOMD or fast-QMMD (ZEROSCF = .TRUE.)
        CFORCE = .FALSE.
        MESH        = (/ 100, 100, 100 /)
        LIMITS      = (/ 5.0d0, 5.0d0, 5.0d0 /)
        RELAXN      = .FALSE. ! If true the nuclear positions are relaxed with respect to hartree-fock forces and energy
        NSTEPS      =   500   ! Number of maximum relaxation steps used when relaxing the nuclear positions
        DR          = 0.5     ! Initial intervall length of line search used in nuclear relaxation.
        NLSPOINTS   = 10      ! Number of points used in the polynomial fit used in the line search of the nuclear relaxation
        RELALGO = 1           ! RELALGO = 1, then the relaxation will be done by searching for the configuration that gives FORCES = 0 (relaxf.f90)
                              ! RELALGO = 2, then the relaxation will be done by minimizing the total energy (relax.f90)
        PORDER      = 7       ! PORDER-1 = the order of the polynomial used in the least square fit used in the line search 
        FTol        = 1.0E-04 ! Tolerance for the forces
        LEXCSP      = .FALSE. ! If true any CISD calculation will be limmited to
                              ! to basis sets consisting of slaterdeterminants
                              ! corresponding to the NEEXC highest occupated
                              ! HF-levels and only HF-oneparticle energy levels
                              ! of energy < ENEXCM will be allowed
        SPINCONSERVE = .TRUE.
        RESTRICT = .FALSE.
        APPROXEE = .TRUE.     ! If true all the electron-electron tensor elements (i,j|k,l) <= sqrt( (i,j|i,j)*(k,l|kl) ) < EETOL will not be calculated
        EETOL = 1.0E-10       ! Tolerance for the calculation of (i,j|k,l)
        MOLDYN = .FALSE.      ! If true perform molecular dynamics calculation.
        WRITEONFLY = .FALSE.  ! Movie frames and energies are saved on the fly if true, otherwise they are saved after molecular dynamics calculation has finished
        MOVIE = .FALSE.       ! If true, a xfs-file (xcrysden) is saved for animating the molecular dynamics calculation
        TEMPERATURE = 300.0d0 ! Temperature for molecular dynamics calculation
        ZEROSCF = .FALSE.     ! No self consistency is done in molecular dynamics
        XLBOMD = .FALSE.      ! If true then A.M. Niklassons method for time-reversible propagation of the density matrix will be used. See 
                              ! A.M Niklasson et al in Phys. Rev. Lett. 100, 123004 (2008) and A.M Niklasson et al  in Phys. Rev. B 86, 174308 (2012)
        PULAY = 4             ! Option for selecting how the Pulay-force is selected. For more details see the subroutine "forces.f90"
        FIXNSCF = -1          ! If FIXNSCF > 0 the number of SCF cykles in a Molecular dynamics run is kept fixed to FIXNSCF
        !==================================================================================
        ! Here the XLBOMD default parameters are given: ( Phys. Rev. B 86, 174308 (2012) )
        !==================================================================================
        kappa = -1.0d0
        alpha = -1.0d0
        DORDER = 8
        !CN    =  (/ -6.0d0, 14.0d0, -8.0d0, -3.0d0, 4.0d0, -1.0d0 /)
        !===========================================================
        Ne = 0
        SAMPLERATE = 10
        NREPLICAS = 2000
        TIMESTEP = 0.00250d0
        TEND = 10.0d0
        TSTART = TIMESTEP
        BETA = 0.0d0  ! If BETA .NE. 0.0d0, then the population control energy ER will be continiously recalculated
        BJASTROW = 1.0d0
        CJASTROW = 0.50d0
        NPERSIST = 50
        NRECALC = INT(1.0d0/TIMESTEP)
        CUTTOFFFACTOR = 1.0d0
        MIX = 0.0d0     ! Linear mixing of the density matrices in scf-calculations.
        DIISORD = 3     ! The order ( number of previous density matrices ) to
                        ! be used in direct inversion in the iterative subspace
                        ! method DIIS. The integere m in Eqn (3) in 
                        ! Chem. Phys. Lett. 73, 393 (1980)
        DIISSTART = 50  ! The maximum number of iterations allowed without the
                        ! DIIS-mixing sheme starts. After DIISSTART electronic 
                        ! iterations the DIIS starts regrdless the size of the 
                        ! density matrix difference vector deltaP ( Eqn. (4) in 
                        ! Chem. Phys. Lett. 73, 393 (1980)
        CUSPCORR = .TRUE.       ! If true the basis functions consisting
                                ! of contracted primitive gaussians are 
                                ! corrected at the nuclear cusps. ( see: J.Chem.Phys. 115, 5362 (2001) )
        rc = 0.10d0             ! The cusp correction comes into play for r <= rc
        CORRALCUSP = .FALSE.    ! If true all basis functions will be cusp
                                ! corrected. If false only basis functions
                                ! constructed from contractions of more than 1
                                ! primitive gaussiam will be cusp corrected.
        NVMC = 100000           ! Number of metropolis moves used in the
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
                IF ( VARNAME(I) .EQ. 'NSTEPS' ) THEN
                        READ(10,*)DUMMY,NSTEPS
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'DR' ) THEN
                        READ(10,*)DUMMY,DR
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( VARNAME(I) .EQ. 'RELAXN' ) THEN
                        READ(10,*)DUMMY,RELAXN
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
                IF ( VARNAME(I) .EQ. 'SOFTSTART') THEN
                        READ(10,*)DUMMY,SOFTSTART
                        LINEHASBEANREAD = .TRUE.
                ENDIF
                IF ( .not. LINEHASBEANREAD )READ(10,*)DUMMY
        ENDDO
        
        CLOSE(10)
      
      ELSE
              print*,'ABORTING SINCE FILE NAMED INPUTFILE IS MISSING'
              STOP
      ENDIF
END SUBROUTINE readin
