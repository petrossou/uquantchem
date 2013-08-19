PROGRAM uquantchem
      USE datatypemodule
      IMPLICIT NONE
      INTEGER  :: NATOMS
      INTEGER, ALLOCATABLE :: ATOMICNUMBERS(:)
      DOUBLE PRECISION, ALLOCATABLE :: APOS(:,:)
      DOUBLE PRECISION :: Tol,ETOT,NucE,Rn(3)
      CHARACTER(LEN=20) :: CORRLEVEL
      INTEGER :: BASISMAP(120)
      TYPE(ATOM), ALLOCATABLE :: ATOMS(:)
      TYPE(BASIS) :: BAS
      CHARACTER(LEN=6) :: DUMMY
      CHARACTER(Len=20) :: date,time,zone
      INTEGER :: NLINES,I,J,K,L,M,NB,Ne,NRED,Lmax,MESH(3),NVMC,IOSA,FNATOMS,FNRED,NSTEPS,NLSPOINTS,PORDER
      LOGICAL :: finns,WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,HYBRID,USEGTO,CORRALCUSP,VMCCALC,HFORBWRITE,CFORCE,RELAXN,WRITEONFLY,MOVIE,MOLDYN
      DOUBLE PRECISION, ALLOCATABLE :: S(:,:),T(:,:),V(:,:),H0(:,:),Intsv(:),EIGENVECT(:,:),Ints(:,:,:,:),gradIntsv(:,:,:),DMAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gradS(:,:,:,:),gradT(:,:,:,:),gradV(:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigenup(:),EHFeigendown(:),Cup(:,:),Cdown(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigen(:),P(:,:),Pup(:,:),Pdown(:,:),C1(:,:),C2(:,:),force(:,:),CN(:),LQ(:,:),CGQ(:,:),ZMUL(:)
      DOUBLE PRECISION  :: PRYSR(25,25),PRYSW(25,25),rts(25),wts(25),TOTALTIME,ECISD,LIMITS(3),ENEXCM,EMP2,TIMESTEP,TEND,TSTART,DR,FTol,CNSTART(4),alphastart,kappastart
      DOUBLE PRECISION :: ESIGMA,EDQMC,a,b,BETA,BJASTROW,CJASTROW,NPERSIST,CUTTOFFFACTOR,MIX,Aexpan,LAMDA,rc,GAMA,BETAA,POLY(6),deltar,TEMPERATURE,kappa,alpha,EETOL
      INTEGER :: STARTTIME(8),STOPTIME(8),RUNTIME(8),NEEXC,SAMPLERATE,NREPLICAS,NRECALC,DIISORD,DIISSTART,NTIMESTEPS,NSCF,PULAY,FIXNSCF,DORDER
      INTEGER :: NLEBEDEV,NCHEBGAUSS,LEBPOINTS(24),CGORDER,LORDER,RELALGO
      INTEGER, EXTERNAL :: ijkl
      DOUBLE PRECISION, EXTERNAL :: massa
      LOGICAL :: RESTART,ZEROSCF,XLBOMD,DFTC,SOFTSTART

      !-------------------
      !STARTING THE CLOCK: 
      !-------------------
     
      call DATE_AND_TIME(date, time, zone,STARTTIME)
     
      DFTC = .FALSE. 

      !=======================================================================
      ! (1) HERE THE DATA FROM THE INPUTFILES: INPUTFILE & BASISFILE IS LOADED
      !=======================================================================

      ! preparing the readin of data by aquiring the number of atoms, NATOMS,
      ! and the number of lines, NLINES, appearing before the entry of the
      ! characters 'NATOMS' in the file 'INPUTFILE'
        
      CALL checkinputfile(NLINES,NATOMS)
      
      ALLOCATE(ATOMICNUMBERS(NATOMS),APOS(NATOMS,3),ATOMS(NATOMS),force(NATOMS,3),ZMUL(NATOMS))
      
      CALL readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBERS,APOS,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, &
      & SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,HYBRID,rc,CORRALCUSP,NVMC,HFORBWRITE,IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol, &
      & NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS,EETOL,RELALGO,SOFTSTART)

      IF ( DORDER .LT. 4 )  DORDER = 4
      ALLOCATE(CN(DORDER))
        
      IF ( RELALGO .NE. 1 .AND. RELALGO .NE. 2 ) RELALGO = 1

      ! The array containing the number of Lebedev quadrature points is assigned here:
      LEBPOINTS = (/ 110, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)

      ! Here the arrays containg the Lebedev and Chebyshev-Gauss quadrature
      ! points and abssiscae are allocated and created in the case that one has
      ! choosen to perform a DFT calculation
      IF ( CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'B3LYP') THEN
              DFTC = .TRUE.
              CGORDER = NCHEBGAUSS
              LORDER  = LEBPOINTS(NLEBEDEV)
              ALLOCATE(LQ(LORDER,3),CGQ(CGORDER,2))
              IF ( NLEBEDEV .GE. 1  .AND. NLEBEDEV .LE. 5   ) CALL lebedev1(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 6  .AND. NLEBEDEV .LE. 10  ) CALL lebedev2(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 11 .AND. NLEBEDEV .LE. 15  ) CALL lebedev3(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 16 .AND. NLEBEDEV .LE. 20  ) CALL lebedev4(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 21 .AND. NLEBEDEV .LE. 24  ) CALL lebedev5(NLEBEDEV,LORDER,LQ)
              CALL chebgauss(CGORDER,CGQ)
      ELSE
              CGORDER = 1
              LORDER = 1
              ALLOCATE(LQ(LORDER,3),CGQ(CGORDER,2))
              LQ = 0.0d0
              CGQ = 0.0d0
      ENDIF

      !================================================================================================
      ! assigning the disipative coeficients depending of the order of
      ! disipation in the XL-BOMD time propagation of the density matrix 
      ! J. Chem. Phys. 130, 214109 (2009)
      !================================================================================================
      alphastart = 0.150d0
      kappastart = 1.690d0
      CNSTART = (/ -2.0d0, 3.0d0, 0.0d0, -1.0d0 /)
      
      IF ( DORDER .EQ. 4 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.150d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.690d0
              CN = (/ -2.0d0, 3.0d0, 0.0d0, -1.0d0 /)
      ENDIF
      IF ( DORDER .EQ. 5 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.0570d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.750d0
              CN = (/ -3.0d0, 6.0d0, -2.0d0, -2.0d0, 1.0d0 /)
      ENDIF
      IF ( DORDER .EQ. 6 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.0180d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.820d0
              CN = (/ -6.0d0, 14.0d0, -8.0d0, -3.0d0, 4.0d0, -1.0d0 /)
      ENDIF
      IF ( DORDER .EQ. 7 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.00550d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.840d0
              CN = (/ -14.0d0, 36.0d0, -27.0d0, -2.0d0, 12.0d0, -6.0d0, 1.0d0/)
      ENDIF
      IF ( DORDER .EQ. 8 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.00160d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.860d0
              CN = (/ -36.0d0, 99.0d0, -88.0d0, 11.0d0, 32.0d0, -25.0d0, 8.0d0, -1.0d0 /)
      ENDIF
      IF ( DORDER .EQ. 9 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.000440d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.880d0
              CN = (/ -99.0d0, 286.0d0, -286.0d0, 78.0d0, 78.0d0, -90.0d0, 42.0d0, -10.0d0, 1.0d0 /)
      ENDIF
      IF ( DORDER .EQ. 10 ) THEN
              IF ( alpha .LT. 0.0d0 ) alpha = 0.000120d0
              IF ( kappa .LT. 0.0d0 ) kappa = 1.890d0
              CN = (/ -286.0d0, 858.0d0, -936.0d0, 364.0d0, 168.0d0, -300.0d0, 184.0d0, -63.0d0, 12.0d0, -1.0d0 /)
      ENDIF
      !================================================================================================

      IF ( FIXNSCF .EQ. 0 ) FIXNSCF = -1

      IF ( CORRLEVEL .NE. 'DQMC' .OR. CORRLEVEL .NE. 'VMC' ) HYBRID = .FALSE.
      
      VMCCALC = .FALSE.
      IF ( CORRLEVEL .EQ. 'VMC' ) VMCCALC = .TRUE.
      IF ( NATOMS .EQ. 1 ) CFORCE = .FALSE.
      
      IF ( XLBOMD ) MOLDYN = .TRUE.

      IF ( RELAXN .OR. MOLDYN ) CFORCE = .TRUE.
      ! If the user accidently forgets to set the number of electrons in the
      ! input file.
      IF ( Ne .EQ. 0 ) THEN
               DO I=1,NATOMS
                        Ne = Ne + ATOMS(I)%Z
               ENDDO
       ENDIF

      IF ( .not. LEXCSP ) NEEXC = Ne

      !======================================================================
      ! If we are running a continued MD-run we need to udate the atomic
      ! positions from the MOLDYNRESTART.dat-file
      !======================================================================
      RESTART = .FALSE.
      inquire(file='MOLDYNRESTART.dat',exist=RESTART)
      IF ( MOLDYN .AND. RESTART ) THEN
              OPEN(300,FILE='MOLDYNRESTART.dat',STATUS='OLD',ACTION='READ')
              READ(300,'(I10)')I
              DO I=1,NATOMS
                        READ(300,'(3(E30.20))')APOS(I,:)
              ENDDO
              CLOSE(300)
      ENDIF
      !======================================================================
      DO I=1,NATOMS
        ATOMS(I)%Z = ATOMICNUMBERS(I)
        ATOMS(I)%M = massa(ATOMICNUMBERS(I))
        ATOMS(I)%R = APOS(I,:)
      ENDDO

      ! Here the nuclear-nuclear repulsion energy is calcluated:
      
      NucE = 0.0d0
      
      DO I=1,NATOMS
        DO J=I+1,NATOMS
                Rn = ATOMS(I)%R - ATOMS(J)%R
                NucE = NucE + ATOMS(I)%Z*ATOMS(J)%Z/sqrt(DOT_PRODUCT(Rn,Rn))
        ENDDO
      ENDDO

      ! Loading the basis set data such as the gaussian contraction coefficients
      ! and exponents.

      inquire(file='BASISFILE',exist=finns)
      IF ( finns ) THEN
             BASISMAP(:) = 0
             ! Here the number of basis functions of all the atoms in the 
             ! file "BASISFILE" are loaded into the array BASISMAP
             
             CALL readbasis(NATOMS,ATOMS,BASISMAP,.TRUE.)
             
             DO I=1,NATOMS
                ALLOCATE(ATOMS(I)%PSI(BASISMAP(ATOMS(I)%Z)))
             ENDDO
             
             ! Here the local basis for each atom is loaded
             
             CALL readbasis(NATOMS,ATOMS,BASISMAP,.FALSE.)

             ! CALCULATING THE STO-EXPONENTS AND NORMS
             ! SAVING AND SAVING IT IN THE HYBRID BASIS
             
             deltar = 0.0010d0

             DO I= 1,NATOMS
                DO J=1,ATOMS(I)%NBAS
                        call stoexponent(ATOMS(I)%PSI(J),deltar,LAMDA,Aexpan,GAMA,BETAA,POLY,rc,USEGTO,CORRALCUSP)
                        ATOMS(I)%PSI(J)%STOEXPON = LAMDA
                        ATOMS(I)%PSI(J)%A = Aexpan
                        ATOMS(I)%PSI(J)%RC = rc
                        ATOMS(I)%PSI(J)%P = POLY
                        ATOMS(I)%PSI(J)%DR = deltar
                        ATOMS(I)%PSI(J)%GAMA = GAMA
                        ATOMS(I)%PSI(J)%BETA = BETAA

                        IF ( .not. USEGTO .AND. HYBRID ) THEN
                                ATOMS(I)%PSI(J)%HYBRID = .TRUE.
                        ELSE
                                ATOMS(I)%PSI(J)%HYBRID = .FALSE.
                        ENDIF
                ENDDO
             ENDDO

             
      ELSE
             print*,'ABORTING SINCE THE FILE  "BASISFILE" IS MISSING'
      ENDIF
        !===========================================================================
        ! (2) HERE THE MOLECULAR BASIS SET IS CONSTRUCTED FROM ALL THE LOCAL ATOMIC
        !     BASISSETS, BY ADDING ALL THE L-SUBSHELLS AND CONCATANATING THE
        !     LOCAL ATOMIC BASES.
        !===========================================================================
        
        ! Here the total number of basis functions is counted
        CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.TRUE.)
        
        BAS%NBAS = NB
        ALLOCATE(BAS%PSI(NB))
        
        IF (HFORBWRITE .AND. IOSA .GT. NB ) THEN
                WRITE(*,*)'                    WARNING!  '
                WRITE(*,*)'      THE INDEX OF THE HF-EIGEN FUNCTION TO BE SAVED '
                WRITE(*,*)'EXCEEDS THE TOTAL NUMBER OF BASIS FUNCTIONS NB =',NB
                WRITE(*,*)'THE INDEX HF-EIGENFUNCTION SAVED WILL THERFORE BE SET TO NB'
                IOSA = NB
        ENDIF
        IF (HFORBWRITE .AND. IOSA .EQ. 0 ) THEN
                IOSA = (Ne-MOD(Ne,2))/2 + MOD(Ne,2)
        ENDIF

        ! Here the total molecular basis set is put together
        CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
        
        !=======================================================================
        ! (3) HERE THE OVERLAP-, KINETIC- AND POTENTIAL ENERGY- MATRICES ARE
        !     CALCULATED.
        !=======================================================================
        
        ! Here the normalization of the basis functions is performed

        CALL normalize(BAS)

        ! Here the overlap matrix is calculated 

        ALLOCATE(S(NB,NB),gradS(NATOMS,3,NB,NB),DMAT(NB,NB))
        
        CALL overlap(NATOMS,BAS,S,gradS)
        
        ! Here the kinetic energy matrix is calculated

        ALLOCATE(T(NB,NB),gradT(NATOMS,3,NB,NB))

        CALL kinetic(NATOMS,BAS,T,gradT)
        
        ! Here the potential energy matrix is calculated

        ALLOCATE(V(NB,NB),gradV(NATOMS,3,NB,NB))

        CALL potential(BAS,NATOMS,ATOMS,V,gradV)
       
        ! If a DFT type of calculation is performed the 
        ! exchange correlation matrix elements and their
        ! nuclear gradients are allocated 

        !IF ( DFTC ) ALLOCATE(Vxc(2,NB,NB),gVxc(NATOMS,2,3,NB,NB))

        !===========================================================================
        ! (4) Here the electron-electron repulsion 'tensor', Intsv(:), is calulated
        !===========================================================================
        
        ! The number of symmetry reduced electron-electron tensor elements

        NRED = ( (((NB+1)*NB)/2)**2 + ((NB+1)*NB)/2 )/2
        
        ALLOCATE(Intsv(NRED))
        IF ( CFORCE ) THEN
                FNRED = NRED
                FNATOMS = NATOMS
        ELSE
                FNRED =1 
                FNATOMS =1 
        ENDIF
        
        ALLOCATE(gradIntsv(FNATOMS,3,FNRED))

        ! Here we pre calculate the roots and weights of all the rys polynomials 
        ! up to 2*Lmax + 1, (Lmax = maximum L guanum number of basisset). Since 
        ! these roots and weights correspnd to the calculation of (ij|kl) when 
        ! all four orbitals i,j,k,l are located on the same cite, i.e X = 0, 
        ! this will save a lot of time in the calculation of (ij|kl) especially 
        ! for small molecules

        ! Determining Lmax
        Lmax = 0
        DO I=1,NATOMS
                DO J=1,ATOMS(I)%NBAS
                        IF ( Lmax .LT. ATOMS(I)%PSI(J)%L(1) ) Lmax = ATOMS(I)%PSI(J)%L(1)
                ENDDO
        ENDDO
        
        ! Pre calculating the Rys weights and roots:
        IF ( .not. RELAXN .AND. .not. MOLDYN ) THEN
                print*,'    ========================================================'
                print*,'       Pre calculation of Rys-polynomial weights and roots  '
                print*,'    ========================================================'
        ENDIF

        PRYSR = 0.0d0
        PRYSW = 0.0d0

        DO I=1,2*Lmax+1

                rts = 0.0d0
                wts = 0.0d0
                
               call RysQuad(I,0.0d0,rts,wts)
               
               PRYSR(I,:) = rts
               PRYSW(I,:) = wts

        ENDDO
        
        IF ( .not. RELAXN .AND. .not. MOLDYN ) THEN
                print*,'    ========================================================'
                print*,'        Calculation of electron repulsion tensor (ij|kl)    '
                print*,'    ========================================================'
        ENDIF
        
        DMAT = 1.0d0

        CALL eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,EETOL,CFORCE,DMAT)

        ALLOCATE(H0(NB,NB))
        
        H0 = T+V

        IF ( RELAXN .AND. .not. MOLDYN ) THEN
                ! Relxing by searching for minimum energy
                IF ( RELALGO .EQ. 2 ) THEN
                        CALL relax(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS,& 
                                        & APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
                ! Relaxing by searching for FORCES = 0
                IF (RELALGO .EQ. 1 ) THEN
                        CALL relaxf(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS,& 
                                        & APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
        ENDIF

        IF ( MOLDYN ) THEN
                print*,' ========================================================'
                print*,'           Comencing with molecular dynamics calculation '
                print*,' ========================================================'
                
                NTIMESTEPS = INT(TEND/TIMESTEP)
                IF ( SOFTSTART ) THEN
                  CALL moleculardynamicssoft(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,TIMESTEP,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW,&
                  & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,LORDER,CGORDER,LQ,CGQ,EETOL,CNSTART,alphastart,kappastart)
                ELSE
                  CALL moleculardynamics(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,TIMESTEP,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW,&
                  & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
        ENDIF

        !==============================================================
        ! (5) Here the Unrestricted or Restricted Hartree-FocK
        !     equations are selfconsistently solved
        !==============================================================
        
        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                
                ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB),P(NB,NB),C1(NB,NB))
                P = 0.0d0
                CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1, .TRUE., .TRUE.,.FALSE. )
                
                IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)
                !====================
                ! CALCULATING FORCES
                !================================================================================================================== 
                IF ( CFORCE ) THEN
                        P = P/2.0d0
                        CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,LORDER,CGORDER,LQ,CGQ,force)
                        
                        WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                        DO I=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                        ENDDO
                ENDIF
                !================================================================================================================== 

                IF ( HFORBWRITE ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,EIGENVECT,EIGENVECT)
                
                ! If specified, the charge density is calculated and saved here:
                IF ( WRITEDENS ) THEN 
                        C1 = 0.0d0
                        DO I=1,Ne/2
                                C1(:,I) = EIGENVECT(:,I)
                        ENDDO
                        CALL makedens(C1,NB,P)
                        P = 2*P
                        CALL chargedenssave(LIMITS,MESH,BAS,P,Ne)
                ENDIF

                IF ( WHOMOLUMO ) call homolumosave(LIMITS,MESH,BAS,EIGENVECT,EIGENVECT,EHFeigen,EHFeigen,Ne,NATOMS,ATOMS,DFTC)
        ENDIF

        IF ( CORRLEVEL .EQ. 'URHF' .OR. CORRLEVEL .EQ. 'CISD' .OR. CORRLEVEL .EQ. 'MP2' .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC ) THEN
                
                ALLOCATE(EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),P(NB,NB),Pup(NB,NB),Pdown(NB,NB),C1(NB,NB),C2(NB,NB))
                
                IF ( .not. RESTRICT .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC ) THEN
                        Pup = 0.0d0
                        Pdown = 0.0d0
                        IF ( .not. DFTC ) CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE., .TRUE.,.FALSE. )
                        
                        
                        IF ( DFTC ) CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE.,.TRUE.,.FALSE.)

                        IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,Pup+Pdown,S,ZMUL)
                        
                        !====================
                        ! CALCULATING FORCES
                        !==================================================================================================================
                        IF ( CFORCE ) THEN
                                
                                CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,LORDER,CGORDER,LQ,CGQ,force)
                                
                                WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                                WRITE(*,*)' '
                                WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                                DO I=1,NATOMS
                                        WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                                ENDDO
                                        
                        ENDIF
                        !==================================================================================================================

                        IF ( HFORBWRITE ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,Cup,Cdown)
                ENDIF

                IF ( (RESTRICT .AND. CORRLEVEL .EQ. 'CISD') .OR. (RESTRICT .AND.  CORRLEVEL .EQ. 'MP2')  ) THEN
                        
                        ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB))
                        P = 0.0d0
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE., .TRUE.,.FALSE. )
                        
                        IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)

                        !====================
                        ! CALCULATING FORCES
                        !==================================================================================================================
                        IF ( CFORCE ) THEN
                                P = P/2.0d0
                                CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,LORDER,CGORDER,LQ,CGQ,force)
                        
                                WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                                WRITE(*,*)' '
                                WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                                DO I=1,NATOMS
                                        WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                                ENDDO

                        ENDIF
                        !==================================================================================================================

                        IF ( HFORBWRITE ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,EIGENVECT,EIGENVECT)
                        
                        Cup   = EIGENVECT
                        Cdown = EIGENVECT
                        EHFeigenup   = EHFeigen
                        EHFeigendown = EHFeigen
                ENDIF
                
                ! If specified, the charge density is calculated and saved here:
                IF ( WRITEDENS .AND. ( CORRLEVEL .EQ. 'RHF'  .OR. CORRLEVEL .EQ. 'URHF' .OR. DFTC ) ) THEN 
                        C1 = 0.0d0
                        C2 = 0.0d0
                        DO I=1,( Ne - MOD(Ne,2) )/2
                                IF ( .not. DFTC ) THEN
                                        C1(:,I) = Cup(:,I)
                                ELSE
                                        C1(:,I) = Cdown(:,I)
                                ENDIF
                        ENDDO
                        DO I=1,( Ne + MOD(Ne,2) )/2
                                IF (  .not. DFTC  ) THEN
                                        C2(:,I) = Cdown(:,I)
                                ELSE
                                        C2(:,I) = Cup(:,I)
                                ENDIF
                        ENDDO
                        CALL makedens(C1,NB,Pup)
                        CALL makedens(C2,NB,Pdown)
                        P = Pup + Pdown
                        CALL chargedenssave(LIMITS,MESH,BAS,P,Ne)
                ENDIF

                IF ( WHOMOLUMO ) call homolumosave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,DFTC)
                
                
                IF ( CORRLEVEL .EQ. 'CISD' .OR. CORRLEVEL .EQ. 'MP2' ) THEN
                        !print*,'========================================================'
                        !print*,'  Tranfering (ij|kl) from vector to tensor form         '
                        !print*,'========================================================'
                
                        ALLOCATE(Ints(NB,NB,NB,NB))
               
                        DO I=1,NB
                                DO J=1,NB
                                        DO K=1,NB
                                                DO L=1,NB
                                                        Ints(I,J,K,L) = Intsv(ijkl(I,J,K,L))
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDDO
                 ENDIF


                DEALLOCATE(Intsv)
                
                IF ( CORRLEVEL .EQ. 'MP2'  ) CALL MP2(Cup,Cdown,Ints,NB,Ne,ETOT-nucE,nuce,EMP2,EHFeigenup,EHFeigendown,SPINCONSERVE)
                 
                ! Here the CISD calculation starts
                IF ( CORRLEVEL .EQ. 'CISD' ) THEN       
                        print*,' '
                        print*,'  ========================================================'
                        print*,'    Starting Configuration interaction calculation (CISD) '
                        print*,'  ========================================================'
                        
                        CALL CISD(Cup,Cdown,Ints,H0,S,NB,Ne,ETOT-nucE,nuce,ECISD,WRITECICOEF,EHFeigenup,EHFeigendown,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE)
                        
                ENDIF
                
                ! Here the DQMC or the VMC calculation starts
                IF ( CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' ) THEN       
                        IF ( CORRLEVEL .EQ. 'DQMC' ) THEN
                                print*,' '
                                print*,'  =============================================================='
                                print*,'    Starting the Diffusion Quantum Montecarlo Calulation (DQMC) '
                                print*,'  =============================================================='
                        ELSE
                                print*,' '
                                print*,'  =============================================================='
                                print*,'   Starting the Variational Quantum Montecarlo Calulation (VMC) '
                                print*,'  =============================================================='
                        ENDIF
                        
                        CALL dqmc(Ne,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,BETA,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,ETOT-nucE,nucE,ESIGMA,EDQMC,NPERSIST,NRECALC, & 
                        &       CUTTOFFFACTOR,NVMC,VMCCALC,WRITEDENS,MESH,LIMITS)
                        
                ENDIF
        ENDIF
       
       !-----------------------------------------------------
       ! THE TOTAL EXECUTION TIME IS CALCULATED AND PRINTED
       !-----------------------------------------------------

       call DATE_AND_TIME(date, time, zone,STOPTIME)
       print*,' '
       print*,'      ================ TOTAL RUNTIME: ================'
       print*,' '
       RUNTIME = STOPTIME-STARTTIME
       TOTALTIME = 24*3600*RUNTIME(3)+RUNTIME(5)*3600+RUNTIME(6)*60+RUNTIME(7)
       RUNTIME(3) = INT(TOTALTIME/(24*3600))
       TOTALTIME = MOD(TOTALTIME,24*3600*1.0)
       RUNTIME(5) = INT(TOTALTIME/3600)
       TOTALTIME = MOD(TOTALTIME,3600*1.0)
       RUNTIME(6) = INT(TOTALTIME/60)
       TOTALTIME = MOD(TOTALTIME,60*1.0)
       RUNTIME(7) = INT(TOTALTIME)
       WRITE(*,'(A15,I5,A3,I5,A3,I5,A4,I5,A2,I5,A3)')'       RUNTIME=',RUNTIME(3),' d ',RUNTIME(5),' h ',RUNTIME(6),' min ',RUNTIME(7),' s',RUNTIME(8),' ms'
       print*,' '
       print*,'      ================================================'
END PROGRAM uquantchem
