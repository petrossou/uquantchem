PROGRAM uquantchem
      USE datatypemodule
      USE exchcorrmodule
      USE omp_lib
      IMPLICIT NONE
      INTEGER  :: NATOMS
      INTEGER, ALLOCATABLE :: ATOMICNUMBERS(:)
      DOUBLE PRECISION, ALLOCATABLE :: APOS(:,:)
      DOUBLE PRECISION :: Tol,ETOT,NucE,Rn(3)
      CHARACTER(LEN=20) :: CORRLEVEL,EPROFILE
      INTEGER :: BASISMAP(120),BASISMAPAUX(120)
      TYPE(ATOM), ALLOCATABLE :: ATOMS(:),ATOMSAUX(:)
      TYPE(BASIS) :: BAS,BASAUX
      CHARACTER(LEN=6) :: DUMMY
      CHARACTER(Len=20) :: date,time,zone
      INTEGER :: NLINES,I,J,K,L,M,NB,NBAUX,Ne,Lmax,LmaxAUX,MESH(3),NVMC,IOSA,FNATOMS,NSTEPS,NLSPOINTS,PORDER,ZEROSCFTYPE
      INTEGER*8 :: NRED,FNRED
      LOGICAL :: finns,WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,HYBRID,USEGTO,CORRALCUSP,VMCCALC,HFORBWRITE,CFORCE,RELAXN,WRITEONFLY,MOVIE,MOLDYN
      DOUBLE PRECISION, ALLOCATABLE :: S(:,:),T(:,:),V(:,:),H0(:,:),Intsv(:),EIGENVECT(:,:),Ints(:,:,:,:),gradIntsv(:,:,:),SINV(:,:),DMAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gradS(:,:,:,:),gradT(:,:,:,:),gradV(:,:,:,:),PEXu(:,:,:),PEXd(:,:,:),PEXuu(:,:,:),PEXdd(:,:,:),Pupp(:,:),Pdownn(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigenup(:),EHFeigendown(:),Cup(:,:),Cdown(:,:),H00(:,:),DPTENSOR(:,:,:),DIPOLET(:,:),WRI(:,:,:),gradWRI(:,:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: VRI(:,:),gradVRI(:,:,:,:),VRIINV(:,:)
      COMPLEX*16, ALLOCATABLE ::  Cupc(:,:),Cdownc(:,:),Pupc(:,:),Pdownc(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigen(:),P(:,:),Pup(:,:),Pdown(:,:),C1(:,:),C2(:,:),force(:,:),CN(:),LQ(:,:),CGQ(:,:),ZMUL(:),PHOLE(:,:),PEXCITED(:,:),Pgup(:,:),Pgdown(:,:)
      INTEGER, ALLOCATABLE :: Q1(:),Q2(:),Q3(:)
      DOUBLE PRECISION  :: PRYSR(25,25),PRYSW(25,25),rts(25),wts(25),TOTALTIME,ECISD,LIMITS(3),ENEXCM,EMP2,TIMESTEP,TEND,TSTART,DR,FTol,EETOL,CNSTART(4),ETEMP,ENTROPY,mu,ETEMPE,autoev
      DOUBLE PRECISION :: ESIGMA,EDQMC,a,b,BETA,BJASTROW,CJASTROW,CUTTOFFFACTOR,MIX,Aexpan,LAMDA,rc,radius,GAMA,BETAA,POLY(6),deltar,TEMPERATURE,kappa,alpha,alphastart,kappastart,Egs
      DOUBLE PRECISION :: OMEGA,EFIELDMAX,TID,DRF,NORM,DAMPING
      INTEGER :: STARTTIME(8),STOPTIME(8),RUNTIME(8),NEEXC,SAMPLERATE,NREPLICAS,NRECALC,NPERSIST,DIISORD,DIISSTART,NTIMESTEPS,NSCF,PULAY,FIXNSCF,DORDER,TIMEINDEX
      INTEGER :: NLEBEDEV,NCHEBGAUSS,LEBPOINTS(24),CGORDER,LORDER,NTOTALQUAD,RELALGO,IORBNR(2),AORBS,DIISORDEX,EDIR,NEPERIOD,II
      INTEGER, EXTERNAL :: ijkl
      DOUBLE PRECISION, EXTERNAL :: gto,massa
      LOGICAL :: RESTART,ZEROSCF,XLBOMD,DFTC,DOTDFT,ADEF,DOABSSPECTRUM,DIFFDENS,AFORCE,RIAPPROX,SCRATCH,DIAGDG,FIELDREAD
      INTEGER :: NSCCORR,LIMPRECALC
      DOUBLE PRECISION :: MIXTDDFT,SCERR,FIELDDIR(3),TRANSCOORD(3,3),cosTH,sinTH,cosFI,sinFI,RIKTNING(3)
      !-------------------
      !STARTING THE CLOCK: 
      !-------------------
      
      SCRATCH = .TRUE.
      AFORCE = .TRUE.
      DRF = 0.00010d0

      call DATE_AND_TIME(date, time, zone,STARTTIME)
     
      DFTC = .FALSE.

      autoev = 27.211383860d0
      !=======================================================================
      ! (1) HERE THE DATA FROM THE INPUTFILES: INPUTFILE & BASISFILE IS LOADED
      !=======================================================================

      ! preparing the readin of data by aquiring the number of atoms, NATOMS,
      ! and the number of lines, NLINES, appearing before the entry of the
      ! characters 'NATOMS' in the file 'INPUTFILE'
        
      CALL checkinputfile(NLINES,NATOMS)
      
      ALLOCATE(ATOMICNUMBERS(NATOMS),APOS(NATOMS,3),ATOMS(NATOMS),ATOMSAUX(NATOMS),force(NATOMS,3),ZMUL(NATOMS))
      
      CALL readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBERS,APOS,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, &
      & SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,HYBRID,rc,CORRALCUSP,NVMC,HFORBWRITE, &
      & IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol,NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS, &
      & EETOL,RELALGO,ZEROSCFTYPE,ETEMP,IORBNR,AORBS,DIISORDEX,DOTDFT,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF,DOABSSPECTRUM,DIFFDENS,AFORCE,DRF,&
      & NSCCORR,MIXTDDFT,SCERR,RIAPPROX,LIMPRECALC,DIAGDG,FIELDDIR,FIELDREAD,DAMPING)

      !=======================================================================
      ! If the direction of the electric field has been set explicitly by the
      ! user and a TDDFT/THF or ADEF-calculation is to be performed, a change 
      ! of coordinate system is performed so that the z-axis of the new
      ! coordinate system is paralell with the E-field. The x-direction of the 
      ! new coordinate system will be ( cos(FI), -sin(FI) ) and the y-direction
      ! ( cos(TH)*sin(FI), cos(TH)*cos/(FI), -sin(TH) ) where:
      ! TH = arccos(FIELDDIR(3)) and FI = arctan(FIELDDIR(1)/FIELDDIR(2))
      !========================================================================

      IF ( FIELDREAD .AND. ( DOTDFT .OR. ADEF ) ) THEN
         IF ( FIELDDIR(1) .NE. 0.0d0 .OR. FIELDDIR(2) .NE. 0.0d0 .OR.  FIELDDIR(3) .NE. 0.0d0 ) THEN
                RIKTNING(1) = FIELDDIR(1)/sqrt(DOT_PRODUCT(FIELDDIR,FIELDDIR))
                RIKTNING(2) = FIELDDIR(2)/sqrt(DOT_PRODUCT(FIELDDIR,FIELDDIR))
                RIKTNING(3) = FIELDDIR(3)/sqrt(DOT_PRODUCT(FIELDDIR,FIELDDIR))
                cosTH = RIKTNING(3)
                sinTH = sqrt( 1.0d0 - cosTH**2)
                IF ( sinTH .NE. 0.0d0 ) THEN
                        sinFI = RIKTNING(1)/sinTH
                        cosFI = RIKTNING(2)/sinTH
                        ! constructing the coordinate transformation matrix:
                        TRANSCOORD(1,:) =(/     cosFI,       -sinFI,    0.0d0 /)
                        TRANSCOORD(2,:) =(/ cosTH*sinFI,  cosTH*cosFI, -sinTH /)
                        TRANSCOORD(3,:) =(/ sinTH*sinFI,  sinTH*cosFI,  cosTH /)
                        ! Performing the change of coordinate system:
                        DO I=1,NATOMS
                                APOS(I,:) = MATMUL(TRANSCOORD,APOS(I,:))
                        ENDDO
                ELSE
                        EDIR = 2
                ENDIF
         ELSE 
                WRITE(*,*)' '
                WRITE(*,*)'****************************************************************************'
                WRITE(*,*)'***                          WARNING!                                    ***'
                WRITE(*,*)'      The direction of the electric field has been set to (0, 0, 0 ).'
                WRITE(*,*)'     Since this do not make sense the field direction is set to (0,0,1)'
                WRITE(*,*)'and in the case of an AC-field the direction of propagation is set to (0,1,0)'
                WRITE(*,*)'****************************************************************************'
                WRITE(*,*)' '
                EDIR = 2
         ENDIF
      ENDIF

      IF ( .not. AFORCE .AND. ( MOLDYN .OR. RELAXN .OR. XLBOMD ) ) THEN
                AFORCE = .TRUE.
                WRITE(*,*)' '
                WRITE(*,*)'*******************************************************************************************************************************'
                WRITE(*,*)'                                                    WARNING!                                                                   '
                WRITE(*,*)' AFORCE = .FALSE is not allowed for molecular dynamics calculations (MOLDYN=.TRUE.) or relaxation calculations (RELAXN=.TRUE.) '
                WRITE(*,*)'                          AFORCE is thus set to TRUE, and analytical forces are used.                                          '
                WRITE(*,*)'*******************************************************************************************************************************'
                WRITE(*,*)' '
      ENDIF

      IF ( ADEF .AND. DOTDFT ) THEN
              WRITE(*,*)' '
              WRITE(*,*)'************************************************************************************************************************'
              WRITE(*,*)'                                    WARNING CALCULATION STOPPED, REASON:                                                '
              WRITE(*,*)'Time dependent DFT/HF (DOTDFT=.TRUE.) option cannot be combined with adiabatic electric field calculations (ADEF=.TRUE.)'
              WRITE(*,*)'      You can only choose one of the two above options to be true, i.e either DOTDFT=.TRUE. excusive or ADEF=.TRUE.     '
              WRITE(*,*)'************************************************************************************************************************'
              WRITE(*,*)' '
              STOP
      ENDIF

      IF ( ADEF ) THEN
              WRITE(*,*)' '
              WRITE(*,*)'=========================================================================================================================='
              WRITE(*,*)'                             CALCULATION PERFORMED WITH THE PRESENCE OF ELECTRIC FIELD                                    '
              IF ( XLBOMD .OR. MOLDYN .AND. EPROFILE .NE. 'ST' ) THEN
                        WRITE(*,*)'                    THE ELECTRONS WILL BE ASSUMED TO FOLLOW THE TIME-DEPENDENT FIELD ADIABATICALLY                        '
              ENDIF
              WRITE(*,*)'=========================================================================================================================='
              WRITE(*,*)' ' 
      ENDIF

      IF ( DORDER .LT. 4 )  DORDER = 4
      ALLOCATE(CN(DORDER))
       
      IF ( RELALGO .NE. 1 .AND. RELALGO .NE. 2 ) RELALGO = 1

      IF ( CORRLEVEL .EQ. 'RHF' .AND. DOTDFT ) THEN
               WRITE(*,*)'=================================================================='
               WRITE(*,*)'Restricted time dependent Hartree-Fock caclulation not implemented'
               WRITE(*,*)'Changing CORRLEVEL from "RHF" to "URHF"'
               WRITE(*,*)'=================================================================='
               CORRLEVEL = 'URHF'
      ENDIF

      ! The array containing the number of Lebedev quadrature points is assigned here:
      LEBPOINTS = (/ 110, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)

      ! Here the arrays containg the Lebedev and Chebyshev-Gauss quadrature
      ! points and abssiscae are allocated and created in the case that one has
      ! choosen to perform a DFT calculation
      IF ( CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'B3LYP') THEN
              DFTC = .TRUE.
              HYBRID = .FALSE.
              CGORDER = NCHEBGAUSS 
              LORDER  = LEBPOINTS(NLEBEDEV)
              ALLOCATE(LQ(LORDER,3),CGQ(CGORDER,2))
              IF ( NLEBEDEV .GE. 1  .AND. NLEBEDEV .LE. 5   ) CALL lebedev1(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 6  .AND. NLEBEDEV .LE. 10  ) CALL lebedev2(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 11 .AND. NLEBEDEV .LE. 15  ) CALL lebedev3(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 16 .AND. NLEBEDEV .LE. 20  ) CALL lebedev4(NLEBEDEV,LORDER,LQ)
              IF ( NLEBEDEV .GE. 21 .AND. NLEBEDEV .LE. 24  ) CALL lebedev5(NLEBEDEV,LORDER,LQ) 
              CALL chebgauss(CGORDER,CGQ)
              NTOTALQUAD = CGORDER*LORDER*NATOMS
              ALLOCATE(Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD))
              L = 0
              DO I=1,CGORDER
                DO J=1,LORDER
                   DO K=1,NATOMS
                        L = L+1
                        Q1(L) = I
                        Q2(L) = J
                        Q3(L) = K
                   ENDDO
                ENDDO
              ENDDO
      ELSE
              CGORDER = 1
              LORDER = 1
              ALLOCATE(LQ(LORDER,3),CGQ(CGORDER,2))
              NTOTALQUAD = 1
              ALLOCATE(Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD))
              Q1 = 0.0d0
              Q2 = 0.0d0
              Q3 = 0.0d0
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
        
      IF ( XLBOMD ) MOLDYN = .TRUE.

      IF ( EPROFILE .NE. 'ST' .AND. .not. DOTDFT .AND. .not. MOLDYN .AND. ADEF ) THEN
             EPROFILE = 'ST'
             WRITE(*,*)' '
             WRITE(*,*)'TIME-INDEPENDENT CALCULATIONS ONLY PERMITTED WITH STATIC HOMOGENEOUS FIELD'
             WRITE(*,*)'EPROFILE has been changed to ST'
             WRITE(*,*)' '
      ENDIF

      !IF ( (RELAXN .OR. MOLDYN ) .AND. CORRLEVEL .NE. 'RHF' .AND. CORRLEVEL .NE. 'URHF'  ) THEN
      !        WRITE(*,*)'Relaxation of nuclea/MD only possible for RHF or URHF'
      !        WRITE(*,*)'Set CORRLEVEL to RHF/URHF'
      !        STOP
      !ENDIF
      
      
      ! Ionization hole calculations only work at zero-electronic temperature
      ! and with the ZEROSCFTYPE = 1 alternative when calculating forces and 
      ! energies in a md-calculation
      IF (IORBNR(1) .NE. 0  ) THEN
              ZEROSCFTYPE = 1
              !ETEMP = -1000.0d0
      ENDIF

      IF ( RELAXN .OR. MOLDYN ) CFORCE = .TRUE.
              
      ! Ensuring that idempotency is not used when calculating the pulay
      ! contribution to the force when non-integere occupation of electronic
      ! states is used.
      IF ( ETEMP .GT. 0.0d0 .AND. CFORCE ) PULAY = 2
      IF ( ZEROSCFTYPE .EQ. 2  ) PULAY = 2
      
      IF ( NATOMS .EQ. 1 ) CFORCE = .FALSE.
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
        !------------------------------------
        ! same thing for the auxiliary basis:
        !------------------------------------
        ATOMSAUX(I)%Z = ATOMICNUMBERS(I)
        ATOMSAUX(I)%M = massa(ATOMICNUMBERS(I))
        ATOMSAUX(I)%R = APOS(I,:)
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
                        ATOMS(I)%PSI(J)%HYBRID = .FALSE.
                        call stoexponent(ATOMS(I)%PSI(J),deltar,LAMDA,Aexpan,GAMA,BETAA,POLY,rc,USEGTO,CORRALCUSP)
                        ATOMS(I)%PSI(J)%STOEXPON = LAMDA
                        ATOMS(I)%PSI(J)%A = Aexpan
                        ATOMS(I)%PSI(J)%RC = rc
                        ATOMS(I)%PSI(J)%P = POLY
                        ATOMS(I)%PSI(J)%DR = deltar
                        ATOMS(I)%PSI(J)%GAMA = GAMA
                        ATOMS(I)%PSI(J)%BETA = BETAA
                        !print*,I,J,Aexpan,LAMDA,BETAA,GAMA,rc
                        IF ( .not. USEGTO .AND. HYBRID ) THEN
                                ATOMS(I)%PSI(J)%HYBRID = .TRUE.
                        ELSE
                                ATOMS(I)%PSI(J)%HYBRID = .FALSE.
                        ENDIF
                ENDDO
            ENDDO
            !-------------------------------------------
            ! DEBUGG PURPOSES:
            !-------------------------------------------
            !OPEN(9999,FILE='STOCOMPAIRE.dat',ACTION='WRITE')
            !DO I=1,10000
            !    radius = 0.10d0*I/10000
            !    WRITE(9999,*)radius,ATOMS(1)%PSI(1)%A*EXP(-ATOMS(1)%PSI(1)%STOEXPON*radius)+ATOMS(1)%PSI(1)%BETA*radius+ATOMS(1)%PSI(1)%GAMA,gto(ATOMS(1)%PSI(1),radius)
            !ENDDO
            !close(9999)

             
      ELSE
             print*,'ABORTING SINCE THE FILE  "BASISFILE" IS MISSING'
             STOP
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
                WRITE(*,*)'                    WARNING!                               '
                WRITE(*,*)'      THE INDEX OF THE HF-EIGEN FUNCTION TO BE SAVED       '
                WRITE(*,*)'EXCEEDS THE TOTAL NUMBER OF BASIS FUNCTIONS NB =',NB
                WRITE(*,*)'THE INDEX HF-EIGENFUNCTION SAVED WILL THERFORE BE SET TO NB'
                IOSA = NB
        ENDIF
        IF (HFORBWRITE .AND. IOSA .EQ. 0 ) THEN
                IOSA = (Ne-MOD(Ne,2))/2 + MOD(Ne,2)
        ENDIF
        WRITE(*,*)' '
        WRITE(*,'(A50,I4)')'        Total number of basis-functions used NB = ',NB

        ! Here the total molecular basis set is put together
        CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
        
        !=======================================================================
        ! (3) HERE THE OVERLAP-, KINETIC- AND POTENTIAL ENERGY- MATRICES ARE
        !     CALCULATED.
        !=======================================================================
        
        ! Here the normalization of the basis functions is performed

        CALL normalize(BAS)

        ! Here the overlap matrix is calculated 

        ALLOCATE(S(NB,NB),SINV(NB,NB),gradS(NATOMS,3,NB,NB),DMAT(NB,NB),PEXu(24,NB,NB),PEXd(24,NB,NB))
        ALLOCATE(PEXuu(24,NB,NB),PEXdd(24,NB,NB))

        CALL overlap(NATOMS,BAS,S,gradS)
        
        ! Here the kinetic energy matrix is calculated
        ALLOCATE(T(NB,NB),gradT(NATOMS,3,NB,NB))

        CALL kinetic(NATOMS,BAS,T,gradT)
        
        ! Here the potential energy matrix is calculated

        ALLOCATE(V(NB,NB),gradV(NATOMS,3,NB,NB))

        CALL potential(BAS,NATOMS,ATOMS,V,gradV)

        !================================================================================================
        ! If the auxiliary basis-file, BASISFILEAUX, is present in the run-directory and the user has not
        ! explicitly stated that he/she does not want to utilize the Resoulution ! of Identity approximation 
        ! for the (ij|kl)-integrals, the auxiliary basis will be set up here in the exact same manner as 
        ! the standard basis-set
        !================================================================================================
        
      NBAUX = 1

      inquire(file='BASISFILEAUX',exist=finns)
      IF ( finns .AND.  RIAPPROX ) THEN
             WRITE(*,*)' '
             WRITE(*,'(A61)')'    ========================================================='
             WRITE(*,'(A60)')'     Using the Resolution of the Identity (RI) approximation'
             WRITE(*,'(A61)')'    ========================================================='
             WRITE(*,*)' '
             BASISMAPAUX(:) = 0
             ! Here the number of basis functions of all the atoms in the 
             ! file "BASISFILE" are loaded into the array BASISMAP
             
             CALL readbasisaux(NATOMS,ATOMSAUX,BASISMAPAUX,.TRUE.)

             DO I=1,NATOMS
                ALLOCATE(ATOMSAUX(I)%PSI(BASISMAPAUX(ATOMS(I)%Z)))
             ENDDO

             ! Here the local basis for each atom is loaded
             
             CALL readbasisaux(NATOMS,ATOMSAUX,BASISMAPAUX,.FALSE.)

             ! CALCULATING THE STO-EXPONENTS AND NORMS
             ! SAVING AND SAVING IT IN THE HYBRID BASIS
             
             deltar = 0.0010d0

             DO I= 1,NATOMS
                DO J=1,ATOMSAUX(I)%NBAS
                        ATOMSAUX(I)%PSI(J)%HYBRID = .FALSE.
                        call stoexponent(ATOMSAUX(I)%PSI(J),deltar,LAMDA,Aexpan,GAMA,BETAA,POLY,rc,USEGTO,CORRALCUSP)
                        ATOMSAUX(I)%PSI(J)%STOEXPON = LAMDA
                        ATOMSAUX(I)%PSI(J)%A = Aexpan
                        ATOMSAUX(I)%PSI(J)%RC = rc
                        ATOMSAUX(I)%PSI(J)%P = POLY
                        ATOMSAUX(I)%PSI(J)%DR = deltar
                        ATOMSAUX(I)%PSI(J)%GAMA = GAMA
                        ATOMSAUX(I)%PSI(J)%BETA = BETAA
                ENDDO
            ENDDO

           !===========================================================================
           ! (2) HERE THE MOLECULAR BASIS SET IS CONSTRUCTED FROM ALL THE LOCAL ATOMIC
           !     BASISSETS, BY ADDING ALL THE L-SUBSHELLS AND CONCATANATING THE
           !     LOCAL ATOMIC BASES.
           !===========================================================================
        
           ! Here the total number of aux-basis functions is counted
           CALL gettotalbasis(NATOMS,ATOMSAUX,NBAUX,BASAUX,.TRUE.)

           BASAUX%NBAS = NBAUX
           ALLOCATE(BASAUX%PSI(NBAUX))
           WRITE(*,*)' '
           WRITE(*,'(A59,I4)')'   Total number of auxiliary basis-functions used NBAUX = ',NBAUX
           WRITE(*,*)' '

           ! Here the total aux-molecular basis set is put together
           CALL gettotalbasis(NATOMS,ATOMSAUX,NBAUX,BASAUX,.FALSE.)

           !=======================================================================
           ! (3) HERE THE OVERLAP-, KINETIC- AND POTENTIAL ENERGY- MATRICES ARE
           !     CALCULATED.
           !=======================================================================
        
           ! Here the normalization of the basis functions is performed

           CALL normalize(BASAUX)

        ELSE 
                IF ( RIAPPROX ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A72)')'************************************************************************'
                        WRITE(*,'(A72)')'*****                        WARNING !                             *****'
                        WRITE(*,'(A72)')'************************************************************************'
                        WRITE(*,'(A72)')' You have explicitly stated that the RI-approximation should be used.   '
                        WRITE(*,'(A72)')'    However, since the the BASISFILEAUX-file is not present in the      '
                        WRITE(*,'(A72)')' run-directory, the RI-approximation will not be used. This will result '
                        WRITE(*,'(A72)')' in a much slower calculation than expected and the over-use of memory. '
                        WRITE(*,'(A72)')'************************************************************************'
                        WRITE(*,*)' '
                ENDIF
                RIAPPROX = .FALSE.
        ENDIF
123 CONTINUE
        IF ( CORRLEVEL .EQ. 'MP2' .AND. RIAPPROX ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A62)')'**************************************************************'
                WRITE(*,'(A62)')'***                       WARNING!                         ***'
                WRITE(*,'(A62)')' RI-approximation in MP2 calculations are not yet implemented.'
                WRITE(*,'(A62)')'Switching off RI-approximation, i.e setting RIAPPROX = .FALSE.'
                WRITE(*,'(A62)')'**************************************************************'
                WRITE(*,*)' '
                RIAPPROX = .FALSE.
        ENDIF
        IF ( CORRLEVEL .EQ. 'CISD' .AND. RIAPPROX ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A61)')'**************************************************************'
                WRITE(*,'(A61)')'***                      WARNING!                          ***'
                WRITE(*,'(A61)')'RI-approximation in CISD calculations are not yet implemented.'
                WRITE(*,'(A61)')'Switching off RI-approximation, i.e setting RIAPPROX = .FALSE.'
                WRITE(*,'(A61)')'**************************************************************'
                WRITE(*,*)' '
                RIAPPROX = .FALSE.
        ENDIF
        ! If a DFT type of calculation is performed the 
        ! exchange correlation matrix elements and their 
        ! nuclear gradients are allocated 
        
        !IF ( DFTC ) ALLOCATE(Vxc(2,NB,NB),gVxc(NATOMS,2,3,NB,NB))

        !===========================================================================
        ! (4) Here the electron-electron repulsion 'tensor', Intsv(:), is calulated
        !===========================================================================
      
        ! The number of symmetry reduced electron-electron tensor elements

        !NRED = ( (((NB+1)*NB)/2)**2 + ((NB+1)*NB)/2 )/2
        IF ( .not. RIAPPROX .OR. ( RIAPPROX .AND. NB .LE. LIMPRECALC ) ) THEN
                NRED = (NB+1)*NB
                NRED = NRED/2
                NRED = (NRED*NRED + ((NB+1)*NB)/2 )/2
        ELSE
                NRED = 1
        ENDIF

        ALLOCATE(Intsv(NRED))

        IF ( CFORCE .AND. AFORCE ) THEN
                FNRED = NRED
                FNATOMS = NATOMS
        ELSE
                FNRED =1
                FNATOMS =1
        ENDIF
       
        IF ( .not. RIAPPROX ) THEN 
                ALLOCATE(gradIntsv(FNATOMS,3,FNRED))
        ELSE
                IF ( NB .LE. LIMPRECALC ) THEN
                        ALLOCATE(gradIntsv(FNATOMS,3,FNRED))
                ELSE
                        ALLOCATE(gradIntsv(1,3,1))
                ENDIF
        ENDIF
                        
                        

        
        ! Here we pre calculate the roots and weights of all the rys polynomials 
        ! up to 2*Lmax + 1, (Lmax = maximum L guanum number of basisset). Since 
        ! these roots and weights correspnd to the calculation of (ij|kl) when 
        ! all four orbitals i,j,k,l are located on the same cite, i.e X = 0, 
        ! this will save a lot of time in the calculation of (ij|kl) especially 
        ! for small molecules

        ! Determining Lmax for the standard basis
        Lmax = 0
        
        DO I=1,NATOMS
                DO J=1,ATOMS(I)%NBAS
                        IF ( Lmax .LT. ATOMS(I)%PSI(J)%L(1) ) Lmax = ATOMS(I)%PSI(J)%L(1)
                ENDDO
        ENDDO

        ! Determining Lmax for the auxiliary basis
        IF ( RIAPPROX ) THEN
                LmaxAUX = 0
                DO I=1,NATOMS
                        DO J=1,ATOMSAUX(I)%NBAS
                                IF ( Lmax .LT. ATOMSAUX(I)%PSI(J)%L(1) ) LmaxAUX = ATOMSAUX(I)%PSI(J)%L(1)
                        ENDDO
                ENDDO
                IF ( LmaxAUX .GT. Lmax ) Lmax = LmaxAUX
        ENDIF

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

        ALLOCATE(WRI(BAS%NBAS,BAS%NBAS,NBAUX),gradWRI(FNATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX))
        ALLOCATE(VRI(NBAUX,NBAUX),gradVRI(FNATOMS,3,NBAUX,NBAUX),VRIINV(NBAUX,NBAUX))
        IF ( RIAPPROX ) THEN
                CALL calcWRI(NATOMS,FNATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,CFORCE)
                CALL calcVRI(NATOMS,FNATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,CFORCE)
                CALL invert(VRI,VRIINV,BASAUX%NBAS)
                IF ( NB .LE. LIMPRECALC ) CALL eeintsRI(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,EETOL,CFORCE,DMAT,NBAUX,VRI,WRI,gradVRI,gradWRI)
        ELSE
                IF ( AFORCE ) THEN 
                        CALL eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,EETOL,CFORCE,DMAT)
                ELSE
                        CALL eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,EETOL,.FALSE.,DMAT)
                ENDIF
        ENDIF
        
        IF ( NB .LE. LIMPRECALC .AND. .not. MOLDYN .AND. .not. RELAXN ) RIAPPROX = .FALSE.

        ALLOCATE(H0(NB,NB),PHOLE(NB,NB),PEXCITED(NB,NB),H00(NB,NB),DPTENSOR(3,NB,NB),DIPOLET(NB,NB))
        
        H0 = T+V
        
        IF ( ADEF .AND. EPROFILE .EQ. 'DP' ) THEN
                WRITE(*,*)' '
                WRITE(*,*)'WARNING ADIABATIC ELECRIC FIELD CALCULATIONS ARE NOT PERMITTED WITH DIRAC PULSE ENVELOPE' 
                WRITE(*,*)'CHANGING EPROFILE FROM DP TO ST'
                WRITE(*,*)' '
                EPROFILE = 'ST'
        ENDIF

        IF ( ADEF .AND. .not. DOTDFT ) THEN
                IF ( RELAXN .OR. EPROFILE .EQ. 'ST' ) THEN
                        CALL dipoletensor(BAS,DPTENSOR)
                        IF ( EDIR .EQ. 1 ) DIPOLET = DPTENSOR(2,:,:)
                        IF ( EDIR .EQ. 2 ) DIPOLET = DPTENSOR(3,:,:)
                        IF ( EDIR .EQ. 3 ) DIPOLET = DPTENSOR(1,:,:)
             
                        H0 = H0 + EFIELDMAX*DIPOLET

                        DO II=1,NATOMS
                                IF ( EDIR .EQ. 1 ) THEN
                                        NucE = NucE - EFIELDMAX*ATOMS(II)%Z*ATOMS(II)%R(2)
                                ELSE IF ( EDIR .EQ. 2 ) THEN
                                        NucE = NucE - EFIELDMAX*ATOMS(II)%Z*ATOMS(II)%R(3)
                                ELSE IF ( EDIR .EQ. 3 ) THEN
                                        NucE = NucE - EFIELDMAX*ATOMS(II)%Z*ATOMS(II)%R(1)
                                ENDIF
                        ENDDO
                ENDIF
        ENDIF
        
        !==============================================================================
        ! Here the starting-guess density matrices are put together from
        ! density matrices obtained from previous atomic calculations. The
        ! atomic density matrices are read from the file 'DENSMATSTARTGUESS.dat'
        !==============================================================================
        ALLOCATE(Pup(NB,NB),Pdown(NB,NB))
        IF ( NATOMS .GT. 1 ) THEN
                CALL Rdensmatstartguess(ATOMS,NATOMS,NB,Pup,Pdown,SCRATCH)
                NORM = SUM(S*(Pup+Pdown))
                IF ( NORM .NE. 0.0d0 .AND. .not. SCRATCH ) THEN
                        Pup = 0.50d0*(Ne/NORM)*(Pup+Pdown)
                        Pdown = Pup
                ENDIF
        ENDIF

        !==========================================================================
        ! Here a slightly naive starting-guess for the SCF-calculations is
        ! created, by assuming a diagonal density matrix for the first iteration.
        !==========================================================================
        inquire(file='DENSMATSTARTGUESS.dat',exist=finns)
        IF ( (DIAGDG .AND. .not. finns) .OR. (DIAGDG .AND. NATOMS .EQ. 1) ) THEN
                Pup = 0.0d0
                Pdown = 0.0d0
                DO I=1,NB
                        Pup(I,I) = 1.0d0
                ENDDO
                NORM = SUM(S*Pup)
                Pup = 0.50d0*(Ne/NORM)*Pup
                Pdown = Pup
                SCRATCH = .FALSE.
        ENDIF

        IF ( RELAXN .AND. .not. MOLDYN ) THEN
                IF ( EPROFILE .NE. 'ST' .AND. ADEF ) THEN
                        EPROFILE = 'ST'
                        WRITE(*,*)' '
                        WRITE(*,*)'WARNING RELAXATION ONLY PERMITTED FOR STATIC HOMOGENEOUS FIELD'
                        WRITE(*,*)'EPROFILE has been changed to ST'
                        WRITE(*,*)' '
                ENDIF
                ! Relxing by searching for minimum energy
                IF ( RELALGO .EQ. 2 ) THEN
                        CALL relax(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS, &
                                              & APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY,NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,&
                                              & ETEMP,EDIR,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX,LIMPRECALC)
                ENDIF
                ! Relaxing by searching for FORCES = 0
                IF ( RELALGO .EQ. 1 ) THEN 
                        CALL relaxf(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS, &
                                              & APPROXEE,CORRLEVEL,PRYSR,PRYSW,PULAY,NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,&
                                              & ETEMP,EDIR,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX,LIMPRECALC)
                ENDIF
                GOTO 4000
        ENDIF

        IF ( MOLDYN ) THEN
                print*,'     ========================================================'
                print*,'           Comencing with molecular dynamics calculation     '
                print*,'     ========================================================'
                
                NTIMESTEPS = INT(TEND/TIMESTEP)
                CALL moleculardynamics(gradS,gradT,gradV,gradIntsv,S,H0,Intsv,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISORDEX,DIISSTART,NATOMS,NTIMESTEPS,TIMESTEP,BAS,ATOMS,APPROXEE,&
                & CORRLEVEL,PRYSR,PRYSW,WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,&
                & ZEROSCFTYPE,ETEMP,IORBNR,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX,LIMPRECALC,DAMPING)
        ENDIF
        
        !==============================================================
        ! (5) Here the Unrestricted or Restricted Hartree-FocK
        !     equations are selfconsistently solved
        !==============================================================
        
        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                
                ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB),P(NB,NB),C1(NB,NB))
                P = 0.0d0
                IF ( .not. SCRATCH ) P = Pup
                CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE., SCRATCH,.FALSE. ,NBAUX,VRI,WRI,RIAPPROX)
                
                IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)
                !====================
                ! CALCULATING FORCES
                !================================================================================================================== 
                IF ( CFORCE ) THEN
                        P = P/2.0d0
                        IF ( AFORCE ) THEN
                                CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY, & 
                                & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
                        ELSE
                                CALL numforce(S,H0,Intsv,NB,NRED,Ne,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,DRF,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                                & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,NEPERIOD,EDIR,EPROFILE,EFIELDMAX,&
                                & OMEGA,0.0d0,ADEF,Pup,Pdown,.TRUE.,NBAUX,ATOMSAUX,BASAUX,VRI,WRI,RIAPPROX,force)
                        ENDIF
                        WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                        DO I=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                        ENDDO
                ENDIF
                
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
                
                ALLOCATE(EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),P(NB,NB),C1(NB,NB),C2(NB,NB),Pgup(NB,NB),Pgdown(NB,NB))
                ALLOCATE(Cupc(NB,NB),Cdownc(NB,NB),Pupc(NB,NB),Pdownc(NB,NB),Pupp(NB,NB),Pdownn(NB,NB))
                
                IF ( .not. RESTRICT .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC ) THEN
                        IF ( SCRATCH ) THEN
                                Pup = 0.0d0
                                Pdown = 0.0d0
                        ENDIF
                        IF ( .not. DFTC ) THEN
                                CALL URHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE.,SCRATCH,.FALSE. &
                                &,ETEMP,ENTROPY,NBAUX,VRI,WRI,RIAPPROX)
                                IF ( DOTDFT ) THEN
                                        ETEMPE = ETEMP
                                        Pupc = Pup
                                        Pdownc = Pdown
                                        Cupc = Cup
                                        Cdownc = Cdown
                                        NTIMESTEPS = INT(TEND/TIMESTEP)
                                        CALL TDFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                        & ETOT,Cupc,Cdownc,Pupc,Pdownc,.TRUE.,ETEMPE,mu,OMEGA,EDIR,NEPERIOD,EPROFILE,NTIMESTEPS,TIMESTEP,EFIELDMAX,&
                                        & IORBNR,PEXu,PEXd,PEXuu,PEXdd,DOABSSPECTRUM,DIFFDENS,NSCCORR,MIXTDDFT,SCERR,NBAUX,VRI,WRI,RIAPPROX)
                                        !==========================================================================================
                                        ! Saving the chargensity at time index = TIMESTEP/10,2TIMESTEP/10,3TIMESTEP/10,...,TIMESTEP
                                        !==========================================================================================
                                        IF (  WRITEDENS ) THEN
                                                DO TIMEINDEX=1,24
                                                        Pup = PEXu(TIMEINDEX,:,:)
                                                        Pdown = PEXd(TIMEINDEX,:,:)
                                                        Pupp = PEXuu(TIMEINDEX,:,:)
                                                        Pdownn = PEXdd(TIMEINDEX,:,:)
                                                        CALL exciteddenssave(LIMITS,MESH,BAS,Pup,Pdown,Pupp,Pdownn,NATOMS,ATOMS,TIMEINDEX,DIFFDENS)
                                                ENDDO
                                        ENDIF
                                        IF ( AORBS .NE. 0 ) CALL exitedarbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS,NTIMESTEPS)
                                        GOTO 4000
                                ENDIF
                        ENDIF
                        !IF ( IORBNR .NE. 0 ) THEN
                        !        ETEMPE = -1000.0d0
                        !ELSE
                                ETEMPE = ETEMP
                        !ENDIF

                        IF ( DFTC ) THEN    
                                CALL DFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE.,SCRATCH,.FALSE.,ETEMPE,mu,ENTROPY,NBAUX,VRI,WRI,RIAPPROX)
                                IF ( DOTDFT ) THEN
                                        Pupc = Pup
                                        Pdownc = Pdown
                                        Cupc = Cup
                                        Cdownc = Cdown
                                        NTIMESTEPS = INT(TEND/TIMESTEP)
                                        CALL TDFT(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                        & ETOT,Cupc,Cdownc,Pupc,Pdownc,.TRUE.,ETEMPE,mu,OMEGA,EDIR,NEPERIOD,EPROFILE,NTIMESTEPS,TIMESTEP,EFIELDMAX,&
                                        & IORBNR,PEXu,PEXd,PEXuu,PEXdd,DOABSSPECTRUM,DIFFDENS,NSCCORR,MIXTDDFT,SCERR,NBAUX,VRI,WRI,RIAPPROX)
                                        IF (  WRITEDENS ) THEN
                                                DO TIMEINDEX=1,24
                                                        Pup = PEXu(TIMEINDEX,:,:)
                                                        Pdown = PEXd(TIMEINDEX,:,:)
                                                        Pupp = PEXuu(TIMEINDEX,:,:)
                                                        Pdownn = PEXdd(TIMEINDEX,:,:)
                                                        CALL exciteddenssave(LIMITS,MESH,BAS,Pup,Pdown,Pupp,Pdownn,NATOMS,ATOMS,TIMEINDEX,DIFFDENS)
                                                ENDDO
                                        ENDIF
                                        IF ( AORBS .NE. 0 ) CALL exitedarbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS,NTIMESTEPS)
                                        GOTO 4000
                                ENDIF
                                CALL PRINTENERGYEIGENVAL(BAS,EHFeigenup,EHFeigendown,Cup,Cdown,.FALSE.)
                                CALL PRINTDIPOLETENSOR(BAS,Cup,Cdown)
                        ENDIF
                        
                        IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,Pup+Pdown,S,ZMUL)

                        IF ( IORBNR(1) .NE. 0 .AND. DFTC ) THEN
                                IF ( ABS(IORBNR(1)) .GT. NB ) THEN
                                        IORBNR(1) = NB
                                        WRITE(*,'(A70)')'**********************************************************************'
                                        WRITE(*,'(A27,I4,A7,I4,A40)')'        WARNING! IORBNR(1)=', IORBNR(1), '> NB = ',NB,' (NB = Number of basis functions)       '
                                        WRITE(*,'(A70)')'                    Setting IORBNR(1) = NB                            '
                                        WRITE(*,'(A70)')'**********************************************************************'
                                ENDIF
                                IF ( ABS(IORBNR(2)) .GT. NB ) THEN
                                        IORBNR(2) = NB
                                        WRITE(*,'(A70)')'**********************************************************************'
                                        WRITE(*,'(A27,I4,A7,I4,A40)')'        WARNING! IORBNR(2)=', IORBNR(2), '> NB = ',NB,' (NB = Number of basis functions)       '
                                        WRITE(*,'(A70)')'                    Setting IORBNR(2) = NB                            '
                                        WRITE(*,'(A70)')'**********************************************************************'
                                ENDIF
                                WRITE(*,'(A70)')'======================================================================'
                                WRITE(*,'(A70)')'             Performing Core-hole excitation calculation              '
                                WRITE(*,'(A70)')'======================================================================'
                                IF ( IORBNR(2) .EQ. 0 ) THEN
                                        WRITE(*,'(A60,I4)')' Orbital being ionized and treated as a hole is orbital nr =',IORBNR
                                ELSE
                                         WRITE(*,'(A46,I4)')'         Orbital being excited is orbital nr =',IORBNR(1)
                                         WRITE(*,'(A46,I4)')'     Orbital is being promoted to orbital nr =',IORBNR(2)
                                ENDIF
                                
                                Egs = ETOT
                                Pgup = Pup
                                Pgdown = Pdown
                                CALL DFTCHI(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
                                & ETOT,Cup,Cdown,Pup,Pdown,Pgup,Pgdown,MIX,DIISORDEX,DIISSTART,NSCF,-1,.TRUE.,SCRATCH,.FALSE., &
                                & ETEMP,mu,ENTROPY,IORBNR,PHOLE,PEXCITED,NBAUX,VRI,WRI,RIAPPROX)
                                
                                CALL PRINTENERGYEIGENVAL(BAS,EHFeigenup,EHFeigendown,Cup,Cdown,.TRUE.)
                                WRITE(*,*)'  '
                                WRITE(*,'(A24,I4,A4,A7,A8,E27.20,A8,E27.20,A3)')'The Exitation energy E (',IORBNR(1),'  ->',' EXITED','  ) = ',ETOT-Egs,' a.u. = ',autoev*(ETOT-Egs),' eV'
                                WRITE(*,*)'  '
                        ENDIF
                        !====================
                        ! CALCULATING FORCES
                        !================================================================================================================== 
                        IF ( CFORCE ) THEN
                                IF ( AFORCE ) THEN
                                        CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,&
                                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
                                ELSE
                                        CALL numforce(S,H0,Intsv,NB,NRED,Ne,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,DRF,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,NEPERIOD,EDIR,EPROFILE,EFIELDMAX,&
                                        & OMEGA,0.0d0,ADEF,Pup,Pdown,.TRUE.,NBAUX,ATOMSAUX,BASAUX,VRI,WRI,RIAPPROX,force)
                                ENDIF
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
                        CALL RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,.TRUE.,SCRATCH,.FALSE.,NBAUX,VRI,WRI,RIAPPROX)
                        
                        IF ( NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)

                        !====================
                        ! CALCULATING FORCES
                        !================================================================================================================== 
                        IF ( CFORCE ) THEN
                                P = P/2.0d0
                                IF ( AFORCE ) THEN
                                        CALL forces(CORRLEVEL,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,&
                                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
                                ELSE        
                                        CALL numforce(S,H0,Intsv,NB,NRED,Ne,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,DRF,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                                        & NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,NEPERIOD,EDIR,EPROFILE,EFIELDMAX,&
                                        & OMEGA,0.0d0,ADEF,P,P,.TRUE.,NBAUX,ATOMSAUX,BASAUX,VRI,WRI,RIAPPROX,force)
                                ENDIF
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
                IF ( WRITEDENS .AND. ( CORRLEVEL .EQ. 'RHF'  .OR. CORRLEVEL .EQ.  'URHF' .OR. DFTC  ) ) THEN 
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
                                IF ( .not. DFTC ) THEN
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
               
                IF (   AORBS .NE. 0 ) THEN
                        IF ( ABS(AORBS) .GT. NB ) AORBS = NB*(AORBS/ABS(AORBS))
                        call arbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS)
                ENDIF
                
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
                        CALL dqmc(Ne,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,BETA,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,ETOT-nucE,nucE,ESIGMA,EDQMC,NPERSIST,NRECALC,CUTTOFFFACTOR, &
                        & NVMC,VMCCALC,WRITEDENS,MESH,LIMITS)
                        
                ENDIF
        ENDIF
       
       !-----------------------------------------------------
       ! THE TOTAL EXECUTION TIME IS CALCULATED AND PRINTED
       !-----------------------------------------------------
4000 CONTINUE
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

       !-----------------------------------------------------------------------
       ! Saving the density matrice(s) in the case of a single atom calculation.
       ! To be used as a starting guess for molecular type of calculations with
       ! convergence problems.
       !------------------------------------------------------------------------
       IF ( NATOMS .EQ. 1 ) THEN
                IF ( CORRLEVEL .EQ. 'URHF' .OR. CORRLEVEL .EQ. 'LDA' .OR. CORRLEVEL .EQ.  'PBE' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL Wdensmatstartguess(ATOMS(1)%Z,NB,Pup,Pdown)
                ENDIF
                IF ( CORRLEVEL .EQ. 'RHF' ) CALL Wdensmatstartguess(ATOMS(1)%Z,NB,P,P)
       ENDIF
       !STOP
END PROGRAM uquantchem
