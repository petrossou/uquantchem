PROGRAM uquantchem
      USE OMP_LIB
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER  :: id, ierr, numprocessors, NATOMS
      INTEGER, ALLOCATABLE :: ATOMICNUMBERS(:)
      DOUBLE PRECISION, ALLOCATABLE :: APOS(:,:)
      DOUBLE PRECISION :: Tol,ETOT,NucE,Rn(3)
      CHARACTER(LEN=20) :: CORRLEVEL,EPROFILE
      INTEGER :: BASISMAP(120),BASISMAPAUX(120)
      TYPE(ATOM), ALLOCATABLE :: ATOMS(:),ATOMSAUX(:)
      TYPE(BASIS) :: BAS,BASAUX
      CHARACTER(LEN=6) :: DUMMY
      CHARACTER(Len=20) :: date,time,zone
      INTEGER :: NLINES,I,J,K,L,M,NB,Ne,Lmax,MESH(3),REDISTRIBUTIONFREQ,NPERSIST,NVMC,IOSA,FNATOMS,NSTEPS,NLSPOINTS,PORDER,NSI,NBAUX,LmaxAUX
      INTEGER*8 ::NRED,FNRED,TOTALNONZERO,NONZERO,NONZEROO,Istart,Iend,Istartg,Iendg,Istarts,Iends,NDIAG,MODEIGHT,N0
      LOGICAL :: finns,WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,NOREDIST,HYBRID,USEGTO,CORRALCUSP,VMCCALC,HFORBWRITE,CFORCE,RELAXN
      LOGICAL :: WRITEONFLY,MOVIE,MOLDYN
      DOUBLE PRECISION, ALLOCATABLE :: S(:,:),T(:,:),V(:,:),H0(:,:),MASSCORR(:,:),NUCDIRAC(:,:),Intsv(:),IntsvR(:),IntsDirac(:),EIGENVECT(:,:),Ints(:,:,:,:),gradIntsvR(:,:,:),HUCKELH(:,:),CHUCKEL(:,:),EHUCKEL(:)
      DOUBLE PRECISION, ALLOCATABLE :: VRI(:,:),gradVRI(:,:,:,:),WRI(:,:,:),gradWRI(:,:,:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gradS(:,:,:,:),gradT(:,:,:,:),gradV(:,:,:,:),DMAT(:,:),CHUCKEL2(:,:),dInts(:),Intso(:),Intsoo(:),Pupp(:,:),Pdownn(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigenup(:),EHFeigendown(:),Cup(:,:),Cdown(:,:),PEXCITED(:,:),H00(:,:),DPTENSOR(:,:,:),DIPOLET(:,:),PEXu(:,:,:),PEXd(:,:,:),PEXuu(:,:,:),PEXdd(:,:,:)
      COMPLEX*16, ALLOCATABLE :: Cupc(:,:),Cdownc(:,:),Pupc(:,:),Pdownc(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigen(:),P(:,:),Pup(:,:),Pdown(:,:),C1(:,:),C2(:,:),force(:,:),CN(:),LQ(:,:),CGQ(:,:),ZMUL(:),PHOLE(:,:),Pgup(:,:),Pgdown(:,:)
      DOUBLE PRECISION  :: PRYSR(25,25),PRYSW(25,25),rts(25),wts(25),TOTALTIME,ECISD,LIMITS(3),ENEXCM,EMP2,TIMESTEP,TEND,TSTART,CUTTOFFFACTOR,MIX,Aexpan,LAMDA,rc,INFOH,HK,ETEMPE,OMEGA
      DOUBLE PRECISION :: ESIGMA,EDQMC,a,b,BETA,BJASTROW,CJASTROW,GAMA,BETAA,POLY(6),deltar,DR,FTol,TEMPERATURE,kappa,alpha,EETOL,CNSTART(4),alphastart,kappastart,ETEMP,ENTROPY,mu,Egs,autoev
      DOUBLE PRECISION :: EFIELDMAX,TID,ETEMPEX,MIXEX,NORM,DAMPING
      INTEGER :: STARTTIME(8),STOPTIME(8),RUNTIME(8),NEEXC,SAMPLERATE,NREPLICAS,NRECALC,scount,DIISORD,DIISSTART,NTIMESTEPS,NSCF,PULAY,FIXNSCF,DORDER,IORBNR(2),AORBS,EDIR,NEPERIOD
      INTEGER :: NLEBEDEV,NCHEBGAUSS,LEBPOINTS(24),CGORDER,LORDER,NTOTALQUAD,Qstart,Qend,RELALGO,KK,ZEROSCFTYPE,DIISORDEX,TIMEINDEX,II
      INTEGER, EXTERNAL :: ijkl
      INTEGER, ALLOCATABLE :: rcounts(:),displs(:),Istart2(:),Istart3(:),IND1(:),IND2(:),IND3(:),IND4(:),Q1(:),Q2(:),Q3(:),Q11(:),Q22(:),Q33(:),EEMAP(:)
      INTEGER*8, ALLOCATABLE :: N0p(:)
      DOUBLE PRECISION, EXTERNAL :: massa
      LOGICAL :: RESTART,ZEROSCF,XLBOMD,DFTC,SOFTSTART,HUCKEL,SCRATCH,DOTDFT,ADEF,OPTH,DOABSSPECTRUM,DIFFDENS,RIAPPROX,DIAGDG,FIELDREAD,SCALARRELC,WRITEEEINTS
      INTEGER :: NSCCORR,LIMPRECALC,eeindex,NUMBEROFEEINTS
      DOUBLE PRECISION :: MIXTDDFT,SCERR,FIELDDIR(3),TRANSCOORD(3,3),cosTH,sinTH,cosFI,sinFI,RIKTNING(3),integraal
      DOUBLE PRECISION, PARAMETER :: lightspeed = 137.035999084
      character (len=1000) :: ee_file_name
      LOGICAL :: URHFTORHF

      call MPI_Init ( ierr )
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
      call MPI_Comm_size ( MPI_COMM_WORLD, numprocessors, ierr )

      !-------------------
      !STARTING THE CLOCK: 
      !-------------------
      call DATE_AND_TIME(date, time, zone,STARTTIME)

      WRITEEEINTS = .TRUE. ! Hard coded (There is no input that can be given. If true the (ij | kl) electron-electron tesor will be written to disk)

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
      ALLOCATE(N0p(numprocessors),rcounts(numprocessors),displs(numprocessors),Istart2(numprocessors),Istart3(numprocessors))
      
      CALL readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBERS,APOS,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, &
      & SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,REDISTRIBUTIONFREQ,NOREDIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,HYBRID,rc,CORRALCUSP,NVMC, &
      & HFORBWRITE,IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol,NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS,&
      & EETOL,RELALGO,SOFTSTART,HUCKEL,ZEROSCFTYPE,ETEMP,IORBNR,AORBS,DIISORDEX,DOTDFT,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF,OPTH,MIXEX,DOABSSPECTRUM,& 
      & DIFFDENS,NSCCORR,MIXTDDFT,SCERR,RIAPPROX,LIMPRECALC,DIAGDG,FIELDDIR,FIELDREAD,DAMPING,SCALARRELC,URHFTORHF)

      IF ( MOD(Ne,2) .NE. 0 .AND. CORRLEVEL .EQ. 'RHF' ) THEN
                IF ( id .EQ. 0 ) THEN
                        print*,'Number of electrons is odd and you have chosen CORRLEVEL = RHF'
                        print*,'changing CORRLEVEL = URHF enabaling an open shell calculation '
                ENDIF
                CORRLEVEL = 'URHF'
      ENDIF
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
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,*)'****************************************************************************'
                        WRITE(*,*)'***                          WARNING!                                    ***'
                        WRITE(*,*)'      The direction of the electric field has been set to (0, 0, 0 ).'
                        WRITE(*,*)'     Since this do not make sense the field direction is set to (0,0,1)'
                        WRITE(*,*)'and in the case of an AC-field the direction of propagation is set to (0,1,0)'
                        WRITE(*,*)'****************************************************************************'
                        WRITE(*,*)' '
                ENDIF
                EDIR = 2
         ENDIF
      ENDIF
     
      IF ( ADEF .AND. DOTDFT ) THEN
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A120)')'************************************************************************************************************************'
                        WRITE(*,'(A120)')'                                    WARNING CALCULATION STOPPED, REASON:                                                '
                        WRITE(*,'(A120)')'Time dependent DFT/HF (DOTDFT=.TRUE.) option cannot be combined with adiabatic electric field calculations (ADEF=.TRUE.)'
                        WRITE(*,'(A120)')'      You can only choose one of the two above options to be true, i.e either DOTDFT=.TRUE. excusive or ADEF=.TRUE.     '
                        WRITE(*,'(A120)')'************************************************************************************************************************'
                        WRITE(*,*)' '
                ENDIF
                STOP
      ENDIF


      IF ( ADEF ) THEN
              IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A120)')'=========================================================================================================================='
                        WRITE(*,'(A120)')'                             CALCULATION PERFORMED WITH THE PRESENCE OF ELECTRIC FIELD                                    '
                        IF ( XLBOMD .OR. MOLDYN .AND. EPROFILE .NE. 'ST' ) THEN
                                WRITE(*,'(A120)')'                    THE ELECTRONS WILL BE ASSUMED TO FOLLOW THE TIME-DEPENDENT FIELD ADIABATICALLY                        '
                        ENDIF
                        WRITE(*,'(A120)')'=========================================================================================================================='
                        WRITE(*,*)' '
                ENDIF
      ENDIF

      IF ( DORDER .LT. 4 )  DORDER = 4
      ALLOCATE(CN(DORDER))

      IF ( RELALGO .NE. 1 .AND. RELALGO .NE. 2 ) RELALGO = 1

      IF ( CORRLEVEL .EQ. 'RHF' .AND. DOTDFT ) THEN
              IF ( id .EQ. 0 ) THEN
                      WRITE(*,*)'=================================================================='
                      WRITE(*,*)'Restricted time dependent Hartree-Fock caclulation not implemented'
                      WRITE(*,*)'Changing CORRLEVEL from "RHF" to "URHF"'
                      WRITE(*,*)'=================================================================='
                      CORRLEVEL = 'URHF'
              ENDIF
      ENDIF
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
              NTOTALQUAD = CGORDER*LORDER*NATOMS
              
              ALLOCATE(Q11(NTOTALQUAD),Q22(NTOTALQUAD),Q33(NTOTALQUAD))
              L = 0
              DO I=1,CGORDER
                 DO J=1,LORDER
                    DO K=1,NATOMS
                       L=L+1
                       Q11(L) = I
                       Q22(L) = J
                       Q33(L) = K
                    ENDDO
                 ENDDO
             ENDDO
      ELSE
              CGORDER = 1
              LORDER = 1
              NTOTALQUAD = 1
              ALLOCATE(Q11(NTOTALQUAD),Q22(NTOTALQUAD),Q33(NTOTALQUAD))
              Q11 = 0.0d0
              Q22 = 0.0d0
              Q33 = 0.0d0
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
              CN = (/ -14.0d0, 36.0d0, -27.0d0, -2.0d0, 12.0d0, -6.0d0, 1.0d0 /)
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

      IF ( EPROFILE .NE. 'ST' .AND. .not. DOTDFT .AND. .not. MOLDYN .AND. ADEF ) THEN
              EPROFILE = 'ST'
              IF ( id .EQ. 0 ) THEN
                WRITE(*,*)' '
                WRITE(*,*)'TIME-INDEPENDENT CALCULATIONS ONLY PERMITTED WITH STATIC HOMOGENEOUS FIELD'
                WRITE(*,*)'EPROFILE has been changed to ST'
                WRITE(*,*)' '
              ENDIF
      ENDIF
      !IF ( ( RELAXN .OR. MOLDYN ) .AND. CORRLEVEL .NE. 'RHF' .AND. CORRLEVEL .NE. 'URHF'  ) THEN
      !        IF ( id .EQ. 0 ) THEN
      !                WRITE(*,*)'Relaxation/molecular-dynamics of nuclea only possible for RHF or URHF'
      !                WRITE(*,*)'Change CORRLEVEL to RHF or URHF'
      !        ENDIF
      !        STOP
      !ENDIF
    
      ! Ionization hole calculations only work at zero-electronic temperature
      ! and with the ZEROSCFTYPE = 1 alternative when calculating forces and 
      ! energies in a md-calculation
      IF (IORBNR(1) .NE. 0  ) THEN
              ZEROSCFTYPE = 1
              SOFTSTART = .FALSE.
              !ETEMP = -1000.0d0
      ENDIF


      IF ( RELAXN .OR. MOLDYN ) CFORCE = .TRUE.
     
      ! Ensuring that idempotency is not used when calculating the pulay
      ! contribution to the force when non-integere occupation of
      ! electronic
      ! states is used.

      IF ( ETEMP .GT. 0.0d0 .AND. CFORCE ) PULAY = 2
      IF ( ZEROSCFTYPE .EQ. 2 .AND. DFTC ) PULAY = 2

      IF ( MOLDYN ) RELAXN = .FALSE.

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
                        call stoexponent(ATOMS(I)%PSI(J),deltar,LAMDA,Aexpan,GAMA,BETAA,POLY,rc,USEGTO,CORRALCUSP,id)
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
             IF ( id .EQ. 0 ) print*,'ABORTING SINCE THE FILE  "BASISFILE" IS MISSING'
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
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)'                    WARNING!                        '
                        WRITE(*,*)'      THE INDEX OF THE HF-EIGEN FUNCTION TO BE SAVED '
                        WRITE(*,*)'EXCEEDS THE TOTAL NUMBER OF BASIS FUNCTIONS NB =',NB
                        WRITE(*,*)'THE INDEX HF-EIGENFUNCTION SAVED WILL THERFORE BE SET TO NB'
                ENDIF
               IOSA = NB
       ENDIF
       IF (HFORBWRITE .AND. IOSA .EQ. 0 ) THEN
               IOSA = (Ne-MOD(Ne,2))/2 + MOD(Ne,2)
       ENDIF
     
        IF ( id .EQ. 0 ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A50,I4)')'        Total number of basis-functions used NB = ',NB
        ENDIF

        ! Here the total molecular basis set is put together
        CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)
     
        !=======================================================================
        ! (3) HERE THE OVERLAP-, KINETIC- AND POTENTIAL ENERGY- MATRICES ARE
        !     CALCULATED.
        !=======================================================================
        
        ! Here the normalization of the basis functions is performed

        CALL normalize(BAS)
        
        NSI = 1
        NDIAG = ( NB*(NB+1) )/2

        IF ( CFORCE ) NSI = NB

        ! Here the overlap matrix is calculated 

        ALLOCATE(S(NB,NB),gradS(NATOMS,3,NSI,NSI),DMAT(NB,NB),HUCKELH(NB,NB),CHUCKEL(NB,NB),CHUCKEL2(NB,NB),EHUCKEL(NB),dInts(NDIAG),PEXu(24,NB,NB),PEXd(24,NB,NB))
        ALLOCATE(PEXuu(24,NB,NB),PEXdd(24,NB,NB))

        CALL overlap(NATOMS,BAS,S,gradS,NSI,CFORCE)
       
        ! Here the kinetic energy matrix is calculated

        ALLOCATE(T(NB,NB),gradT(NATOMS,3,NSI,NSI))

        CALL kinetic(NATOMS,BAS,T,gradT,NSI,CFORCE)
        
        ! Here the potential energy matrix is calculated
        
        ALLOCATE(V(NB,NB),gradV(NATOMS,3,NSI,NSI))
        
        CALL potential(BAS,NATOMS,ATOMS,V,gradV,NSI,CFORCE,id,numprocessors)
        
        ! Here the scalar-relatevistic correction terms: the p**4 mass-correction term and the nuclear Dirac term are calculated:
        ALLOCATE(MASSCORR(NB,NB),NUCDIRAC(NB,NB))
        IF ( SCALARRELC ) THEN
                CALL diraccorrelnucpot(BAS,NATOMS,ATOMS,NUCDIRAC,id,numprocessors)
                MASSCORR = -(1.0d0/(2.0d0*lightspeed**2))*MATMUL(T,T)
        ELSE
                NUCDIRAC(:,:) = 0.0d0
                MASSCORR(:,:) = 0.0d0
        ENDIF

        !================================================================================================
        ! If the auxiliary basis-file, BASISFILEAUX, is present in the run-directory and the user has not
        ! explicitly stated that he/she does not want to utilize the Resoulution ! of Identity approximation 
        ! for the (ij|kl)-integrals, the auxiliary basis will be set up here in the exact same manner as 
        ! the standard basis-set
        !================================================================================================

      NBAUX = 1

      inquire(file='BASISFILEAUX',exist=finns)
      IF ( finns .AND.  RIAPPROX ) THEN
             IF ( id .EQ. 0 ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A61)')'    ========================================================='
                WRITE(*,'(A60)')'     Using the Resolution of the Identity (RI) approximation'
                WRITE(*,'(A61)')'    ========================================================='
                WRITE(*,*)' '
             ENDIF
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
           IF ( id .EQ. 0 ) THEN
                WRITE(*,*)' '
                WRITE(*,'(A59,I4)')'   Total number of auxiliary basis-functions used NBAUX = ',NBAUX
                WRITE(*,*)' '
           ENDIF
           ! Here the total aux-molecular basis set is put together
           CALL gettotalbasis(NATOMS,ATOMSAUX,NBAUX,BASAUX,.FALSE.)

           !=======================================================================
           ! (3) HERE THE OVERLAP-, KINETIC- AND POTENTIAL ENERGY- MATRICES ARE
           !     CALCULATED.
           !=======================================================================

           ! Here the normalization of the basis functions is performed

           CALL normalize(BASAUX)

        ELSE
                IF ( RIAPPROX .AND. id .EQ. 0 ) THEN
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

        IF ( CORRLEVEL .EQ. 'MP2' .AND. RIAPPROX ) THEN
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A62)')'**************************************************************'
                        WRITE(*,'(A62)')'***                       WARNING!                         ***'
                        WRITE(*,'(A62)')' RI-approximation in MP2 calculations are not yet implemented.'
                        WRITE(*,'(A62)')'Switching off RI-approximation, i.e setting RIAPPROX = .FALSE.'
                        WRITE(*,'(A62)')'**************************************************************'
                        WRITE(*,*)' '
                ENDIF
                RIAPPROX = .FALSE.
        ENDIF
        IF ( CORRLEVEL .EQ. 'CISD' .AND. RIAPPROX ) THEN
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A61)')'**************************************************************'
                        WRITE(*,'(A61)')'***                      WARNING!                          ***'
                        WRITE(*,'(A61)')'RI-approximation in CISD calculations are not yet implemented.'
                        WRITE(*,'(A61)')'Switching off RI-approximation, i.e setting RIAPPROX = .FALSE.'
                        WRITE(*,'(A61)')'**************************************************************'
                        WRITE(*,*)' '
                ENDIF
                RIAPPROX = .FALSE.
        ENDIF

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        ! If a DFT type of calculation is performed the 
        ! exchange correlation matrix elements and their 
        ! nuclear gradients are allocated 

        !IF ( DFTC ) ALLOCATE(Vxc(2,NB,NB),gVxc(NATOMS,2,3,NB,NB))
        
        !===========================================================================
        ! (4) Here the electron-electron repulsion 'tensor', Intsv(:), is calulated
        !===========================================================================
        
        ! The number of symmetry reduced electron-electron tensor elements
        !  In the MPI-version the full symetry is not employed. Here the number
        ! of symmetry reduced elements are twice as many as in the openmp and serial version.

         NRED =  NB*(NB+1)
         NRED = NRED/2
         NRED = NRED*NRED

        
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
        IF ( id .EQ. 0 .AND. .not. RELAXN .AND. .not. MOLDYN  ) THEN
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
        
        IF ( id .EQ. 0 .AND. .not. RELAXN .AND. .not. MOLDYN  ) THEN
         print*,'    ========================================================'
         print*,'        Calculation of electron repulsion tensor (ij|kl)    '
         print*,'    ========================================================'
        ENDIF
        !===================================================================================
        ! PREPARATION FOR THE PARALELL COMPUTATION OF THE EXCHANGE-CORR ENERGY AND POTENTIAL
        !===================================================================================
        IF ( DFTC ) THEN
           N0 = INT(NTOTALQUAD/numprocessors)

           Qstart = 0
           Qend = 0
           N0p(:) = 0
           DO I=1,numprocessors
             N0p(I) = N0
           ENDDO

           I = 1
           DO WHILE( I .LE. MOD(NTOTALQUAD,numprocessors) )
               N0p(I) = N0p(I) + 1
               I = I+1
           ENDDO

       
           IF ( id .EQ. 0 ) THEN
              Qstart = 1
           ELSE
              Qstart = 1
               DO I=1,id
                  Qstart = Qstart+N0p(I)
               ENDDO
           ENDIF
           Qend = 0
           DO I=1,id+1
               Qend = Qend + N0p(I)
           ENDDO
                
           ALLOCATE(Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend))

           DO I=Qstart,Qend
                Q1(I) = Q11(I)
                Q2(I) = Q22(I)
                Q3(I) = Q33(I)
           ENDDO
        
           DEALLOCATE(Q11,Q22,Q33)
        ELSE
                Qstart = 1
                Qend = 1
                ALLOCATE(Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend))
                DEALLOCATE(Q11,Q22,Q33)
        ENDIF

        !==========================================================================================
        ! END OF PREPARATION FOR THE PARALELL COMPUTATION OF THE EXCHANGE-CORR ENERGY AND POTENTIAL
        !==========================================================================================
        !============================================================================
        ! Here we distribute the calculation of the columns of the ee-integral array
        ! Intsv over the (numprocessors) number of mpi threads.
        !----------------------------------------------------------------------------

        N0 = INT(NRED/numprocessors,8)

        Istart = 0
        Iend = 0
        N0p(:) = 0
        DO I=1,numprocessors
          N0p(I) = N0
        ENDDO

        !--------------------------------------------------------------------------------------------------------------
        ! If the number of elements (NRED) is not divisible by the number of processors, the remaining 
        ! number of elements of Intsv(I) to be calculated are distributed over the (numprocessors) number of mpi threads
        !--------------------------------------------------------------------------------------------------------------
        I = 1
        MODEIGHT = NRED - INT(NRED/numprocessors,8)*numprocessors
        !DO WHILE( I .LE. MOD(NRED,numprocessors) )
        DO WHILE( I .LE. MODEIGHT ) 
               N0p(I) = N0p(I) + 1
               I = I+1
        ENDDO


        IF ( id .EQ. 0 ) THEN
             Istart = 1
        ELSE
            Istart = 1
             DO I=1,id
               Istart = Istart+N0p(I)
             ENDDO
        ENDIF
        Iend = 0
        DO I=1,id+1
            Iend = Iend + N0p(I)
        ENDDO
        
        !=================================================================================
        ! END OF PREPARATION FOR THE PARALELL COMPUTATION OF THE ee-repulsion-tensor Intsv
        !==================================================================================
        IF ( RIAPPROX .AND. CFORCE ) THEN
                FNATOMS = NATOMS
        ELSE
                FNATOMS = 1
        ENDIF
        
        
        ALLOCATE(VRI(NBAUX,NBAUX),gradVRI(FNATOMS,3,NBAUX,NBAUX))
        ALLOCATE(WRI(BAS%NBAS,BAS%NBAS,NBAUX),gradWRI(FNATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX))
        IF ( RIAPPROX ) THEN
                CALL calcWRI(NATOMS,FNATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,CFORCE,id,numprocessors)
                CALL calcVRI(NATOMS,FNATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,CFORCE,id,numprocessors)
        ENDIF
        
        !==============================================================================
        ! Here the starting-guess density matrices are put together from
        ! density matrices obtained from previous atomic calculations. The
        ! atomic density matrices are read from the file 'DENSMATSTARTGUESS.dat'
        !==============================================================================
        SCRATCH = .TRUE.
        
        DMAT = 1.0d0

        Istarts =  Istart
        Iends   =  Iend
        NONZERO = 1
        IF ( .not. MOLDYN .AND. .not. RELAXN ) THEN
        ALLOCATE(IND1(1),IND2(1),IND3(1),IND4(1))
        
        ! Counting the number of "non-zero" elements og the ee-tensor (ij|kl), elements satisfying |(ij|kl)| < EETOL 
        IF ( RIAPPROX ) THEN
                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
        ELSE
                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts)
        ENDIF
        
       DEALLOCATE(IND1,IND2,IND3,IND4)
        
        IF ( NONZERO .EQ. 0 ) THEN
            print*,'Number of (ij|kl) to be calculated on thread,',id,', is zero!'
            ! We cannot allocate IND1(1:0),IND1(1:0),....
            ! Therfore we must create a dummy entry:
            NONZEROO = NONZERO +1
        ELSE
            NONZEROO = NONZERO
        ENDIF

        ALLOCATE(IND1(NONZEROO),IND2(NONZEROO),IND3(NONZEROO),IND4(NONZEROO))
        
        ! creating a mapping from the non-zero index running from 1 to NONZERO on each thread to the contracted index ( see the routine ijkl.f90 for definition )
        IF ( RIAPPROX ) THEN
                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
        ELSE
                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts)
        ENDIF
        
        Istart = 1
        Iend = NONZEROO
        
        IF ( CFORCE ) THEN
                Istartg = Istart
                Iendg = Iend
                FNATOMS = NATOMS
        ELSE
                Istartg = 1
                Iendg = 1
                FNATOMS = 1
        ENDIF
        
        ALLOCATE(IntsvR(Istart:Iend),gradIntsvR(FNATOMS,3,Istartg:Iendg))
        ALLOCATE(IntsDirac(Istart:Iend))
        
        write (ee_file_name,"('EETENSOR_',i0,'.dat')") id
        inquire(file=trim(ee_file_name),exist=finns)
        IF ( id .EQ. 0 .AND. finns .AND. WRITEEEINTS ) THEN
                print*,'    ========================================================='
                print*,'     Precalulated tensor elements (ij|lm) are read from disc '
                print*,'    ========================================================='
        ENDIF
        
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
        IF ( .not. finns .OR. CFORCE .OR. .not. WRITEEEINTS ) THEN
                IF ( RIAPPROX ) THEN
                        CALL eeintsRI(NATOMS,FNATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istartg,Iendg,&
                                & PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,CFORCE,id,NONZERO,DMAT,dInts,NBAUX,VRI,WRI,gradVRI,gradWRI)
                ELSE
                        CALL eeints(NATOMS,FNATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istartg,Iendg,&
                                & PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,CFORCE,id,NONZERO,DMAT,dInts)
                ENDIF

                IF ( WRITEEEINTS ) THEN 
                        write (ee_file_name,"('EETENSOR_',i0,'.dat')") id
                        OPEN(id,FILE=trim(ee_file_name),ACTION='WRITE')
                        WRITE(id,*)Iend
                        DO I=1,Iend
                                WRITE(id,*)I,IntsvR(I)
                        ENDDO
                        CLOSE(id)
                ENDIF

        ELSE
                write (ee_file_name,"('EETENSOR_',i0,'.dat')") id
                OPEN(id,FILE=trim(ee_file_name),STATUS='OLD',ACTION='READ')
                READ(id,*)NUMBEROFEEINTS
                DO I=1,NUMBEROFEEINTS
                        READ(id,*)eeindex,integraal
                        IntsvR(eeindex) = integraal
                ENDDO
                CLOSE(id)
        ENDIF
   
                
               
        !IF ( SCALARRELC ) THEN
        !        write (ee_file_name,"('DOUBLEOVERLAP_',i0,'.dat')") id
        !        inquire(file=trim(ee_file_name),exist=finns)
        !        IF ( id .EQ. 0 .AND. finns .AND. WRITEEEINTS ) THEN
        !                print*,'    =========================================================='
        !                print*,'     Precalulated double overlap integrals are read from disc '
        !                print*,'    =========================================================='
        !        ENDIF
        !      
        !
        !        IF ( .not. finns .OR. .not. WRITEEEINTS ) THEN
        !                CALL doubleoverlap(NATOMS,BAS,IntsDirac,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,NONZERO)
        !                
        !                IF ( WRITEEEINTS ) THEN 
        !                        write (ee_file_name,"('DOUBLEOVERLAP_',i0,'.dat')") id
        !                        OPEN((id+1)*10,FILE=trim(ee_file_name),ACTION='WRITE')
        !                        WRITE((id+1)*10,*)Iend
        !                        DO I=1,Iend
        !                                WRITE((id+1)*10,*)I,IntsDirac(I)
        !                        ENDDO
        !                        CLOSE((id+1)*10)
        !                ENDIF
        !        ELSE
        !                write (ee_file_name,"('DOUBLEOVERLAP_',i0,'.dat')") id
        !                OPEN((id+1)*10,FILE=trim(ee_file_name),STATUS='OLD',ACTION='READ')
        !                READ((id+1)*10,*)NUMBEROFEEINTS
        !                DO I=1,NUMBEROFEEINTS
        !                        READ((id+1)*10,*)eeindex,integraal
        !                        IntsDirac(eeindex) = integraal
        !                ENDDO
        !                CLOSE((id+1)*10)
        !        ENDIF
        !ENDIF

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

        IF ( id .EQ. 0 ) THEN
         print*,' '
         WRITE(*,'(A62)')'    =========================================================='
         WRITE(*,'(A47,I12)')'    Number of non-reduced ee-tensor elements = ',NRED
         WRITE(*,'(A47,I12)')'    Number of  reduced    ee-tensor elements = ',TOTALNONZERO
         WRITE(*,'(A62)')'    =========================================================='
         print*,' '
        ENDIF
        
        !-------------------------------------------------------------------------------------
        ! preparation for gathering all the different parts of the ee-tensor, IntsvR, that are 
        ! spread out on all the computational threds, so that they can be put into the single 
        ! tensor Intsv located at thread id=0
        !--------------------------------------------------------------------------------------
        Istart2(1+id) = Istart
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
         DO I=1,numprocessors
             rcounts(I) = 1
             displs(I) = I -1
         ENDDO
        
         CALL MPI_GATHERV(Istart2(1+id),1,MPI_INTEGER,Istart3,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
         CALL MPI_BCAST(Istart3,numprocessors,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        
         DO I=1,numprocessors
             rcounts(I) = N0p(I)
             displs(I) = Istart3(I)-1
         ENDDO
        
        scount = N0p(id+1)
        ENDIF

        ALLOCATE(H0(NB,NB),PHOLE(NB,NB),PEXCITED(NB,NB),H00(NB,NB),DPTENSOR(3,NB,NB),DIPOLET(NB,NB))
        
        H0 = T+V

        IF ( SCALARRELC ) THEN
                H0 = H0 + NUCDIRAC + MASSCORR
        ENDIF


        IF ( ADEF .AND. EPROFILE .EQ. 'DP' ) THEN
                IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,*)'WARNING ADIABATIC ELECRIC FIELD CALCULATIONS ARE NOT PERMITTED WITH DIRAC PULSE ENVELOPE'
                        WRITE(*,*)'CHANGING EPROFILE FROM DP TO ST'
                        WRITE(*,*)' '
                ENDIF
                EPROFILE = 'HO'
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
        IF ( NATOMS .GT. 0 ) THEN
                CALL Rdensmatstartguess(ATOMS,NATOMS,NB,Pup,Pdown,SCRATCH,id)
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
        IF ( DIAGDG .AND. .not. finns .OR. (DIAGDG .AND. NATOMS .EQ. 1) ) THEN
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

        IF ( RELAXN ) THEN
                IF ( EPROFILE .NE. 'ST' .AND. ADEF ) THEN
                        EPROFILE = 'ST'
                        IF ( id .EQ. 0 ) THEN
                                WRITE(*,*)' '
                                WRITE(*,*)'WARNING RELAXATION ONLY PERMITTED FOR STATIC HOMOGENEOUS FIELD'
                                WRITE(*,*)'EPROFILE has been changed to ST'
                                WRITE(*,*)' '
                        ENDIF
                ENDIF
                ! Relxing by searching for minimum energy
                IF ( RELALGO .EQ. 2 ) THEN
                        CALL relax(gradS,gradT,gradV,S,H0,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS, &
                                & APPROXEE,CORRLEVEL,PRYSR,PRYSW,Istarts,Iends,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,&
                                & LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,EDIR,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
                        GOTO 4000
                ENDIF
                ! Relaxing by searching for FORCES = 0
                IF ( RELALGO .EQ. 1 ) THEN
                        CALL relaxf(gradS,gradT,gradV,S,H0,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS, &
                                & APPROXEE,CORRLEVEL,PRYSR,PRYSW,Istarts,Iends,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,&
                                & LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,EDIR,EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
                        GOTO 4000
                ENDIF
        ENDIF

        IF ( MOLDYN ) THEN
                IF ( id .EQ. 0 ) THEN
                        print*,' ========================================================'
                        print*,'           Comencing with molecular dynamics calculation '
                        print*,' ========================================================'
                ENDIF
                NTIMESTEPS = INT(TEND/TIMESTEP)
                CALL moleculardynamics(gradS,gradT,gradV,S,H0,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISORDEX,DIISSTART,NATOMS,NTIMESTEPS,&
                   & TIMESTEP,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                   & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,Istarts,Iends,numprocessors,id,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,&
                   & Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL,ZEROSCFTYPE,ETEMP,IORBNR,OMEGA,EDIR,NEPERIOD,EPROFILE,EFIELDMAX,ADEF,OPTH,MIXEX, &
                   & BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX,DAMPING)
                GOTO 4000
                                                                                                             
        ENDIF

        !==============================================================
        ! (5) Here the Unrestricted or Restricted Hartree-FocK
        !     equations are selfconsistently solved
        !==============================================================
        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                
                ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB),P(NB,NB),C1(NB,NB))
                IF ( HUCKEL ) THEN
                   P = DMAT
                ENDIF 
                IF ( SCRATCH ) THEN
                        P = 0.0d0
                ELSE
                        P = 0.50d0*(Pup+Pdown)
                ENDIF
                CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART, &
                & NSCF,-1,numprocessors,id,.TRUE.,SCRATCH,.FALSE., SCALARRELC,IntsDirac,ETEMP,ENTROPY)
               
                IF ( id .EQ. 0 .AND. NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)
                !====================
                ! CALCULATING FORCES
                !================================================================================================================== 
                IF ( CFORCE ) THEN
                        P = P/2.0d0
                        CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,&
                             & S,H0,IntsvR,PULAY,force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF)
                        IF ( id .EQ. 0 ) THEN
                           WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                           WRITE(*,*)' '
                           WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                           DO I=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                           ENDDO
                        ENDIF
                ENDIF
                !================================================================================================================== 

                IF ( HFORBWRITE .AND. id .EQ. 0 ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,EIGENVECT,EIGENVECT)

                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_BCAST(ETOT,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(EHFeigen,NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(EIGENVECT,NB**2,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
               
                ! If specified, the charge density is calculated and saved here:
                IF ( WRITEDENS .AND. id .EQ. 0 .AND. ( CORRLEVEL .EQ. 'RHF' .OR. CORRLEVEL .EQ.  'URHF' .OR. DFTC )  ) THEN 
                        C1 = 0.0d0
                        DO I=1,Ne/2
                                C1(:,I) = EIGENVECT(:,I)
                        ENDDO
                        CALL makedens(C1,NB,P)
                        P = 2*P
                        CALL chargedenssave(LIMITS,MESH,BAS,P,Ne)
                ENDIF

                IF ( WHOMOLUMO .AND. id .EQ. 0 ) call homolumosave(LIMITS,MESH,BAS,EIGENVECT,EIGENVECT,EHFeigen,EHFeigen,Ne,id,NATOMS,ATOMS,DFTC)
        ENDIF

        IF ( CORRLEVEL .EQ. 'URHF' .OR. CORRLEVEL .EQ. 'CISD' .OR. CORRLEVEL .EQ. 'MP2' .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC ) THEN
                
                ALLOCATE(EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),P(NB,NB),C1(NB,NB),C2(NB,NB),Pgup(NB,NB),Pgdown(NB,NB))
                ALLOCATE(Cupc(NB,NB),Cdownc(NB,NB),Pupc(NB,NB),Pdownc(NB,NB),Pupp(NB,NB),Pdownn(NB,NB))
                IF ( HUCKEL ) THEN
                   Pup = DMAT/2.0d0
                   Pdown = DMAT/2.0d0
                ENDIF
                IF ( SCRATCH ) THEN
                   Pup = 0.0d0
                   Pdown = 0.0d0
                ENDIF
                IF ( .not. RESTRICT .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC )  THEN
                   IF ( .not. DFTC ) THEN
                           CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol, & 
                                & EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1, &
                                & numprocessors,id,.TRUE.,SCRATCH,.FALSE.,ETEMP,ENTROPY,SCALARRELC,IntsDirac,URHFTORHF)
                           IF ( DOTDFT ) THEN
                                ETEMPE = ETEMP
                                Pupc = Pup
                                Pdownc = Pdown
                                Cupc = Cup
                                Cdownc = Cdown
                                NTIMESTEPS = INT(TEND/TIMESTEP)
                                CALL TDFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD, &
                                & NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cupc,Cdownc,Pupc,Pdownc,&
                                & numprocessors,id,.TRUE.,mu,OMEGA,EDIR,NEPERIOD,EPROFILE,NTIMESTEPS,TIMESTEP,EFIELDMAX,IORBNR,PEXu,PEXd,PEXuu,PEXdd,DOABSSPECTRUM,DIFFDENS,NSCCORR,MIXTDDFT,SCERR)
                                !==========================================================================================
                                ! Saving the chargensity at time index = TIMESTEP/10,2TIMESTEP/10,3TIMESTEP/10,...,TIMESTEP
                                !==========================================================================================
                                IF (  WRITEDENS ) THEN
                                        DO TIMEINDEX=1,24
                                                Pup = PEXu(TIMEINDEX,:,:)
                                                Pdown = PEXd(TIMEINDEX,:,:)
                                                Pupp = PEXuu(TIMEINDEX,:,:)
                                                Pdownn = PEXdd(TIMEINDEX,:,:)
                                                CALL exciteddenssave(LIMITS,MESH,BAS,Pup,Pdown,Pupp,Pdownn,NATOMS,ATOMS,TIMEINDEX,id,DIFFDENS)
                                        ENDDO
                                ENDIF
                                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                                GOTO 4000
                          ENDIF
                   ENDIF 
                  
                   ETEMPE = ETEMP
                   
                   IF ( DFTC ) THEN
                           CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD, &
                           & NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,&
                           &NSCF,-1,numprocessors,id,.TRUE.,SCRATCH,.FALSE.,ETEMPE,ENTROPY,mu)
                           IF ( DOTDFT ) THEN
                                Pupc = Pup
                                Pdownc = Pdown
                                Cupc = Cup
                                Cdownc = Cdown
                                NTIMESTEPS = INT(TEND/TIMESTEP)
                                CALL TDFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD, &
                                & NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cupc,Cdownc,Pupc,Pdownc,&
                                & numprocessors,id,.TRUE.,mu,OMEGA,EDIR,NEPERIOD,EPROFILE,NTIMESTEPS,TIMESTEP,EFIELDMAX,IORBNR,PEXu,PEXd,PEXuu,PEXdd,DOABSSPECTRUM,DIFFDENS,NSCCORR,MIXTDDFT,SCERR)
                                !==========================================================================================
                                ! Saving the chargensity at time index = TIMESTEP/10,2TIMESTEP/10,3TIMESTEP/10,...,TIMESTEP
                                !==========================================================================================
                                IF (  WRITEDENS ) THEN
                                        DO TIMEINDEX=1,24
                                                Pup = PEXu(TIMEINDEX,:,:)
                                                Pdown = PEXd(TIMEINDEX,:,:)
                                                Pupp = PEXuu(TIMEINDEX,:,:)
                                                Pdownn = PEXdd(TIMEINDEX,:,:)
                                                CALL exciteddenssave(LIMITS,MESH,BAS,Pup,Pdown,Pupp,Pdownn,NATOMS,ATOMS,TIMEINDEX,id,DIFFDENS)
                                        ENDDO
                                ENDIF
                                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                                STOP
                           ENDIF
                           IF ( id .EQ. 0 ) THEN 
                                   CALL PRINTENERGYEIGENVAL(BAS,EHFeigenup,EHFeigendown,Cup,Cdown,.FALSE.)
                                   CALL PRINTDIPOLETENSOR(BAS,Cup,Cdown)
                           ENDIF
                   ENDIF
                   IF ( id .EQ. 0 .AND. NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,Pup+Pdown,S,ZMUL)
                   
                        IF ( IORBNR(1) .NE. 0 .AND. DFTC ) THEN
                                IF ( ABS(IORBNR(1)) .GT. NB ) THEN
                                        IORBNR(1) = NB
                                        IF ( id .EQ. 0 ) THEN
                                                WRITE(*,'(A70)')'**********************************************************************'
                                                WRITE(*,'(A27,I4,A7,I4,A40)')'        WARNING! IORBNR(1)=', IORBNR(1), '> NB = ',NB,' (NB = Number of basis functions)       '
                                                WRITE(*,'(A70)')'                    Setting IORBNR(1) = NB                            '
                                                WRITE(*,'(A70)')'**********************************************************************'
                                        ENDIF
                                ENDIF
                                IF ( ABS(IORBNR(2)) .GT. NB ) THEN
                                        IORBNR(2) = NB
                                        IF ( id .EQ. 0 ) THEN
                                                WRITE(*,'(A70)')'**********************************************************************'
                                                WRITE(*,'(A27,I4,A7,I4,A40)')'        WARNING! IORBNR(2)=', IORBNR(2), '> NB = ',NB,' (NB = Number of basis functions)       '
                                                WRITE(*,'(A70)')'                    Setting IORBNR(2) = NB                            '
                                                WRITE(*,'(A70)')'**********************************************************************'
                                        ENDIF
                                ENDIF
                                IF ( id .EQ. 0 ) THEN
                                        WRITE(*,'(A70)')'======================================================================'
                                        WRITE(*,'(A70)')'             Performing Core-hole excitation calculation              '
                                        WRITE(*,'(A70)')'======================================================================'
                                        IF ( IORBNR(2) .EQ. 0 ) THEN
                                                WRITE(*,'(A60,I4)')' Orbital being ionized and treated as a hole is orbital nr =',IORBNR
                                        ELSE
                                                WRITE(*,'(A46,I4)')'         Orbital being excited is orbital nr =',IORBNR(1)
                                                WRITE(*,'(A46,I4)')'     Orbital is being promoted to orbital nr =',IORBNR(2)
                                        ENDIF
                                ENDIF
                                
                                Egs = ETOT
                                Pgup = Pup
                                Pgdown = Pdown
                                ETEMPEX = -1.0d0
                                CALL DFTCHI(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol, &
                                & EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,Pgup,Pgdown,MIX,DIISORDEX,DIISSTART,NSCF,-1,numprocessors,id,.TRUE.,.TRUE.,.FALSE.,ETEMPEX,ENTROPY,mu,&
                                & IORBNR,PHOLE,OPTH)
                                IF ( id .EQ. 0 ) THEN
                                   CALL PRINTENERGYEIGENVAL(BAS,EHFeigenup,EHFeigendown,Cup,Cdown,.TRUE.)
                                   WRITE(*,*)'  '
                                   WRITE(*,'(A24,I4,A4,A7,A8,E27.20,A8,E27.20,A3)')'The Exitation energy E (',IORBNR(1),'  ->',' EXITED','  ) = ',ETOT-Egs,' a.u. = ',autoev*(ETOT-Egs),' eV'
                                   WRITE(*,*)'  '
                                ENDIF
                        ENDIF

                   !====================
                   ! CALCULATING FORCES
                   !================================================================================================================== 
                   IF ( CFORCE ) THEN
                          CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY,&
                          & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF)
                          
                          IF (id .EQ. 0 ) THEN
                             WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                             WRITE(*,*)' '
                             WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                             DO I=1,NATOMS
                                  WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                             ENDDO
                          ENDIF
                   ENDIF
                   !================================================================================================================== 
                  
                IF ( HFORBWRITE .AND. id .EQ. 0 ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,Cup,Cdown)

                   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                   CALL MPI_BCAST(EHFeigenup,NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   CALL MPI_BCAST(EHFeigendown,NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   CALL MPI_BCAST(ETOT,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   CALL MPI_BCAST(Cup,NB**2,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   CALL MPI_BCAST(Cdown,NB**2,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                ENDIF

                IF ( (RESTRICT .AND. CORRLEVEL .EQ. 'CISD') .OR. (RESTRICT .AND.  CORRLEVEL .EQ. 'MP2')  ) THEN
                        ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB))
                        P = 0.0d0
                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD, &
                        & DIISSTART,NSCF,-1,numprocessors,id,.TRUE.,SCRATCH,.FALSE., SCALARRELC,IntsDirac,ETEMP,ENTROPY)

                        IF ( id .EQ. 0 .AND. NATOMS .GT. 1 ) CALL mulliken(NATOMS,ATOMS,BAS,P,S,ZMUL)
                        !====================
                        ! CALCULATING FORCES
                        !================================================================================================================== 
                        IF ( CFORCE ) THEN
                           p = P/2.0d0
                           CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY, &
                                & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL,EDIR,EPROFILE,EFIELDMAX,OMEGA,TID,ADEF)
                        
                           IF (id .EQ. 0 ) THEN
                               WRITE(*,'(A29)')'CALCULATED INTERATOMC FORCES:'
                               WRITE(*,*)' '
                               WRITE(*,'(A84)')'  ATOM           Fx [au]                       Fy [au]                       Fz [au]'
                               DO I=1,NATOMS
                                    WRITE(*,'(I4,E30.20,E30.20,E30.20)')I,force(I,1),force(I,2),force(I,3)
                               ENDDO
                           ENDIF
                       ENDIF
                       !================================================================================================================== 
                       
                        IF ( HFORBWRITE .AND. id .EQ. 0 ) CALL hforbitalsave(IOSA,LIMITS,MESH,BAS,EIGENVECT,EIGENVECT)

                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        CALL MPI_BCAST(ETOT,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                        CALL MPI_BCAST(EHFeigen,NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                        CALL MPI_BCAST(EIGENVECT,NB**2,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                        Cup   = EIGENVECT
                        Cdown = EIGENVECT
                        EHFeigenup   = EHFeigen
                        EHFeigendown = EHFeigen
                ENDIF
                
                ! If specified, the charge density is calculated and saved here:
                IF ( WRITEDENS .AND. id .EQ. 0 .AND. ( CORRLEVEL .EQ. 'RHF' .OR. CORRLEVEL .EQ.  'URHF' .OR. DFTC  ) ) THEN 
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

                IF ( WHOMOLUMO .AND. id .EQ. 0 ) call homolumosave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,id,NATOMS,ATOMS,DFTC)
               
                IF ( AORBS .NE. 0  .AND. id .EQ. 0  ) THEN
                        IF ( ABS(AORBS) .GT. NB ) AORBS = NB*(AORBS/ABS(AORBS))
                        call arbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS)
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'CISD' .OR. CORRLEVEL .EQ. 'MP2' ) THEN
                        !print*,'========================================================'
                        !print*,'  Tranfering (ij|kl) from vector to tensor form         '
                        !print*,'========================================================'
                       
                       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                       
                       ALLOCATE(Intso(NRED),Intsoo(NRED))
                       Intso = 0.0d0
                       Intsoo = 0.0d0
                       
                       
                       DO KK=1,NONZERO
                           I = IND1(KK)
                           J = IND2(KK)
                           K = IND3(KK)
                           L = IND4(KK)
                           Intso(ijkl(I,J,K,L,NB)) = IntsvR(KK)
                       ENDDO

                       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                       DEALLOCATE(IntsvR)

                       CALL MPI_REDUCE(Intso,Intsoo,NRED,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
                       
                       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                       DEALLOCATE(Intso)

                       CALL MPI_BCAST(Intsoo,NRED,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                       ALLOCATE(Ints(NB,NB,NB,NB))
                       Ints  = 0.0d0
                       DO I=1,NB
                        DO J=1,NB
                         DO K=1,NB
                          DO L=1,NB
                           Ints(I,J,K,L) = Intsoo(ijkl(I,J,K,L,NB))
                          ENDDO
                         ENDDO
                        ENDDO
                       ENDDO
                       DEALLOCATE(Intsoo)
                 ENDIF


                
                IF ( CORRLEVEL .EQ. 'MP2'  ) CALL MP2(Cup,Cdown,Ints,NB,Ne,ETOT-nucE,nuce,EMP2,EHFeigenup,EHFeigendown,SPINCONSERVE,numprocessors,id)
                 
                ! Here the CISD calculation starts
                IF ( CORRLEVEL .EQ. 'CISD' ) THEN       
                        IF ( id .EQ. 0 ) THEN
                                print*,' '
                                print*,'  ========================================================'
                                print*,'    Starting Configuration interaction calculation (CISD) '
                                print*,'  ========================================================'
                        ENDIF
                        !CALL CISD(Cup,Cdown,Ints,H0,S,NB,Ne,ETOT-nucE,nuce,ECISD,WRITECICOEF,EHFeigenup,EHFeigendown,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,id)
                        !CALL CISD(NB,Ne,ETOT-nucE,nucE,Cup,Cdown,H0,Ints,id,SPINCONSERVE)
                        CALL CISD(NB,Ne,ETOT-nucE,nucE,Cup,Cdown,H0,Ints,SPINCONSERVE,WRITECICOEF,EHFeigenup,EHFeigendown,ENEXCM,LEXCSP,NEEXC,ECISD,numprocessors,id)
                ENDIF

                ! Here the DQMC or VMV calculation starts
                IF ( CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' ) THEN
                        IF ( id .EQ. 0 .AND. CORRLEVEL .EQ. 'DQMC' ) THEN
                                print*,' '
                                print*,' =============================================================='
                                print*,'   Starting the Diffusion Quantum Montecarlo Calulation (DQMC) '
                                print*,' =============================================================='
                        ENDIF
                        IF ( id .EQ. 0 .AND. CORRLEVEL .EQ. 'VMC' ) THEN
                                print*,' '
                                print*,' =============================================================='
                                print*,'  Starting the Variational Quantum Montecarlo Calulation (VMC) '
                                print*,' =============================================================='
                        ENDIF
                        
                        CALL dqmc(Ne,BJASTROW,CJASTROW,NATOMS,Cup,Cdown,BAS,ATOMS,BETA,SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,ETOT-nucE,nucE,ESIGMA,EDQMC,NPERSIST,NRECALC,NOREDIST, &
                                  & REDISTRIBUTIONFREQ,CUTTOFFFACTOR,NVMC,VMCCALC,WRITEDENS,MESH,LIMITS,numprocessors,id)

                ENDIF
        ENDIF
       
       !-----------------------------------------------------
       ! THE TOTAL EXECUTION TIME IS CALCULATED AND PRINTED
       !-----------------------------------------------------
4000   CONTINUE
       call DATE_AND_TIME(date, time, zone,STOPTIME)
       IF ( id .EQ. 0 ) THEN
       print*,' '
       print*,'        ============== TOTAL RUNTIME: ================'
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
        WRITE(*,'(A20,I5,A3,I5,A3,I5,A4,I5,A2)')'          RUNTIME=',RUNTIME(3),' d ',RUNTIME(5),' h ',RUNTIME(6),' min ',RUNTIME(7),' s'
        print*,' '
        print*,'        =============================================='
       ENDIF
       IF ( NATOMS .EQ. 1 .AND. id .EQ. 0 ) THEN
                IF ( CORRLEVEL .EQ. 'URHF' .OR. CORRLEVEL .EQ. 'LDA' .OR.  CORRLEVEL .EQ.  'PBE' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL Wdensmatstartguess(ATOMS(1)%Z,NB,Pup,Pdown)
                ENDIF
                IF ( CORRLEVEL .EQ. 'RHF' ) CALL Wdensmatstartguess(ATOMS(1)%Z,NB,P,P)
       ENDIF
       CALL MPI_FINALIZE(ierr)
END PROGRAM uquantchem
