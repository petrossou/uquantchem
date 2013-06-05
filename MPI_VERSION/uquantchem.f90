PROGRAM uquantchem
      USE OMP_LIB
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER  :: id, ierr, numprocessors, NATOMS
      INTEGER, ALLOCATABLE :: ATOMICNUMBERS(:)
      DOUBLE PRECISION, ALLOCATABLE :: APOS(:,:)
      DOUBLE PRECISION :: Tol,ETOT,NucE,Rn(3)
      CHARACTER(LEN=20) :: CORRLEVEL
      INTEGER :: BASISMAP(120)
      TYPE(ATOM), ALLOCATABLE :: ATOMS(:)
      TYPE(BASIS) :: BAS
      CHARACTER(LEN=6) :: DUMMY
      CHARACTER(Len=20) :: date,time,zone
      INTEGER :: NLINES,I,J,K,L,M,NB,Ne,NRED,Lmax,MESH(3),REDISTRIBUTIONFREQ,NPERSIST,N0,Istart,Iend,Istartg,Iendg,NVMC,IOSA,FNATOMS,FNRED,NSTEPS,NLSPOINTS,PORDER
      LOGICAL :: finns,WRITECICOEF,WRITEDENS,WHOMOLUMO,LEXCSP,SPINCONSERVE,RESTRICT,APPROXEE,NOREDIST,HYBRID,USEGTO,CORRALCUSP,VMCCALC,HFORBWRITE,CFORCE,RELAXN
      LOGICAL :: WRITEONFLY,MOVIE,MOLDYN
      DOUBLE PRECISION, ALLOCATABLE :: S(:,:),T(:,:),V(:,:),H0(:,:),Intsv(:),IntsvR(:),EIGENVECT(:,:),Ints(:,:,:,:),gradIntsvR(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE :: gradS(:,:,:,:),gradT(:,:,:,:),gradV(:,:,:,:),DMAT(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigenup(:),EHFeigendown(:),Cup(:,:),Cdown(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: EHFeigen(:),P(:,:),Pup(:,:),Pdown(:,:),C1(:,:),C2(:,:),force(:,:),CN(:),LQ(:,:),CGQ(:,:)
      DOUBLE PRECISION  :: PRYSR(25,25),PRYSW(25,25),rts(25),wts(25),TOTALTIME,ECISD,LIMITS(3),ENEXCM,EMP2,TIMESTEP,TEND,TSTART,CUTTOFFFACTOR,MIX,Aexpan,LAMDA,rc
      DOUBLE PRECISION :: ESIGMA,EDQMC,a,b,BETA,BJASTROW,CJASTROW,GAMA,BETAA,POLY(6),deltar,DR,FTol,TEMPERATURE,kappa,alpha,EETOL,CNSTART(4),alphastart,kappastart
      INTEGER :: STARTTIME(8),STOPTIME(8),RUNTIME(8),NEEXC,SAMPLERATE,NREPLICAS,NRECALC,scount,DIISORD,DIISSTART,NTIMESTEPS,NSCF,PULAY,FIXNSCF,DORDER
      INTEGER :: NLEBEDEV,NCHEBGAUSS,LEBPOINTS(24),CGORDER,LORDER,NTOTALQUAD,Qstart,Qend,RELALGO
      INTEGER, EXTERNAL :: ijkl
      INTEGER, ALLOCATABLE :: N0p(:),rcounts(:),displs(:),Istart2(:),Istart3(:),IND1(:),IND2(:),IND3(:),IND4(:),Q1(:),Q2(:),Q3(:),Q11(:),Q22(:),Q33(:)
      DOUBLE PRECISION, EXTERNAL :: massa
      LOGICAL :: RESTART,ZEROSCF,XLBOMD,DFTC,SOFTSTART

      call MPI_Init ( ierr )
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
      call MPI_Comm_size ( MPI_COMM_WORLD, numprocessors, ierr )

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
      
      ALLOCATE(ATOMICNUMBERS(NATOMS),APOS(NATOMS,3),ATOMS(NATOMS),force(NATOMS,3))
      ALLOCATE(N0p(numprocessors),rcounts(numprocessors),displs(numprocessors),Istart2(numprocessors),Istart3(numprocessors))
      
      CALL readin(CORRLEVEL,NATOMS,NLINES,ATOMICNUMBERS,APOS,Ne,Tol,WRITECICOEF,WRITEDENS,WHOMOLUMO,MESH,LIMITS,LEXCSP,NEEXC,ENEXCM,SPINCONSERVE,RESTRICT,APPROXEE, &
      & SAMPLERATE,NREPLICAS,TIMESTEP,TEND,TSTART,BETA,BJASTROW,CJASTROW,NPERSIST,REDISTRIBUTIONFREQ,NOREDIST,NRECALC,CUTTOFFFACTOR,MIX,DIISORD,DIISSTART,HYBRID,rc,CORRALCUSP,NVMC, &
      & HFORBWRITE,IOSA,CFORCE,RELAXN,NSTEPS,DR,FTol,NLSPOINTS,PORDER,WRITEONFLY,MOVIE,MOLDYN,TEMPERATURE,ZEROSCF,XLBOMD,kappa,alpha,DORDER,PULAY,FIXNSCF,NLEBEDEV,NCHEBGAUSS,EETOL,RELALGO,SOFTSTART)
      
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
              CALL lebedev(NLEBEDEV,LORDER,LQ)
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

      !IF ( ( RELAXN .OR. MOLDYN ) .AND. CORRLEVEL .NE. 'RHF' .AND. CORRLEVEL .NE. 'URHF'  ) THEN
      !        IF ( id .EQ. 0 ) THEN
      !                WRITE(*,*)'Relaxation/molecular-dynamics of nuclea only possible for RHF or URHF'
      !                WRITE(*,*)'Change CORRLEVEL to RHF or URHF'
      !        ENDIF
      !        STOP
      !ENDIF
      
      IF ( RELAXN .OR. MOLDYN ) CFORCE = .TRUE.
      
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
        !  In the MPI-version the full symetry is not employed. Here the number
        ! of symmetry reduced elements are twice as many as in the openmp and serial version.

        NRED =  ( (NB*(NB+1))/2)**2
        
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
     
        !===========================================================================
        ! PREPARATION FOR THE PARALELL COMPUTATION OF THE ee-repulsion-tensor Intsv
        !============================================================================
        ! Here we distribute the calculation of the columns of the ee-integral array
        ! Intsv over the (numprocessors) number of mpi threads.
        !----------------------------------------------------------------------------

        N0 = INT(NRED/numprocessors)

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
        DO WHILE( I .LE. MOD(NRED,numprocessors) )
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
        
        IF ( CFORCE ) THEN
                Istartg = Istart
                Iendg = Iend
                FNATOMS = NATOMS
        ELSE
                Istartg = 1
                Iendg = 1
                FNATOMS = 1
        ENDIF

        ALLOCATE(IntsvR(Istart:Iend),gradIntsvR(FNATOMS,3,Istartg:Iendg),IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend))
        
        DMAT = 1.0d0

        CALL eeints(NATOMS,FNATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istartg,Iendg,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,CFORCE,id,DMAT)
        
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
        

        ALLOCATE(H0(NB,NB))
        
        H0 = T+V

        IF ( RELAXN ) THEN
                ! Relxing by searching for minimum energy
                IF ( RELALGO .EQ. 2 ) THEN
                        CALL relax(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS, &
                                & APPROXEE,CORRLEVEL,PRYSR,PRYSW,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
                ! Relaxing by searching for FORCES = 0
                IF ( RELALGO .EQ. 1 ) THEN
                        CALL relaxf(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,DR,BAS,ATOMS, &
                                & APPROXEE,CORRLEVEL,PRYSR,PRYSW,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
        ENDIF

        IF ( MOLDYN ) THEN
                IF ( id .EQ. 0 ) THEN
                        print*,' ========================================================'
                        print*,'           Comencing with molecular dynamics calculation '
                        print*,' ========================================================'
                ENDIF
                NTIMESTEPS = INT(TEND/TIMESTEP)
                IF ( SOFTSTART ) THEN
                   CALL moleculardynamicssoft(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,TIMESTEP,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                   & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,&
                   & Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL,CNSTART,alphastart,kappastart)
                ELSE
                   CALL moleculardynamics(gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,NB,NRED,Ne,nucE,Tol,MIX,DIISORD,DIISSTART,NATOMS,NTIMESTEPS,TIMESTEP,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
                   & WRITEONFLY,MOVIE,SAMPLERATE,TEMPERATURE,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,ZEROSCF,XLBOMD,kappa,alpha,CN,PULAY,FIXNSCF,DORDER,&
                   & Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL)
                ENDIF
                                                                                                             
        ENDIF

        !==============================================================
        ! (5) Here the Unrestricted or Restricted Hartree-FocK
        !     equations are selfconsistently solved
        !==============================================================
        IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                
                ALLOCATE(EIGENVECT(NB,NB),EHFeigen(NB),P(NB,NB),C1(NB,NB))
                P = 0.0d0
                CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.TRUE.,.TRUE.,.FALSE.)

                !====================
                ! CALCULATING FORCES
                !================================================================================================================== 
                IF ( CFORCE ) THEN
                        P = P/2.0d0
                        CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,&
                             & S,H0,IntsvR,PULAY,force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL)
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
                
                ALLOCATE(EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),P(NB,NB),Pup(NB,NB),Pdown(NB,NB),C1(NB,NB),C2(NB,NB))
                
                IF ( .not. RESTRICT .OR. CORRLEVEL .EQ. 'DQMC' .OR. CORRLEVEL .EQ. 'VMC' .OR. DFTC )  THEN
                   Pup = 0.0d0
                   Pdown = 0.0d0
                   IF ( .not. DFTC ) CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1, &
                             & numprocessors,id,.TRUE.,.TRUE.,.FALSE.)


                   IF ( DFTC ) CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD, &
                   & NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.TRUE.,.TRUE.,.FALSE.)

                   !====================
                   ! CALCULATING FORCES
                   !================================================================================================================== 
                   IF ( CFORCE ) THEN

                          CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY,&
                               & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL)
                          
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
                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.TRUE.,.TRUE.,.FALSE.)

                        !====================
                        ! CALCULATING FORCES
                        !================================================================================================================== 
                        IF ( CFORCE ) THEN
                           p = P/2.0d0
                           CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,EIGENVECT,EIGENVECT,P,P,EHFeigen,EHFeigen,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY, &
                                & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL)
                        
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
                
                
                IF ( CORRLEVEL .EQ. 'CISD' .OR. CORRLEVEL .EQ. 'MP2' ) THEN
                        !print*,'========================================================'
                        !print*,'  Tranfering (ij|kl) from vector to tensor form         '
                        !print*,'========================================================'
                        ALLOCATE(Intsv(NRED))
                        Intsv(:) = 0.0d0

                       ! Here the Intsv is gathered together so that the entire tensor is available at all threads.
                       ! In practice the doubles IntsvR(Istart:Iend) is sent from thread id to Intsv at thread = 0
                       call MPI_GATHERV(IntsvR(Istart:Iend),scount,MPI_DOUBLE_PRECISION,Intsv,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        DEALLOCATE(IntsvR)
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        
                        CALL MPI_BCAST(Intsv,NRED,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        
                        ALLOCATE(Ints(NB,NB,NB,NB))
               
                        DO I=1,NB
                                DO J=1,NB
                                        DO K=1,NB
                                                DO L=1,NB
                                                        Ints(I,J,K,L) = Intsv(ijkl(I,J,K,L,NB))
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDDO
                       
                       DEALLOCATE(Intsv)

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
       CALL MPI_FINALIZE(ierr)
END PROGRAM uquantchem
