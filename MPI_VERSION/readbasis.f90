SUBROUTINE readbasis(NATOMS,ATOMS,BASISMAP,MAP)
      ! If MAP = .TRUE. this routine opens the file "BASISFILE" and 
      ! loads the information about the number of basis functions used 
      ! for each atom in the periodic table

      ! If MAP = .FALSE. this routine loads the basis set information
      ! for the all the atoms stored in the "array" ATOMS

      USE datatypemodule
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS
      TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS)
      INTEGER, INTENT(INOUT) :: BASISMAP(120)
      LOGICAL, INTENT(IN) :: MAP
      LOGICAL :: STORE
      INTEGER :: LAS,NLINES,I,J,K,Z,TEMP
      CHARACTER(LEN=2) :: ATYPE,DUMMY
      CHARACTER(LEN=1) :: LMOM
      CHARACTER(LEN=200) :: LINE,DUMMY2
      INTEGER, ALLOCATABLE :: PRIMS(:)
     
      1 FORMAT(A200)
      IF ( MAP ) THEN
              OPEN(20,FILE='BASISFILE',STATUS='OLD',ACTION='READ')
              NLINES = 0
              
              ! Two first lines in the BASISFILE should alwas contain
              ! Name of basis and the string "BASIS:", which is not
              ! any information used here

              READ(20,*)DUMMY
              READ(20,*)DUMMY
              
              DO WHILE (.TRUE. )
                        IF ( NLINES .EQ. 0 ) THEN
                                READ(20,FMT=1,end=1000)LINE
                                READ(LINE,*)ATYPE,Z,TEMP
                                ALLOCATE(PRIMS(TEMP))
                                READ(LINE,*)ATYPE,Z,BASISMAP(Z),PRIMS
                                !WRITE(*,FMT=1)LINE
                                NLINES = 1 + BASISMAP(Z) + SUM(PRIMS)
                                DEALLOCATE(PRIMS)
                        ELSE
                                READ(20,FMT=1,end=1000)LINE
                                NLINES = NLINES -1
                        ENDIF
              ENDDO
              1000 CONTINUE
              CLOSE(20)
      ELSE
                OPEN(20,FILE='BASISFILE',STATUS='OLD',ACTION='READ')
                DO K=1,NATOMS
                        NLINES = 0

                        ! Two first lines in the BASISFILE should alwas contain
                        ! Name of basis and the string "BASIS:", which is not
                        ! any information used here

                        READ(20,*)DUMMY
                        READ(20,*)DUMMY
              
                        DO WHILE (.TRUE.)
                                IF ( NLINES .EQ. 0 ) THEN
                                        READ(20,FMT=1,end=2000)LINE
                                        READ(LINE,*)ATYPE,Z,TEMP
                                        ALLOCATE(PRIMS(TEMP))
                                        READ(LINE,*)ATYPE,Z,BASISMAP(Z),PRIMS
                                        
                                        IF ( ATOMS(K)%Z  .EQ. Z ) THEN
                                                ATOMS(K)%NBAS = BASISMAP(Z)
                                                DO I=1,BASISMAP(Z)
                                                        ATOMS(K)%PSI(I)%NPRIM = PRIMS(I)
                                                ENDDO
                                        ENDIF
                                        
                                        NLINES = 1 + BASISMAP(Z) + SUM(PRIMS)
                                ELSE
                                        STORE = .FALSE.
                                        IF ( ATOMS(K)%Z  .EQ. Z ) STORE = .TRUE.

                                        READ(20,FMT=1,end=2000)LINE
                                
                                        DO I=1,BASISMAP(Z)
                                                        DO J=1,PRIMS(I)
                                                        IF ( STORE ) READ(20,*,end=2000)LMOM,ATOMS(K)%PSI(I)%EXPON(J),ATOMS(K)%PSI(I)%CONTRCOEFF(J)
                                                                IF ( .not. STORE)READ(20,FMT=1,end=2000)LINE
                                                        ENDDO
                                                        IF ( STORE ) THEN
                                                                ATOMS(K)%PSI(I)%L(:) = 0
                                                                IF ( LMOM .EQ. 's' ) ATOMS(K)%PSI(I)%L(1) = 0
                                                                IF ( LMOM .EQ. 'p' ) ATOMS(K)%PSI(I)%L(1) = 1
                                                                IF ( LMOM .EQ. 'd' ) ATOMS(K)%PSI(I)%L(1) = 2
                                                                IF ( LMOM .EQ. 'f' ) ATOMS(K)%PSI(I)%L(1) = 3
                                                                IF ( LMOM .EQ. 'g' ) ATOMS(K)%PSI(I)%L(1) = 4
                                                                IF ( LMOM .EQ. 'h' ) ATOMS(K)%PSI(I)%L(1) = 5
                                                                IF ( LMOM .EQ. 'i' ) ATOMS(K)%PSI(I)%L(1) = 6
                                                                IF ( LMOM .EQ. 'j' ) ATOMS(K)%PSI(I)%L(1) = 7
                                                                IF ( LMOM .EQ. 'k' ) ATOMS(K)%PSI(I)%L(1) = 8
                                                        ENDIF
                                                        READ(20,FMT=1,end=2000)LINE
                                        ENDDO
                                        NLINES = 0
                                        DEALLOCATE(PRIMS)
                                ENDIF
                        ENDDO
                        2000 CONTINUE
                        REWIND(20)
                ENDDO
                CLOSE(20)
      ENDIF

END SUBROUTINE readbasis
