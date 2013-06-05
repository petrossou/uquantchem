SUBROUTINE checkinputfile(NLINES,NATOMS)
      ! This subroutine checks wheater or not the file 'INPUTFILE'
      ! exists and if so, how many lines, NLINES, it contains before the entry of
      ! the number of atoms,NATOMS.
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: NLINES,NATOMS
      INTEGER :: I
      LOGICAL :: finns
      CHARACTER(LEN=6) :: DUMMY

      inquire(file='INPUTFILE',exist=finns)
      IF ( finns ) THEN
              OPEN(10,FILE='INPUTFILE',STATUS='OLD',ACTION='READ')
              NLINES = 0
              DO WHILE( DUMMY .NE. 'NATOMS' )
                        NLINES = NLINES+1
                        READ(10,*)DUMMY
              ENDDO
              rewind(10)
              DO I=1,NLINES
                        IF ( I .LT. NLINES ) READ(10,*)DUMMY
                        IF ( I .EQ. NLINES ) READ(10,*)DUMMY,NATOMS
              ENDDO
              CLOSE(10)
        ELSE
                print*,'ABORTING SINCE THE FILE: "INPUTFILE" IS MISSING'
                STOP
        ENDIF
END SUBROUTINE checkinputfile
