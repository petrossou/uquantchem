      SUBROUTINE blockset( NB,MB, NGLOBAL, NPROW, NPCOL )
!
!     This subroutine chooses a block size
!     for the distributd matrix.
!
!     Written by Petros Souvatzis Uppsala University
!
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: NB,MB
      INTEGER, INTENT(IN) :: NPROW, NPCOL, NGLOBAL

      NB = NGLOBAL/NPROW
      MB = NGLOBAL/NPCOL

      END SUBROUTINE blockset
