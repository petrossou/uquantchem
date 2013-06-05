FUNCTION exciteornot(Ne,Nel,I,J)
      IMPLICIT NONE
      LOGICAL :: exciteornot
      INTEGER :: I,J,Ne,Nel,DN

      exciteornot = .FALSE.

      DN = Ne-Nel

      IF ( ( I .GT. ( DN-MOD(DN,2))/2 .AND. I .LE. ( Ne-MOD(Ne,2))/2 ) .OR. ( I .GT. ( Ne-MOD(Ne,2))/2+( DN+MOD(DN,2))/2 .AND. I .LE. Ne )) THEN
        IF ( ( J .GT. ( DN-MOD(DN,2))/2 .AND. J .LE. ( Ne-MOD(Ne,2))/2 ) .OR. ( J .GT. ( Ne-MOD(Ne,2))/2+( DN+MOD(DN,2))/2 .AND. J .LE. Ne )) exciteornot = .TRUE.
      ENDIF

END FUNCTION exciteornot

