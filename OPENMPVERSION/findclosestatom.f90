SUBROUTINE findclosestatom(NATOMS,ATOMS,rvect,z,zvect,atomnr,rnuc)
      USE datatypemodule
      ! All the numbers in parenthesis to the right are in reference 
      ! to the numbering of the steps in Urigars algorithm,
      ! of the so called improved algorithm on pages 2886-2887 in
      ! J. Chem. Phys. 99, 2865-2890, (1993)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: rvect(3)
      DOUBLE PRECISION, INTENT(OUT) :: z,zvect(3),rnuc(3)
      INTEGER, INTENT(OUT) :: atomnr
      INTEGER :: KK
      DOUBLE PRECISION :: distance,vect(3)

      !---------------------------------------------------
      ! Finding the nucleus closesest to the electron (11)
      !---------------------------------------------------
      
        zvect = 0.0d0

        DO KK=1,NATOMS
                vect = rvect-ATOMS(KK)%R
                distance = sqrt(DOT_PRODUCT(vect,vect))
                IF ( KK .EQ. 1 ) THEN
                        z = distance
                        IF ( distance .NE. 0.0d0 ) zvect = vect/distance
                        atomnr = ATOMS(KK)%Z
                        rnuc = ATOMS(KK)%R
                ENDIF
                IF ( KK .GT. 1 .AND. z .LT. distance ) THEN
                        z = distance
                        IF ( distance .NE. 0.0d0 ) zvect = vect/distance
                        atomnr = ATOMS(KK)%Z
                        rnuc = ATOMS(KK)%R
                ENDIF
        ENDDO
      
        END SUBROUTINE findclosestatom
