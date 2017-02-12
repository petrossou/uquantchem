SUBROUTINE mulliken(NATOMS,ATOMS,BAS,P,S,ZM)
        ! Calculates Mulliken charge of the different atoms
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: NATOMS
        DOUBLE PRECISION, INTENT(OUT) :: ZM(NATOMS) ! Mulliken charges
        INTEGER :: NBAS
        TYPE(BASIS), INTENT(IN) :: BAS
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),S(BAS%NBAS,BAS%NBAS)
        INTEGER :: I,J,K

        DO I=1,NATOMS
                ZM(I) = 1.0d0*ATOMS(I)%Z
                DO J=1,BAS%NBAS
                        IF ( BAS%PSI(J)%ATYPE .EQ. I ) THEN
                                DO K=1,BAS%NBAS
                                        ZM(I) = ZM(I) - P(J,K)*S(J,K)
                                ENDDO
                        ENDIF
                ENDDO
         ENDDO
         WRITE(*,*)
         WRITE(*,'(A55)')'      ================================================='
         WRITE(*,'(A49)')'          The Mulliken charges are the following:'
         WRITE(*,'(A55)')'      ================================================='
         DO I=1,NATOMS
                WRITE(*,'(A21,I3,A8,f12.8)')'                ATOM ',I,':  Zm = ',ZM(I)
         ENDDO
         WRITE(*,'(A55)')'      ================================================='
         WRITE(*,*)
        
END SUBROUTINE mulliken

