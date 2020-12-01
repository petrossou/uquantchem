SUBROUTINE  Wdensmatstartguess(ATOMICNUMBER,NB,Pup,Pdown)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ATOMICNUMBER,NB
        DOUBLE PRECISION, INTENT(IN) :: Pup(NB,NB),Pdown(NB,NB)
        INTEGER :: I,J
        LOGICAL :: finns
        CHARACTER(LEN=200) :: LINE
        
        1 FORMAT(A200)
        inquire(file='DENSMATSTARTGUESS.dat',exist=finns)
        
        IF ( finns ) THEN
                OPEN(25,FILE='DENSMATSTARTGUESS_LATEST.dat',ACTION='READWRITE')
                WRITE(25,'(I4,I6)')ATOMICNUMBER,NB
                DO I=1,NB
                        DO J=I,NB
                                WRITE(25,'(E30.20,E30.20,I6,I6)')Pup(I,J),Pdown(I,J),I,J
                        ENDDO
                ENDDO
        ELSE
                OPEN(25,FILE='DENSMATSTARTGUESS.dat',STATUS='NEW',ACTION='READWRITE')
                WRITE(25,'(I4,I6)')ATOMICNUMBER,NB
                DO I=1,NB
                        DO J=I,NB
                                WRITE(25,'(E30.20,E30.20,I6,I6)')Pup(I,J),Pdown(I,J),I,J
                        ENDDO
                ENDDO
        ENDIF
END SUBROUTINE  Wdensmatstartguess
                
                
                
