SUBROUTINE  Rdensmatstartguess(ATOMS,NATOMS,NB,Pup,Pdown,SCRATCH)
        USE datatypemodule
        IMPLICIT NONE
        TYPE(ATOM), INTENT(IN)  :: ATOMS(NATOMS)
        INTEGER, INTENT(IN) :: NATOMS,NB
        DOUBLE PRECISION, INTENT(INOUT) :: Pup(NB,NB),Pdown(NB,NB)
        LOGICAL, INTENT(INOUT) :: SCRATCH
        DOUBLE PRECISION  :: Pu(120,NB,NB),Pd(120,NB,NB)
        INTEGER :: I,J,K,L,M,Z,NBB,ATOMICNUMBER,NBFMAP(120),OFFSET,BASSUM
        LOGICAL :: finns
        CHARACTER(LEN=200) :: LINE
     
        
        1 FORMAT(A200)

        BASSUM = 0
        NBFMAP = 0

        Pup = 0.0d0
        Pdown = 0.0d0
        
        Pd = 0.0d0
        Pu = 0.0d0 
        inquire(file='DENSMATSTARTGUESS.dat',exist=finns)
        IF ( finns ) THEN
                SCRATCH = .FALSE.
                OPEN(25,FILE='DENSMATSTARTGUESS.dat',STATUS='OLD',ACTION='READWRITE')
                DO WHILE (.TRUE. ) 
                        READ(25,FMT=1,end=1000)LINE
                        READ(LINE,'(I4,I6)')ATOMICNUMBER,NBB
                        NBFMAP(ATOMICNUMBER) = NBB
                        DO I=1,NBB
                                DO J=I,NBB
                                        READ(25,FMT=1,end=1000)LINE
                                        IF ( I .LE. NB .AND. J .LE. NB ) THEN
                                                READ(LINE,'(E30.20,E30.20,I6,I6)')Pu(ATOMICNUMBER,I,J),Pd(ATOMICNUMBER,I,J),L,M
                                        ENDIF
                                ENDDO
                        ENDDO
                ENDDO
        1000    CONTINUE
                CLOSE(25)
               DO K=1,NATOMS
                        BASSUM = BASSUM + NBFMAP(ATOMS(K)%Z)
               ENDDO
               
               IF ( BASSUM .NE. NB ) THEN
                                WRITE(*,*)' '
                                WRITE(*,*)'***********************************************************'
                                WRITE(*,*)'***                      WARNING!                       ***'
                                WRITE(*,*)'     The sum of the number of atomic basis functions       '
                                WRITE(*,*)'     used to create the start-guess density matrices       ' 
                                WRITE(*,*)'do not equal the total number of molecular basis functions.'
                                WRITE(*,*)'Please check that you have generated the start density file'
                                WRITE(*,*)'DENSMATSTARTGUESS.dat using the same basis-set as in the   '
                                WRITE(*,*)'current calculation.                                       '
                                WRITE(*,*)'             THIS CALCULATION WILL BE DONE FROM SCRATCH!   '
                                WRITE(*,*)'***********************************************************'
                                WRITE(*,*)' '
                                SCRATCH = .TRUE.
                                RETURN
                ENDIF
                OFFSET = 0
                DO K=1,NATOMS
                        IF ( K .GT. 1 ) OFFSET = OFFSET + NBFMAP(ATOMS(K-1)%Z)
                        NBB = NBFMAP(ATOMS(K)%Z)
                        DO I=1,NBB
                                DO J=I,NBB
                                        IF ( J + OFFSET .LE. NB .AND. I + OFFSET .LE. NB ) THEN
                                                Pup(I + OFFSET, J + OFFSET)   = Pu(ATOMS(K)%Z,I,J)
                                                Pdown(I + OFFSET, J + OFFSET) = Pd(ATOMS(K)%Z,I,J)
                                                IF ( I .NE. J ) THEN
                                                        Pup(J + OFFSET, I + OFFSET)   = Pup(I + OFFSET, J + OFFSET)
                                                        Pdown(J + OFFSET, I + OFFSET) = Pdown(I + OFFSET, J + OFFSET)
                                                ENDIF
                                        ENDIF
                                ENDDO
                        ENDDO
                ENDDO
        ENDIF
        
END SUBROUTINE  Rdensmatstartguess
                
                
                
