SUBROUTINE hforbitalsave(N,LIMITS,MESH,BAS,Cup,Cdown)
      ! This subroutine saves the hartree-fock orbital with idex, N, on disk.
      ! Input that is needed is the MESH, i.e how many grid points are
      ! going to be used in the x, y and z directions to discretizise
      ! the density. The LIMITS limits in space for the calculated
      ! charge density: x in [-LIMITS(1),LIMITS(1)]  ,  
      ! y in [-LIMITS(1),LIMITS(1)] , z in [-LIMITS(1),LIMITS(1)].
      USE datatypemodule
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      INTEGER, INTENT(IN) :: MESH(3)
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,K
      DOUBLE PRECISION :: R(3),wfnktup,wfnktdown
      DOUBLE PRECISION, EXTERNAL :: hforbitalval

      2 FORMAT(I6,I6,I6,F50.20,F50.20)
      OPEN(23,FILE='ORBITAL.dat',ACTION='WRITE')
      DO I=1,MESH(1)
                R(1) = -LIMITS(1) + 2.0d0*LIMITS(1)*(I-1)/(1.0d0*(MESH(1)-1))
                DO J=1,MESH(2)
                        R(2) = -LIMITS(2) + 2.0d0*LIMITS(2)*(J-1)/(1.0d0*(MESH(2)-1))
                        DO K=1,MESH(3)
                                R(3) = -LIMITS(3) + 2.0d0*LIMITS(3)*(K-1)/(1.0d0*(MESH(3)-1))
                                wfnktup = hforbitalval(BAS,Cup,N,R)
                                wfnktdown = hforbitalval(BAS,Cdown,N,R)
                                WRITE(23,FMT=2)I,J,K,wfnktup,wfnktdown
                        ENDDO
                ENDDO
       ENDDO
       close(23)
END SUBROUTINE hforbitalsave

      
