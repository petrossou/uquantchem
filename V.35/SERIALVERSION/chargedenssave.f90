SUBROUTINE chargedenssave(LIMITS,MESH,BAS,P,Ne)
      ! This subroutine saves the hartree-fock charge-density on disk.
      ! Input that is needed is the MESH, i.e how many grid points are
      ! going to be used in the x, y and z directions to discretizise
      ! the density. The LIMITS limits in space for the calculated
      ! charge density: x in [-LIMITS(1),LIMITS(1)]  ,  
      ! y in [-LIMITS(1),LIMITS(1)] , z in [-LIMITS(1),LIMITS(1)].
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      INTEGER, INTENT(IN) :: MESH(3),Ne
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, EXTERNAL :: rho
      INTEGER :: I,J,K
      DOUBLE PRECISION :: R(3),NORM,MAXRHO,dens

      NORM = 1.0d0
      MAXRHO = -1000000.0d0

      2 FORMAT(I6,I6,I6,F30.20)
      OPEN(9,FILE='CHARGEDENS.dat',ACTION='WRITE')
      DO I=1,MESH(1)
                R(1) = -LIMITS(1) + 2.0d0*LIMITS(1)*(I-1)/(1.0d0*(MESH(1)-1))
                DO J=1,MESH(2)
                        R(2) = -LIMITS(2) + 2.0d0*LIMITS(2)*(J-1)/(1.0d0*(MESH(2)-1))
                        DO K=1,MESH(3)
                                R(3) = -LIMITS(3) + 2.0d0*LIMITS(3)*(K-1)/(1.0d0*(MESH(3)-1))
                                dens = NORM*rho(BAS,P,R)
                                IF ( MAXRHO .LT. dens ) MAXRHO = dens
                                WRITE(9,FMT=2)I,J,K,dens
                        ENDDO
                ENDDO
       ENDDO
       close(9)
       print*,' '
       print*,'   ========================================================'
       print*,'     Charge density has been calculated and saved to file. '
       print*,'   ========================================================'
       print*,' '
       WRITE(*,'(A37,F25.20)')'  Max charge density value, rho_max =',MAXRHO
       print*,' '
       print*,'   ========================================================'
       print*,' '
END SUBROUTINE chargedenssave

      
