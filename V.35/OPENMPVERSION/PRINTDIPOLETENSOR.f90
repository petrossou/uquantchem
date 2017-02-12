SUBROUTINE PRINTDIPOLETENSOR(BAS,Cup,Cdown)
      ! This subroutine prints the energy eigenvalues and corresponding
      ! angular-momentum character to file
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASIS), INTENT(IN)  :: BAS
      DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION :: BDPTENSOR(3,BAS%NBAS,BAS%NBAS),DPTESORu(3,BAS%NBAS,BAS%NBAS),DPTESORd(3,BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,NB

      NB = BAS%NBAS
        
     CALL dipoletensor(BAS,BDPTENSOR)
     CALL makeltesor(Cup,BDPTENSOR,NB,DPTESORu)
     CALL makeltesor(Cdown,BDPTENSOR,NB,DPTESORd)
     
     OPEN(21,FILE='DIPOLETENSOR.dat',ACTION='WRITE')
     
     WRITE(21,'(A164)')'===================================================================================================================================================================='
     WRITE(21,'(A164)')'                                              ANGULAR MOMENTUM (UP)                     |                          DIPOLE MOMENT (down)                            |'
     WRITE(21,'(A164)')'===================================================================================================================================================================='
     WRITE(21,'(A164)')'   (I,J)                 <Psi(I)|x|Psi(J)>    <Psi(I)|y|Psi(J)>    <Psi(I)|z|Psi(J)>    |      <Psi(I)|x|Psi(J)>    <Psi(I)|y|Psi(J)>    <Psi(I)|z|Psi(J)>          '
     WRITE(21,'(A164)')'===================================================================================================================================================================='
     DO I=1,NB
        DO J=I,NB
                WRITE(21,'(I4,I4,6(F25.15))')I,J,DPTESORu(1,I,J),DPTESORu(2,I,J),DPTESORu(3,I,J),DPTESORd(1,I,J),DPTESORd(2,I,J),DPTESORd(3,I,J)
        ENDDO
     ENDDO
     CLOSE(21)
END SUBROUTINE PRINTDIPOLETENSOR

