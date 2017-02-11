SUBROUTINE PRINTENERGYEIGENVAL(BAS,EHFeigenup,EHFeigendown,Cup,Cdown,EXITED)
      ! This subroutine prints the energy eigenvalues and corresponding
      ! angular-momentum character to file
      USE datatypemodule
      IMPLICIT NONE
      TYPE(BASIS), INTENT(IN)  :: BAS
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(BAS%NBAS),EHFeigendown(BAS%NBAS),Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      LOGICAL, INTENT(IN) :: EXITED
      DOUBLE PRECISION :: BLTENSOR(3,BAS%NBAS,BAS%NBAS),LTESORu(3,BAS%NBAS,BAS%NBAS),LTESORd(3,BAS%NBAS,BAS%NBAS),L2u(BAS%NBAS,BAS%NBAS),L2d(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION :: Lx2u(BAS%NBAS,BAS%NBAS),Ly2u(BAS%NBAS,BAS%NBAS),Lz2u(BAS%NBAS,BAS%NBAS),Lx2d(BAS%NBAS,BAS%NBAS),Ly2d(BAS%NBAS,BAS%NBAS),Lz2d(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,NB

      NB = BAS%NBAS
        
     CALL orbitalmomentumtensor(BAS,BLTENSOR)
     CALL makeltesor(Cup,BLTENSOR,NB,LTESORu)
     CALL makeltesor(Cdown,BLTENSOR,NB,LTESORd)
     Lx2u = -MATMUL(LTESORu(1,:,:),LTESORu(1,:,:))
     Ly2u = -MATMUL(LTESORu(2,:,:),LTESORu(2,:,:)) 
     Lz2u = -MATMUL(LTESORu(3,:,:),LTESORu(3,:,:))
     Lx2d = -MATMUL(LTESORd(1,:,:),LTESORd(1,:,:)) 
     Ly2d = -MATMUL(LTESORd(2,:,:),LTESORd(2,:,:))
     Lz2d = -MATMUL(LTESORd(3,:,:),LTESORd(3,:,:))
     L2u = Lx2u + Ly2u + Lz2u
     L2d = Lx2d + Ly2d + Lz2d
     
     IF (.not. EXITED ) OPEN(21,FILE='EIGENVALUES.dat',ACTION='WRITE')
     IF ( EXITED ) OPEN(21,FILE='EXCITEDEIGENVALUES.dat',ACTION='WRITE')
     
     WRITE(21,'(A164)')'==================================================================================================================================================================='
     WRITE(21,'(A164)')'                       ENERGY EIGENVALUES:              |                 ANGULAR MOMENTUM (UP)               |            ANGULAR MOMENTUM (down)                 |'
     WRITE(21,'(A164)')'===================================================================================================================================================================='
     WRITE(21,'(A164)')'   N            Eup [a.u.]              Edown [a.u.]    |      Lx**2          Ly**2          Lz**2      L**2  |    Lx**2          Ly**2          Lz**2      L**2   |'
     WRITE(21,'(A164)')'===================================================================================================================================================================='
     DO I=1,NB
        WRITE(21,'(I4,2(F25.15),3(F15.6),I7,3(F15.6),I7)')I,EHFeigenup(I),EHFeigendown(I),Lx2u(I,I),Ly2u(I,I),Lz2u(I,I),NINT(L2u(I,I)),Lx2d(I,I),Ly2d(I,I),Lz2d(I,I),NINT(L2d(I,I))
     ENDDO
     CLOSE(21)
END SUBROUTINE PRINTENERGYEIGENVAL

