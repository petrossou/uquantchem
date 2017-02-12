SUBROUTINE arbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS)
      ! This subroutine saves the hartree-fock charge-density on disk.
      ! Input that is needed is the MESH, i.e how many grid points are
      ! going to be used in the x, y and z directions to discretizise
      ! the density. The LIMITS limits in space for the calculated
      ! charge density: x in [-LIMITS(1),LIMITS(1)]  ,  
      ! y in [-LIMITS(1),LIMITS(1)] , z in [-LIMITS(1),LIMITS(1)].
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      INTEGER, INTENT(IN) :: MESH(3),Ne,NATOMS,AORBS
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(BAS%NBAS),EHFeigendown(BAS%NBAS)
      DOUBLE PRECISION :: rho
      DOUBLE PRECISION :: Phomo(BAS%NBAS,BAS%NBAS),Plumo(BAS%NBAS,BAS%NBAS),Ch(BAS%NBAS,BAS%NBAS),Cl(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,K,Nup,Ndown,Ih,Il
      DOUBLE PRECISION :: R(3),NORM,MAXHOMO,MAXLUMO,denshomo,denslumo,homo,lumo
      DOUBLE PRECISION, EXTERNAL :: hforbitalval
       print*,' '
       print*,'============================================================='
       print*,'          ORBITAL WITH INDEX = AORBS IS SAVED TO FILE        '
       print*,'============================================================='
       print*,' '
       
       MAXHOMO = -1000000.0d0
       MAXLUMO = -1000000.0d0
       Phomo = 0.0d0
       Plumo = 0.0d0
       Ch = 0.0d0
       Cl = 0.0d0
      
      IF ( AORBS .LT. 0  ) THEN
             Ch(:,-AORBS)  = Cdown(:,-AORBS)
      ELSE
             Ch(:,AORBS)  = Cup(:,AORBS)
      ENDIF
      
      2 FORMAT(I6,I6,I6,F30.20)
      
      OPEN(100,FILE='ARBORB.xsf',ACTION='WRITE')
      
      WRITE(100,'(A5)')'ATOMS'
      
      DO I=1,NATOMS
                WRITE(100,'(I4,3(F15.10))')ATOMS(I)%Z,ATOMS(I)%R(1)*0.52917720859,ATOMS(I)%R(2)*0.52917720859,ATOMS(I)%R(3)*0.52917720859
      ENDDO
      
      WRITE(100,'(A22)')'BEGIN_BLOCK_DATAGRID3D'
      WRITE(100,'(A15)')'3D_ORBITAL_#001'
      WRITE(100,'(A15)')'DATAGRID_3D_G98CUBE'

      WRITE(100,'(3(I6))')MESH

      WRITE(100,'(3(F15.10))') -LIMITS(1)*0.52917720859,-LIMITS(2)*0.52917720859,-LIMITS(3)*0.52917720859

      
      WRITE(100,'(3(F15.10))') 2*LIMITS(1)*0.52917720859,0.0d0,0.0d0
      WRITE(100,'(3(F15.10))') 0.0d0,2*LIMITS(2)*0.52917720859,0.0d0
      WRITE(100,'(3(F15.10))') 0.0d0, 0.0d0,2*LIMITS(3)*0.52917720859
      
      DO K=1,MESH(3)
                R(3) = -LIMITS(3) + 2.0d0*LIMITS(3)*(K-1)/(1.0d0*(MESH(3)-1))
                DO J=1,MESH(2)
                        R(2) = -LIMITS(2) + 2.0d0*LIMITS(2)*(J-1)/(1.0d0*(MESH(2)-1))
                        DO I=1,MESH(1)
                                R(1) = -LIMITS(1) + 2.0d0*LIMITS(1)*(I-1)/(1.0d0*(MESH(1)-1))
                                homo = hforbitalval(BAS,Ch,ABS(AORBS),R)
                                WRITE(100,'(F18.15)')homo
                        ENDDO
                ENDDO
       ENDDO
       WRITE(100,'(A15)')'END_DATAGRID_3D'
       WRITE(100,'(A20)')'END_BLOCK_DATAGRID3D'
       close(100)
       print*,' '
       WRITE(*,'(A77)')'============================================================================='
       IF ( AORBS .LT. 0 ) THEN
               WRITE(*,'(A34,I4,A41)')'The spin-down orbital with index =',ABS(AORBS),' have been written to the file ARBORB.xsf'
       ELSE
               WRITE(*,'(A32,I4,A41)')'The spin-up orbital with index =',ABS(AORBS),' have been written to the file ARBORB.xsf'
       ENDIF
       WRITE(*,'(A77)')'============================================================================='
       print*,' '
       print*,' '
       WRITE(*,'(A77)')'============================================================================='
       print*,' '
END SUBROUTINE arbitaryorbitalsave

      
