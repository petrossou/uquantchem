SUBROUTINE homolumosave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,DFTC)
      ! This subroutine saves the hartree-fock charge-density on disk.
      ! Input that is needed is the MESH, i.e how many grid points are
      ! going to be used in the x, y and z directions to discretizise
      ! the density. The LIMITS limits in space for the calculated
      ! charge density: x in [-LIMITS(1),LIMITS(1)]  ,  
      ! y in [-LIMITS(1),LIMITS(1)] , z in [-LIMITS(1),LIMITS(1)].
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      INTEGER, INTENT(IN) :: MESH(3),Ne,NATOMS
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(BAS%NBAS),EHFeigendown(BAS%NBAS)
      LOGICAL, INTENT(IN) :: DFTC
      DOUBLE PRECISION :: rho
      DOUBLE PRECISION :: Phomo(BAS%NBAS,BAS%NBAS),Plumo(BAS%NBAS,BAS%NBAS),Ch(BAS%NBAS,BAS%NBAS),Cl(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,K,Nup,Ndown,Ih,Il
      DOUBLE PRECISION :: R(3),NORM,MAXHOMO,MAXLUMO,denshomo,denslumo,homo,lumo
      DOUBLE PRECISION, EXTERNAL :: hforbitalval
       print*,' '
       print*,'============================================================='
       print*,'   HOMO and LUMO densities are being calculated and saved.   '
       print*,'============================================================='
       print*,' '
       
       MAXHOMO = -1000000.0d0
       MAXLUMO = -1000000.0d0
       Phomo = 0.0d0
       Plumo = 0.0d0
       Ch = 0.0d0
       Cl = 0.0d0
       
       IF ( .not. DFTC ) THEN
                Ndown = (Ne + MOD(Ne,2))/2
                Nup   = (Ne - MOD(Ne,2))/2
       ELSE
                Ndown = (Ne - MOD(Ne,2))/2
                Nup   = (Ne + MOD(Ne,2))/2
      ENDIF


      IF ( EHFeigendown(Ndown) .GT. EHFeigenup(Nup) ) THEN
             Ch(:,Ndown)  = Cdown(:,Ndown)
             Ih = Ndown
             CALL makedens(Ch,BAS%NBAS,Phomo)
      ELSE
             Ch(:,Nup)  = Cup(:,Nup)
             Ih = Nup
             CALL makedens(Ch,BAS%NBAS,Phomo)
      ENDIF

      IF ( EHFeigendown(Ndown+1) .LT. EHFeigenup(Nup+1) ) THEN
             Cl(:,Ndown+1)  = Cdown(:,Ndown+1)
             Il = Ndown+1
             CALL makedens(Cl,BAS%NBAS,Plumo)
      ELSE
             Cl(:,Nup+1)  = Cup(:,Nup+1)
             Il = Nup+1
             CALL makedens(Cl,BAS%NBAS,Plumo)
      ENDIF
      
      2 FORMAT(I6,I6,I6,F30.20)
      
      OPEN(10,FILE='HOMODENS.dat',ACTION='WRITE')
      OPEN(11,FILE='LUMODENS.dat',ACTION='WRITE')
      
      OPEN(100,FILE='HOMO.xsf',ACTION='WRITE')
      OPEN(110,FILE='LUMO.xsf',ACTION='WRITE')
      
      WRITE(100,'(A5)')'ATOMS'
      WRITE(110,'(A5)')'ATOMS'
      
      DO I=1,NATOMS
                WRITE(100,'(I4,3(F15.10))')ATOMS(I)%Z,ATOMS(I)%R(1)*0.52917720859,ATOMS(I)%R(2)*0.52917720859,ATOMS(I)%R(3)*0.52917720859
                WRITE(110,'(I4,3(F15.10))')ATOMS(I)%Z,ATOMS(I)%R(1)*0.52917720859,ATOMS(I)%R(2)*0.52917720859,ATOMS(I)%R(3)*0.52917720859
      ENDDO
      
      WRITE(100,'(A22)')'BEGIN_BLOCK_DATAGRID3D'
      WRITE(110,'(A22)')'BEGIN_BLOCK_DATAGRID3D'
      WRITE(100,'(A15)')'3D_ORBITAL_#001'
      WRITE(110,'(A15)')'3D_ORBITAL_#001'
      WRITE(100,'(A15)')'DATAGRID_3D_G98CUBE'
      WRITE(110,'(A15)')'DATAGRID_3D_G98CUBE'

      WRITE(100,'(3(I6))')MESH
      WRITE(110,'(3(I6))')MESH

      WRITE(100,'(3(F15.10))') -LIMITS(1)*0.52917720859,-LIMITS(2)*0.52917720859,-LIMITS(3)*0.52917720859
      WRITE(110,'(3(F15.10))') -LIMITS(1)*0.52917720859,-LIMITS(2)*0.52917720859,-LIMITS(3)*0.52917720859

      !WRITE(100,'(3(F15.10))') 0.0d0,0.0d0,0.0d0
      !WRITE(110,'(3(F15.10))') 0.0d0,0.0d0,0.0d0
      
      WRITE(100,'(3(F15.10))') 2*LIMITS(1)*0.52917720859,0.0d0,0.0d0
      WRITE(100,'(3(F15.10))') 0.0d0,2*LIMITS(2)*0.52917720859,0.0d0
      WRITE(100,'(3(F15.10))') 0.0d0, 0.0d0,2*LIMITS(3)*0.52917720859
      WRITE(110,'(3(F15.10))') 2*LIMITS(1)*0.52917720859,0.0d0,0.0d0
      WRITE(110,'(3(F15.10))') 0.0d0,2*LIMITS(2)*0.52917720859,0.0d0
      WRITE(110,'(3(F15.10))') 0.0d0, 0.0d0,2*LIMITS(3)*0.52917720859
      
      DO K=1,MESH(3)
                R(3) = -LIMITS(3) + 2.0d0*LIMITS(3)*(K-1)/(1.0d0*(MESH(3)-1))
                DO J=1,MESH(2)
                        R(2) = -LIMITS(2) + 2.0d0*LIMITS(2)*(J-1)/(1.0d0*(MESH(2)-1))
                        DO I=1,MESH(1)
                                R(1) = -LIMITS(1) + 2.0d0*LIMITS(1)*(I-1)/(1.0d0*(MESH(1)-1))
                                denshomo = rho(BAS,Phomo,R)
                                denslumo = rho(BAS,Plumo,R)
                                homo = hforbitalval(BAS,Ch,Ih,R)
                                lumo = hforbitalval(BAS,Cl,Il,R)
                                IF ( MAXHOMO .LT. denshomo ) MAXHOMO = denshomo
                                IF ( MAXLUMO .LT. denslumo ) MAXLUMO = denslumo
                                WRITE(10,FMT=2)I,J,K,denshomo
                                WRITE(11,FMT=2)I,J,K,denslumo
                                WRITE(100,'(F18.15)')homo
                                WRITE(110,'(F18.15)')lumo
                        ENDDO
                ENDDO
       ENDDO
       WRITE(100,'(A15)')'END_DATAGRID_3D'
       WRITE(100,'(A20)')'END_BLOCK_DATAGRID3D'
       WRITE(110,'(A15)')'END_DATAGRID_3D'
       WRITE(110,'(A20)')'END_BLOCK_DATAGRID3D'
       close(10)
       close(11)
       close(100)
       close(110)
       print*,' '
       print*,'============================================================='
       print*,'HOMO and LUMO densies has been calculated and saved to files.'
       print*,'============================================================='
       print*,' '
       WRITE(*,'(A37,F25.20)')'  Max HOMO-density value, rho_max =',MAXHOMO
       WRITE(*,'(A37,F25.20)')'  Max LUMO-density value, rho_max =',MAXLUMO
       print*,' '
       print*,'   ========================================================'
       print*,' '
END SUBROUTINE homolumosave

      
