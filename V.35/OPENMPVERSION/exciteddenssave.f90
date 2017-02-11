SUBROUTINE exciteddenssave(LIMITS,MESH,BAS,Pup,Pdown,Pupp,Pdownn,NATOMS,ATOMS,TIMEINDEX,DIFFDENS)
      ! This subroutine saves the hartree-fock charge-density on disk.
      ! Input that is needed is the MESH, i.e how many grid points are
      ! going to be used in the x, y and z directions to discretizise
      ! the density. The LIMITS limits in space for the calculated
      ! charge density: x in [-LIMITS(1),LIMITS(1)]  ,  
      ! y in [-LIMITS(1),LIMITS(1)] , z in [-LIMITS(1),LIMITS(1)].
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      LOGICAL, INTENT(IN) :: DIFFDENS
      INTEGER, INTENT(IN) :: MESH(3),NATOMS,TIMEINDEX
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS),Pdown(BAS%NBAS,BAS%NBAS),Pupp(BAS%NBAS,BAS%NBAS),Pdownn(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION :: rho
      INTEGER :: I,J,K,TID
      DOUBLE PRECISION :: R(3),NORM,MAXUP,MAXDOWN,densup,densdown,dens
      DOUBLE PRECISION, EXTERNAL :: hforbitalval

       TID = TIMEINDEX - 1
       IF ( TIMEINDEX .EQ. 1 ) THEN
                print*,' '
                print*,'==================================================================='
                print*,'   Exited State Charge-densities are being calculated and saved.   '
                print*,'==================================================================='
                print*,' '
                WRITE(*,'(A3,I4,A14)')'t =',TID,'x Timesteps/24'
        ELSE
                WRITE(*,'(A3,I4,A14)')'t =',TID,'x Timesteps/24'
        ENDIF
       
       MAXUP = -1000000.0d0
       MAXDOWN = -1000000.0d0
      
      IF ( TIMEINDEX .EQ. 1 ) OPEN(100,FILE='DENSEXITED_0.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 2 ) OPEN(100,FILE='DENSEXITED_1.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 3 ) OPEN(100,FILE='DENSEXITED_2.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 4 ) OPEN(100,FILE='DENSEXITED_3.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 5 ) OPEN(100,FILE='DENSEXITED_4.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 6 ) OPEN(100,FILE='DENSEXITED_5.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 7 ) OPEN(100,FILE='DENSEXITED_6.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 8 ) OPEN(100,FILE='DENSEXITED_7.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 9 ) OPEN(100,FILE='DENSEXITED_8.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 10 ) OPEN(100,FILE='DENSEXITED_9.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 11 ) OPEN(100,FILE='DENSEXITED_10.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 12 ) OPEN(100,FILE='DENSEXITED_11.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 13 ) OPEN(100,FILE='DENSEXITED_12.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 14 ) OPEN(100,FILE='DENSEXITED_13.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 15 ) OPEN(100,FILE='DENSEXITED_14.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 16 ) OPEN(100,FILE='DENSEXITED_15.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 17 ) OPEN(100,FILE='DENSEXITED_16.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 18 ) OPEN(100,FILE='DENSEXITED_17.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 19 ) OPEN(100,FILE='DENSEXITED_18.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 20 ) OPEN(100,FILE='DENSEXITED_19.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 21 ) OPEN(100,FILE='DENSEXITED_20.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 22 ) OPEN(100,FILE='DENSEXITED_21.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 23 ) OPEN(100,FILE='DENSEXITED_22.xsf',ACTION='WRITE')
      
      IF ( TIMEINDEX .EQ. 24 ) OPEN(100,FILE='DENSEXITED_23.xsf',ACTION='WRITE')
      
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
                                IF ( DIFFDENS ) THEN
                                        densup = rho(BAS,Pup-Pupp,R)
                                        densdown = rho(BAS,Pdown-Pdownn,R)
                                        dens = densup + densdown
                                ELSE
                                        densup = rho(BAS,Pup,R)
                                        densdown = rho(BAS,Pdown,R)
                                        dens = densup + densdown
                                ENDIF

                                WRITE(100,'(F18.15)')dens
                        ENDDO
                ENDDO
       ENDDO
       WRITE(100,'(A15)')'END_DATAGRID_3D'
       WRITE(100,'(A20)')'END_BLOCK_DATAGRID3D'
       close(100)
END SUBROUTINE exciteddenssave

      
