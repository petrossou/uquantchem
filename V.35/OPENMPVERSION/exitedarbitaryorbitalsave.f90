SUBROUTINE exitedarbitaryorbitalsave(LIMITS,MESH,BAS,Cup,Cdown,EHFeigenup,EHFeigendown,Ne,NATOMS,ATOMS,AORBS,NTIMESTEPS)
      ! This subroutine saves the square modulus of an orbital weighted with the
      ! occupation-number obtained from a TDDFT/THF calculation at 24 different times of the TDDFT/THF calculation
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: LIMITS(3)
      INTEGER, INTENT(IN) :: MESH(3),Ne,NATOMS,AORBS,NTIMESTEPS
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(BAS%NBAS),EHFeigendown(BAS%NBAS)
      DOUBLE PRECISION :: rho
      DOUBLE PRECISION :: Phomo(BAS%NBAS,BAS%NBAS),Plumo(BAS%NBAS,BAS%NBAS),Ch(BAS%NBAS,BAS%NBAS),Cl(BAS%NBAS,BAS%NBAS)
      INTEGER :: I,J,K,Nup,Ndown,Ih,Il,STEP
      DOUBLE PRECISION :: R(3),NORM,MAXHOMO,MAXLUMO,denshomo,denslumo,homo,lumo,LAS(BAS%NBAS),TIME
      DOUBLE PRECISION, EXTERNAL :: hforbitalval
       print*,' '
       print*,'============================================================='
       print*,'          ORBITAL WITH INDEX = AORBS IS SAVED TO FILE        '
       print*,'============================================================='
       print*,' '
      
       IF ( AORBS .GT. 0 ) OPEN(222,FILE='OCCUPATIONUP.dat',ACTION='READ')
       IF ( AORBS .LT. 0 ) OPEN(222,FILE='OCCUPATIONDOWN.dat',ACTION='READ')
       MAXHOMO = -1000000.0d0
       MAXLUMO = -1000000.0d0
       Phomo = 0.0d0
       Plumo = 0.0d0
       Ch = 0.0d0
       Cl = 0.0d0
      
      2 FORMAT(I6,I6,I6,F30.20)
      
      DO STEP=1,NTIMESTEPS
        READ(222,*)TIME,LAS
        IF ( MOD(STEP,INT(NTIMESTEPS/24)) .EQ. 0 ) THEN
                IF ( INT(NTIMESTEPS/24)  .EQ. 1 ) THEN
                        print*,' '
                        print*,'=================================================================================='
                        print*,'  The orbital**2 with orbital-index = AORBS is being saved at different time-steps'
                        print*,'=================================================================================='
                        print*,' '
                        WRITE(*,'(A3,I4,A14)')'t =',INT(24.0d0*STEP/(1.0d0*NTIMESTEPS)) ,'x Timesteps/24'
                ELSE
                        WRITE(*,'(A3,I4,A14)')'t =',INT(24.0d0*STEP/(1.0d0*NTIMESTEPS)) ,'x Timesteps/24'
                ENDIF 

                IF ( AORBS .LT. 0  ) THEN
                        Ch(:,-AORBS)  = sqrt(DABS(LAS(-AORBS)))*Cdown(:,-AORBS)
                ELSE
                        Ch(:,AORBS)  = sqrt(DABS(LAS(AORBS)))*Cup(:,AORBS)
                ENDIF
                
                IF ( STEP .EQ. 1*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_1.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 2*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_2.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 3*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_3.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 4*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_4.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 5*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_5.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 6*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_6.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 7*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_7.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 8*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_8.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 9*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_9.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 10*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_10.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 11*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_11.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 12*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_12.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 13*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_13.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 14*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_14.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 15*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_15.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 16*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_16.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 17*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_17.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 18*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_18.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 19*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_19.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 20*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_20.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 21*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_21.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 22*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_22.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 23*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_23.xsf',ACTION='WRITE')
                IF ( STEP .EQ. 24*INT(NTIMESTEPS/24) ) OPEN(100,FILE='ARBORB_24.xsf',ACTION='WRITE')
      
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
                                        WRITE(100,'(F18.15)')homo**2
                                ENDDO
                        ENDDO
                ENDDO
                WRITE(100,'(A15)')'END_DATAGRID_3D'
                WRITE(100,'(A20)')'END_BLOCK_DATAGRID3D'
                close(100)
        ENDIF
      ENDDO
END SUBROUTINE exitedarbitaryorbitalsave

      
