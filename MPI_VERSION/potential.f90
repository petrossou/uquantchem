SUBROUTINE potential(BAS,NATOMS,ATOMS,V,gradV,NSI,CFORCE,id,numprocessors)
      ! This subroutine calculates the potential energy matrix
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primpotential
      INTEGER, INTENT(IN) :: NATOMS,NSI,id,numprocessors
      LOGICAL, INTENT(IN) :: CFORCE
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: V(BAS%NBAS,BAS%NBAS),gradV(NATOMS,3,NSI,NSI)
      DOUBLE PRECISION :: Vl(BAS%NBAS,BAS%NBAS),gradVl(NATOMS,3,NSI,NSI)
      INTEGER :: NB,I,J,N,K,M,L1,M1,N1,L2,M2,N2,N0,NPP,NTOT,Istart,Istop,NCONTRACT,ierr
      INTEGER, ALLOCATABLE :: IND1(:), IND2(:)
      DOUBLE PRECISION :: NO1,NO2,NP1,NP2,gradient(3,3)
      DOUBLE PRECISION :: A(3),B(3),alpha1,alpha2,coeff1,coeff2
      
      !print*,'V_ij='
      NB = BAS%NBAS

      gradVl(:,:,:,:) = 0.0d0
      Vl(:,:) = 0.0d0 

      !=============================================================================================
      ! Preparation for calculating potential energy matrix elements and their gradients in paralell
      !=============================================================================================

      NTOT = ( (BAS%NBAS+1)*BAS%NBAS )/2
      ALLOCATE(IND1(NTOT),IND2(NTOT))

      IF ( numprocessors > NTOT ) THEN
         Istart = id+1
         Istop  = id+1
      ELSE
         NPP = ( NTOT - MOD(NTOT,numprocessors) )/numprocessors
         Istart = id*NPP+1
         IF ( id < numprocessors -1 ) THEN
            Istop  = (id+1)*NPP
         ELSE
            Istop  = (id+1)*NPP + MOD(NTOT,numprocessors)
         ENDIF
      ENDIF
      
      !==========================================================
      ! Creating the mapping between the contracted index and the
      ! regular indexes: (I,J)
      !========================================================== 

      DO J=1,BAS%NBAS
         DO I=1,J
            NCONTRACT= ( J*(J-1) )/2 + I
            IND1(NCONTRACT) = I
            IND2(NCONTRACT) = J
         ENDDO
      ENDDO
      !print*,id,NPP,Istart,Istop
      !STOP
            
IF ( ( numprocessors .LT. NTOT ) .OR. ( numprocessors .GE. NTOT .AND. id + 1 .LE. NTOT ) ) THEN
     
    DO NCONTRACT=Istart,Istop
                I = IND1(NCONTRACT)
                J = IND2(NCONTRACT)

                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
              
                        Vl(I,J) = 0.0d0
                        
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R
                        NO2 = BAS%PSI(J)%NORM
                        DO N=1,BAS%PSI(I)%NPRIM
                                
                                alpha1 = BAS%PSI(I)%EXPON(N)
                                coeff1 = BAS%PSI(I)%CONTRCOEFF(N)
                                NP1 = BAS%PSI(I)%PRIMNORM(N)
                                DO K=1,BAS%PSI(J)%NPRIM
                                        alpha2 = BAS%PSI(J)%EXPON(K)
                                        coeff2 = BAS%PSI(J)%CONTRCOEFF(K)
                                        NP2 = BAS%PSI(J)%PRIMNORM(K)
                                        DO M=1,NATOMS
                                                !----------------------------------------------
                                                ! Here we calculate the potential energy matrix
                                                !----------------------------------------------
                                                
                                                Vl(I,J) = Vl(I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*primpotential(L1,M1,N1,A,alpha1,ATOMS(M)%R,L2,M2,N2,B,alpha2)
                                                
                                                !---------------------------------------------------------------
                                                ! Here we calculate the gradient  of the potential energy matrix
                                                !---------------------------------------------------------------
                                                IF ( CFORCE ) THEN
                                                  CALL gradprimpotential(L1,M1,N1,A,alpha1,ATOMS(M)%R,L2,M2,N2,B,alpha2,gradient)

                                                  gradVl(BAS%PSI(I)%ATYPE,:,I,J) = gradVl(BAS%PSI(I)%ATYPE,:,I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(1,:)
                                                  gradVl(BAS%PSI(J)%ATYPE,:,I,J) = gradVl(BAS%PSI(J)%ATYPE,:,I,J) - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(2,:)
                                                  gradVl(M,:,I,J)                = gradVl(M,:,I,J)                - ATOMS(M)%Z*NO1*NO2*coeff1*coeff2*NP1*NP2*gradient(3,:)
                                                ENDIF
                                        ENDDO
                                ENDDO
                        ENDDO

                        IF ( I .NE. J ) THEN 
                                Vl(J,I) = Vl(I,J)
                                IF ( CFORCE ) THEN
                                  DO M=1,NATOMS
                                        gradVl(M,:,J,I) = gradVl(M,:,I,J) 
                                  ENDDO
                                ENDIF
                        ENDIF
                        !print*,I,J,V(I,J)

      ENDDO
ENDIF
       CALL MPI_REDUCE(Vl,V,NB*NB,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
       CALL MPI_REDUCE(gradVl,gradV,3*NATOMS*NSI*NSI,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)

       CALL MPI_BCAST(V,NB*NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       CALL MPI_BCAST(gradV,3*NATOMS*NSI*NSI,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE potential
