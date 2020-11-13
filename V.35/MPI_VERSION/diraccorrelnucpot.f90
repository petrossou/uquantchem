SUBROUTINE diraccorrelnucpot(BAS,NATOMS,ATOMS,VDIRAC,id,numprocessors)
      ! This subroutine calculates the potential energy matrix
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: basfunkval
      INTEGER, INTENT(IN) :: NATOMS,id,numprocessors
      TYPE(BASIS), INTENT(IN) :: BAS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: VDIRAC(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION :: Vl(BAS%NBAS,BAS%NBAS)
      INTEGER :: NB,I,J,N,K,M,N0,NPP,NTOT,Istart,Istop,NCONTRACT,ierr
      INTEGER, ALLOCATABLE :: IND1(:), IND2(:)
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, PARAMETER :: lightspeed = 137.035999084
      
      !print*,'V_ij='
      NB = BAS%NBAS

      Vl(:,:) = 0.0d0 

      !==========================================================================================================
      ! Preparation for calculating the relatevistic Dirac correction comming from the nuclear-electron potential
      !==========================================================================================================

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

                        Vl(I,J) = 0.0d0
                        
                        DO M=1,NATOMS
                                !--------------------------------------------------------------------------
                                ! Here we calculate the Dirac electron nuclear potential correction matrix
                                !--------------------------------------------------------------------------
                                                
                                Vl(I,J) = Vl(I,J) - (pi/(2.0*(lightspeed**2)))*ATOMS(M)%Z*basfunkval(BAS%PSI(I),ATOMS(M)%R)*basfunkval(BAS%PSI(J),ATOMS(M)%R)

                        ENDDO

                        IF ( I .NE. J ) THEN 
                                Vl(J,I) = Vl(I,J)
                        ENDIF

                ENDDO
      ENDIF
      CALL MPI_REDUCE(Vl,VDIRAC,NB*NB,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)

      CALL MPI_BCAST(VDIRAC,NB*NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE diraccorrelnucpot
