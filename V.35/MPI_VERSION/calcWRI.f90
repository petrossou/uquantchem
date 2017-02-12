SUBROUTINE calcWRI(NATOMS,FNATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,CFORCE,id,numprocessors)
! This subroutine calculates the 
! electron-electron repulsion tensor.
USE datatypemodule
USE OMP_LIB
IMPLICIT NONE
INCLUDE "mpif.h"
DOUBLE PRECISION, EXTERNAL :: primeeintegral
TYPE(BASIS), INTENT(IN) :: BAS,BASAUX
INTEGER, INTENT(IN) :: NATOMS,FNATOMS,id,numprocessors
DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25)
LOGICAL, INTENT(IN) :: CFORCE
DOUBLE PRECISION, INTENT(OUT) :: WRI(BAS%NBAS,BAS%NBAS,BASAUX%NBAS),gradWRI(FNATOMS,3,BAS%NBAS,BAS%NBAS,BASAUX%NBAS)
DOUBLE PRECISION, ALLOCATABLE :: TWRI(:,:,:),TgradWRI(:,:,:,:,:)
INTEGER*8,ALLOCATABLE :: IMAP(:,:)
INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,GG,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,NRED
INTEGER*8 :: N0p(numprocessors),Istart,Iend,N0,MODEIGHT
DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,DMAX
DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,DMAXIMUM(6),NOLL(3)
DOUBLE PRECISION, ALLOCATABLE :: gradient(:,:)
INTEGER :: ierr
INTEGER, EXTERNAL :: ijkl
EXTERNAL :: gradprimeeintegral
LOGICAL :: SAMECITE

NOLL = 0.0d0

NRED = ((BAS%NBAS+1)*BAS%NBAS)/2
ALLOCATE(IMAP(NRED,2))

DO J=1,BAS%NBAS
        DO I=1,J
                M = (J*(J-1))/2 + I
                        IMAP(M,1) = I
                        IMAP(M,2) = J
        ENDDO
ENDDO

!============================================================================
! Here we distribute the calculation of the WRI-array.
!----------------------------------------------------------------------------

N0 = INT(NRED/numprocessors,8)

Istart = 0
Iend = 0
N0p(:) = 0
DO I=1,numprocessors
   N0p(I) = N0
ENDDO

I = 1
MODEIGHT = NRED - INT(NRED/numprocessors,8)*numprocessors
DO WHILE( I .LE. MODEIGHT )
       N0p(I) = N0p(I) + 1
       I = I+1
ENDDO


IF ( id .EQ. 0 ) THEN
       Istart = 1
ELSE
    Istart = 1
    DO I=1,id
       Istart = Istart+N0p(I)
    ENDDO
ENDIF

Iend = 0
DO I=1,id+1
            Iend = Iend + N0p(I)
ENDDO

WRI(:,:,:) = 0.0d0
gradWRI(:,:,:,:,:) = 0.0d0


 ALLOCATE(gradient(4,3))
 DO G=Istart,Iend
        DO Q=1,BASAUX%NBAS
                I = IMAP(G,1) 
                J = IMAP(G,2) 
        
                L1= BAS%PSI(I)%L(1)
                M1= BAS%PSI(I)%L(2)
                N1= BAS%PSI(I)%L(3)
                A = BAS%PSI(I)%R
                NO1 = BAS%PSI(I)%NORM 
        
                L2= BAS%PSI(J)%L(1)
                M2= BAS%PSI(J)%L(2)
                N2= BAS%PSI(J)%L(3)
                B = BAS%PSI(J)%R
                NO2 = BAS%PSI(J)%NORM

                L3= BASAUX%PSI(Q)%L(1)
                M3= BASAUX%PSI(Q)%L(2)
                N3= BASAUX%PSI(Q)%L(3)
                C = BASAUX%PSI(Q)%R
                NO3 = BASAUX%PSI(Q)%NORM 
                
                L4= 0
                M4= 0
                N4= 0
                !D = NOLL
                D = BASAUX%PSI(Q)%R
                NO4 = 1.0d0
                
                 ! SUM OVER PRIMITIVE GAUSSIANS:
                 DO M=1,BAS%PSI(I)%NPRIM
                         al1 = BAS%PSI(I)%EXPON(M)
                         co1 = BAS%PSI(I)%CONTRCOEFF(M)
                         NP1 = BAS%PSI(I)%PRIMNORM(M)
                         DO N=1,BAS%PSI(J)%NPRIM
                                 al2 = BAS%PSI(J)%EXPON(N)
                                 co2 = BAS%PSI(J)%CONTRCOEFF(N)
                                 NP2 = BAS%PSI(J)%PRIMNORM(N)
                                 DO K=1,BASAUX%PSI(Q)%NPRIM
                                         al3 = BASAUX%PSI(Q)%EXPON(K)
                                         co3 = BASAUX%PSI(Q)%CONTRCOEFF(K)
                                         NP3 = BASAUX%PSI(Q)%PRIMNORM(K)
                                                    
                                         al4 = 0.0d0
                                         co4 = 1.0d0
                                         NP4 = 1.0d0

                                         TERM = NO1*NO2*NO3*NO4*NP1*NP2*NP3*NP4*co1*co2*co3*co4
                                         ! Here we check if all the orbitals i,j,Q  are centered on the
                                         ! same atom. If so we can use the pre calculated weights and
                                         ! roots of the Rys polynomals in order to calculate (ij|k)
                                         IF ( DOT_PRODUCT(A-B,A-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-B,C-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-D,C-D) .EQ. 0.0d0 ) THEN
                                                  SAMECITE = .TRUE.
                                         ELSE
                                                  SAMECITE = .FALSE.
                                         ENDIF

                                        !-----------------------------------------------
                                        ! Here we  calculate the  integrals (ij|Q)
                                        !-----------------------------------------------

                                        WRI(I,J,Q) = WRI(I,J,Q) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)

                                        !--------------------------------------------------------------
                                        ! Here we caclulate the gradient of  the two electron integrals
                                        !--------------------------------------------------------------

                                        IF ( CFORCE ) THEN
                                                  CALL gradprimeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW,gradient)
                                                  gradWRI(BAS%PSI(I)%ATYPE,:,I,J,Q)    = gradWRI(BAS%PSI(I)%ATYPE,:,I,J,Q)    + TERM*gradient(1,:)
                                                  gradWRI(BAS%PSI(J)%ATYPE,:,I,J,Q)    = gradWRI(BAS%PSI(J)%ATYPE,:,I,J,Q)    + TERM*gradient(2,:)
                                                  gradWRI(BASAUX%PSI(Q)%ATYPE,:,I,J,Q) = gradWRI(BASAUX%PSI(Q)%ATYPE,:,I,J,Q) + TERM*gradient(3,:)
                                        ENDIF
                                ENDDO
                        ENDDO
                ENDDO
        ENDDO
        IF ( I .NE. J ) THEN 
                WRI(J,I,:)  = WRI(I,J,:)
                gradWRI(:,:,J,I,:) = gradWRI(:,:,I,J,:)
        ENDIF
ENDDO  
DEALLOCATE(gradient)

IF ( id .EQ. 0 ) THEN
        ALLOCATE(TWRI(BAS%NBAS,BAS%NBAS,BASAUX%NBAS),TgradWRI(FNATOMS,3,BAS%NBAS,BAS%NBAS,BASAUX%NBAS))
       TWRI = 0.0d0
       TgradWRI = 0.0d0
ENDIF

CALL MPI_REDUCE(WRI,TWRI,BAS%NBAS*BAS%NBAS*BASAUX%NBAS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
CALL MPI_REDUCE(gradWRI,TgradWRI,BAS%NBAS*BAS%NBAS*BASAUX%NBAS*3*FNATOMS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)

IF ( id .EQ. 0 ) THEN
      WRI = TWRI
      gradWRI = TgradWRI
      DEALLOCATE(TWRI,TgradwRI)
ENDIF

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

CALL MPI_BCAST(gradWRI,BAS%NBAS*BAS%NBAS*BASAUX%NBAS*3*FNATOMS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(WRI,BAS%NBAS*BAS%NBAS*BASAUX%NBAS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE calcWRI
