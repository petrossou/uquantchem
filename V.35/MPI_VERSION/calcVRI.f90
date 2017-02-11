SUBROUTINE calcVRI(NATOMS,FNATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,CFORCE,id,numprocessors)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BASAUX
      INTEGER, INTENT(IN) :: NATOMS,FNATOMS,id,numprocessors
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25)
      LOGICAL, INTENT(IN) :: CFORCE
      DOUBLE PRECISION, INTENT(OUT) :: VRI(BASAUX%NBAS,BASAUX%NBAS),gradVRI(FNATOMS,3,BASAUX%NBAS,BASAUX%NBAS)
      DOUBLE PRECISION, ALLOCATABLE :: TVRI(:,:),TgradVRI(:,:,:,:)
      INTEGER*8,ALLOCATABLE :: IMAP(:,:)
      INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,GG,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART
      INTEGER*8 :: NRED,N0p(numprocessors),Istart,Iend,N0,MODEIGHT
      INTEGER :: ierr
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,DMAX
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,DMAXIMUM(6),NOLL(3)
      DOUBLE PRECISION, ALLOCATABLE :: gradient(:,:)
      INTEGER, EXTERNAL :: ijkl
      EXTERNAL :: gradprimeeintegral
      LOGICAL :: SAMECITE

      NOLL = 0.0d0
    
      NRED = ((BASAUX%NBAS+1)*BASAUX%NBAS)/2
      ALLOCATE(IMAP(NRED,2))

      DO J=1,BASAUX%NBAS
        DO I=1,J
                M = (J*(J-1))/2 + I
                        IMAP(M,1) = I
                        IMAP(M,2) = J
        ENDDO
     ENDDO

        !============================================================================
        ! Here we distribute the calculation of the columns of the ee-integral array
        ! Intsv over the (numprocessors) number of mpi threads.
        !----------------------------------------------------------------------------

        N0 = INT(NRED/numprocessors,8)

        Istart = 0
        Iend = 0
        N0p(:) = 0
        DO I=1,numprocessors
          N0p(I) = N0
        ENDDO

        !--------------------------------------------------------------------------------------------------------------
        ! If the number of elements (NRED) is not divisible by the number of processors, the remaining 
        ! number of elements of Intsv(I) to be calculated are distributed over the (numprocessors) number of mpi threads
        !--------------------------------------------------------------------------------------------------------------
     
        I = 1
        MODEIGHT = NRED - INT(NRED/numprocessors,8)*numprocessors
        !DO WHILE( I .LE. MOD(NRED,numprocessors) )
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
    
     VRI(:,:) = 0.0d0
     gradVRI(:,:,:,:) = 0.0d0
     
      ALLOCATE(gradient(4,3))
         
        DO G=Istart,Iend
                 I = IMAP(G,1) 
                 J = IMAP(G,2) 
                
                 L1= BASAUX%PSI(I)%L(1)
                 M1= BASAUX%PSI(I)%L(2)
                 N1= BASAUX%PSI(I)%L(3)
                 A = BASAUX%PSI(I)%R
                 NO1 = BASAUX%PSI(I)%NORM 
                
                 L2= 0
                 M2= 0
                 N2= 0
                 !B = NOLL
                 B = BASAUX%PSI(I)%R
                 NO2 = 1.0d0

                 L3= BASAUX%PSI(J)%L(1)
                 M3= BASAUX%PSI(J)%L(2)
                 N3= BASAUX%PSI(J)%L(3)
                 C = BASAUX%PSI(J)%R
                 NO3 = BASAUX%PSI(J)%NORM 
                 
                 L4= 0
                 M4= 0
                 N4= 0
                 !D = NOLL
                 D = BASAUX%PSI(J)%R
                 NO4 = 1.0d0
              
                 
                 ! SUM OVER PRIMITIVE GAUSSIANS:
                 DO M=1,BASAUX%PSI(I)%NPRIM
                          al1 = BASAUX%PSI(I)%EXPON(M)
                          co1 = BASAUX%PSI(I)%CONTRCOEFF(M)
                          NP1 = BASAUX%PSI(I)%PRIMNORM(M)
                          DO N=1,BASAUX%PSI(J)%NPRIM
                                  al3 = BASAUX%PSI(J)%EXPON(N)
                                  co3 = BASAUX%PSI(J)%CONTRCOEFF(N)
                                  NP3 = BASAUX%PSI(J)%PRIMNORM(N)
                         
                                  al2 = 0.0d0
                                  co2 = 1.0d0
                                  NP2 = 1.0d0
                                                       
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

                                   VRI(I,J) = VRI(I,J) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)

                                  !--------------------------------------------------------------
                                  ! Here we caclulate the gradient of  the two electron integrals
                                  !--------------------------------------------------------------

                                  IF ( CFORCE ) THEN
                                           CALL gradprimeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW,gradient)
                                            gradVRI(BASAUX%PSI(I)%ATYPE,:,I,J) = gradVRI(BASAUX%PSI(I)%ATYPE,:,I,J) + TERM*gradient(1,:)
                                            gradVRI(BASAUX%PSI(J)%ATYPE,:,I,J) = gradVRI(BASAUX%PSI(J)%ATYPE,:,I,J) + TERM*gradient(3,:)
                                  ENDIF
                         ENDDO
                ENDDO
                IF ( I .NE. J ) THEN
                        VRI(J,I) = VRI(I,J)
                        gradVRI(:,:,J,I) = gradVRI(:,:,I,J)
                ENDIF
                                
        ENDDO  
        DEALLOCATE(gradient)

       IF ( id .EQ. 0 ) THEN
                ALLOCATE(TVRI(BASAUX%NBAS,BASAUX%NBAS),TgradVRI(FNATOMS,3,BASAUX%NBAS,BASAUX%NBAS))
                TVRI = 0.0d0
                TgradVRI = 0.0d0
        ENDIF

       CALL MPI_REDUCE(VRI,TVRI,BASAUX%NBAS*BASAUX%NBAS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
       CALL MPI_REDUCE(gradVRI,TgradVRI,BASAUX%NBAS*BASAUX%NBAS*3*FNATOMS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      
       IF ( id .EQ. 0 ) THEN
                VRI = TVRI 
                gradVRI = TgradVRI
                DEALLOCATE(TVRI,TgradVRI)
       ENDIF

       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

       CALL MPI_BCAST(gradVRI,BASAUX%NBAS*BASAUX%NBAS*3*FNATOMS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       CALL MPI_BCAST(VRI,BASAUX%NBAS*BASAUX%NBAS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE calcVRI
