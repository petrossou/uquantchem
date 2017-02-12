SUBROUTINE counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,id,DMAT,CREATEMAPPING,NDIAG,NONZERO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: id,numprocessors,NBAUX
      INTEGER*8, INTENT(IN) :: Istart,Iend
      INTEGER, INTENT(INOUT) :: IND1(NONZERO),IND2(NONZERO),IND3(NONZERO),IND4(NONZERO)
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, INTENT(IN) :: VRI(NBAUX,NBAUX),WRI(BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION, INTENT(INOUT) :: dInts(NDIAG)
      LOGICAL, INTENT(IN) :: APPROXEE,CREATEMAPPING
      INTEGER*8, INTENT(INOUT) :: NONZERO,NDIAG
      INTEGER*8, INTENT (OUT) :: TOTALNONZERO
      DOUBLE PRECISION, ALLOCATABLE :: tempintsl(:),tempintsi(:)
      INTEGER*8, ALLOCATABLE :: indexus(:),indexusl(:)
      INTEGER*8 :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,GG,STRL,KK,PP,PPP
      INTEGER*8 :: rcounts(numprocessors),displs(numprocessors),ierr,scount,SI,SIC(numprocessors)
      INTEGER*8 :: NONZEROSAVE
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL,DMAX,INVRI(NBAUX,NBAUX),VEC(NBAUX)
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,gradient(4,3),DMAXIMUM(6)
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 

      call invert(VRI,INVRI,NBAUX)

      NONZEROSAVE = NONZERO
      STRL = ((BAS%NBAS+1)*BAS%NBAS)/2
      ALLOCATE(tempintsl(STRL),tempintsi(STRL),indexus(STRL),indexusl(STRL))
      tempintsl = 0.0d0
      tempintsi = 0.0d0
      IF ( .not. CREATEMAPPING ) dInts = 0.0d0
      indexus = 0
      SI = 0
      NONZERO = 0

      P = 0
      PP = Istart-1

      DO J=1,BAS%NBAS
        DO I=1,J
            DO L=1,BAS%NBAS
              DO K=1,L
                M = (J*(J-1))/2 + I
                N = (L*(L-1))/2 + K
                P = P + 1
                IF ( P .EQ. id+1 ) PP = PP + 1
                IF (  PP .LE. Iend .AND. P .EQ. id+1 ) THEN
                             IMAP(PP,1) = I
                             IMAP(PP,2) = J
                             IMAP(PP,3) = K
                             IMAP(PP,4) = L
                             IMAP(PP,5) = M
                             IMAP(PP,6) = N
                             PPP = PP
                ENDIF
                IF ( P .EQ. numprocessors ) P = 0
              ENDDO
           ENDDO
        ENDDO
     ENDDO

DO GG=1,2
        DO NN=Istart,Iend
                        G = NN
                        I = IMAP(G,1) 
                        J = IMAP(G,2) 
                        K = IMAP(G,3) 
                        L = IMAP(G,4) 
                        IJ = IMAP(G,5)
                        KL = IMAP(G,6) 
                
                        !------------------------
                        ! Swartz inequality part
                        !------------------------
                        ! Here we use Schwartz inequality |(i,j,k,l)| <= ! [ |(i,j,i,j)|*(k,l,k,l) ]**(1/2). 
                        ! See Eqn 9.12.25 p.404 in T. Helgaker et al's Molecular Electronic-Structure Theory
                        ! We use this together with Eqn. (56) and (57) page. 86,
                        ! in Haettig, Multiscale Simulation Methods in Molecular
                        ! Sciences, J. Grotendorst, N. Attig, S. Blugel, D. Marx (Eds.),
                        ! Institute for Advanced Simulation, Forschungszentrum Julich, NIC Series, Vol. 42, ISBN 978-3-9810843-8-2,pp.  77-120, 2009.
                        
                        DMAXIMUM = (/ 2*DABS(DMAT(I,J)),2*DABS(DMAT(K,L)),DABS(DMAT(J,K)),DABS(DMAT(I,K)),DABS(DMAT(J,L)),DABS(DMAT(I,L)) /)

                        DMAX = MAXVAL(DMAXIMUM)

                        DOSUMOVERPRIMS = .FALSE.

                        IF (  GG .EQ. 1 .AND. IJ .EQ. KL ) THEN
                                  DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        
                        IF ( ( GG .EQ. 2 .AND. sqrt(DABS(dInts(IJ)*dInts(KL)))*DMAX .GE. Tol ) .OR. ( GG .EQ. 2  .AND. .not. APPROXEE)  ) THEN
                                ! counting the non-zero ee-integrals.
                                    NONZERO = NONZERO + 1
                                    IF (  CREATEMAPPING ) THEN
                                       IND1(NONZERO) = IMAP(G,1)
                                       IND2(NONZERO) = IMAP(G,2)
                                       IND3(NONZERO) = IMAP(G,3)
                                       IND4(NONZERO) = IMAP(G,4)
                                     ENDIF
                        ENDIF
                       
                        IF ( CREATEMAPPING  ) DOSUMOVERPRIMS = .FALSE.
                        !-----------------------
                        ! End of inequality part
                        !-----------------------
                       INTEGRAAL = 0.0d0
                        IF ( DOSUMOVERPRIMS .AND. GG .EQ. 1 ) THEN
                                !-----------------------------------------------
                                ! Here we calculate the two electron integrals
                                !-----------------------------------------------

                                IF ( GG .EQ. 1 ) THEN
                                        CALL DGEMM ( 'N', 'N', NBAUX, 1, NBAUX, 1.0d0, INVRI, NBAUX, WRI(K,L,:), NBAUX, 0.0d0,VEC,NBAUX )
                                        INTEGRAAL = DOT_PRODUCT(WRI(I,J,:),VEC)
                                ENDIF

                        ENDIF
                        IF ( GG .EQ. 1 .AND. IJ .EQ. KL ) THEN
                              tempintsl(IJ) = INTEGRAAL
                        ENDIF
       ENDDO
       IF ( GG .EQ. 1 .AND. .not. CREATEMAPPING ) THEN
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          CALL MPI_REDUCE(tempintsl,dInts,NDIAG,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
          CALL MPI_BCAST(dInts,NDIAG,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
       ENDIF
       
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              
       
ENDDO
!CALL MPI_REDUCE(NONZERO,TOTALNONZERO,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD, ierr)
!CALL MPI_BCAST(TOTALNONZERO,1,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
CALL MPI_REDUCE(NONZERO,TOTALNONZERO,1,MPI_LONG_LONG_INT,MPI_SUM,0,MPI_COMM_WORLD, ierr)
CALL MPI_BCAST(TOTALNONZERO,1,MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD, ierr)
!--------------------------------------------------------------
! Restoring the NONZERO variable, since after 
! the run NONZERO migth have been set to zero, when it infact 
! should have the dummy value 1.
!--------------------------------------------------------------
IF ( CREATEMAPPING ) NONZERO = NONZEROSAVE
END SUBROUTINE counteeRI
