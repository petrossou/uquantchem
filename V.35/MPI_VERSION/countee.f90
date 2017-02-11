SUBROUTINE countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,id,DMAT,CREATEMAPPING,NDIAG,NONZERO,TOTALNONZERO,dInts)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: id,numprocessors
      INTEGER*8, INTENT(IN) :: Istart,Iend
      INTEGER, INTENT(INOUT) :: IND1(NONZERO),IND2(NONZERO),IND3(NONZERO),IND4(NONZERO)
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS)
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
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL,DMAX
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,gradient(4,3),DMAXIMUM(6)
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
     
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

                        L3= BAS%PSI(K)%L(1)
                        M3= BAS%PSI(K)%L(2)
                        N3= BAS%PSI(K)%L(3)
                        C = BAS%PSI(K)%R
                        NO3 = BAS%PSI(K)%NORM 
                
                        L4= BAS%PSI(L)%L(1)
                        M4= BAS%PSI(L)%L(2)
                        N4= BAS%PSI(L)%L(3)
                        D = BAS%PSI(L)%R
                        NO4 = BAS%PSI(L)%NORM
                     
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

                        IF ( DOSUMOVERPRIMS .AND. GG .EQ. 1 ) THEN
                                INTEGRAAL = 0.0d0
                                ! SUM OVER PRIMITIVE GAUSSIANS:
                                DO M=1,BAS%PSI(I)%NPRIM
                                        al1 = BAS%PSI(I)%EXPON(M)
                                        co1 = BAS%PSI(I)%CONTRCOEFF(M)
                                        NP1 = BAS%PSI(I)%PRIMNORM(M)
                                        DO N=1,BAS%PSI(J)%NPRIM
                                                al2 = BAS%PSI(J)%EXPON(N)
                                                co2 = BAS%PSI(J)%CONTRCOEFF(N)
                                                NP2 = BAS%PSI(J)%PRIMNORM(N)
                                                DO P=1,BAS%PSI(K)%NPRIM
                                                        al3 = BAS%PSI(K)%EXPON(P)
                                                        co3 = BAS%PSI(K)%CONTRCOEFF(P)
                                                        NP3 = BAS%PSI(K)%PRIMNORM(P)
                                                        DO Q=1,BAS%PSI(L)%NPRIM
                                                                al4 = BAS%PSI(L)%EXPON(Q)
                                                                co4 = BAS%PSI(L)%CONTRCOEFF(Q)
                                                                NP4 = BAS%PSI(L)%PRIMNORM(Q)
                                                                TERM = NO1*NO2*NO3*NO4*NP1*NP2*NP3*NP4*co1*co2*co3*co4
                                                                ! Here we check if all the orbitals i,j,k,l are centered on the
                                                                ! same atom. If so we can use the pre calculated weights and
                                                                ! roots of the Rys polynomals in order to calculate (ij|kl)
                                                                IF ( DOT_PRODUCT(A-B,A-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-B,C-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-D,C-D) .EQ. 0.0d0 ) THEN
                                                                        SAMECITE = .TRUE.
                                                                ELSE
                                                                        SAMECITE = .FALSE.
                                                                ENDIF
                                                                
                                                                !-----------------------------------------------
                                                                ! Here we calculate the two electron integrals
                                                                !-----------------------------------------------

                                                                IF ( GG .EQ. 1 ) THEN
                                                                    INTEGRAAL = INTEGRAAL + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)
                                                                ENDIF

                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                ENDDO
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
END SUBROUTINE countee
