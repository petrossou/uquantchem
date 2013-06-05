SUBROUTINE eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,IND1,IND2,IND3,IND4,NRED,Istart,Iend,Istartg,Iendg,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,CFORCE,id,DMAT)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: NRED,id,numprocessors,Istart,Iend,Istartg,Iendg,NATOMS,FNATOMS
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS)
      LOGICAL, INTENT(IN) :: APPROXEE,CFORCE
      INTEGER, INTENT(INOUT) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(INOUT) :: Intsv(Istart:Iend),gradIntsv(FNATOMS,3,Istartg:Iendg)
      DOUBLE PRECISION, ALLOCATABLE :: tempintsl(:),tempintsi(:),tempints(:)
      INTEGER, ALLOCATABLE :: indexus(:),indexusl(:)
      INTEGER :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,GG,STRL,KK
      INTEGER :: rcounts(numprocessors),displs(numprocessors),ierr,scount,SI,SIC(numprocessors)
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL,DMAX
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,gradient(4,3),DMAXIMUM(6)
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
      
      STRL = ((BAS%NBAS+1)*BAS%NBAS)/2
      ALLOCATE(tempintsl(STRL),tempintsi(STRL),tempints(STRL),indexus(STRL),indexusl(STRL))
      tempintsl = 0.0d0
      tempintsi = 0.0d0
      tempints = 0.0d0
      indexus = 0
      SI = 0

      P = 0
      DO J=1,BAS%NBAS
        DO I=1,J
            DO L=1,BAS%NBAS
              DO K=1,L
                M = (J*(J-1))/2 + I
                N = (L*(L-1))/2 + K
                P = P + 1
                IF ( P .LE. Iend .AND. P .GE. Istart ) THEN
                             IMAP(P,1) = I
                             IND1(P) = I
                             IMAP(P,2) = J
                             IND2(P) = J
                             IMAP(P,3) = K
                             IND3(P) = K
                             IMAP(P,4) = L
                             IND4(P) = L
                             IMAP(P,5) = M
                             IMAP(P,6) = N
                ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
      
     Intsv(:) = 0.0d0
     gradIntsv(:,:,:) = 0.0d0

!==================================================
!==================================================
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

                        IF ( ( GG .EQ. 1 .AND. IJ .EQ. KL ) .OR. ( GG .EQ.  2 .AND. IJ .NE. KL .AND. sqrt(DABS(tempints(IJ)*tempints(KL)))*DMAX .NE. 0.0d0 ) ) THEN
                                  DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF (  APPROXEE .AND. GG .EQ.  2 .AND. IJ .NE. KL .AND. sqrt(DABS(tempints(IJ)*tempints(KL)))*DMAX .LT. Tol  ) THEN
                                DOSUMOVERPRIMS = .FALSE.
                        ENDIF
                        !-----------------------
                        ! End of inequality part
                        !-----------------------

                        IF ( DOSUMOVERPRIMS ) THEN
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

                                                                Intsv(G) = Intsv(G) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)

                                                                !-------------------------------------------------------------
                                                                ! Here we calculate the gradient of the two electron integrals
                                                                !-------------------------------------------------------------

                                                                IF ( CFORCE ) THEN
                                                                        CALL gradprimeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW,gradient)
                                                                        gradIntsv(BAS%PSI(I)%ATYPE,:,G) = gradIntsv(BAS%PSI(I)%ATYPE,:,G) + TERM*gradient(1,:) 
                                                                        gradIntsv(BAS%PSI(J)%ATYPE,:,G) = gradIntsv(BAS%PSI(J)%ATYPE,:,G) + TERM*gradient(2,:)
                                                                        gradIntsv(BAS%PSI(K)%ATYPE,:,G) = gradIntsv(BAS%PSI(K)%ATYPE,:,G) + TERM*gradient(3,:)
                                                                        gradIntsv(BAS%PSI(L)%ATYPE,:,G) = gradIntsv(BAS%PSI(L)%ATYPE,:,G) + TERM*gradient(4,:)
                                                                ENDIF
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDIF
                        IF ( GG .EQ. 1 .AND. IJ .EQ. KL ) THEN
                              SI = SI+1  
                              tempintsl(SI) = Intsv(G)
                              indexusl(SI) = IJ
                        ENDIF
       ENDDO
       !=============================================================
       ! The gathering and broadcasting below is nessescairy in order 
       ! for the Swartz inequality part (see code above) to work.
       !--------------------------------------------------------------
       ! Here we gather all the double electron integrals (i,j,i,j)
       ! into the array tempintsi at thread 0 and broad cast them to 
       ! all other threads. 
       ! Furhermore we collect the indexes IJ = J*(J-1) + I int the 
       ! array indexus and reaindex the array tempintsi into the array
       ! tempints.
       !-------------------------------------------------------------
       IF ( GG .EQ. 1 ) THEN
          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          !---------------------------------------------------------------------------------------
          ! preparation for gathering the data in the local array SI int the array SIC at thread 0
          !---------------------------------------------------------------------------------------
          DO KK=1,numprocessors
             displs(KK) = KK - 1
             rcounts(KK) = 1
          ENDDO
          !--------------------------------------------------------------------------------------
          ! gathering the local SI data ( the number of times IJ=IK ) of the threads into the 
          ! array SIC located at thred 0
          !--------------------------------------------------------------------------------------
          CALL MPI_GATHERV(SI,1,MPI_INTEGER,SIC,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
          !--------------------------------------------------------------------
          ! broadcasting the array SIC located at thread 0 to all other threads
          !--------------------------------------------------------------------
          CALL MPI_BCAST(SIC,numprocessors,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         
          !---------------------------------------------------------------------------------
          ! Preparing for the gathering of the arrays indexusl and tempintsl into the arrays
          !  indexus and tempintsi respectively 
          !---------------------------------------------------------------------------------
          displs(1) = 0
          DO KK=1,numprocessors
             rcounts(KK) = SIC(KK)
             IF ( KK .GE. 2 ) displs(KK) = displs(KK-1) + SIC(KK-1)
          ENDDO
          
          IF ( SIC(id+1) .NE. 0 ) THEN
             !------------------------------------------------------------------------------------------------------------------------
             ! gathering the local thread arrays indexusl and tempintsl into the arrays indexus and tempintsi located at thread 0
             !------------------------------------------------------------------------------------------------------------------------
             CALL MPI_GATHERV(indexusl(1:SIC(id+1)),rcounts(id+1),MPI_INTEGER,indexus,rcounts,displs,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
             CALL MPI_GATHERV(tempintsl(1:SIC(id+1)),rcounts(id+1),MPI_DOUBLE_PRECISION,tempintsi,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          ENDIF

          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

          CALL MPI_BCAST(indexus,STRL,MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
          CALL MPI_BCAST(tempintsi,STRL,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
          !-----------------------------------------
          ! reindixing according to IJ = J*(J-1) + I
          !----------------------------------------- 
          DO KK=1,STRL
             IF( indexus(KK) .NE. 0 ) THEN
               tempints(indexus(KK)) = tempintsi(KK)
             ENDIF
          ENDDO
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              
       
ENDDO
END SUBROUTINE eeints
