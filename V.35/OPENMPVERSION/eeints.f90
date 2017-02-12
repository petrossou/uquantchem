SUBROUTINE eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,Tol,CFORCE,DMAT)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: NATOMS,FNATOMS
      INTEGER*8, INTENT(IN) :: NRED,FNRED
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS)
      LOGICAL, INTENT(IN) :: APPROXEE,CFORCE
      DOUBLE PRECISION, INTENT(OUT) :: Intsv(NRED),gradIntsv(FNATOMS,3,FNRED)
      INTEGER :: IMAP(NRED,6)
      INTEGER :: I,J,K,L,M,N,NN,MM,P,Q,G,GG,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,DMAX
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4
      DOUBLE PRECISION, ALLOCATABLE :: gradient(:,:),DMAXIMUM(:)
      INTEGER, EXTERNAL :: ijkl
      EXTERNAL :: gradprimeeintegral
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
      
      DO J=1,BAS%NBAS
        DO I=1,J
            DO L=1,BAS%NBAS
              DO K=1,L
                M = (J*(J-1))/2 + I
                N = (L*(L-1))/2 + K
                IF ( M .LE. N ) THEN
                        P = (N*(N-1))/2 + M
                        IMAP(P,1) = I
                        IMAP(P,2) = J
                        IMAP(P,3) = K
                        IMAP(P,4) = L
                        IMAP(P,5) = M
                        IMAP(P,6) = N
                ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     
     Intsv(:) = 0.0d0
     IF ( FNATOMS .NE. 1 .AND. FNRED .NE. 1 ) gradIntsv(:,:,:) = 0.0d0

     DO GG=1,2
         !$OMP PARALLEL SHARED(Intsv,gradIntsv,IMAP,BAS,CFORCE) &
         !$OMP PRIVATE(NN,MM,NSTART,M,N,P,Q,I,J,K,L,G,IJ,KL,L1,L2,L3,L4,M1,M2,M3,M4,N1,N2,N3,N4,A,B,C,D,NO1,NO2,NO3,NO4,gradient) &
         !$OMP PRIVATE(DOSUMOVERPRIMS,SAMECITE,al1,al2,al3,al4,co1,co2,co3,co4,TERM,NP1,NP2,NP3,NP4,DMAXIMUM,DMAX)
         ALLOCATE(gradient(4,3),DMAXIMUM(6))
         !$OMP DO
        DO NN=1,((BAS%NBAS+1)*BAS%NBAS)/2
                IF ( GG .EQ. 1 ) NSTART = NN
                IF ( GG .EQ. 2 ) NSTART = 1
                DO MM=NSTART,NN
                        G = (NN*(NN-1))/2 + MM
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
               
                        DOSUMOVERPRIMS = .FALSE.
                        ! Here we use Schwartz inequality |(i,j,k,l)| <= ! [ |(i,j,i,j)|*(k,l,k,l) ]**(1/2). 
                        !See Eqn 9.12.25 p.404 in T. Helgaker et al's Molecular Electronic-Structure Theory
                        ! We use this together with Eqn. (56) and (57) page. 86,
                        ! in Haettig, Multiscale Simulation Methods in Molecular
                        ! Sciences, J. Grotendorst, N. Attig, S. Blugel, D. Marx (Eds.),
                        !Institute for Advanced Simulation, Forschungszentrum Julich, NIC Series, Vol. 42, ISBN 978-3-9810843-8-2, pp.  77-120, 2009.
                        
                        
                        DMAXIMUM = (/ 2*DABS(DMAT(I,J)),2*DABS(DMAT(K,L)),DABS(DMAT(J,K)),DABS(DMAT(I,K)),DABS(DMAT(J,L)),DABS(DMAT(I,L)) /)
                       
                        DMAX = MAXVAL(DMAXIMUM)

                        IF ( ( GG .EQ. 1 .AND. IJ .EQ. KL ) .OR. ( GG .EQ.  2 .AND. IJ .NE. KL .AND. sqrt(DABS(Intsv(ijkl(I,J,I,J))*Intsv(ijkl(K,L,K,L))))*DMAX .NE. 0.0d0 ) ) THEN
                                DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF (  APPROXEE .AND. GG .EQ.  2 .AND. IJ .NE. KL .AND. sqrt(DABS(Intsv(ijkl(I,J,I,J))*Intsv(ijkl(K,L,K,L))))*DMAX .LT. Tol  ) THEN
                                DOSUMOVERPRIMS = .FALSE.
                        ENDIF
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
                                                                ! Here we ! calculate the two electron integrals
                                                                !-----------------------------------------------

                                                                Intsv(G) = Intsv(G) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)

                                                                !--------------------------------------------------------------
                                                                ! Here we caclulate the gradient of  the two electron integrals
                                                                !--------------------------------------------------------------

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
                        !Ints(J,I,K,L) = Ints(I,J,K,L)
                        !Ints(I,J,L,K) = Ints(I,J,K,L)
                        !Ints(J,I,L,K) = Ints(I,J,K,L)

                        !Ints(K,L,I,J) = Ints(I,J,K,L)
                        !Ints(L,K,I,J) = Ints(I,J,K,L)
                        !Ints(K,L,J,I) = Ints(I,J,K,L)
                        !Ints(L,K,J,I) = Ints(I,J,K,L)

                                        
                ENDDO 
        ENDDO  
        !$OMP END DO
        DEALLOCATE(gradient,DMAXIMUM)
        !$OMP END PARALLEL
ENDDO
END SUBROUTINE eeints
