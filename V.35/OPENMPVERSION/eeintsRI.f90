SUBROUTINE eeintsRI(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,NRED,FNRED,PRYSR,PRYSW,APPROXEE,Tol,CFORCE,DMAT,NBAUX,VRI,WRI,gradVRI,gradWRI)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE functionsmodule
      IMPLICIT NONE
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: NATOMS,FNATOMS,NBAUX
      INTEGER*8, INTENT(IN) :: NRED,FNRED
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS)
      LOGICAL, INTENT(IN) :: APPROXEE,CFORCE
      DOUBLE PRECISION, INTENT(IN) :: VRI(NBAUX,NBAUX),WRI(BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION, INTENT(IN) :: gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: Intsv(NRED),gradIntsv(FNATOMS,3,FNRED)
      INTEGER :: IMAP(NRED,6),Imem,Jmem,Kmem,Lmem
      INTEGER :: I,J,K,L,M,N,NN,MM,P,Q,G,GG,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,LL
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,DMAX,MAT1(BAS%NBAS,BAS%NBAS,NBAUX),MAT2(FNATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,INVRI(NBAUX,NBAUX),MAT(FNATOMS,3,NBAUX,NBAUX)
      DOUBLE PRECISION, ALLOCATABLE :: gradient(:,:),VEC(:),VEC1(:),VEC2(:),DMAXIMUM(:)
      INTEGER, EXTERNAL :: ijkl
      EXTERNAL :: gradprimeeintegral
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
        
      call invert(VRI,INVRI,NBAUX)
      !INVRI = inv(VRI)  ! See functionsmodule.f90 

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
     gradIntsv(:,:,:) = 0.0d0

     IF ( CFORCE ) THEN
             DO LL=1,NATOMS
                DO K =1,3
                        MAT(LL,K,:,:) = MATMUL(INVRI,MATMUL(gradVRI(LL,K,:,:),INVRI))
                ENDDO
             ENDDO

                DO I=1,BAS%NBAS
                        DO J=I,BAS%NBAS
                                MAT1(I,J,:) = MATMUL(INVRI,WRI(I,J,:))
                                IF ( I .NE. J ) MAT1(J,I,:) = MAT1(I,J,:)
                                
                                DO LL=1,NATOMS
                                        DO K=1,3
                                                MAT2(LL,K,I,J,:) = MATMUL(MAT(LL,K,:,:),WRI(I,J,:))
                                                IF ( I .NE. J ) MAT2(LL,K,J,I,:) = MAT2(LL,K,I,J,:)
                                        ENDDO
                                ENDDO
                        ENDDO
                ENDDO

     ENDIF
                       
     DO GG=1,2
         !$OMP PARALLEL SHARED(Intsv,gradIntsv,IMAP,BAS,CFORCE) &
         !$OMP PRIVATE(NN,MM,NSTART,M,N,P,Q,I,J,K,L,G,IJ,KL,L1,L2,L3,L4,M1,M2,M3,M4,N1,N2,N3,N4,A,B,C,D,NO1,NO2,NO3,NO4,gradient,VEC,VEC1,VEC2) &
         !$OMP PRIVATE(DOSUMOVERPRIMS,SAMECITE,al1,al2,al3,al4,co1,co2,co3,co4,TERM,NP1,NP2,NP3,NP4,LL,Imem,Jmem,Kmem,Lmem,DMAXIMUM,DMAX)
         ALLOCATE(gradient(4,3),VEC(NBAUX),VEC1(NBAUX),VEC2(NBAUX),DMAXIMUM(6))
         Imem = 0
         Jmem = 0
         Kmem = 0
         Lmem = 0
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
                                !Intsv(G) = DOT_PRODUCT(WRI(I,J,:),MATMUL(INVRI,WRI(K,L,:)))

                                !=======================================================================
                                ! Here we use lapack-blas matrix multiplication routine instead of 
                                ! matmul since it is much faster. C = ALPHA*op(A)*op(B) + beta*C
                                ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
                                ! op[A](M,K), op[B](K,N), C(M,N) , C(LDC,*) , A(LDA,*)  B(LDB,*)
                                !=======================================================================
                                
                                CALL DGEMM ( 'N', 'N', NBAUX, 1, NBAUX, 1.0d0, INVRI, NBAUX, WRI(K,L,:), NBAUX, 0.0d0,VEC,NBAUX )
                                Intsv(G) = DOT_PRODUCT(WRI(I,J,:),VEC)
                                !--------------------------------------------------------------
                                ! Here we caclulate the gradient of  the two electron integrals
                                !--------------------------------------------------------------

                                IF ( CFORCE ) THEN
                                        DO LL=1,NATOMS
                                                
                                                gradIntsv(LL,:,G) = MATMUL(gradWRI(LL,:,I,J,:),MAT1(K,L,:))+MATMUL(gradWRI(LL,:,K,L,:),MAT1(I,J,:))
                                                
                                                !gradIntsv(LL,1,G) = gradIntsv(LL,1,G) - DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,1,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,1,G) = gradIntsv(LL,1,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,1,K,L,:))
                                               
                                                !gradIntsv(LL,2,G) = gradIntsv(LL,2,G) - DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,2,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,2,G) = gradIntsv(LL,2,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,2,K,L,:))
                                                
                                                !gradIntsv(LL,3,G) = gradIntsv(LL,3,G) - DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,3,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,3,G) = gradIntsv(LL,3,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,3,K,L,:))
                                        ENDDO
                                ENDIF
                        ENDIF
                ENDDO
        ENDDO  
        !$OMP END DO
        DEALLOCATE(gradient,VEC,VEC1,VEC2,DMAXIMUM)
        !$OMP END PARALLEL
ENDDO
END SUBROUTINE eeintsRI

