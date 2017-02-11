SUBROUTINE eeintsRI(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istartg,Iendg,PRYSR,PRYSW,numprocessors,&
                & APPROXEE,Tol,CFORCE,id,NONZERO,DMAT,dInts,NBAUX,VRI,WRI,gradVRI,gradWRI)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: id,numprocessors,NATOMS,FNATOMS,NBAUX
      INTEGER*8, INTENT(IN) :: Istart,Iend,Istartg,Iendg,NONZERO,NDIAG
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS),dInts(NDIAG)
      DOUBLE PRECISION, INTENT(IN) :: VRI(NBAUX,NBAUX),WRI(BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION, INTENT(IN) :: gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX)
      LOGICAL, INTENT(IN) :: APPROXEE,CFORCE
      INTEGER, INTENT(INOUT) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(INOUT) :: Intsv(Istart:Iend),gradIntsv(FNATOMS,3,Istartg:Iendg)
      INTEGER*8 :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,GG,STRL,KK,PPP,II,LL
      INTEGER*8 :: rcounts(numprocessors),displs(numprocessors),ierr,scount,SI,SIC(numprocessors)
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL,DMAX
      DOUBLE PRECISION :: MAT1(BAS%NBAS,BAS%NBAS,NBAUX),MAT2(FNATOMS,3,BAS%NBAS,BAS%NBAS,NBAUX)
      DOUBLE PRECISION :: VEC(NBAUX),VEC1(NBAUX),VEC2(NBAUX),INVRI(NBAUX,NBAUX),MAT(FNATOMS,3,NBAUX,NBAUX)
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,gradient(4,3),DMAXIMUM(6)
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE

        
      CALL invert(VRI,INVRI,NBAUX)

      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
     
      IF ( NONZERO .EQ. 0 ) THEN
        ! These are dummy entries:
           IND1(1) = 1
           IND2(1) = 1
           IND3(1) = 1
           IND4(1) = 1
      ENDIF      

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

DO GG=1,1
        IF ( NONZERO .NE. 0 ) THEN
        DO NN=Istart,Iend
                        G = NN
                        I = IND1(G) 
                        J = IND2(G) 
                        K = IND3(G) 
                        L = IND4(G) 
                        IJ = (J*(J-1))/2 + I
                        KL = (L*(L-1))/2 + K
                
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

                        IF ( ( IJ .EQ. KL .AND. sqrt(DABS(dInts(IJ)*dInts(KL)))*DMAX .GE. Tol ) .OR. ( IJ .EQ. KL .AND. .not. APPROXEE ) )  THEN
                                  Intsv(G) = dInts(IJ)
                                  IF ( CFORCE ) DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF (  APPROXEE .AND. IJ .NE. KL .AND. sqrt(DABS(dInts(IJ)*dInts(KL)))*DMAX .GE. Tol  ) THEN
                                DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF ( .not. APPROXEE .AND. IJ .NE. KL ) DOSUMOVERPRIMS = .TRUE.
                        !-----------------------
                        ! End of inequality part
                        !-----------------------

                        IF ( DOSUMOVERPRIMS ) THEN
                                !-----------------------------------------------
                                ! Here we calculate the two electron integrals
                                !-----------------------------------------------
                                IF ( KL .NE. IJ ) THEN
                                        !=======================================================================
                                        ! Here we use lapack-blas matrix multiplication routine instead of 
                                        ! matmul since it is much faster. C = ALPHA*op(A)*op(B) + beta*C
                                        ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC ) 
                                        ! op[A](M,K), op[B](K,N), C(M,N) , C(LDC,*) , A(LDA,*)  B(LDB,*)
                                        !=======================================================================

                                        CALL DGEMM ( 'N', 'N', NBAUX, 1, NBAUX, 1.0d0, INVRI, NBAUX, WRI(K,L,:), NBAUX, 0.0d0,VEC,NBAUX )
                                        Intsv(G) = DOT_PRODUCT(WRI(I,J,:),VEC)

                                ENDIF
                                !-------------------------------------------------------------
                                ! Here we calculate the gradient of the two electron integrals
                                !-------------------------------------------------------------

                                IF ( CFORCE ) THEN
                                        DO LL=1,NATOMS

                                                gradIntsv(LL,:,G) = MATMUL(gradWRI(LL,:,I,J,:),MAT1(K,L,:))+MATMUL(gradWRI(LL,:,K,L,:),MAT1(I,J,:))

                                                !gradIntsv(LL,1,G) = gradIntsv(LL,1,G) - DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,1,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,1,G) = gradIntsv(LL,1,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,1,K,L,:))

                                                !gradIntsv(LL,2,G) = gradIntsv(LL,2,G) -DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,2,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,2,G) = gradIntsv(LL,2,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,2,K,L,:))

                                                !gradIntsv(LL,3,G) = gradIntsv(LL,3,G) - DOT_PRODUCT(WRI(I,J,:),MATMUL(MAT(LL,3,:,:),WRI(K,L,:)))
                                                gradIntsv(LL,3,G) = gradIntsv(LL,3,G) - DOT_PRODUCT(WRI(I,J,:),MAT2(LL,3,K,L,:))
                                        ENDDO

                                ENDIF
                        ENDIF
       ENDDO
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              
       
ENDDO
END SUBROUTINE eeintsRI
