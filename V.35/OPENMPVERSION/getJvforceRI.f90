SUBROUTINE getJvforceRI(P,NB,NBAUX,NATOMS,COORD,ATOMINDEX,VRI,WRI,gradVRI,gradWRI,Jout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,NBAUX,NATOMS,COORD,ATOMINDEX
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(IN) :: gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: Jout(NB,NB)
      DOUBLE PRECISION  :: INVRI(NBAUX,NBAUX),MAT(NBAUX,NBAUX)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),WEC1(:),WEC2(:),dWEC1(:),dWEC2(:),PVEC(:),MAT1(:,:),MAT2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: M1(:,:),M2(:,:),M3(:,:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl

      
      N = NB*NB
      ALLOCATE(WEC1(NBAUX),WEC2(NBAUX),dWEC1(NBAUX),dWEC2(NBAUX),MAT1(N,NBAUX),MAT2(N,NBAUX))
      call invert(VRI,INVRI,NBAUX)
      
      MAT = MATMUL(INVRI,MATMUL(gradVRI(ATOMINDEX,COORD,:,:),INVRI))
      
      MAT1 = reshape(WRI,(/N,NBAUX/))
      MAT2 = reshape(gradWRI(ATOMINDEX,COORD,:,:,:),(/N,NBAUX/))
     
      ALLOCATE(M1(NBAUX,N),M2(NBAUX,N),M3(NBAUX,N))
      ALLOCATE(PVEC(N))
       
      M1 = MATMUL(INVRI,TRANSPOSE(MAT1))
      M2 = MATMUL(INVRI,TRANSPOSE(MAT2))
      M3 = MATMUL(MAT,TRANSPOSE(MAT1))
      
      PVEC = reshape(P,(/N/) )
       
      Jout = 0.0d0

      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC)
      ALLOCATE(VEC(N))
      !$OMP DO 
      DO I=1,NB
                DO J=I,NB
                        !=======================================================================
                        ! Here we use lapack-blas matrix multiplication  routine instead of 
                        ! matmul since it is much faster. C =  ALPHA*op(A)*op(B) + beta*C
                        ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A,  LDA, B, LDB, BETA, C, LDC )
                        ! op[A](M,K), op[B](K,N), C(M,N) , C(LDC,*) , A(LDA,*)  B(LDB,*)
                        !=======================================================================
                        !VEC = 0.0d0
        
                        !CALL DGEMM ( 'T', 'N', N, 1, NBAUX, 1.0d0, M1,  NBAUX, gradWRI(ATOMINDEX,COORD,I,J,:), NBAUX, 1.0d0, VEC, N )
                        !CALL DGEMM ( 'T', 'N', N, 1, NBAUX, 1.0d0, M2,  NBAUX, WRI(I,J,:), NBAUX, 1.0d0, VEC, N )
                        !CALL DGEMM ( 'T', 'N', N, 1, NBAUX, -1.0d0, M3,  NBAUX, WRI(I,J,:), NBAUX, 1.0d0, VEC, N )
                        
                        VEC = MATMUL(gradWRI(ATOMINDEX,COORD,I,J,:),M1)
                        VEC = VEC + MATMUL(WRI(I,J,:),M2)
                        VEC = VEC - MATMUL(WRI(I,J,:),M3)
                        Jout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                        !Jout(I,J) = SUM(PVEC*VEC) 
                        IF ( I .NE. J ) Jout(J,I) = Jout(I,J)
                ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(VEC)
      !$OMP END PARALLEL
END SUBROUTINE getJvforceRI
           


