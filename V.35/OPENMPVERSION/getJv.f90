SUBROUTINE getJv(P,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jout)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: RIAPPROX
      INTEGER, INTENT(IN) :: NB,NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(NRED),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: Jout(NB,NB)
      DOUBLE PRECISION  :: INVRI(NBAUX,NBAUX)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),WEC1(:),WEC2(:),PVEC(:),MAT(:,:),MA(:,:),VEC2(:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl

      N = NB*NB
      IF ( RIAPPROX ) THEN
                ALLOCATE(WEC1(NBAUX),WEC2(NBAUX),MAT(N,NBAUX))
                call invert(VRI,INVRI,NBAUX)
                MAT = reshape(WRI,(/N,NBAUX/))
      ENDIF

      ALLOCATE(PVEC(N),MA(NBAUX,N))
      PVEC = reshape(P,(/N/) )
      
      !=======================================================================
      ! Here we use lapack-blas matrix multiplication routine instead of 
      ! matmul since it is much faster. C = ALPHA*op(A)*op(B) + beta*C
      ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
      ! op[A](M,K), op[B](K,N), C(M,N) , C(LDC,*) , A(LDA,*)  B(LDB,*)
      !=======================================================================
                                
      IF ( RIAPPROX ) CALL DGEMM( 'N', 'T', NBAUX, N, NBAUX, 1.0d0, INVRI, NBAUX, MAT, N, 0.0d0, MA, NBAUX )

      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC)
      ALLOCATE(VEC(N))
      !$OMP DO
      DO I=1,NB
                DO J=I,NB
                        IF ( RIAPPROX ) THEN
                                !=======================================================================
                                ! Here we use lapack-blas matrix multiplication routine instead of 
                                ! matmul since it is much faster. C = ALPHA*op(A)*op(B) + beta*C
                                ! DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
                                ! op[A](M,K), op[B](K,N), C(M,N) , C(LDC,*) , A(LDA,*)  B(LDB,*)
                                !=======================================================================
                                
                                
                                CALL DGEMM ( 'T', 'N', N, 1, NBAUX, 1.0d0, MA, NBAUX, WRI(I,J,:), NBAUX, 0.0d0, VEC,N )
                                
                                Jout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                                !Jout(I,J) = DOT_PRODUCT(PVEC,MATMUL(WRI(I,J,:),MA)) 
                        ELSE
                                M = 0
                                DO L=1,NB
                                        DO K =1,NB
                                                M = M+1
                                                VEC(M) = Intsv(ijkl(I,J,K,L))
                                        ENDDO
                                ENDDO
                                Jout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                        ENDIF
                        IF ( I .NE. J ) Jout(J,I) = Jout(I,J)
                ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(VEC)
      !$OMP END PARALLEL
END SUBROUTINE getJv
           


