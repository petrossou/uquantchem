SUBROUTINE getKv(P,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kout)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: RIAPPROX
      INTEGER, INTENT(IN) :: NB,NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(NRED),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION :: INVRI(NBAUX,NBAUX),MAT1(NB,NBAUX),MAT2(NB,NBAUX),MAT3(NB,NB)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),PVEC(:),WEC1(:),WEC2(:),MA(:,:),MA2(:,:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl

      IF ( RIAPPROX ) THEN
                ALLOCATE(WEC1(NBAUX),WEC2(NBAUX))
                call invert(VRI,INVRI,NBAUX)
      ENDIF 
 
      N = NB*NB
      ALLOCATE(PVEC(N))
      
      PVEC = reshape(P,(/N/) )
      
      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC,MA,MA2)
      ALLOCATE(VEC(N),MA(NBAUX,NB),MA2(NB,NB))
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

                                CALL DGEMM ( 'N', 'T', NBAUX, NB, NBAUX, 1.0d0, INVRI, NBAUX, WRI(:,J,:), NB, 0.0d0, MA, NBAUX )
                                
                                CALL  DGEMM ( 'N', 'N', NB, NB, NBAUX, 1.0d0, WRI(I,:,:), NB, MA, NBAUX, 0.0d0, MA2, NB )
                                
                                Kout(I,J) = DOT_PRODUCT(PVEC,reshape(MA2,(/N/)))
                                
                        ELSE
                                M = 0
                                DO L=1,NB
                                        DO K =1,NB
                                                M = M+1
                                                VEC(M) = Intsv(ijkl(I,K,L,J))
                                        ENDDO
                                ENDDO
                                Kout(I,J) = DOT_PRODUCT(PVEC,VEC)
                        ENDIF
                        IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
                ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(VEC,MA,MA2)
      !$OMP END PARALLEL
END SUBROUTINE getKv
           


