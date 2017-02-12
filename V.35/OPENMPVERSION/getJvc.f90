SUBROUTINE getJvc(P,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jout)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: RIAPPROX
      INTEGER, INTENT(IN) :: NB,NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: Intsv(NRED),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION  :: INVRI(NBAUX,NBAUX)
      COMPLEX*16, INTENT(IN) :: P(NB,NB)
      COMPLEX*16, INTENT(OUT) :: Jout(NB,NB)
      COMPLEX*16, ALLOCATABLE :: VEC(:),PVEC(:)
      DOUBLE PRECISION,ALLOCATABLE  :: WEC1(:),WEC2(:),MAT(:,:),MA(:,:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl
      COMPLEX*16 :: RE
     
      N = NB*NB
      RE = (1.0d0,0.0d0)
      IF ( RIAPPROX ) THEN
                ALLOCATE(WEC1(NBAUX),WEC2(NBAUX),MAT(N,NBAUX))
                call invert(VRI,INVRI,NBAUX)
                MAT = reshape(WRI,(/N,NBAUX/))
      ENDIF
 
      ALLOCATE(PVEC(N),MA(NBAUX,N))
      Jout = (0.0d0,0.0d0)

      PVEC = reshape(P,(/N/) )
    
      IF ( RIAPPROX ) MA = MATMUL(INVRI,TRANSPOSE(MAT))

      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC)
      ALLOCATE(VEC(N))
      !$OMP DO
      DO I=1,NB
                DO J=I,NB
                        IF ( RIAPPROX ) THEN
                                Jout(I,J) = DOT_PRODUCT(PVEC,RE*MATMUL(WRI(I,J,:),MA))
                        ELSE
                                M = 0
                                DO L=1,NB
                                        DO K =1,NB
                                                M = M+1
                                                VEC(M) = Intsv(ijkl(I,J,K,L))*RE
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
END SUBROUTINE getJvc
           


