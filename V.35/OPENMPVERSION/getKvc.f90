SUBROUTINE getKvc(P,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kout)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: RIAPPROX
      INTEGER, INTENT(IN) :: NB,NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: Intsv(NRED),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION  :: INVRI(NBAUX,NBAUX)
      COMPLEX*16, INTENT(IN) :: P(NB,NB)
      COMPLEX*16, INTENT(OUT) :: Kout(NB,NB)
      COMPLEX*16, ALLOCATABLE :: VEC(:),PVEC(:)
      DOUBLE PRECISION, ALLOCATABLE :: WEC1(:),WEC2(:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl
      COMPLEX*16 :: RE 
      RE = (1.0d0,0.0d0)
     
      IF ( RIAPPROX ) THEN
                ALLOCATE(WEC1(NBAUX),WEC2(NBAUX))
                call invert(VRI,INVRI,NBAUX)
      ENDIF
 
      N = NB*NB
      ALLOCATE(PVEC(N))
      Kout = (0.0d0,0.0d0)

      PVEC = reshape(P,(/N/) )
      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC)
      ALLOCATE(VEC(N))
      !$OMP DO    
      DO I=1,NB
                DO J=I,NB
                        IF ( RIAPPROX ) THEN
                                Kout(I,J) = DOT_PRODUCT(PVEC,RE*reshape(MATMUL(WRI(I,:,:),MATMUL(INVRI,TRANSPOSE(WRI(:,J,:)))),(/N/)))
                        ELSE
                                M = 0
                                DO L=1,NB
                                        DO K =1,NB
                                                M = M+1
                                                VEC(M) = Intsv(ijkl(I,K,L,J))*RE
                                        ENDDO
                                ENDDO
                                Kout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                        ENDIF
                        IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
                ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(VEC)
      !$OMP END PARALLEL
END SUBROUTINE getKvc
           


