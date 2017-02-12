SUBROUTINE getKv(P,NB,NRED,Intsv,Kout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(NRED)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),PVEC(:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl

      
      N = NB*NB
      ALLOCATE(VEC(N),PVEC(N))
      
      PVEC = reshape(P,(/N/) )
      
      DO I=1,NB
                DO J=I,NB
                        M = 0
                        DO L=1,NB
                                DO K =1,NB
                                        M = M+1
                                        VEC(M) = Intsv(ijkl(I,K,L,J))
                                ENDDO
                        ENDDO
                        Kout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                        IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
                ENDDO
      ENDDO
END SUBROUTINE getKv
           


