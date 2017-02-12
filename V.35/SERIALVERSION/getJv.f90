SUBROUTINE getJv(P,NB,NRED,Intsv,Jout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(NRED)
      DOUBLE PRECISION, INTENT(OUT) :: Jout(NB,NB)
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
                                        VEC(M) = Intsv(ijkl(I,J,K,L))
                                ENDDO
                        ENDDO
                        Jout(I,J) = DOT_PRODUCT(PVEC,VEC) 
                        IF ( I .NE. J ) Jout(J,I) = Jout(I,J)
                ENDDO
      ENDDO
END SUBROUTINE getJv
           


