SUBROUTINE getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kout)
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER, INTENT(IN) :: NB,numprocessors,id
      INTEGER*8, INTENT(IN) :: Istart,Iend,NRED
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(Istart:Iend)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),PVEC(:)
      INTEGER*8 :: I,J,N,M,MM,K,L,G,ierr,displs(numprocessors),rcounts(numprocessors)
      INTEGER, EXTERNAL :: ijkl
      DOUBLE PRECISION :: Ksend(NB,NB),Krecieve(numprocessors),t1,t2,t3,t4,fact
      N = NB*NB
      ALLOCATE(VEC(N),PVEC(N))
      
      !PVEC = reshape(P,(/N/) )
 
      Kout = 0.0d0
      Krecieve = 0.0d0
      Ksend = 0.0d0
      
      DO I=1,numprocessors
         displs(I) = I-1
         rcounts = 1
      ENDDO

      DO G=Istart,Iend
                I = IND1(G)
                K = IND2(G)
                J = IND3(G)
                L = IND4(G)
                
               IF ( I .LE. J ) THEN
                  Ksend(I,J) = Ksend(I,J) + Intsv(G)*P(K,L)
                  Ksend(J,I) =  Ksend(I,J)
               ENDIF

               IF ( K .LE. L .AND. J .LT. L ) THEN
                  Ksend(K,L) = Ksend(K,L) + Intsv(G)*P(I,J)
                  Ksend(L,K) = Ksend(K,L) 
               ENDIF

               IF ( K .LE. J .AND. I .LT. K ) THEN
                  Ksend(K,J) = Ksend(K,J) + Intsv(G)*P(I,L)
                  Ksend(J,K) = Ksend(K,J) 
               ENDIF
               
               IF ( I .LE. L .AND. I .LT. K .AND. J .LT. L ) THEN
                  Ksend(I,L) = Ksend(I,L) + Intsv(G)*P(K,J)
                  Ksend(L,I) = Ksend(I,L) 
               ENDIF
               
      ENDDO
      CALL MPI_REDUCE(Ksend,Kout,NB*NB,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(Kout,NB*NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE getKv
