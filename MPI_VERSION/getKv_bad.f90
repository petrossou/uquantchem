SUBROUTINE getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kout)
      IMPLICIT NONE
      INCLUDE "mpif.h"
      INTEGER, INTENT(IN) :: NB,NRED,Istart,Iend,numprocessors,id
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),Intsv(Istart:Iend)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),PVEC(:)
      INTEGER :: I,J,N,M,MM,K,L,G,ierr,displs(numprocessors),rcounts(numprocessors)
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
                
                IF (  I .LE. J .AND. K .LE. L ) THEN
                   Ksend(I,J) = Ksend(I,J) + Intsv(G)*P(K,L)
                   IF ( K .NE. L .AND. I .EQ. J ) Ksend(I,J) = Ksend(I,J) + Intsv(G)*P(K,L)
                ENDIF
                IF ( K .LE. L .AND. J .LT. L .AND. I .LE. J ) THEN
                   Ksend(K,L) = Ksend(K,L) + Intsv(G)*P(I,J)
                   IF ( I .NE. J .AND. K .EQ. L ) Ksend(K,L) = Ksend(K,L) + Intsv(G)*P(I,J)
                ENDIF
                IF ( K .LE. J .AND. I .LT. K  .AND.  I .LE. L ) THEN
                   Ksend(K,J) = Ksend(K,J) + Intsv(G)*P(I,L)
                   IF ( I .NE. L .AND. K .EQ. J ) Ksend(K,J) = Ksend(K,J) + Intsv(G)*P(I,L)
                ENDIF
                IF ( I .LE. L  .AND. I .LT. K .AND. J .LT. L .AND.  K .LE. J) THEN
                   Ksend(I,L) = Ksend(I,L) + Intsv(G)*P(K,J)
                   IF ( K .NE. J .AND. I .EQ. L ) Ksend(I,L) = Ksend(I,L) + Intsv(G)*P(K,J)
                ENDIF
                
               IF (  I .LT. J .AND. K .GT. L ) THEN
                   Ksend(I,J) = Ksend(I,J) + Intsv(G)*P(K,L)
                ENDIF
                IF ( K .LT. L  .AND. J .LT. L .AND.  I .GT. J ) THEN
                   Ksend(K,L) = Ksend(K,L) + Intsv(G)*P(I,J)
                ENDIF
                IF ( K .LT. J .AND. I .LT. K .AND. I .GT. L ) THEN
                   Ksend(K,J) = Ksend(K,J) + Intsv(G)*P(I,L)
                ENDIF
                IF ( I .LT. L .AND. I .LT. K .AND. J .LT. L .AND. K .GT. J) THEN
                   Ksend(I,L) = Ksend(I,L) + Intsv(G)*P(K,J)
                ENDIF
      ENDDO
      DO J=1,NB
         DO I = 1,J
             Krecieve = 0.0d0
             CALL MPI_GATHER(Ksend(I,J), 1, MPI_DOUBLE_PRECISION, Krecieve(1+id), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
             !CALL MPI_GATHERV(Ksend(I,J),1,MPI_DOUBLE_PRECISION,Krecieve,rcounts,displs,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
             CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
             Kout(I,J) = SUM(Krecieve)  ! The true Kout(I,J) only exists on thread 0 now !!
             IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
             IF ( id .EQ. 0 ) THEN 
                print*,I,J,Kout(J,I)
             ENDIF
         ENDDO
      ENDDO
      STOP
      CALL MPI_BCAST(Kout,NB*NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE getKv
