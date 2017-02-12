SUBROUTINE getKvforceRI(P,NB,NBAUX,NATOMS,COORD,ATOMINDEX,VRI,WRI,gradVRI,gradWRI,Kout)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,NBAUX,NATOMS,COORD,ATOMINDEX
      DOUBLE PRECISION, INTENT(IN) :: P(NB,NB),VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(IN) :: gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: Kout(NB,NB)
      DOUBLE PRECISION  :: INVRI(NBAUX,NBAUX),MAT(NBAUX,NBAUX)
      DOUBLE PRECISION, ALLOCATABLE :: VEC(:),WEC1(:),WEC2(:),dWEC1(:),dWEC2(:),PVEC(:),MA(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MAT1(:,:),MAT2(:,:),dMAT1(:,:),dMAT2(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: M1(:,:,:),M2(:,:,:),M3(:,:,:)
      INTEGER :: I,J,N,M,K,L
      INTEGER, EXTERNAL :: ijkl

      ALLOCATE(WEC1(NBAUX),WEC2(NBAUX),dWEC1(NBAUX),dWEC2(NBAUX))
      ALLOCATE(MAT1(NB,NBAUX),MAT2(NB,NBAUX),dMAT1(NB,NBAUX),dMAT2(NB,NBAUX))
      call invert(VRI,INVRI,NBAUX)
       
      MAT = MATMUL(INVRI,MATMUL(gradVRI(ATOMINDEX,COORD,:,:),INVRI))
      N = NB*NB
      ALLOCATE(VEC(N),PVEC(N))
      ALLOCATE(M1(NB,NBAUX,NB),M2(NB,NBAUX,NB),M3(NB,NBAUX,NB))

      DO J=1,NB
                M1(J,:,:) = MATMUL(INVRI,TRANSPOSE(WRI(:,J,:)))
                M2(J,:,:) = MATMUL(INVRI,TRANSPOSE(gradWRI(ATOMINDEX,COORD,:,J,:)))
                M3(J,:,:) = MATMUL(MAT,TRANSPOSE(WRI(:,J,:)))
      ENDDO
      PVEC = reshape(P,(/N/) )

      Kout = 0.0d0
      VEC = 0.0d0
      !$OMP PARALLEL &
      !$OMP PRIVATE(M,I,J,L,K,VEC,MA)
      ALLOCATE(MA(NB,NB))
      !$OMP DO
      DO I=1,NB
                DO J=I,NB
                        !MA = MATMUL(gradWRI(ATOMINDEX,COORD,I,:,:),MATMUL(INVRI,TRANSPOSE(WRI(:,J,:))))
                        !MA = MA + MATMUL(WRI(I,:,:),MATMUL(INVRI,TRANSPOSE(gradWRI(ATOMINDEX,COORD,:,J,:))))
                        !MA = MA - MATMUL(WRI(I,:,:),MATMUL(MAT,TRANSPOSE(WRI(:,J,:))))
                        MA = MATMUL(gradWRI(ATOMINDEX,COORD,I,:,:),M1(J,:,:))
                        MA = MA + MATMUL(WRI(I,:,:),M2(J,:,:))
                        MA = MA - MATMUL(WRI(I,:,:),M3(J,:,:))
                        Kout(I,J) = DOT_PRODUCT(PVEC,reshape(MA,(/N/))) 
                        !Kout(I,J) = SUM(PVEC*reshape(MA,(/N/))) 
                        IF ( I .NE. J ) Kout(J,I) = Kout(I,J)
                ENDDO
      ENDDO
      !$OMP END DO
      DEALLOCATE(MA)
      !$OMP END PARALLEL
END SUBROUTINE getKvforceRI
           


