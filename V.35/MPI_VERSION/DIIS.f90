SUBROUTINE DIIS(NB,N,ERR,Pprevious,PNEW,LAMDA,INFO)
     ! This is the direct inversion in the iterative subspace method DIIS,
     ! by P. Puley in CHEM. Phys. Lett. 73, 393 (1984). The routine is 
     ! used to estimate the self-consistent density-matrix in the region 
     ! where the density dependence of the energy can be approximated well
     ! up to second order.
     ! Also, see T. Halgaker, P. Jorgensen and J. Olsen 
     ! "Molecular Electronic Structure Theory" p.460-463
     IMPLICIT NONE
     INTEGER, INTENT(IN) :: NB,N
     DOUBLE PRECISION, INTENT(IN) :: Pprevious(50,NB,NB),ERR(50,NB,NB)
     DOUBLE PRECISION, INTENT(OUT) :: PNEW(NB,NB),LAMDA
     INTEGER, INTENT(OUT) :: INFO
     DOUBLE PRECISION, ALLOCATABLE :: DP(:,:),B(:,:),BINV(:,:),ALPHA(:),C(:),EIGENVAL(:),EIGENVECT(:,:)
     DOUBLE PRECISION :: CONDNUMBER 
     INTEGER :: I,J,K,NN,NH,LDA,LDB,NRHS
     INTEGER, ALLOCATABLE :: IPIV(:)
     EXTERNAL :: DGESV

     IF ( N .GT. 50 ) THEN
        WRITE(*,*)'Error in DIIS.f90, maximum dimension of subspace is 50!'
        STOP
     ENDIF

     NN = NB*NB
     NH = N

     ALLOCATE(DP(NH,NN),B(NH+1,NH+1),ALPHA(NH+1),BINV(NH+1,NH+1),C(NH+1),IPIV(NH+1),EIGENVAL(NH+1),EIGENVECT(NH+1,NH+1))
    
     
     ! Calculating the matrix B ( Eq. (6) in CHEM. Phys. Lett. 73, 393 (1984) )
     DO I =1,NH
      DO J=1,NH
         B(I,J) = SUM(ERR(I,:,:)*ERR(J,:,:))
      ENDDO
     ENDDO

     B(NH+1,:) = -1.0d0
     B(:,NH+1) = -1.0d0
     B(NH+1,NH+1) = 0.0d0
     C = 0.0d0 
     C(NH+1) = -1.0d0

     BINV = B
     ALPHA = C
     
     LDA = NH+1
     LDB = NH+1
     NRHS = 1
     
     CALL diagh( B,NH+1,EIGENVAL,EIGENVECT,INFO)
     IF ( INFO .NE. 0 ) THEN
        CONDNUMBER = DABS(MAXVAL(EIGENVAL)/MINVAL(EIGENVAL))
     ELSE
        CONDNUMBER = 1.0E6
     ENDIF

     ! Solving B*ALPHA = C 
     IF ( INFO .EQ. 0  ) THEN
        CALL DGESV( NH+1, NRHS, BINV, LDA, IPIV, ALPHA, LDB, INFO )
     ELSE
             INFO = -1
     ENDIF
     
     ! Saving the Lagrange multilyer for the ouput
     LAMDA = ALPHA(NH+1)
     
     ! Calculating the new density matrix:
     PNEW = 0.0d0
     DO I=1,NH
         PNEW = PNEW + Pprevious(I,:,:)*ALPHA(I)
     ENDDO
     DEALLOCATE(DP,B,ALPHA,BINV,C,IPIV)

END SUBROUTINE DIIS
