SUBROUTINE MP2(Cup,Cdown,Ints,NB,Ne,E0,nuce,EMP2,EHFeigenup,EHFeigendown)
      ! This subroutine constructs the many-electron basis set in the form
      ! of slater determinants, these slater determinants are in turn
      ! constructed from the single and double excitations of the
      ! unrestricted solutions to the Hartree-Fock equations. With this
      ! basis-set the hamiltonian matrix in the CI-basis is calculated
      ! from the slater rules (Cook Book page.71-74), and diagonalized.
      ! The electron repulsion tensor, nuclear attraction matrix and
      ! kinetic energy matrix are expressed in terms of primitive
      ! gaussian orbitals. The matrices are preconstructed together with
      ! the corresponding gaussian basis set expansion coefficients by
      ! the python program CISD.py, and consequently is needed as input
      ! to this program. By Petros Souvatzis Thu Nov 24 10:25:47 CET
      ! 2011

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Ne
      DOUBLE PRECISION, INTENT(IN) :: Cup(NB,NB),Cdown(NB,NB),Ints(NB,NB,NB,NB),E0,nuce
      DOUBLE PRECISION, INTENT(IN) :: EHFeigenup(NB),EHFeigendown(NB)
      DOUBLE PRECISION, INTENT(OUT) :: EMP2
      INTEGER :: I,J,N,M,O,P,N2,K,L,K1,II,JJ,DIFFER,Q,A(2),B(2),Nel,NBl
      INTEGER, ALLOCATABLE :: EXC(:,:),V1(:)
      DOUBLE PRECISION, ALLOCATABLE :: Calpha(:),Ca(:),Cbeta(:),Cb(:),EMP2TEMP(:)
      DOUBLE PRECISION :: TEMP1,DEIJKL
      LOGICAL :: EXCITE,SPINCONSERVE,ALLPARALELL

      SPINCONSERVE = .TRUE.
      NBl = NB
      Nel = Ne
      
      ! Calculating the number of possible double excitations
      ! used in the MP2 calculation
      N2 =  (Ne-MOD(Ne,2))*(NB-(Ne-MOD(Ne,2))/2)*(Ne+MOD(Ne,2))*(NB-(Ne+MOD(Ne,2))/2)/4+1
      N2 = N2 + (Ne+MOD(Ne,2))*((Ne+MOD(Ne,2))/2-1)*(NB - (Ne+MOD(Ne,2))/2 )*(NB - (Ne+MOD(Ne,2))/2 - 1 )/8
      N2 = N2 + (Ne-MOD(Ne,2))*((Ne-MOD(Ne,2))/2-1)*(NB - (Ne-MOD(Ne,2))/2 )*(NB - (Ne-MOD(Ne,2))/2 - 1 )/8
      
      print*,' ==========================================================='
      print*,'          NUMBER OF MP2 EXCITATIONS USED=',N2
      print*,' ==========================================================='

      
      ALLOCATE(EXC(N2,2*NB))
      
      ! Setting the HF-ground state in the excitations matrix EXC:
      
      EXC(:,:) = 0
      
      M = (Ne-MOD(Ne,2))/2

      DO I=1,M
        EXC(1,I) = 1
        EXC(1,I+NB) = 1
      ENDDO
      EXC(1,M+NB+MOD(Ne,2)) = 1

      
      K1 = 1

       !---------------------------------------------------------------------
       ! (1) Storing all the different double excitations in the matrix EXC:
       !---------------------------------------------------------------------

       DO I=1,Ne
        DO J=I+1,Ne
          DO K=1,2*NB
           DO L=K+1,2*NB
            IF ( K .GT. M  .AND.  ( K .LE. NB .OR. K .GT. NB+M+MOD(Ne,2)) .AND. L .GT. M .AND. ( L .LE. NB .OR. L .GT. NB+M+MOD(Ne,2)) ) THEN
                EXCITE = .FALSE.
                II = I
                JJ = J
                IF ( I .GT. (Ne-MOD(Ne,2))/2 ) II = I - (Ne-MOD(Ne,2))/2 + NB
                IF ( J .GT. (Ne-MOD(Ne,2))/2 ) JJ = J - (Ne-MOD(Ne,2))/2 + NB
                IF ( .NOT. SPINCONSERVE ) THEN
                        EXCITE = .TRUE.
                ELSE
                        IF ( (II .LE. NB .AND. K .LE. NB ) .OR. (II .GT. NB .AND. K .GT. NB ) ) THEN
                                IF ( (JJ .LE. NB .AND. L .LE. NB ) .OR. (JJ .GT. NB .AND. L .GT. NB ) ) EXCITE = .TRUE.
                        ENDIF
                        IF ( (II .LE. NB .AND. L .LE. NB ) .OR. (II .GT. NB .AND. L .GT. NB ) ) THEN
                                IF ( (JJ .LE. NB .AND. K .LE. NB ) .OR. (JJ .GT. NB .AND. K .GT. NB ) ) EXCITE = .TRUE.
                        ENDIF
                ENDIF
                IF ( EXCITE ) THEN
                        K1 = K1+1
                        DO N=1,M
                                EXC(K1,N) = 1
                                EXC(K1,NB+N) =1
                        ENDDO
                        EXC(K1,NB+M+MOD(Ne,2)) = 1
                        EXC(K1,II) = 0
                        EXC(K1,JJ) = 0
                        EXC(K1,K) = 1
                        EXC(K1,L) = 1
                ENDIF
            ENDIF
           ENDDO
          ENDDO
         ENDDO
        ENDDO
        !------------------------------------------------
        ! (2) Calculating the MP2 energy correction
        !------------------------------------------------
        ALLOCATE(EMP2TEMP(N2))
        EMP2TEMP = 0.0d0
        !$OMP PARALLEL SHARED(Ints,Cup,Cdown,EHFeigenup,EHFeigendown,E0,NB,NBl,N2,EXC,EMP2TEMP) &
        !$OMP & PRIVATE(M,N,K,L,J,A,B,Q,DIFFER,V1,Ca,Calpha,Cb,Cbeta,TEMP1,DEIJKL,ALLPARALELL)
        ALLOCATE(Calpha(NB),Ca(NB),Cbeta(NB),Cb(NB),V1(2*NB))
        !$OMP DO
        DO I=1,N2
            DEIJKL = 0.0d0
            TEMP1  = 0.0d0
            V1 = EXC(1,:)-EXC(I,:)
            DIFFER = Int(DOT_PRODUCT(V1,V1))/2
            IF ( DIFFER .EQ. 2 ) THEN
                     L = 1
                     Q = 1
                     A(:) = 0
                     B(:) = 0
                     DO K=1,2*NB
                        IF  ( V1(K) .EQ. 1 ) THEN
                                A(L) = K
                                L = L+1
                        ENDIF
                        IF  ( V1(K) .EQ. -1 ) THEN
                                B(Q) = K
                                Q = Q+1
                        ENDIF
                     ENDDO
                     
                     ALLPARALELL = .FALSE.
                     ! Checking wheater or not all the spins point in the same
                     ! direction:
                     IF ( A(1) .LE. NB .AND. A(2) .LE. NB .AND. B(1) .LE. NB .AND. B(2) .LE. NB ) ALLPARALELL = .TRUE.
                     IF ( A(1) .GT. NB .AND. A(2) .GT. NB .AND. B(1) .GT. NB .AND. B(2) .GT. NB ) ALLPARALELL = .TRUE.
                        
                     IF ( ALLPARALELL ) THEN
                                     IF ( A(1) .LE. NB ) THEN
                                             Calpha = Cup(:,A(1))
                                             Cbeta = Cup(:,B(1))
                                             Ca = Cup(:,A(2))
                                             Cb = Cup(:,B(2))
                                             DEIJKL = EHFeigenup(A(1))+EHFeigenup(A(2))-EHFeigenup(B(1))-EHFeigenup(B(2))
                                     ENDIF
                                     IF ( A(1) .GT. NB ) THEN
                                             Calpha = Cdown(:,A(1)-NB)
                                             Cbeta = Cdown(:,B(1)-NB)
                                             a = Cdown(:,A(2)-NB)
                                             Cb = Cdown(:,B(2)-NB)
                                             DEIJKL = EHFeigendown(A(1)-NB)+EHFeigendown(A(2)-NB)-EHFeigendown(B(1)-NB)-EHFeigendown(B(2)-NB)
                                     ENDIF
                                     DO M=1,NB
                                        DO N=1,NB
                                                DO K=1,NB
                                                        DO L=1,NB
                                                                TEMP1 = TEMP1 + Calpha(M)*Cbeta(K)*Ca(N)*Cb(L)*Ints(M,K,N,L) - Calpha(M)*Cb(L)*Ca(N)*Cbeta(K)*Ints(M,L,N,K)
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                     ENDDO
                      ENDIF
                     ! if the spins of the differing orbitals are parallell, i.e S(a(1)) = S(b(1))
                     ! and S(a(2)) = S(b(2)):
                     IF ( (( A(1) .LE. NB .AND. B(1) .LE. NB ) .OR. ( A(1) .GT. NB .AND. B(1) .GT. NB ) ) .AND. .not. ALLPARALELL ) THEN
                             IF ( ( A(2) .LE. NB .AND. B(2) .LE. NB ) .OR. ( A(2) .GT. NB .AND. B(2) .GT. NB ) ) THEN
                                     
                                     IF ( A(1) .LE. NB ) THEN
                                             Calpha = Cup(:,A(1))
                                             Cbeta = Cup(:,B(1))
                                             DEIJKL = EHFeigenup(A(1))-EHFeigenup(B(1))
                                     ENDIF
                                     IF ( A(1) .GT. NB ) THEN
                                             Calpha = Cdown(:,A(1)-NB)
                                             Cbeta = Cdown(:,B(1)-NB)
                                             DEIJKL = EHFeigendown(A(1)-NB)-EHFeigendown(B(1)-NB)
                                     ENDIF
                                     IF ( A(2) .LE. NB ) THEN
                                             Ca = Cup(:,A(2))
                                             Cb = Cup(:,B(2))
                                             DEIJKL = DEIJKL + EHFeigenup(A(2))-EHFeigenup(B(2))
                                     ENDIF
                                     IF ( A(2) .GT. NB ) THEN
                                             Ca = Cdown(:,A(2)-NB)
                                             Cb = Cdown(:,B(2)-NB)
                                             DEIJKL = DEIJKL + EHFeigendown(A(2)-NB)-EHFeigendown(B(2)-NB)
                                     ENDIF
                                     DO M=1,NB
                                        DO N=1,NB
                                                DO K=1,NB
                                                        DO L=1,NB
                                                                TEMP1 = TEMP1 + Calpha(M)*Cbeta(K)*Ca(N)*Cb(L)*Ints(M,K,N,L)
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                     ENDDO
                               ENDIF
                      ENDIF
                      ! if the spinns of the differing orbitals are parallell, i.e S(a[1]) == S(b[2])
                      ! and S(a[2]) = S(b[1])
                      IF ( (( A(1) .LE. NB .AND. B(2) .LE. NB ) .OR. ( A(1) .GT. NB .AND. B(2) .GT. NB )) .AND. .not. ALLPARALELL ) THEN
                             IF ( ( A(2) .LE. NB .AND. B(1) .LE. NB ) .OR. ( A(2) .GT. NB .AND. B(1) .GT. NB ) ) THEN
                                     
                                     IF ( A(1) .LE. NB ) THEN
                                             Calpha = Cup(:,A(1))
                                             Cbeta = Cup(:,B(2))
                                             DEIJKL = EHFeigenup(A(1))-EHFeigenup(B(2))
                                     ENDIF
                                     IF ( A(1) .GT. NB ) THEN
                                             Calpha = Cdown(:,A(1)-NB)
                                             Cbeta = Cdown(:,B(2)-NB)
                                             DEIJKL = EHFeigendown(A(1)-NB)-EHFeigendown(B(2)-NB)
                                     ENDIF
                                     IF ( A(2) .LE. NB ) THEN
                                             Ca = Cup(:,A(2))
                                             Cb = Cup(:,B(1))
                                             DEIJKL = DEIJKL + EHFeigenup(A(2))-EHFeigenup(B(1))
                                     ENDIF
                                     IF ( A(2) .GT. NB ) THEN
                                             Ca = Cdown(:,A(2)-NB)
                                             Cb = Cdown(:,B(1)-NB)
                                             DEIJKL = DEIJKL + EHFeigendown(A(2)-NB)-EHFeigendown(B(1)-NB)
                                     ENDIF
                                     
                                     DO M=1,NB
                                        DO N=1,NB
                                                DO K=1,NB
                                                        DO L=1,NB
                                                                TEMP1 = TEMP1 - Calpha(M)*Cbeta(L)*Ca(N)*Cb(K)*Ints(M,L,N,K)
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                     ENDDO
                              ENDIF
                      ENDIF
                      IF ( DEIJKL .NE. 0.0d0 ) EMP2TEMP(I)  = EMP2TEMP(I) + (TEMP1**2)/DEIJKL
             ENDIF
        ENDDO
        !$OMP END DO
        !$OMP END PARALLEL
        EMP2 = SUM(EMP2TEMP)
        print*,' '
        print*,' '
        print*,'==========================================================='
        print*,'            Results from the MP2 calculation:              '
        print*,'==========================================================='
        print*,' '              
        print*,' '
        print*,' '
        !-----------------------------
        ! HERE THE OUTPUT IS GENERATED
        !-----------------------------
        WRITE(*,'(A27,F30.20,A3)')'      Correlation Energy =',EMP2,' au'
        WRITE(*,'(A27,F30.20,A3)')'     Hartree-Fock Energy =',E0,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')' Nuclear repulsio Energy =',nuce,' au'
        IF ( nuce .NE. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2+nuce,' au'
        IF ( nuce .EQ. 0) WRITE(*,'(A27,F30.20,A3)')'            Total Energy =',E0+EMP2,' au'
        print*,' '
        print*,' '
       
        END SUBROUTINE MP2
