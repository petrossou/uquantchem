SUBROUTINE RHF(S,H0,Intsv,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,C,P,MIX,DIISORD,DIISSTART,NSCF,FIXNSCF,POUT,SCRATCH,ZEROSCF)
      ! This subroutine calculates the self consistent Restricted
      ! Hartree-Fock solution
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF
      INTEGER, INTENT(IN) :: NB,Ne,FIXNSCF
      INTEGER*8, INTENT(IN) :: NRED
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),H0(NB,NB),Intsv(NRED),Tol,nucE,MIX
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigen(NB),ETOT
      DOUBLE PRECISION, INTENT(INOUT) :: C(NB,NB),P(NB,NB)
      DOUBLE PRECISION :: J(NB,NB),K(NB,NB),F(NB,NB),G(NB,NB),C2(NB,NB),C3(NB,NB),DE,EOLD
      DOUBLE PRECISION :: Pold(NB,NB),Ps(50,NB,NB),Pt(NB,NB),LAMDA,DELTAP,MIXING
      DOUBLE PRECISION :: ERRS(50,NB,NB),ERR(NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB)
      INTEGER :: I,II,III,L,M,N,INFO
      INTEGER :: MAXITER
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO

      MIXING = MIX

      IF ( DIISORD .GT. 25 ) THEN
              WRITE(*,*)'The order of DIIS mixing set =',DIISORD,' Is too Big.'
              WRITE(*,*)'Maximum allowed value is 25. This value has been used instead.'
              DIISORD = 25
      ENDIF
     
      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIGENVECT)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
      SH = MATMUL(EIGENVECT,MATMUL(SL,TRANSPOSE(EIGENVECT)))

      IF ( SCRATCH ) CALL diaghHF( H0,S,NB,EHFeigen,C)
        
      IF ( MOD(Ne,2) .EQ. 0 ) THEN
                N = Ne/2
      ELSE
                print*,'----------------------------------------------------------------------------'
                print*,'Attempting to run Restricted Hartree-Fock for open-shelled system. ABORTING!'
                print*,'              Change input parameter CORRLEVEL to URHF                      '
                print*,'----------------------------------------------------------------------------'
                STOP
      ENDIF
        

      DE = 2.0d0*Tol
      I = 0
      II = 0
      Pold = 0.0d0
      STARTPRINTDIISIFO = .FALSE.
      IF ( FIXNSCF .GT. 0 ) THEN
              MAXITER =  FIXNSCF-1
      ELSE
              MAXITER = 400
      ENDIF
                
      !=======================================================
      ! This is used for the zero scf cycle (FAST-MD)
      !=======================================================
      IF ( ZEROSCF ) THEN
                CALL getJv(P,NB,NRED,Intsv,J)
                CALL getKv(P,NB,NRED,Intsv,K)
                F = H0 + J - 0.50d0*K
                G = J - 0.50d0*K
                CALL diaghHF( F,S,NB,EHFeigen,C)
                ETOT = SUM(F*P)-0.5*SUM(G*P) + nucE
                C2(:,:) = 0.0d0
                DO M=1,N
                        C2(:,M) = C(:,M)
                ENDDO
                CALL makedens(C2,NB,P)
                P = 2.0d0*P
                NSCF = 0
                RETURN
      ENDIF
      !=======================================================

      DO WHILE ( (DABS(DE) .GT. Tol .AND. I .LE. MAXITER) .OR. ( DELTAP .GT.  sqrt(Tol) .AND. I .LE. MAXITER) )

                C2(:,:) = 0.0d0
                DO M=1,N
                        C2(:,M) = C(:,M)
                ENDDO
                
                 !==============================================================
                 ! In the case of XL-BOMD, it is crusial that the scf starts 
                 ! with the density matrices Pup and Pdown provided by the input, 
                 ! since it is these matrices that have been propagated by the 
                 ! XL-BOMD scheme and not the coefficient matrices Cup and Cdown
                 ! Therefor we have the if-construct here.
                 !==============================================================
                 IF ( I .GT. 0 ) CALL makedens(C2,NB,P)

                ! Calculating the change in the density matrix
                DELTAP = sqrt(DOT_PRODUCT(reshape(2.0d0*P-Pold,(/NB**2/)),reshape(2.0d0*P-Pold,(/NB**2/))))

                ! Linear mixing step
                IF ( I .GT. 0 ) THEN
                        P = 2.0d0*P*(1.0d0 - MIXING) + POLD*MIXING
                ELSE
                        P = 2.0d0*P
                ENDIF

                POLD = P
               
                ! Saving the density matrices from the previous
                ! iterations, to be used with the DIIS-method.
                IF ( II .LT. 2*DIISORD ) THEN
                        II = II + 1
                        Ps(II,:,:) = P
                        ERRS(II,:,:) = ERR
                ELSE
                        ! refreshing the density matrices used by the DIIS
                        IF ( DIISORD .NE. 0 ) THEN
                                DO III=1,II-1
                                        Ps(III,:,:) = Ps(III+1,:,:)
                                        ERRS(III,:,:) = ERRS(III+1,:,:)
                                ENDDO
                                Ps(II,:,:) = P
                                ERRS(II,:,:) = ERR
                        ENDIF
                ENDIF

                ! Here we use the DIIS method ( CHEM. Phys. Lett. 73, 393 (1980) )
                ! in order to estimate the self-consistent density matrices:
                IF ( II .GE. 4 ) THEN
                        
                        CALL DIIS(NB,II,ERRS,Ps,Pt,LAMDA,INFO)
                
                        IF ( (INFO .EQ. 0 .AND. LAMDA .LT. 1.0d0) .OR. ( INFO .EQ. 0 .AND. I .GT. DIISSTART ) ) THEN
                                IF ( I .GE. DIISSTART ) THEN
                                        P = Pt
                                        STARTPRINTDIISIFO = .TRUE.
                                        MIXING = 0.0d0
                                ENDIF
                        ENDIF
                ENDIF

                CALL getJv(P,NB,NRED,Intsv,J)
                CALL getKv(P,NB,NRED,Intsv,K)

                F = H0 + J - 0.50d0*K
                
                G = J - 0.50d0*K
                
                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERR = MATMUL(F,MATMUL(P,S)) - MATMUL(S,MATMUL(P,F))
                        ERR = MATMUL(TRANSPOSE(SH),MATMUL(ERR,SH))
                ENDIF

                CALL diaghHF( F,S,NB,EHFeigen,C)

                ETOT = SUM(F*P)-0.50d0*SUM(G*P) + nucE

                IF ( I .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'         Entering the scf restricted Hartree-Fock loop        '
                        print*,'   =========================================================='
                        print*,' '
                        WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LT. 0 ) DE = ETOT-EOLD
        
                IF ( STARTPRINTDIISIFO ) THEN
                        IF ( POUT ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDA
                ELSE
                        IF ( POUT ) WRITE(*,'(I4,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP
                ENDIF
        
                EOLD = ETOT
                I = I+1
     ENDDO
     
     NSCF = I
     C3 = 0.0d0
     DO M=1,N
         C3(:,M) = C(:,M)
     ENDDO
     CALL makedens(C3,NB,P)
     P = 2.0d0*P

     IF ( I .LE. MAXITER .OR. FIXNSCF .GT. 0 ) THEN
             IF ( POUT )  THEN
                print*,' '
                WRITE(*,'(A54)')'                 Convergence reached within tolerance:'
                WRITE(*,'(A22,E8.1,A23)')'Tol=',Tol,' au .Aborting scf loop.'
                print*,' '
                IF ( ETOT .GT. -1.0E03) WRITE(*,'(A33,E27.20,A3)'),' Hartree-Fock energy:   E = ',ETOT,' au'
                IF ( ETOT .LT. -1.0E03) WRITE(*,'(A33,E30.20,A3)'),' Hartree-Fock energy:   E = ',ETOT,' au'
                print*,' '
             ENDIF
     ELSE
             print*,'---------------------------------------------------------'
             print*,'CALCULATION FAILED TO CONVERGE WITHIN ,',Tol,'au'
             print*,'           UQUANTCHEM ABORTS!                    '
             print*,'---------------------------------------------------------'
             STOP
     ENDIF
     OPEN(22,FILE='RHFEIGENVALUES.dat',ACTION='WRITE')
     DO I=1,NB
        WRITE(22,*)I,EHFeigen(I)
     ENDDO
     CLOSE(22)
END SUBROUTINE RHF
