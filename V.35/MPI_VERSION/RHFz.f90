SUBROUTINE RHFz(S,H0,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,C,P,MIX,DIISORD,DIISSTART,NSCF,FIXNSCF,numprocessors,id,POUT,SCRATCH,ZEROSCF,Jup,Jdown,Kup,Kdown)
      ! This subroutine calculates the self consistent Restricted
      ! Hartree-Fock solution
      IMPLICIT NONE
      INCLUDE "mpif.h"
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF
      INTEGER, INTENT(IN) :: NB,NRED,Istart,Iend,Ne,numprocessors,id,FIXNSCF
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),H0(NB,NB),Intsv(Istart:Iend),Tol,nucE,MIX
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigen(NB),ETOT
      DOUBLE PRECISION, INTENT(INOUT) :: C(NB,NB),P(NB,NB)
      DOUBLE PRECISION, INTENT(OUT) :: Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB)
      DOUBLE PRECISION :: J(NB,NB),K(NB,NB),F(NB,NB),G(NB,NB),C2(NB,NB),DE,EOLD,DELTAP
      DOUBLE PRECISION :: Pold(NB,NB),Ps(50,NB,NB),Pt(NB,NB),LAMDA,MIXING,C3(NB,NB)
      DOUBLE PRECISION :: ERRS(50,NB,NB),ERR(NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB)
      INTEGER :: I,II,III,L,M,N,ierr,INFO
      INTEGER :: MAXITER
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO

      MIXING = MIX
      
      IF ( DIISORD .GT. 25 ) THEN
              IF ( id .EQ. 0 ) THEN
                      WRITE(*,*)'The order of DIIS mixing set =',DIISORD,' Is too Big.'
                      WRITE(*,*)'Maximum allowed value is 25. This value has been used instead.'
              ENDIF
              DIISORD = 25
      ENDIF
        
      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIGENVECT,INFO)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
      SH = MATMUL(EIGENVECT,MATMUL(SL,TRANSPOSE(EIGENVECT)))

      IF ( SCRATCH ) CALL diaghHF( H0,S,NB,EHFeigen,C,INFO)
        
      IF ( MOD(Ne,2) .EQ. 0 ) THEN
                N = Ne/2
      ELSE
                IF ( id .EQ. 0 ) THEN
                 print*,'----------------------------------------------------------------------------'
                 print*,'Attempting to run Restricted Hartree-Fock for open-shelled system. ABORTING!'
                 print*,'              Change input parameter CORRLEVEL to URHF                      '
                 print*,'----------------------------------------------------------------------------'
                ENDIF
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
                CALL getJv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,J)
                CALL getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,K)
                F = H0 + J - 0.50d0*K
                G = J - 0.50d0*K
                CALL diaghHF( F,S,NB,EHFeigen,C,INFO)
                !ETOT = SUM(F*P)-0.5*SUM(G*P) + nucE
                C2(:,:) = 0.0d0
                DO M=1,N
                        C2(:,M) = C(:,M)
                ENDDO
                CALL makedens(C2,NB,P)
                P = 2.0d0*P
                NSCF = 0
                Jup   = J/2.0d0
                Jdown = Jup
                Kup   = K/2.0d0
                Kdown = Kup
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
                DELTAP = sqrt(DOT_PRODUCT(reshape(P-Pold,(/NB**2/)),reshape(P-Pold,(/NB**2/)) ))

                ! Linear mixing step
                IF ( I .GT. 0 ) THEN 
                    P = 2*P*(1.0d0 - MIXING) + POLD*MIXING
                ELSE
                    P = 2*P
                ENDIF

                POLD = P 
                
                ! Saving the density matrices from the previous
                ! iterations, to be used with the DIIS-method.
                IF ( II .LT. 2*DIISORD ) THEN
                         II = II + 1
                         Ps(II,:,:) = P
                         ERRS(II,:,:) = ERR
                 ELSE
                        IF ( DIISORD .NE. 0 ) THEN
                                ! refreshing the density matrices used by the DIIS
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
                                 P = Pt
                                 STARTPRINTDIISIFO = .TRUE.
                                 MIXING = 0.0d0
                         ENDIF
                 ENDIF
                
                !================================================================================
                ! Here we make sure that the density matrix is consistent over all the threads
                ! the "golden standard" is here set by the master thread.
                !================================================================================
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_BCAST(P,NB*NB,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !================================================================================

                CALL getJv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,J)
                CALL getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,K)
               

                F = H0 + J - 0.50d0*K
                
                G = J - 0.50d0*K
                
                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERR = MATMUL(F,MATMUL(P,S)) - MATMUL(S,MATMUL(P,F))
                        ERR = MATMUL(TRANSPOSE(SH),MATMUL(ERR,SH))
                ENDIF

                CALL diaghHF( F,S,NB,EHFeigen,C,INFO)

                IF ( FIXNSCF .LT. 0 ) ETOT = SUM(F*P)-0.5*SUM(G*P) + nucE

                IF ( I .EQ. 0 .AND. id .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'         Entering the scf restricted Hartree-Fock loop        '
                        print*,'   =========================================================='
                        print*,' '
                        WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LT. 0  ) DE = ETOT-EOLD
        
                IF ( id .EQ. 0 .AND. POUT ) THEN
                        IF ( STARTPRINTDIISIFO ) THEN
                                WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDA
                        ELSE
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP
                        ENDIF
                ENDIF
        
                EOLD = ETOT
                I = I+1
                !================================================================================
                ! Here we make sure that the energy change and density change
                ! are consistent over all threads, so that to avoid deadlock.
                ! the "golden standard" is here set by the master thread.
                !================================================================================
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_BCAST(DE,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(DELTAP,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !=================================================================================

     ENDDO

     C3 = 0.0d0
     DO M=1,N
        C3(:,M) = C(:,M)
     ENDDO
     CALL makedens(C3,NB,P)
     P = P*2.0d0
     NSCF = I
     Jup   = J/2.0d0
     Jdown = Jup
     Kup   = K/2.0d0
     Kdown = Kup
     
     IF ( I .LE. MAXITER  .OR. FIXNSCF .GT. 0  ) THEN
            IF ( id .EQ. 0 .AND. POUT ) THEN
             print*,' '
             WRITE(*,'(A54)')'                 Convergence reached within tolerance:'
             WRITE(*,'(A22,E8.1,A23)')'Tol=',Tol,' au .Aborting scf loop.'
             print*,' '
             IF ( ETOT .GT. -1.0E03) WRITE(*,'(A32,E27.20,A3)'),' Hartree-Fock energy:   E =',ETOT,' au'
             IF ( ETOT .LT. -1.0E03) WRITE(*,'(A32,E30.20,A3)'),' Hartree-Fock energy:   E =',ETOT,' au'
             print*,' '
           ENDIF
     ELSE
            IF (id .EQ. 0 ) THEN
             print*,'---------------------------------------------------------'
             print*,'CALCULATION FAILED TO CONVERGE WITHIN ,',Tol,'au'
             print*,'           UQUANTCHEM ABORTS!                    '
             print*,'---------------------------------------------------------'
            ENDIF
             STOP
     ENDIF
     IF (id .EQ. 0 .AND. POUT ) THEN
             OPEN(22,FILE='RHFEIGENVALUES.dat',ACTION='WRITE')
             DO I=1,NB
                  WRITE(22,*)I,EHFeigen(I)
             ENDDO
             CLOSE(22)
     ENDIF
END SUBROUTINE RHFz
