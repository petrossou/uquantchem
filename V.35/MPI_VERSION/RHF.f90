SUBROUTINE RHF(S,H0,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol,EHFeigen,ETOT,C,P,MIX,DIISORD, & 
      & DIISSTART,NSCF,FIXNSCF,numprocessors,id,POUT,SCRATCH,ZEROSCF,SCALARRELC,IntsDirac,ETEMP,ENTROPY)
      ! This subroutine calculates the self consistent Restricted
      ! Hartree-Fock solution
      IMPLICIT NONE
      INCLUDE "mpif.h"
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF,SCALARRELC
      INTEGER, INTENT(IN) :: NB,Ne,numprocessors,id,FIXNSCF
      INTEGER*8, INTENT(IN) :: NRED,Istart,Iend
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),H0(NB,NB),Intsv(Istart:Iend),Tol,nucE,MIX, IntsDirac(Istart:Iend),ENTROPY,ETEMP
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigen(NB),ETOT
      DOUBLE PRECISION, INTENT(INOUT) :: C(NB,NB),P(NB,NB)
      DOUBLE PRECISION :: J(NB,NB),K(NB,NB),JD(NB,NB),F(NB,NB),G(NB,NB),C2(NB,NB),DE,EOLD,DELTAP,P1(NB,NB),P2(NB,NB)
      DOUBLE PRECISION :: Pold(NB,NB),Ps(50,NB,NB),Pt(NB,NB),LAMDA,MIXING,C3(NB,NB),mu,FTOT,FOLD
      DOUBLE PRECISION :: ERRS(50,NB,NB),ERR(NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB),TOLDNe
      INTEGER :: I,II,III,L,M,N,ierr,INFO
      INTEGER :: MAXITER
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO, CONTINUEWITHSREL
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0      
      DOUBLE PRECISION, PARAMETER :: lightspeed = 137.035999084

      MIXING = MIX
      CONTINUEWITHSREL = .FALSE. 
      ! The tolerance used when calculating
      ! the chemical potential
      TOLDNe = 1.0E-8

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
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL getJv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,J)
                CALL getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,K)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                F = H0 + J - 0.50d0*K
                G = J - 0.50d0*K
                !IF ( SCALARRELC ) THEN
                !        ! In the case of scalar relatevistic corrections. The nuclear attraction Dirac correction 
                !        ! term is allready included in H0 together with the mass-correction term. Here we add the 
                !        ! Dirac correction term emanating from the electro-electron (Hartree) repulsion term.
                !        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !        CALL getJv(P,NB,NRED,Istart,Iend,IntsDirac,IND1,IND2,IND3,IND4,numprocessors,id,JD)
                !        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !
                !        F   = F + JD*(pi/(2.0*(lightspeed**2)))
                !        G   = G + JD*(pi/(2.0*(lightspeed**2)))
                !
                !ENDIF
                CALL diaghHF( F,S,NB,EHFeigen,C,INFO)
                !ETOT = SUM(F*P)-0.5*SUM(G*P) + nucE
                C2(:,:) = 0.0d0
                DO M=1,N
                        C2(:,M) = C(:,M)
                ENDDO
                IF ( ETEMP .LT. 0 ) THEN
                        CALL makedens(C2,NB,P)
                        P = 2.0d0*P
                ELSE
                        CALL makedensT(TOLDNe,C2,C2,EHFeigen,EHFeigen,ETEMP,NB,Ne,P1,P2,mu,ENTROPY)
                        P = P1 + P2
                ENDIF

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
                IF ( I .GT. 0 ) THEN
                        IF ( ETEMP .LT. 0 ) THEN
                                CALL makedens(C2,NB,P)
                                P = 2.0d0*P
                        ELSE
                                CALL makedensT(TOLDNe,C2,C2,EHFeigen,EHFeigen,ETEMP,NB,Ne,P1,P2,mu,ENTROPY)
                                P = P1 + P2
                        ENDIF
                ENDIF
                        
                ! Calculating the change in the density matrix
                DELTAP = sqrt(DOT_PRODUCT(reshape(P-Pold,(/NB**2/)),reshape(P-Pold,(/NB**2/)) ))

                ! Linear mixing step
                IF ( I .GT. 0 ) THEN 
                    P = P*(1.0d0 - MIXING) + POLD*MIXING
                ELSE
                    P = P
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
                                 IF ( I .GE. DIISSTART ) THEN
                                        P = Pt
                                        STARTPRINTDIISIFO = .TRUE.
                                        MIXING = 0.0d0
                                 ENDIF
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

                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL getJv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,J)
                CALL getKv(P,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,K)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                
                F = H0 + J - 0.50d0*K
                
                G = J - 0.50d0*K
                !IF ( SCALARRELC ) THEN
                !        ! In the case of scalar relatevistic corrections. The nuclear attraction Dirac correction 
                !        ! term is allready included in H0 together with the mass-correction term. Here we add the 
                !        ! Dirac correction term emanating from the electro-electron (Hartree) repulsion term.
                !        ! This screened Dirac term is unstable thus it is only calculated at the first step.
                !        ! This also demand that the Non-rel and convered Pup and Pdown are provided at first iteration.
                !        ! This term can only be considered as a Perturbative term and not part of the SCF-cycle.
                !
                !        IF ( I .EQ. 0 ) THEN
                !                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !                CALL getJv(P,NB,NRED,Istart,Iend,IntsDirac,IND1,IND2,IND3,IND4,numprocessors,id,JD)
                !                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !        ENDIF
                !
                !       F   = F + JD*(pi/(2.0*(lightspeed**2)))
                !       G   = G + JD*(pi/(2.0*(lightspeed**2)))
                !       CONTINUEWITHSREL = .TRUE.
                !
                !ENDIF
                
                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERR = MATMUL(F,MATMUL(P,S)) - MATMUL(S,MATMUL(P,F))
                        ERR = MATMUL(TRANSPOSE(SH),MATMUL(ERR,SH))
                ENDIF

                CALL diaghHF( F,S,NB,EHFeigen,C,INFO)

                IF ( FIXNSCF .LT. 0 ) THEN
                        ETOT = SUM(F*P)-0.5*SUM(G*P) + nucE
                        IF ( ETEMP .GT. 0.0d0  ) FTOT = ETOT - ETEMP*ENTROPY
                ENDIF

                IF ( I .EQ. 0 .AND. id .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'         Entering the scf restricted Hartree-Fock loop        '
                        print*,'   =========================================================='
                        print*,' '
                        IF ( ETEMP .LT. 0.0d0  ) WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                        IF ( ETEMP .GT. 0.0d0  ) WRITE(*,'(A4,A20,A29,A31,A27,A31)')'N','E [au]','F [au]','DF [au]',' DP','DIIS'

                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .LT. 0.0d0 ) DE = ETOT-EOLD
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .GT. 0.0d0 ) DE = FTOT-FOLD


                IF ( id .EQ. 0 .AND. POUT ) THEN
                        IF ( STARTPRINTDIISIFO ) THEN
                                IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP,LAMDA
                                IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDA

                        ELSE
                                IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP
                                IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP

                        ENDIF
                ENDIF
        
                EOLD = ETOT
                FOLD = FTOT
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

     IF ( I .LE. MAXITER  .OR. FIXNSCF .GT. 0  ) THEN
            IF ( id .EQ. 0 .AND. POUT ) THEN
             print*,' '
             WRITE(*,'(A54)')'                 Convergence reached within tolerance:'
             WRITE(*,'(A22,E8.1,A23)')'Tol=',Tol,' au .Aborting scf loop.'
             print*,' '
             IF ( ETOT .GT. -1.0E03) WRITE(*,'(A33,E27.20,A3)'),'      Hartree-Fock energy:   E = ',ETOT,' au'
             IF ( ETOT .LT. -1.0E03) WRITE(*,'(A33,E30.20,A3)'),'      Hartree-Fock energy:   E = ',ETOT,' au'
             IF ( FTOT .GT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A7,A6,A20,E27.20,A3)'),'       ','  RHF ',' free-energy: F = ',FTOT,' au'
             IF ( FTOT .LT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A7,A6,A20,E30.20,A3)'),'       ','  RHF ',' free-energy: F = ',FTOT,' au'
             print*,' '
             IF ( mu .GT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A33,E27.20,A3)'),'           Fermi Energy:   E_F = ',mu,' au'
             IF ( mu .LT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A33,E30.20,A3)'),'           Fermi Energy:   E_F = ',mu,' au'
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
END SUBROUTINE RHF
