SUBROUTINE TDFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
           & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,FIXNSCF,numprocessors,id,POUT,SCRATCH,ZEROSCF,ETEMP,mu,ENTROPY,OMEGA,EDIR,NEPERIOD,EPROFILE)
      ! This subroutine calculates the self consistent UN-Restricted
      ! Hartree-Fock solution
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      CHARACTER(LEN=20) :: CORRLEVEL,EPROFILE
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF
      INTEGER, INTENT(IN) :: NEPERIOD
      INTEGER, INTENT(IN) :: NTOTALQUAD,Qstart,Qend,Istart,Iend,NB,NRED,Ne,numprocessors,id,FIXNSCF,NATOMS,LORDER,CGORDER,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend),EDIR
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) ::  S(NB,NB),H0(NB,NB),Intsv(Istart:Iend),Tol,nucE,MIX,LQ(LORDER,3),CGQ(CGORDER,2),ETEMP,OMEGA
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigenup(NB),EHFeigendown(NB),ETOT,ENTROPY
      COMPLEX*16, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB)
      DOUBLE PRECISION, INTENT(INOUT) :: mu
      COMPLEX*16 :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),C3(NB,NB),C4(NB,NB)
      DOUBLE PRECISION :: DE,EOLD,DELTAP,Pupr(NB,NB),Pdownr(NB,NB)
      COMPLEX*16 :: Fup(NB,NB),Fdown(NB,NB),G(NB,NB),C1(NB,NB),C2(NB,NB),PTnomix(NB,NB),PTold(NB,NB)
      DOUBLE PRECISION ::  LAMDAu,LAMDAd,TOLDNe,FTOT,FOLD
      COMPLEX*16 :: Pups(50,NB,NB),Pdowns(50,NB,NB),Pupold(NB,NB),Pdownold(NB,NB),Pupt(NB,NB),Pdownt(NB,NB)
      DOUBLE PRECISION :: SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB)
      COMPLEX*16 :: ERRSU(50,NB,NB),ERRSD(50,NB,NB),ERRU(NB,NB),ERRD(NB,NB),Fupo(NB,NB),Fdowno(NB,NB),IM,RE
      DOUBLE PRECISION :: Vxc(2,NB,NB),Nequad,EIG(NB,NB),CDB(NB,NB),MIXING
      INTEGER :: III,II,I,L,M,N,Neup,Nedown,ierr,INFO1,INFO2,RESET
      INTEGER :: MAXITER
      DOUBLE PRECISION, EXTERNAL :: exc,quadcheck
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO
        
      RE = (1.0d0,0.0d0)
      IM = (0.0d0,1.0d0)
      Fup =   (0.0d0,0.0d0)
      Fdown = (0.0d0,0.0d0)
      Fupo = (0.0d0,0.0d0)
      Fdowno = (0.0d0,0.0d0)
      Pupold = (0.0d0,0.0d0)
      Pdownold = (0.0d0,0.0d0)
      
      MIXING = MIX
      IF ( FIXNSCF .GT. 0 ) THEN
              MAXITER = FIXNSCF-1
      ELSE
              MAXITER = 400
      ENDIF
      RESET = 0
      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIG,INFO1)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
      SH = MATMUL(EIG,MATMUL(SL,TRANSPOSE(EIG)))

      IF ( DIISORD .GT. 25 ) THEN
              IF ( id .EQ. 0 ) THEN
                      WRITE(*,*)'The order of DIIS mixing set =',DIISORD,' Is too Big.'
                      WRITE(*,*)'Maximum allowed value is 25. This value has been used instead.'
              ENDIF
              DIISORD = 25
      ENDIF
      
      IF ( SCRATCH ) THEN
              CALL diaghHF( H0,S,NB,EHFeigenup,CDB,INFO1)
                C2 = CDB*RE
                C1 = C2
                Pup = (0.0d0,0.0d0)
                Pdown = (0.0d0,0.0d0)
      ELSE
                C1 = Cup
                C2 = Cdown
      ENDIF
        
      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2

      Pups = (0.0d0,0.0d0)
      Pdowns = (0.0d0,0.0d0)
      ERRSU = (0.0d0,0.0d0)
      ERRSD = (0.0d0,0.0d0)

      DE = 2.0d0*Tol
      DELTAP = 2.0d0*Tol
      
      I = 0
      II = 0
      PTold = (0.0d0,0.0d0)
      STARTPRINTDIISIFO = .FALSE.
      ! The tolerance used when calculating 
      ! the chemical potential 
      TOLDNe = 1.0E-8
      
      !=======================================================
      ! This is used for the zero scf cycle (FAST-MD)
      !=======================================================
      IF ( ZEROSCF ) THEN
                PT = Pup + Pdown
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL getJvc(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jup)
                CALL getJvc(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jdown)
                IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                    CALL getKvc(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kup)
                    CALL getKvc(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kdown)
                ENDIF
                Pupr = DBLE(Pup)
                Pdownr = DBLE(Pdown)
                CALL getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                        Fup   = H0 + Vxc(1,:,:) + Jdown - 0.20d0*Kup   + Jup
                        Fdown = H0 + Vxc(2,:,:) + Jdown - 0.20d0*Kdown + Jup
                ELSE
                        Fup   = H0 + Vxc(1,:,:) + Jdown + Jup
                        Fdown = H0 + Vxc(2,:,:) + Jdown + Jup
                ENDIF                
                
                ! Transform to orthogonal representation:
                Fupo   = MATMUL(SH,MATMUL(Fup,SH))
                Fdowno = MATMUL(SH,MATMUL(Fdown,SH))

                CALL diaghc( Fupo,NB,EHFeigenup,C1,INFO2)
                CALL diaghc( Fdowno,NB,EHFeigendown,C2,INFO2)

                ! Transform to non-orthogonal representation:
                C1 = MATMUL(SH,C1)
                C2 = MATMUL(SH,C2)

                Fup   = Fup   - Vxc(1,:,:)
                Fdown = Fdown - Vxc(2,:,:)
                  
                !ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + & 
                !& exc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id) + nucE
                IF ( ETEMP .LT. 0.0d0 ) THEN
                        Cup(:,:) = (0.0d0,0.0d0)
                        DO M=1,Neup
                                Cup(:,M) = C1(:,M)
                        ENDDO
                        Cdown(:,:) = (0.0d0,0.0d0)
                        DO M=1,Nedown
                                Cdown(:,M) = C2(:,M)
                        ENDDO
                        CALL makedensc(Cup,NB,Pup)
                        CALL makedensc(Cdown,NB,Pdown)
                ELSE
                        Cup = C1
                        Cdown = C2
                        CALL makedensTc(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                ENDIF
                NSCF = 0
                RETURN
      ENDIF
      !=======================================================

      DO WHILE ( (DABS(DE) .GT. Tol .AND. I .LE. MAXITER) .OR. ( DELTAP .GT.  sqrt(Tol) .AND. I .LE. MAXITER) )
                !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                IF ( ETEMP .LE. 0.0d0 ) THEN
                        Cup(:,:) = (0.0d0,0.0d0)
                        DO M=1,Neup
                                Cup(:,M) = C1(:,M)
                        ENDDO

                        Cdown(:,:) = (0.0d0,0.0d0)
                        DO M=1,Nedown
                                Cdown(:,M) = C2(:,M)
                        ENDDO
                ELSE
                        Cup = C1
                        Cdown = C2
                ENDIF
                !==============================================================
                ! In the case of XL-BOMD, it is crusial that the scf starts 
                ! with the density matrices Pup and Pdown provided by the input
                ! since it is these matrices that have been propagated by the 
                ! XL-BOMD scheme and not the coefficient matrices Cup and Cdown
                ! Therefor we have the if-construct here.
                !==============================================================
                IF ( I .GT. 0 ) THEN
                        IF ( ETEMP .LE. 0.0d0 ) THEN
                                CALL makedensc(Cup,NB,Pup)
                                CALL makedensc(Cdown,NB,Pdown)
                        ELSE
                                CALL makedensTc(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                        ENDIF
                ENDIF
  
                PT = Pup + Pdown
               
                ! Calculating the change of the total density matrix:
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0  ) DELTAP = sqrt(DOT_PRODUCT(reshape(PT-PTold,(/NB**2/)),reshape(PT-PTold,(/NB**2/)) ))
                
                ! Linear mixing:
                IF ( I .GT. 0 ) THEN
                   Pup = Pup*(1.0d0-MIXING) + Pupold*MIXING
                   Pdown = Pdown*(1.0d0-MIXING) + Pdownold*MIXING
                ELSE
                   Pup = Pup
                   Pdown = Pdown
                ENDIF
                
                Pupold = Pup 
                Pdownold = Pdown 
 
                PTold = Pupold + Pdownold
              
                ! Saving the density matrices from the previous 
                ! iterations, to be used with the DIIS-method.
                IF ( II .LT. 2*DIISORD ) THEN
                   II = II + 1
                   Pups(II,:,:) = Pup
                   Pdowns(II,:,:) = Pdown
                   ERRSU(II,:,:) = ERRU
                   ERRSD(II,:,:) = ERRD
                ELSE
                        ! refreshing the density matrices used by the DIIS
                        IF ( DIISORD .NE. 0 ) THEN
                                DO III=1,II-1
                                        Pups(III,:,:) = Pups(III+1,:,:)
                                        Pdowns(III,:,:) = Pdowns(III+1,:,:)
                                        ERRSU(III,:,:) = ERRSU(III+1,:,:)
                                        ERRSD(III,:,:) = ERRSD(III+1,:,:)
                                ENDDO
                                Pups(II,:,:) = Pup
                                Pdowns(II,:,:)  = Pdown
                                ERRSU(II,:,:) = ERRU
                                ERRSD(II,:,:) = ERRD
                        ENDIF
                ENDIF
                
                ! Here we use the DIIS method ( CHEM. Phys. Lett. 73, 393 (1980) )
                ! in order to estimate the self-consistent density matrices:
                IF ( II .GE. 4 ) THEN
                
                    CALL DIISc(NB,II,ERRSU,Pups,Pupt,LAMDAu,INFO1)
                    CALL DIISc(NB,II,ERRSD,Pdowns,Pdownt,LAMDAd,INFO2)

                    IF ( (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. LAMDAu+LAMDAd .LT. 1.00d0) .OR. (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. I .GT. DIISSTART )  ) THEN
                            IF ( DBLE(SUM(Pupt*S)) .GT. 0.750d0*Neup .AND. DBLE(SUM(Pdownt*S)) .GT. 0.750d0*Nedown ) THEN
                                Pup = Pupt
                                Pdown = Pdownt
                                STARTPRINTDIISIFO = .TRUE.
                                MIXING = 0.0d0
                            ENDIF
                    ENDIF
                ENDIF

                !================================================================================
                ! Here we make sure that the density matrices are consistent over all the threads
                ! the "golden standard" is here set by the master thread.
                !================================================================================
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_BCAST(Pup,NB*NB,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(Pdown,NB*NB,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !=================================================================================
10 CONTINUE
                PT = Pup + Pdown
                
                CALL getJvc(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jup)
                
                CALL getJvc(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jdown)
                 
                IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                   CALL getKvc(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kup)
                   CALL getKvc(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kdown)
                ENDIF

                Pupr = DBLE(Pup)
                Pdownr = DBLE(Pdown)
                CALL getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)

                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

                ! First  iteration is done with Hartree-Fock
                IF ( I .LT. 1 .AND. SCRATCH ) THEN
                        Fup   = H0 +  Jdown + Jup - Kup
                        Fdown = H0 +  Jdown + Jup - Kdown
                ELSE
                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                Fup   = H0 +  Jdown + Jup  + Vxc(1,:,:) - 0.20d0*Kup
                                Fdown = H0 +  Jdown + Jup  + Vxc(2,:,:) - 0.20d0*Kdown
                        ELSE
                                Fup   = H0 +  Jdown + Jup + Vxc(1,:,:)
                                Fdown = H0 +  Jdown + Jup + Vxc(2,:,:)
                        ENDIF
                ENDIF

                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERRU = MATMUL(Fup,MATMUL(Pup,S)) - MATMUL(S,MATMUL(Pup,Fup))
                        ERRD = MATMUL(Fdown,MATMUL(Pdown,S)) - MATMUL(S,MATMUL(Pdown,Fdown))
                        ERRU = MATMUL(TRANSPOSE(SH),MATMUL(ERRU,SH))
                        ERRD = MATMUL(TRANSPOSE(SH),MATMUL(ERRD,SH))
                ENDIF
                
                ! Transform to orthogonal representation:
                Fupo = MATMUL(SH,MATMUL(Fup,SH))
                Fdowno = MATMUL(SH,MATMUL(Fdown,SH))
                !print*,'13'
                CALL diaghc( Fupo,NB,EHFeigenup,C1,INFO2)
                CALL diaghc( Fdowno,NB,EHFeigendown,C2,INFO2)

                ! Transform to non-orthogonal representation:
                C1 = MATMUL(SH,C1)
                C2 = MATMUL(SH,C2)


                IF ( (INFO1 .NE. 0 .OR. INFO2 .NE. 0 ) .AND. RESET .EQ. 0 ) THEN
                        Pup   =  Pupold 
                        Pdown = Pdownold
                        RESET = 1
                        GOTO 10 
                ENDIF
                RESET = 0

                IF ( I .LT. 1 .AND. SCRATCH ) THEN
                        ETOT = 0.50d0*DBLE(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + nucE
                ELSE
                        Fup   = Fup   - Vxc(1,:,:)
                        Fdown = Fdown - Vxc(2,:,:)
                        IF ( FIXNSCF .LT. 0 ) THEN
                		Pupr = DBLE(Pup)
                		Pdownr = DBLE(Pdown)
                                ETOT = 0.50d0*DBLE(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + &
                                & exc(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id) + nucE
                                IF ( ETEMP .GT. 0.0d0  ) FTOT = ETOT - ETEMP*ENTROPY
                        ENDIF
                ENDIF 

                IF ( I .EQ. 0 .AND. id .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'                 Entering the scf DFT loop                   '
                        print*,'   =========================================================='
                        print*,' '
                        IF ( ETEMP .LT. 0.0d0  ) WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                        IF ( ETEMP .GT. 0.0d0  ) WRITE(*,'(A4,A20,A29,A31,A27,A31)')'N','E [au]','F [au]','DF [au]',' DP','DIIS'
                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .LT. 0.0d0 ) DE = ETOT-EOLD
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .GT. 0.0d0 ) DE = FTOT-FOLD
        
                IF ( id .EQ. 0 .AND. POUT ) THEN 
                   IF ( STARTPRINTDIISIFO ) THEN
                           IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP,LAMDAd+LAMDAu
                           IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDAd+LAMDAu
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

     NSCF = I
     
     Cup = C1
     Cdown = C2
    
     IF ( ETEMP .LT. 0.0d0 ) THEN
        C3 = (0.0d0,0.0d0)
        DO M=1,Neup
                C3(:,M) = C1(:,M)
        ENDDO
        C4 = (0.0d0,0.0d0)
        DO M=1,Nedown
                C4(:,M) = C2(:,M)
        ENDDO
        CALL makedensc(C3,NB,Pup)
        CALL makedensc(C4,NB,Pdown)
     ELSE
        CALL makedensTc(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
     ENDIF
     Pupr = DBLE(Pup)
     Pdownr = DBLE(Pdown)

     IF ( I .LE. MAXITER  .OR. FIXNSCF .GT. 0  ) THEN
           IF( POUT ) THEN
                   Nequad =  quadcheck(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id)
           ENDIF
           IF ( id .EQ. 0 .AND. POUT ) THEN
             print*,' '
             WRITE(*,'(A54)')'                 Convergence reached within tolerance:'
             WRITE(*,'(A22,E8.1,A23)')'Tol=',Tol,' au .Aborting scf loop.'
             print*,' '
             IF ( CORRLEVEL .NE. 'B3LYP' ) THEN
                        IF ( ETOT .GT. -1.0E03) WRITE(*,'(A6,A3,A20,E27.20,A3)'),'      ',CORRLEVEL,' energy:        E = ',ETOT,' au'
                        IF ( ETOT .LT. -1.0E03) WRITE(*,'(A6,A3,A20,E30.20,A3)'),'      ',CORRLEVEL,' energy:        E = ',ETOT,' au'
                        IF ( FTOT .GT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A6,A3,A20,E27.20,A3)'),'      ',CORRLEVEL,' free-energy:   F = ',FTOT,' au'
                        IF ( FTOT .LT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A6,A3,A20,E30.20,A3)'),'      ',CORRLEVEL,' free-energy:   F = ',FTOT,' au'
                        
             ELSE
                        IF ( ETOT .GT. -1.0E03) WRITE(*,'(A6,A5,A20,E27.20,A3)'),'      ',CORRLEVEL,' energy:        E = ',ETOT,' au'
                        IF ( ETOT .LT. -1.0E03) WRITE(*,'(A6,A5,A20,E30.20,A3)'),'      ',CORRLEVEL,' energy:        E = ',ETOT,' au'
                        IF ( FTOT .GT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A6,A5,A20,E27.20,A3)'),'      ',CORRLEVEL,' free-energy:   F = ',FTOT,' au'
                        IF ( FTOT .LT. -1.0E03 .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(A6,A5,A20,E30.20,A3)'),'      ',CORRLEVEL,' free-energy:   F = ',FTOT,' au'  
             ENDIF 
             print*,' '
             WRITE(*,'(A75,F9.6)')'     The exact number of electrons calculated from the trace of P*S, Ne =  ',SUM(PT*S)
             WRITE(*,'(A75,F9.6)' )'Number of electrons calculated from integrating the charge-density, Ne = ',Nequad
                                                                                                               
           ENDIF
     ELSE
           IF ( id .EQ. 0 ) THEN
             print*,'---------------------------------------------------------'
             print*,'CALCULATION FAILED TO CONVERGE WITHIN ,',Tol,'au'
             print*,'           UQUANTCHEM ABORTS!                    '
             print*,'---------------------------------------------------------'
           ENDIF
           STOP
     ENDIF
     IF ( id .EQ. 0 .AND. POUT ) THEN
             OPEN(21,FILE='UHFEIGENVALUES.dat',ACTION='WRITE')
             DO I=1,NB
                WRITE(21,*)I,EHFeigenup(I),EHFeigendown(I)
             ENDDO
             CLOSE(21)
     ENDIF
     IF ( ETEMP .GT. 0.0d0 ) ETOT = FTOT
END SUBROUTINE TDFT
