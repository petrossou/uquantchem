SUBROUTINE TDFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
& ETOT,Cup,Cdown,Pup,Pdown,numprocessors,id,POUT,mu,OMEGA,EDIR,NEPERIOD,EPROFILE,NTIMESTEPS,TIMESTEP,EFIELDMAX,IORBNR,PEXu,PEXd,PEXuu,PEXdd,DOABSSPECTRUM,DIFFDENS,NSCCORR,MIXTDDFT,SCERR)
      ! This subroutine calculates the self consistent UN-Restricted
      ! Hartree-Fock solution
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      CHARACTER(LEN=20) :: CORRLEVEL,EPROFILE
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      LOGICAL, INTENT(IN) :: POUT,DOABSSPECTRUM,DIFFDENS
      INTEGER, INTENT(IN) :: NEPERIOD,NTIMESTEPS,IORBNR(2),NSCCORR
      INTEGER, INTENT(IN) :: NTOTALQUAD,Qstart,Qend,NB,Ne,numprocessors,id,NATOMS,LORDER,CGORDER,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend),EDIR
      INTEGER*8, INTENT(IN) :: Istart,Iend,NRED
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(IN) ::  S(NB,NB),H0(NB,NB),Intsv(Istart:Iend),Tol,nucE,LQ(LORDER,3),CGQ(CGORDER,2),OMEGA,EFIELDMAX,TIMESTEP,MIXTDDFT,SCERR
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigenup(NB),EHFeigendown(NB),ETOT,PEXu(24,NB,NB),PEXd(24,NB,NB),PEXuu(24,NB,NB),PEXdd(24,NB,NB)
      COMPLEX*16, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB)
      COMPLEX*16 ::  Cupo(NB,NB),Cdowno(NB,NB)
      DOUBLE PRECISION, INTENT(INOUT) :: mu
      COMPLEX*16 :: Puu(2,NB,NB),Pdd(2,NB,NB),Expu(NB,NB),Expd(NB,NB),Pupno(NB,NB),Pdownno(NB,NB),C0U(NB,NB),C0D(NB,NB)
      COMPLEX*16 :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),C3(NB,NB),C4(NB,NB),F3(NB,NB),F4(NB,NB),PHOLE(NB,NB)
      DOUBLE PRECISION :: DE,EOLD,DELTAP,Pupr(NB,NB),Pdownr(NB,NB),PTEMP(NB,NB),CONV
      COMPLEX*16 :: Fup(NB,NB),Fdown(NB,NB),G(NB,NB),C1(NB,NB),C2(NB,NB),PTnomix(NB,NB),PTold(NB,NB)
      DOUBLE PRECISION ::  LAMDAu,LAMDAd,TOLDNe,FTOT,FOLD
      COMPLEX*16 :: Pups(50,NB,NB),Pdowns(50,NB,NB),Pupold(NB,NB),Pdownold(NB,NB),Pupt(NB,NB),Pdownt(NB,NB),Fups(2,NB,NB),Fdowns(2,NB,NB)
      DOUBLE PRECISION :: SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB),DTENS(NB,NB),SHP(NB,NB),Csu(NB,NB),Csd(NB,NB)
      COMPLEX*16 :: ERRSU(50,NB,NB),ERRSD(50,NB,NB),ERRU(NB,NB),ERRD(NB,NB),Fupo(NB,NB),Fdowno(NB,NB),IM,RE,DIPOLE(NB,NB),SHH(NB,NB),SHHP(NB,NB)
      DOUBLE PRECISION :: Vxc(2,NB,NB),Nequad,EIG(NB,NB),CDB(NB,NB),MIXING,DPTENSOR(3,NB,NB),DIPOLET(NB,NB),OCCU(NB),OCCD(NB),DT,EFIELD,TIME,MOMENT(NTIMESTEPS),EFIELDP
      INTEGER :: III,II,I,L,M,N,Neup,Nedown,ierr,INFO1,INFO2,RESET,ISAVE
      INTEGER :: MAXITER,ISTARTA,Ncorr
      DOUBLE PRECISION, EXTERNAL :: exc,quadcheck,trace
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: STARTPRINTDIISIFO,MIDPOINT
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

      MIDPOINT = .FALSE. 
      ISAVE = 0
      RE = (1.0d0,0.0d0)
      IM = (0.0d0,1.0d0)
      Fup =   (0.0d0,0.0d0)
      Fdown = (0.0d0,0.0d0)
      Fupo = (0.0d0,0.0d0)
      Fdowno = (0.0d0,0.0d0)
      Pupold = (0.0d0,0.0d0)
      Pdownold = (0.0d0,0.0d0)
      ISTARTA = 0
      Ncorr = 0

      Cupo = Cup
      Cdowno = Cdown

      ! If the profile of the EM-wave is squre the dipolmatrix elements 
      ! are calculated here
      CALL dipoletensor(BAS,DPTENSOR)
      IF ( EDIR .EQ. 1 ) DIPOLET = DPTENSOR(2,:,:)
      IF ( EDIR .EQ. 2 ) DIPOLET = DPTENSOR(3,:,:)
      IF ( EDIR .EQ. 3 ) DIPOLET = DPTENSOR(1,:,:)
      DTENS = DIPOLET

      IF ( id .EQ. 0 ) THEN
        OPEN(111,FILE='TDFTOUT.dat',ACTION='WRITE')
        OPEN(222,FILE='OCCUPATIONUP.dat',ACTION='WRITE')
        OPEN(333,FILE='OCCUPATIONDOWN.dat',ACTION='WRITE')
      ENDIF

      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIG,INFO1)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
      SH = MATMUL(EIG,MATMUL(SL,TRANSPOSE(EIG)))
      SHH = RE*SH

      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2
      
      ! If we want to run the TDDFT/THF caclulation with 
      ! a hole present at t = 0 and see how it propagates 
      ! for t > 0, with or without field, here is where the 
      ! whole contribution to the density matrix is computed.
      IF ( IORBNR(1) .NE. 0 ) THEN
              C3 = 0.0d0
              IF ( IORBNR(1) .GT. 0 ) C3(:,IORBNR(1)) = Cup(:,IORBNR(1))
              IF ( IORBNR(1) .LT. 0 ) C3(:,-IORBNR(1)) = Cdown(:,-IORBNR(1))
              CALL makedensc(C3,NB,PHOLE)
              IF ( IORBNR(1) .GT. 0 ) Pup   = Pup   - PHOLE
              IF ( IORBNR(1) .LT. 0 ) Pdown = Pdown - PHOLE
              ! Calculating the lowdin S^(1/2) matrix
              SL = 0.0d0
              DO I=1,NB
                        SL(I,I) = sqrt(LAM(I))
              ENDDO
              SHP = MATMUL(EIG,MATMUL(SL,TRANSPOSE(EIG)))
              SHHP = RE*SHP
              ! Calculating the projection vectors (groundstate withouth hole)
              ! in the orthogonal representation.
              C0U = MATMUL(SHHP,Cup)
              C0D = MATMUL(SHHP,Cdown)
      ENDIF

      
      I = 0
      DT = TIMESTEP

      Puu(2,:,:) = Pup
      Pdd(2,:,:) = Pdown
      
      DO WHILE (  I .LT. NTIMESTEPS  )
                !================================================================================
                ! Here we make sure that the density matrices are consistent over all the threads
                ! the "golden standard" is here set by the master thread.
                !================================================================================
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                CALL MPI_BCAST(Pup,NB*NB,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BCAST(Pdown,NB*NB,MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
                CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                !=================================================================================
               
                I = I+1
                TIME = (I-1)*TIMESTEP

                Puu(1,:,:) = Puu(2,:,:)
                Pdd(1,:,:) = Pdd(2,:,:)
                Puu(2,:,:) = Pup
                Pdd(2,:,:) = Pdown

                Fups(1,:,:)   = Fups(2,:,:)
                Fdowns(1,:,:) = Fdowns(2,:,:) 
                Fups(2,:,:)   = Fupo
                Fdowns(2,:,:) = Fdowno


                IF ( I .GT. 1 ) THEN
                        ! ==========================================
                        ! Transform to non-orthogonal representation:
                        !===========================================
                        Pupno   = MATMUL(SHH,MATMUL(Pup,TRANSPOSE(SHH)))
                        Pdownno = MATMUL(SHH,MATMUL(Pdown,TRANSPOSE(SHH)))
                ELSE
                        Pupno   = Pup
                        Pdownno = Pdown
                ENDIF

                Pupno = 0.50d0*( Pupno + TRANSPOSE(CONJG(Pupno)) )
                Pdownno = 0.50d0*( Pdownno + TRANSPOSE(CONJG(Pdownno)) )

                CALL getJvc(Pupno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jup)
                CALL getJvc(Pdownno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jdown)

                Pupr = DBLE(Pupno)
                Pdownr = DBLE(Pdownno)
                 
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                
                        CALL getKvc(Pupno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kup)
                        CALL getKvc(Pdownno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kdown)
             
                        Fup   = H0 +  Jdown + Jup - Kup
                        Fdown = H0 +  Jdown + Jup - Kdown
                ELSE

                        CALL getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)

                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                CALL getKvc(Pupno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kup)
                                CALL getKvc(Pdownno,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kdown)
                                Fup   = H0 +  Jdown + Jup  + Vxc(1,:,:) - 0.20d0*Kup
                                Fdown = H0 +  Jdown + Jup  + Vxc(2,:,:) - 0.20d0*Kdown
                        ELSE
                                Fup   = H0 +  Jdown + Jup + Vxc(1,:,:)
                                Fdown = H0 +  Jdown + Jup + Vxc(2,:,:)
                        ENDIF
                ENDIF

                !======================================================
                ! Here the amplitude modulation of the E-field is done:
                !======================================================
                DIPOLE = (0.0d0,0.0d0)
                IF ( EPROFILE .EQ. 'AC' ) CALL dipoletensorinhom(BAS,EDIR,OMEGA,TIME,DIPOLET)

                IF ( TIME .GE. 0.0d0 .AND. TIME .LE. 2.0d0*pi/OMEGA ) THEN
                               EFIELD = EFIELDMAX*(OMEGA*TIME/(2.0d0*pi))
                ENDIF

                IF ( TIME .GT. 2.0d0*pi/OMEGA .AND. TIME .LE. (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                EFIELD = EFIELDMAX
                ENDIF

                IF ( TIME .GT. (NEPERIOD + 1.0d0)*(2.0d0*pi/OMEGA) .AND.  TIME .LT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) THEN
                                EFIELD = EFIELDMAX*( (NEPERIOD + 2.0d0) - OMEGA*TIME/(2.0d0*pi) )
                ENDIF

                IF ( TIME .GT. (NEPERIOD + 2.0d0)*(2.0d0*pi/OMEGA) ) EFIELD = 0.0d0

                IF ( EPROFILE .EQ. 'AC' ) DIPOLE = EFIELD*DIPOLET*RE
                IF ( EPROFILE .EQ. 'AC' ) EFIELD = EFIELD*sin(OMEGA*TIME)*RE
                
                IF ( EPROFILE .EQ. 'HO' ) EFIELD = EFIELD*sin(OMEGA*TIME)*RE
                IF ( EPROFILE .EQ. 'HO' ) DIPOLE = EFIELD*DIPOLET*RE

                IF ( EPROFILE .EQ. 'CIRC' ) THEN
                        IF ( EDIR .EQ. 1 ) THEN
                                DIPOLE = EFIELD*(DPTENSOR(2,:,:)*sin(OMEGA*TIME) + DPTENSOR(3,:,:)*cos(OMEGA*TIME) )*RE
                        ELSE IF ( EDIR .EQ. 2 ) THEN
                                DIPOLE = EFIELD*(DPTENSOR(3,:,:)*sin(OMEGA*TIME) + DPTENSOR(1,:,:)*cos(OMEGA*TIME) )*RE
                        ELSE IF ( EDIR .EQ. 3 ) THEN
                                DIPOLE = EFIELD*(DPTENSOR(1,:,:)*sin(OMEGA*TIME) + DPTENSOR(2,:,:)*cos(OMEGA*TIME) )*RE
                        ENDIF
                ENDIF

                IF ( EPROFILE .EQ. 'DP' ) THEN
                        IF ( I .EQ. 1 ) THEN
                                F3 = Fup
                                F4 = Fdown
                                Fup   = Fup   + EFIELDMAX*DIPOLET*RE
                                Fdown = Fdown + EFIELDMAX*DIPOLET*RE
                                EFIELD = EFIELDMAX
                        ELSE
                                DIPOLE = (0.0d0,0.0d0)
                                EFIELD = 0.0d0
                        ENDIF
                ENDIF


                IF ( I .GT. 1 ) THEN
                        Fup   = Fup + DIPOLE
                        Fdown = Fdown + DIPOLE
                ENDIF
                
                Fup   = 0.50d0*(Fup + TRANSPOSE(CONJG(Fup)) )
                Fdown = 0.50d0*(Fdown + TRANSPOSE(CONJG(Fdown)) )
                
                !=======================================
                ! Transform to orthogonal representation:
                !=======================================
                
                Fupo = MATMUL(SH,MATMUL(Fup,SH))
                Fdowno = MATMUL(SH,MATMUL(Fdown,SH))

                Fupo = 0.50d0*( Fupo + TRANSPOSE(CONJG(Fupo)) )
                Fdowno = 0.50d0*( Fdowno + TRANSPOSE(CONJG(Fdowno)) )
               
                IF (  EPROFILE .EQ. 'DP' .AND. I .EQ. 1 ) THEN
                        F3 = 0.50d0*(F3 + TRANSPOSE(CONJG(F3)) )
                        F4 = 0.50d0*(F4 + TRANSPOSE(CONJG(F4)) )
                        F3 = MATMUL(SH,MATMUL(F3,SH))
                        F4 = MATMUL(SH,MATMUL(F4,SH))
                        F3 = 0.50d0*(F3 + TRANSPOSE(CONJG(F3)) )
                        F4 = 0.50d0*(F4 + TRANSPOSE(CONJG(F4)) )
                        CALL diaghc( F3,NB,EHFeigenup,C3,INFO1)
                        CALL diaghc( F4,NB,EHFeigendown,C4,INFO2)
                ENDIF
 
                IF ( I .EQ. 1 ) THEN 
                        CALL diaghc( Fupo,NB,EHFeigenup,C1,INFO2)
                        CALL diaghc( Fdowno,NB,EHFeigendown,C2,INFO2)
                ELSE
                        IF ( MIDPOINT .AND. I .GT. 2 ) THEN
                                CALL diaghc( 0.50d0*(2.0d0*Fups(2,:,:)+Fupo-Fups(1,:,:)),NB,EHFeigenup,C1,INFO2)
                                CALL diaghc( 0.50d0*(2.0d0*Fdowns(2,:,:)+Fdowno-Fdowns(1,:,:)),NB,EHFeigendown,C2,INFO2)
                        ELSE
                                ! Change made by me Tue Mar  3 22:46:10 EET 2015
                                ! to try and increase the accuracy of the integrand (Fockian) from Ortho(DT) to Ortho(DT**2)
                                IF ( I .GT. 2  ) THEN
                                        CALL diaghc( 0.25d0*(7.0d0*Fupo-4*Fups(2,:,:)+Fups(1,:,:)),NB,EHFeigenup,C1,INFO2)
                                        CALL diaghc( 0.25d0*(7.0d0*Fdowno-4*Fdowns(2,:,:)+Fdowns(1,:,:)),NB,EHFeigendown,C2,INFO2)
                                ELSE
                                        CALL diaghc( 0.50d0*(3.0d0*Fupo-Fups(2,:,:)),NB,EHFeigenup,C1,INFO2)
                                        CALL diaghc( 0.50d0*(3.0d0*Fdowno-Fdowns(2,:,:)),NB,EHFeigendown,C2,INFO2)
                                ENDIF
                        ENDIF
                ENDIF

                IF ( I .EQ. 1 ) THEN
                     !===========================================================================================
                     ! At time = 0 we transform the density matrices to orthogonal form by using the eiegnvectors
                     ! to the orthogonal Fockians Fupo and Fdowno to construct them (When no hole is present)
                     !============================================================================================
                     Neup   = ( Ne - MOD(Ne,2) )/2
                     Nedown = ( Ne + MOD(Ne,2) )/2
                     Cup(:,:) = 0.0d0
                     DO M=1,Neup
                        Cup(:,M) = C1(:,M)
                     ENDDO
                     Cdown(:,:) = 0.0d0
                     DO M=1,Nedown
                        Cdown(:,M) = C2(:,M)
                     ENDDO
                     
                     IF ( IORBNR(1) .EQ. 0 ) THEN
                                ! When no hole is present
                                CALL makedensc(Cup,NB,Pup)
                                CALL makedensc(Cdown,NB,Pdown)
                                IF (  EPROFILE .EQ. 'DP' ) THEN
                                        C0U = C3
                                        C0D = C4
                                ELSE
                                        C0U = C1
                                        C0D = C2
                                ENDIF
                     ELSE
                                ! ============================================================
                                ! Transform to orthogonal representation when hole is present:
                                !=============================================================
                                Pup   = MATMUL(SHHP,MATMUL(Pup,TRANSPOSE(SHHP)))
                                Pdown = MATMUL(SHHP,MATMUL(Pdown,TRANSPOSE(SHHP)))
                     ENDIF

                     Puu(1,:,:) = 0.50d0*(Pup + TRANSPOSE(CONJG(Pup)) )
                     Pdd(1,:,:) = 0.50d0*(Pdown + TRANSPOSE(CONJG(Pdown)) )

                     Puu(2,:,:) = Puu(1,:,:) 
                     Pdd(2,:,:) = Pdd(1,:,:)
                ENDIF

                PT = Pup + Pdown

                !=======================================================================================
                ! Calculating the dipolemoment:
                !=======================================================================================
                MOMENT(I) = 0.0d0
                DO M=1,NATOMS
                        IF ( EDIR .EQ. 1 ) MOMENT(I) = MOMENT(I) + ATOMS(M)%R(2)*ATOMS(M)%Z
                        IF ( EDIR .EQ. 2 ) MOMENT(I) = MOMENT(I) + ATOMS(M)%R(3)*ATOMS(M)%Z
                        IF ( EDIR .EQ. 3 ) MOMENT(I) = MOMENT(I) + ATOMS(M)%R(1)*ATOMS(M)%Z
                ENDDO

                MOMENT(I) = MOMENT(I) - DBLE(SUM((Pupno+Pdownno)*DTENS))
               
                !=============================================================================
                ! Calculating the occupation numbers, see Eqn  (12) J. Chem. Phys. 128, 114113
                !=============================================================================
                DO M=1,NB
                        OCCU(M) = DBLE(DOT_PRODUCT(CONJG(C0U(:,M)),MATMUL(Puu(1,:,:),C0U(:,M))))
                        OCCD(M) = DBLE(DOT_PRODUCT(CONJG(C0D(:,M)),MATMUL(Pdd(1,:,:),C0D(:,M))))
                ENDDO
                IF ( id .EQ. 0 ) THEN
                        WRITE(222,'(E30.20,1000(F12.8))')TIME,OCCU
                        WRITE(333,'(E30.20,1000(F12.8))')TIME,OCCD
                ENDIF
                
                ! Here samples of the density matrices projected on the excited
                ! orbitals of the ground-state calculation are saved
                IF ( INT(NTIMESTEPS/24) .GT. 0 ) THEN
                        IF ( MOD(I-1,INT(NTIMESTEPS/24)) .EQ. 0 .AND. ISAVE .LT. 24 ) THEN
                                ISAVE = ISAVE + 1
                                Csu = 0.0d0
                                Csd = 0.0d0
                                !DO M=Neup,NB
                                DO M=1,NB
                                        IF ( DIFFDENS ) THEN
                                                IF ( M .LE. Neup ) THEN
                                                        Csu(:,M) = sqrt(DABS(1.0d0-OCCU(M)))*C0U(:,M)
                                                ENDIF
                                        ELSE
                                                Csu(:,M) = sqrt(DABS(OCCU(M)))*C0U(:,M)
                                        ENDIF
                                ENDDO
                                !DO M=Nedown,NB
                                DO M=1,NB
                                        IF ( DIFFDENS ) THEN
                                                IF ( M .LE. Nedown ) THEN
                                                        Csd(:,M) = sqrt(DABS(1.0d0-OCCD(M)))*C0D(:,M)
                                                ENDIF
                                        ELSE
                                                Csd(:,M) = sqrt(DABS(OCCD(M)))*C0D(:,M)
                                        ENDIF
                                ENDDO
                                CALL makedens(Csu,NB,PTEMP)
                                PEXu(ISAVE,:,:) =MATMUL(SHH,MATMUL(PTEMP,TRANSPOSE(SHH)))
                                CALL makedens(Csd,NB,PTEMP)
                                PEXd(ISAVE,:,:) =MATMUL(SHH,MATMUL(PTEMP,TRANSPOSE(SHH)))
                                !PEXu(ISAVE,:,:) = DBLE(Pupno(:,:))
                                !PEXd(ISAVE,:,:) = DBLE(Pdownno(:,:))
                                IF ( DIFFDENS ) THEN
                                        Csu = 0.0d0
                                        Csd = 0.0d0
                                        DO M=Neup+1,NB
                                                Csu(:,M) = sqrt(DABS(OCCU(M)))*C0U(:,M)
                                        ENDDO
                                        DO M=Nedown+1,NB
                                                Csd(:,M) = sqrt(DABS(OCCD(M)))*C0D(:,M)
                                        ENDDO
                                        CALL makedens(Csu,NB,PTEMP)
                                        PEXuu(ISAVE,:,:) =MATMUL(SHH,MATMUL(PTEMP,TRANSPOSE(SHH)))
                                        CALL makedens(Csd,NB,PTEMP)
                                        PEXdd(ISAVE,:,:) =MATMUL(SHH,MATMUL(PTEMP,TRANSPOSE(SHH)))
                                ENDIF
                        ENDIF
                ENDIF
                !==================================================================================
                ! Here we propagate the density matrix. In the case of t = 0:
                ! P(Dt) = exp(-i*DT*F(0))*P(0)*exp(i*DT*Fi(0)), and in the case t > 0
                ! and MIDPOINT = .TRUE. : P(t+Dt) = exp(-i*2*DT*F(t))*P(t-Dt)*exp(i*2*DT*Fi(t)),
                ! according to the modified midpoint algorithm, Eqn (8), J. Chem. Phys. 128, 114113
                !==================================================================================
                IF ( MIDPOINT .AND. I .GT. 2 ) THEN
                    DT = 2.0d0*TIMESTEP
                ELSE
                    DT = TIMESTEP
                ENDIF

                Expu = (0.0d0,0.0d0)
                Expd = (0.0d0,0.0d0)

                DO M=1,NB
                    Expu(M,M) = EXP(-IM*DT*EHFeigenup(M))
                    Expd(M,M) = EXP(-IM*DT*EHFeigendown(M))
                ENDDO

                Expu = MATMUL(C1,MATMUL(Expu,TRANSPOSE(CONJG(C1))))
                Expd = MATMUL(C2,MATMUL(Expd,TRANSPOSE(CONJG(C2))))
                
                IF ( MIDPOINT .AND. I .GT. 2 ) THEN
                        Pup   = MATMUL(Expu,MATMUL(Puu(1,:,:),TRANSPOSE(CONJG(Expu))))
                        Pdown = MATMUL(Expd,MATMUL(Pdd(1,:,:),TRANSPOSE(CONJG(Expd))))
                ELSE
                        Pup   = MATMUL(Expu,MATMUL(Puu(2,:,:),TRANSPOSE(CONJG(Expu))))
                        Pdown = MATMUL(Expd,MATMUL(Pdd(2,:,:),TRANSPOSE(CONJG(Expd))))
                ENDIF

                Pup = 0.50d0*(Pup + TRANSPOSE(CONJG(Pup)) )
                Pdown = 0.50d0*(Pdown + TRANSPOSE(CONJG(Pdown)) )
                
                !===========================================================================================================================
                ! Here we try and stabilize the time-propagation in the case of DFT. Since the accuracy as regards of the total number of
                ! electrons unfortunately depend on the integration grid which (I suspect) create instabilities. Therefore here we make sure
                ! that the total number of up and down electrons are kept konstant throughout the propagation. Wed Mar  4 09:19:35 CET 2015
                !===========================================================================================================================
                
                !Pup   = (Neup/TRACE(DBLE(Pup),NB))*Pup
                !Pdown = (Nedown/TRACE(DBLE(Pdown),NB))*Pdown

                !=============================================================================================================================================
                ! This subroutine solves the equation P(t+Dt) = exp(0.5*i*Dt*(F[P(t)]+F[P(t+Dt)]))*P(t)*exp(-0.5*i*Dt*(F[P(t)]+F[P(t+Dt)]))*P(t),
                ! by iteration
                !=============================================================================================================================================
                
                IF ( NSCCORR .GT. 0 .AND. I .GT. 1 ) THEN
                        CALL TDFTCORR(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE, &
                        & Cupo,Cdowno,Pup,Pdown,numprocessors,id,OMEGA,EDIR,NEPERIOD,EPROFILE,NSCCORR,MIXTDDFT,TIMESTEP,EFIELDMAX,IORBNR, &
                        & Fups(2,:,:),Fdowns(2,:,:),Puu(2,:,:),Pdd(2,:,:),CONV,SCERR,TIME+DT,Ncorr,Fupo,Fdowno)
                ENDIF
     
                !===========================================================================================================================
                ! Calculating the total energy, which is no longer a good quantum number.
                !===============================================================================================
                IF ( CORRLEVEL .NE. 'URHF' ) THEN
                        Fup   = Fup   - Vxc(1,:,:)*RE
                        Fdown = Fdown - Vxc(2,:,:)*RE
                ENDIF

                IF ( MOD(I-1,20) .EQ. 0 ) THEN
                        ETOT = 0.50d0*DBLE(SUM(H0*(Pupno+Pdownno)) + SUM(Fup*Pupno) + SUM(Fdown*Pdownno) )+nucE
                        IF ( CORRLEVEL .NE. 'URHF' ) THEN
                                ETOT = ETOT + exc(CORRLEVEL,NATOMS,ATOMS,BAS,Pupr,Pdownr,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id)
                        ENDIF
                ENDIF

                IF ( EPROFILE .EQ. 'CIRC' ) THEN
                        EFIELDP = EFIELD*sin(OMEGA*TIME) 
                ELSE
                        EFIELDP = EFIELD
                ENDIF

                IF ( I .EQ. 1 .AND.  POUT .AND. id .EQ. 0 ) THEN
                        print*,' '
                        print*,'             =========================================================='
                        print*,'                        Entering the TDFT time propagation             '
                        print*,'             =========================================================='
                        print*,' '
                        IF ( NSCCORR .LE. 0 ) THEN
                                WRITE(*,'(A9,A22,A34,A24,A38,A24)')'N','Time [au]','Dipole moment [au]',' E [au] ','   Number of electrons  ',' E-field [a.u] '
                        ELSE
                                WRITE(*,'(A9,A22,A34,A24,A38,A24,A24,A9)')'N','Time [au]','Dipole moment [au]',' E [au] ','   Number of electrons  ',' E-field [a.u] ',&
                                & '     DP(SC)   ','    Ncorr'
                        ENDIF   
                ENDIF

                IF ( id .EQ. 0 ) THEN
                        IF ( NSCCORR .LE. 0 ) THEN
                                IF ( I .EQ. 1 ) WRITE(111,'(A9,A22,A34,A24,A38,A24)')'N','Time [au]','Dipole moment [au]',' E [au] ','   Number of electrons  ',' E-field [a.u] '

                                IF ( POUT ) WRITE(*,'(I9,E30.20,E30.20,E30.20,E30.20,E30.20,E30.20)'),I-1,TIME,MOMENT(I),ETOT,TRACE(DBLE(PT),NB),EFIELDP
                                WRITE(111,'(I9,E30.20,E30.20,E30.20,E30.20,E30.20,E30.20)'),I-1,TIME,MOMENT(I),ETOT,TRACE(DBLE(PT),NB),EFIELDP
                        ELSE              
                                IF ( I .EQ. 1 ) WRITE(111,'(A9,A22,A34,A24,A38,A24,A24,A9)')'N','Time [au]','Dipole moment [au]',' E [au] ','   Number of electrons  ',' E-field [a.u] ',&
                                & '     DP(SC)   ','    Ncorr'

                                IF ( POUT ) WRITE(*,'(I9,E30.20,E30.20,E30.20,E30.20,E30.20,E14.2,I9)'),I-1,TIME,MOMENT(I),ETOT,TRACE(DBLE(PT),NB),EFIELDP,CONV,Ncorr
                                WRITE(111,'(I9,E30.20,E30.20,E30.20,E30.20,E30.20,E14.2,I9)'),I-1,TIME,MOMENT(I),ETOT,TRACE(DBLE(PT),NB),EFIELDP,CONV,Ncorr
                        ENDIF
                ENDIF

     ENDDO
     IF ( id .EQ. 0 ) THEN
        CLOSE(111)
        CLOSE(222)
        CLOSE(333)
        IF ( DOABSSPECTRUM ) CALL fouriertransf(NTIMESTEPS,TIMESTEP,MOMENT,EFIELDMAX)
     ENDIF
     CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE TDFT
