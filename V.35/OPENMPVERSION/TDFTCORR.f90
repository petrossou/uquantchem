SUBROUTINE TDFTCORR(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE, &
           & Cup,Cdown,Pup,Pdown,ETEMP,OMEGA,EDIR,NEPERIOD,EPROFILE,NSCCORR,MIXTDDFT,TIMESTEP,EFIELDMAX,IORBNR,Fuprevo,&
           & Fdownprevo,Pupprevo,Pdownprevo,SCERR,CONV,TIME,Ncorr,NBAUX,VRI,WRI,RIAPPROX,Fupo,Fdowno)
      ! This subroutine calculates the self consistent DFT solution
      USE datatypemodule
      IMPLICIT NONE
      CHARACTER(LEN=20) :: CORRLEVEL,EPROFILE
      LOGICAL, INTENT(IN) :: RIAPPROX
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      INTEGER, INTENT(OUT) :: Ncorr
      INTEGER, INTENT(IN) :: NB,Ne,LORDER,CGORDER,NATOMS,NTOTALQUAD,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD),NEPERIOD,EDIR,IORBNR(2),NSCCORR,NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      DOUBLE PRECISION, INTENT(IN) :: VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),gradS(NATOMS,3,NB,NB),H0(NB,NB),Intsv(NRED),nucE,LQ(LORDER,3),CGQ(CGORDER,2),ETEMP,OMEGA,EFIELDMAX,TIME
      DOUBLE PRECISION, INTENT(OUT) :: TIMESTEP,CONV
      DOUBLE PRECISION :: EHFeigenup(NB),EHFeigendown(NB)
      COMPLEX*16, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB),Fupo(NB,NB),Fdowno(NB,NB)
      COMPLEX*16, INTENT(IN) :: Fuprevo(NB,NB),Fdownprevo(NB,NB),Pupprevo(NB,NB),Pdownprevo(NB,NB)
      DOUBLE PRECISION, INTENT(IN)  :: MIXTDDFT,SCERR
      DOUBLE PRECISION :: FTOT,FOLD,BLTENSOR(3,NB,NB),LTESORu(3,NB,NB),LTESORd(3,NB,NB)
      COMPLEX*16 :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB)
      DOUBLE PRECISION :: DE,EOLD,DELTAP,LAMDAu,LAMDAd,MIXING,L2u(NB,NB),L2d(NB,NB),EFIELD,PTEMP(NB,NB)
      COMPLEX*16 :: Fup(NB,NB),Fdown(NB,NB),G(NB,NB),C1(NB,NB),C2(NB,NB),C3(NB,NB),C4(NB,NB),C0U(NB,NB),C0D(NB,NB),PHOLE(NB,NB)
      COMPLEX*16 :: Puu(2,NB,NB),Pdd(2,NB,NB),Expu(NB,NB),Expd(NB,NB),Pupno(NB,NB),Pdownno(NB,NB),Fups(2,NB,NB),Fdowns(2,NB,NB),F3(NB,NB),F4(NB,NB)
      DOUBLE PRECISION :: TOLDNe,DPTENSOR(3,NB,NB),DIPOLET(NB,NB),OCCU(NB),OCCD(NB),DTENS(NB,NB),SHP(NB,NB),Csu(NB,NB),Csd(NB,NB)
      COMPLEX*16 :: PTold(NB,NB),Pupold(NB,NB),Pdownold(NB,NB),Pups(50,NB,NB),Pdowns(50,NB,NB),Pupt(NB,NB),Pdownt(NB,NB),EIGENVECT(NB,NB),DIPOLE(NB,NB)
      DOUBLE PRECISION :: Vxc(2,NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),CDB(NB,NB),Pupr(NB,NB),Pdownr(NB,NB),EIG(NB,NB),DT
      COMPLEX*16 :: TESTA(NB,NB),ERRSU(50,NB,NB),ERRSD(50,NB,NB),ERRU(NB,NB),ERRD(NB,NB),RE,IM,SHH(NB,NB),SHHP(NB,NB)
      INTEGER :: I,II,III,L,M,N,Neup,Nedown,INFO1,INFO2,RESET,ISTART
      INTEGER :: MAXITER,ISAVE
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
      Pupold = (0.0d0,0.0d0)
      Pdownold = (0.0d0,0.0d0)

      ISTART = 0
     
      ! If the profile of the EM-wave is squre the dipolmatrix elements 
      ! are calculated here
      CALL dipoletensor(BAS,DPTENSOR)
      EFIELD = EFIELDMAX
      IF ( EDIR .EQ. 1 ) DIPOLET = DPTENSOR(2,:,:)
      IF ( EDIR .EQ. 2 ) DIPOLET = DPTENSOR(3,:,:)
      IF ( EDIR .EQ. 3 ) DIPOLET = DPTENSOR(1,:,:)
      DTENS = DIPOLET

      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIG)
      SL = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
      ENDDO
     
      SH = MATMUL(EIG,MATMUL(SL,TRANSPOSE(EIG)))
      SHH = RE*SH

      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2

      I = 1
      DT = TIMESTEP
      CONV = 2*DABS(SCERR)

      Puu(2,:,:) = Pup
      Pdd(2,:,:) = Pdown
      
      DO WHILE ( I .LE. NSCCORR .AND. CONV .GT. SCERR )
               
                Pupold = Pup
                Pdownold = Pdown
 
                ! ==========================================
                ! Transform to non-orthogonal representation:
                !===========================================
                Pupno   = MATMUL(SHH,MATMUL(Pup,TRANSPOSE(SHH)))
                Pdownno = MATMUL(SHH,MATMUL(Pdown,TRANSPOSE(SHH)))
                
                Pupno = 0.50d0*( Pupno + TRANSPOSE(CONJG(Pupno)) )
                Pdownno = 0.50d0*( Pdownno + TRANSPOSE(CONJG(Pdownno)) )

                CALL getJvc(Pupno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jup)
                CALL getJvc(Pdownno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jdown)
                
                Pupr = DBLE(Pupno)
                Pdownr = DBLE(Pdownno)
                
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                   
                        CALL getKvc(Pupno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                        CALL getKvc(Pdownno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                   
                        Fup   = H0 +  Jdown + Jup - Kup 
                        Fdown = H0 +  Jdown + Jup - Kdown
                ELSE 
                
                        CALL getvxc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pupr,Pdownr,gradS,LORDER,CGORDER,LQ,CGQ,Vxc)
                      
                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                CALL getKvc(Pupno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                                CALL getKvc(Pdownno,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
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
                      DIPOLE = (0.0d0,0.0d0)
                      EFIELD = 0.0d0
                ENDIF

                Fup   = Fup + DIPOLE
                Fdown = Fdown + DIPOLE
                
                Fup   = 0.50d0*(Fup + TRANSPOSE(CONJG(Fup)) )
                Fdown = 0.50d0*(Fdown + TRANSPOSE(CONJG(Fdown)) )
                
                !=======================================
                ! Transform to orthogonal representation:
                !=======================================
                
                Fupo = MATMUL(SH,MATMUL(Fup,SH))
                Fdowno = MATMUL(SH,MATMUL(Fdown,SH))

                Fupo   = 0.50d0*(Fupo + TRANSPOSE(CONJG(Fupo)) )
                Fdowno = 0.50d0*(Fdowno + TRANSPOSE(CONJG(Fdowno)) )

                CALL diaghc( 0.50d0*(Fupo+Fuprevo),NB,EHFeigenup,C1,INFO2)
                CALL diaghc( 0.50d0*(Fdowno+Fdownprevo),NB,EHFeigendown,C2,INFO2)
              
                PT = Pup + Pdown
                
                !==================================================================================
                ! Here we propagate the density matrix. In the case of t = 0:
                ! P(Dt) = exp(-i*DT*F(0))*P(0)*exp(i*DT*Fi(0)), and in the case t > 0
                ! and MIDPOINT = .TRUE. : P(t+Dt) = exp(-i*2*DT*F(t))*P(t-Dt)*exp(i*2*DT*Fi(t)),
                ! according to the modified midpoint algorithm, Eqn (8), J. Chem. Phys. 128, 114113
                !==================================================================================
                Expu = (0.0d0,0.0d0)
                Expd = (0.0d0,0.0d0)
                

                DO M=1,NB
                    Expu(M,M) = EXP(-IM*DT*EHFeigenup(M))
                    Expd(M,M) = EXP(-IM*DT*EHFeigendown(M))
                ENDDO
                
                Expu = MATMUL(C1,MATMUL(Expu,TRANSPOSE(CONJG(C1))))
                Expd = MATMUL(C2,MATMUL(Expd,TRANSPOSE(CONJG(C2))))
                
                Pup   = (1.0d0-MIXTDDFT)*MATMUL(Expu,MATMUL(Pupprevo,TRANSPOSE(CONJG(Expu)))) + MIXTDDFT*Pup
                Pdown = (1.0d0-MIXTDDFT)*MATMUL(Expd,MATMUL(Pdownprevo,TRANSPOSE(CONJG(Expd)))) + MIXTDDFT*Pdown


                Pup = 0.50d0*(Pup + TRANSPOSE(CONJG(Pup)) )
                Pdown = 0.50d0*(Pdown + TRANSPOSE(CONJG(Pdown)) )

                !===========================================================================================================================
                ! Here we try and stabilize the time-propagation in the case of DFT. Since the accuracy as regards of the total number of
                ! electrons unfortunately depend on the integration grid which (I suspect) create instabilities. Therefore here we make sure
                ! that the total number of up and down electrons are kept  konstant throughout the propagation. Wed Mar  4 09:19:35 CET 2015
                !===========================================================================================================================

                !Pup   = (Neup/TRACE(DBLE(Pup),NB))*Pup
                !Pdown = (Nedown/TRACE(DBLE(Pdown),NB))*Pdown

                !====================================================================
                ! Calculating the difference in P between two consecutive iterations
                !====================================================================
                CONV = 0.0d0
                DO M=1,NB
                        DO N=1,NB
                                CONV = CONV + (1.0d0/(2.0d0*NB**2))*(ABS(DBLE(Pup(N,M)-Pupold(N,M))) + ABS(DBLE(Pdown(N,M)-Pdownold(N,M))))
                        ENDDO
                ENDDO

                !===============================================================================================
                ! Calculating the total energy, which is no longer a good quantum number.
                !===============================================================================================
                IF ( CORRLEVEL .NE. 'URHF' ) THEN
                        Fup   = Fup   - Vxc(1,:,:)*RE
                        Fdown = Fdown - Vxc(2,:,:)*RE
                ENDIF
                Ncorr = I
                I = I + 1
     ENDDO

END SUBROUTINE TDFTCORR

