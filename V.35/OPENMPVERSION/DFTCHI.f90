SUBROUTINE DFTCHI(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,S,gradS,H0,Intsv,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown, &
& ETOT,Cup,Cdown,Pup,Pdown,Pgup,Pgdown,MIX,DIISORD,DIISSTART,NSCF,FIXNSCF,POUT,SCRATCH,ZEROSCF,ETEMP,mu,ENTROPY,IORBNR,PHOLE,PEXCITED,NBAUX,VRI,WRI,RIAPPROX)
      ! This subroutine calculates the self consistent DFT solution
      USE datatypemodule
      IMPLICIT NONE
      CHARACTER(LEN=20) :: CORRLEVEL
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      LOGICAL, INTENT(IN) :: POUT,SCRATCH,ZEROSCF,RIAPPROX
      INTEGER, INTENT(IN) :: NB,Ne,FIXNSCF,LORDER,CGORDER,NATOMS,NTOTALQUAD,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD),IORBNR(2),NBAUX
      INTEGER*8, INTENT(IN) :: NRED
      INTEGER, INTENT(INOUT) :: DIISORD,DIISSTART
      INTEGER, INTENT(OUT) :: NSCF
      DOUBLE PRECISION, INTENT(IN) :: S(NB,NB),gradS(NATOMS,3,NB,NB),H0(NB,NB),Intsv(NRED),Tol,nucE,MIX,LQ(LORDER,3),CGQ(CGORDER,2),ETEMP
      DOUBLE PRECISION, INTENT(IN) :: VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX)
      DOUBLE PRECISION, INTENT(OUT) :: EHFeigenup(NB),EHFeigendown(NB),ETOT,ENTROPY,PHOLE(NB,NB),PEXCITED(NB,NB)
      DOUBLE PRECISION, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB),mu,Pgup(NB,NB),Pgdown(NB,NB)
      DOUBLE PRECISION :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),FTOT,FOLD,PHOLEold(NB,NB),PEXCITEDold(NB,NB)
      DOUBLE PRECISION :: Fup(NB,NB),Fdown(NB,NB),G(NB,NB),C1(NB,NB),C2(NB,NB),C3(NB,NB),C4(NB,NB),DE,EOLD,DELTAP,LAMDAu,LAMDAd,LAMDAh,MIXING
      DOUBLE PRECISION :: PTold(NB,NB),Pupold(NB,NB),Pdownold(NB,NB),Pups(50,NB,NB),Pdowns(50,NB,NB),Pupt(NB,NB),Pdownt(NB,NB),TOLDNe
      DOUBLE PRECISION :: Qs(NB,NB),Q(NB,NB),QQ(NB,NB),Ch(NB,NB),Fh(NB,NB),Eholeigen(NB),Fsaved(NB,NB),EN(NB,NB),VEC(NB)
      DOUBLE PRECISION :: Vxc(2,NB,NB),Vxchole(2,NB,NB),TESTA(NB,NB),ERRSU(50,NB,NB),ERRSD(50,NB,NB),ERRU(NB,NB),ERRD(NB,NB),SH(NB,NB),SL(NB,NB),LAM(NB),EIGENVECT(NB,NB)
      DOUBLE PRECISION :: ERRSh(50,NB,NB),ERRh(NB,NB),PHOLEt(NB,NB),PHOLEs(50,NB,NB)
      INTEGER :: I,II,III,JJJ,L,M,N,Neup,Nedown,INFO1,INFO2,INFO3,JJ,RESET
      INTEGER :: MAXITER
      DOUBLE PRECISION, EXTERNAL :: exc,quadcheck
      INTEGER, EXTERNAL :: ijkl,TRACE
      LOGICAL :: STARTPRINTDIISIFO

      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2

      MIXING  = MIX

      IF ( FIXNSCF .GT. 0 ) THEN
              MAXITER = FIXNSCF-1
      ELSE
              MAXITER = 400
      ENDIF
      RESET = 0

      ! Calculating the lowdin S^(-1/2) matrix
      CALL diagh( S,NB,LAM,EIGENVECT)
      SL = 0.0d0
      EN = 0.0d0
      DO I=1,NB
        SL(I,I) = 1.0d0/sqrt(LAM(I))
        EN(I,I) = 1.0d0
      ENDDO
      SH = MATMUL(EIGENVECT,MATMUL(SL,TRANSPOSE(EIGENVECT)))

      IF ( DIISORD .GT. 25 ) THEN
              WRITE(*,*)'The order of DIIS mixing set =',DIISORD,' Is too Big.'
              WRITE(*,*)'Maximum allowed value is 25. This value has been used instead.'
              DIISORD = 25
      ENDIF

       IF ( SCRATCH ) THEN
                ! Calulating the hole and excitation contributions to the density matrix
                C1 = 0.0d0
                IF ( IORBNR(1) .GT. 0 ) THEN
                        C1(:,IORBNR(1)) = Cup(:,IORBNR(1))
                ELSE
                        C1(:,-IORBNR(1)) = Cdown(:,-IORBNR(1))
                ENDIF
                CALL makedens(C1,NB,PHOLE)
       ENDIF
      Pups = 0.0d0
      Pdowns = 0.0d0
      PHOLEs = 0.0d0
      ERRSU = 0.0d0
      ERRSD = 0.0d0
      ERRSh = 0.0d0
      
      DE = 2.0d0*Tol
      DELTAP = 2.0d0*Tol
      I = 0
      II = 0
      PTold = 0.0d0
      STARTPRINTDIISIFO = .FALSE.

      ! The tolerance used when calculating 
      ! the chemical potential 
      TOLDNe = 1.0E-8

      !=======================================================
      ! This is used for the zero scf cycle (FAST-MD)
      !=======================================================
      IF ( ZEROSCF ) THEN
                PT = Pup + Pdown
                CALL getJv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jup)
                CALL getJv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jdown)
                CALL getvxc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,Vxc)
                IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL getKv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                        CALL getKv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                        Fup   = H0 + Vxc(1,:,:) + Jdown - 0.20d0*Kup   + Jup
                        Fdown = H0 + Vxc(2,:,:) + Jdown - 0.20d0*Kdown + Jup
                ELSE
                        Fup   = H0 + Vxc(1,:,:) + Jdown + Jup
                        Fdown = H0 + Vxc(2,:,:) + Jdown + Jup
                ENDIF

                ! Constructing the hole-fockian and projected fockian according to eqn (29) p. 88. in my blue note-book:
                IF ( IORBNR(1) .GT. 0 ) THEN
                        ! The hole-fockian:
                        CALL gethole(Fup,S,NB,Neup,ABS(IORBNR(1)),PHOLE,Pgup,Pup)
                        ! The projected fockian:
                        Fsaved = Fup
                        QQ = EN - MATMUL(S,PHOLE)
                        Fup = MATMUL(MATMUL(QQ,Fup),TRANSPOSE(QQ))
                ELSE
                        ! The hole-fockian:
                        CALL gethole(Fdown,S,NB,Nedown,ABS(IORBNR(1)),PHOLE,Pgdown,Pdown)
                        ! The projected fockian:
                        Fsaved = Fdown
                        QQ = EN - MATMUL(S,PHOLE)
                        Fdown = MATMUL(MATMUL(QQ,Fdown),TRANSPOSE(QQ))
                ENDIF
                
                CALL diaghHF( Fup,S,NB,EHFeigenup,C1,INFO1)
                CALL diaghHF( Fdown,S,NB,EHFeigendown,C2,INFO2)
                
                IF ( IORBNR(1) .GT. 0 ) THEN
                        Fup   = Fsaved   - Vxc(1,:,:)
                        Fdown = Fdown - Vxc(2,:,:)
                ELSE
                        Fup   = Fup   -  Vxc(1,:,:)
                        Fdown = Fsaved - Vxc(2,:,:)
                ENDIF
                !ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + exc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ) + nucE

                IF ( ETEMP .LT. 0.0d0 ) THEN
                        Cup(:,:) = 0.0d0
                        JJJ = 0
                        IF ( IORBNR(1) .GT. 0 ) THEN
                                DO M=1,NB
                                        VEC = C1(:,M)
                                        IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Neup ) THEN
                                                JJJ = JJJ + 1
                                                Cup(:,JJJ) = VEC
                                        ENDIF
                                ENDDO
                        ELSE
                                DO M=1,Neup
                                        Cup(:,M) = C1(:,M)
                                ENDDO
                        ENDIF

                        Cdown(:,:) = 0.0d0
                        JJJ = 0
                        IF ( IORBNR(1) .LT. 0 ) THEN
                                DO M=1,NB
                                        VEC = C2(:,M)
                                        IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Nedown ) THEN
                                                JJJ = JJJ + 1
                                                Cdown(:,JJJ) = VEC
                                        ENDIF
                                ENDDO
                        ELSE
                                DO M=1,Nedown
                                        Cdown(:,M) = C2(:,M)
                                ENDDO
                        ENDIF 
                        CALL makedens(Cup,NB,Pup)
                        CALL makedens(Cdown,NB,Pdown)
                        NSCF = 0
                ELSE
                        Cup = C1
                        Cdown = C2
                        CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                ENDIF
                NSCF = 0
                RETURN
      ENDIF
      !=====================================================================
     
      DO WHILE ( (DABS(DE) .GT. Tol .AND. I .LE. MAXITER) .OR. ( DELTAP .GT. sqrt(Tol) .AND. I .LE. MAXITER) )
                IF ( ETEMP .LE. 0.0d0 ) THEN
                        Cup(:,:) = 0.0d0
                        JJJ = 0
                        IF ( IORBNR(1) .GT. 0 ) THEN
                                DO M=1,NB
                                        VEC = C1(:,M)
                                        IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Neup ) THEN
                                                JJJ = JJJ + 1
                                                Cup(:,JJJ) = VEC
                                        ENDIF
                                ENDDO
                        ELSE
                                DO M=1,Neup
                                        Cup(:,M) = C1(:,M)
                                ENDDO
                        ENDIF
                        
                        Cdown(:,:) = 0.0d0
                        JJJ = 0
                        IF ( IORBNR(1) .LT. 0 ) THEN
                                DO M=1,NB
                                        VEC = C2(:,M)
                                        IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Nedown ) THEN
                                                JJJ = JJJ + 1
                                                Cdown(:,JJJ) = VEC
                                        ENDIF
                                ENDDO
                        ELSE
                                DO M=1,Nedown
                                        Cdown(:,M) = C2(:,M)
                                ENDDO
                        ENDIF
                ELSE
                       Cup = C1
                       Cdown = C2
                ENDIF
               
                !==============================================================
                ! In the case of XL-BOMD, it is crusial that the scf starts 
                ! with the density matrices Pup and Pdown provided by the input,
                ! since it is these matrices that have been propagated by the 
                ! XL-BOMD scheme and not the coefficient matrices Cup and Cdown
                ! Therefor we have the if-construct here.
                !==============================================================
                IF ( I .GT. 0 ) THEN
                        IF ( ETEMP .LE. 0.0d0 ) THEN
                                CALL makedens(Cup,NB,Pup)
                                CALL makedens(Cdown,NB,Pdown)
                        ELSE
                                CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
                        ENDIF
                        C3 = 0.0d0
                        IF ( IORBNR(1) .GT. 0 ) THEN
                                DO M=Neup+1,NB
                                        C3(:,M) = C1(:,M)
                                ENDDO
                         ELSE
                                DO M=Nedown+1,NB
                                        C3(:,M) = C2(:,M)
                                ENDDO
                        ENDIF
                        CALL makedens(C3,NB,Q)
                ENDIF
                PT = Pup + Pdown
               
                ! Calculating the change of the total density matrix:
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0  ) DELTAP = sqrt(DOT_PRODUCT(reshape(PT-PTold,(/NB**2/)),reshape(PT-PTold,(/NB**2/))))
                
                ! Linear mixing:
                IF ( I .GT. 0 ) THEN
                        Pup = Pup*(1.0d0-MIXING) + Pupold*MIXING
                        Pdown = Pdown*(1.0d0-MIXING) + Pdownold*MIXING
                        !PHOLE = PHOLE*(1.0d0-MIXING) + PHOLEold*MIXING
                ELSE
                        Pup = Pup
                        Pdown = Pdown
                ENDIF
                
                Pupold = Pup
                Pdownold = Pdown
                PHOLEold = PHOLE

                PTold = Pupold + Pdownold

                ! Saving the density matrices from the previous
                ! iterations, to be used with the DIIS-method.
                IF ( II .LT. 2*DIISORD ) THEN
                        II = II + 1
                        Pups(II,:,:) = Pup
                        Pdowns(II,:,:) = Pdown
                        PHOLEs(II,:,:) = PHOLE
                        ERRSU(II,:,:) = ERRU
                        ERRSD(II,:,:) = ERRD
                        ERRSh(II,:,:) = ERRh
                ELSE
                        ! refreshing the density matrices used by the DIIS
                        IF ( DIISORD  .NE. 0 ) THEN
                                DO III=1,II-1
                                        Pups(III,:,:) = Pups(III+1,:,:)
                                        Pdowns(III,:,:) = Pdowns(III+1,:,:)
                                        PHOLEs(III,:,:) = PHOLEs(III+1,:,:)
                                        ERRSU(III,:,:) = ERRSU(III+1,:,:)
                                        ERRSD(III,:,:) = ERRSD(III+1,:,:)
                                        ERRSh(III,:,:) = ERRSh(III+1,:,:)
                                ENDDO
                                Pups(II,:,:) = Pup
                                Pdowns(II,:,:)  = Pdown
                                PHOLEs(II,:,:)  = PHOLE
                                ERRSU(II,:,:) = ERRU
                                ERRSD(II,:,:) = ERRD
                                ERRSh(II,:,:) = ERRh
                        ENDIF
                ENDIF

                ! Here we use the DIIS method ( CHEM. Phys. Lett. 73, 393 (1980) )
                ! in order to estimate the self-consistent density matrices:
                IF ( II .GE. 4 ) THEN

                        CALL DIIS(NB,II,ERRSU,Pups,Pupt,LAMDAu,INFO1)
                        CALL DIIS(NB,II,ERRSD,Pdowns,Pdownt,LAMDAd,INFO2)
                        CALL DIIS(NB,II,ERRSh,PHOLEs,PHOLEt,LAMDAh,INFO3)
                        
                        IF ( (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. LAMDAu+LAMDAd .LT. 1.0d0) .OR. (INFO1 .EQ. 0 .AND. INFO2 .EQ. 0 .AND. I .GT. DIISSTART )  ) THEN
                                IF ( SUM(Pupt*S) .GT. 0.750d0*Neup .AND.  SUM(Pdownt*S) .GT. 0.750d0*Nedown  ) THEN
                                        Pup = Pupt
                                        Pdown = Pdownt
                                        !PHOLE = PHOLEt
                                        STARTPRINTDIISIFO = .TRUE.
                                        MIXING = 0.0d0
                                ENDIF
                        ENDIF
                ENDIF

10 CONTINUE     
                PT = Pup + Pdown
                
                CALL getJv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jup)
                CALL getJv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Jdown)
                
                CALL getvxc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,Vxc)
                
                ! First  iteration is done with Hartree-Fock
                IF ( I .LT. 1 .AND. SCRATCH ) THEN
                   
                        CALL getKv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                        CALL getKv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                   
                        Fup   = H0 +  Jdown + Jup - Kup 
                        Fdown = H0 +  Jdown + Jup - Kdown
                ELSE
                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                CALL getKv(Pup,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kup)
                                CALL getKv(Pdown,NB,NBAUX,NRED,Intsv,VRI,WRI,RIAPPROX,Kdown)
                                Fup   = H0 +  Jdown + Jup  + Vxc(1,:,:) - 0.20d0*Kup
                                Fdown = H0 +  Jdown + Jup  + Vxc(2,:,:) - 0.20d0*Kdown
                        ELSE
                                Fup   = H0 +  Jdown + Jup + Vxc(1,:,:)
                                Fdown = H0 +  Jdown + Jup + Vxc(2,:,:)
                        ENDIF
                ENDIF

                ! Constructing the hole-fockian and projected fockian according to eqn (29) p. 88. in my blue note-book:
                IF ( IORBNR(1) .GT. 0 ) THEN
                        ! The hole-fockian:
                        IF ( I .GT. 0 ) CALL gethole(Fup,S,NB,Neup,ABS(IORBNR(1)),PHOLE,Pgup,Pup)
                        !print*,'---',SUM(PHOLE*S)
                        ! The projected fockian:
                        Fsaved = Fup
                        QQ = EN - MATMUL(S,PHOLE)
                        Fup = MATMUL(MATMUL(QQ,Fup),TRANSPOSE(QQ))
                ELSE
                        ! The hole-fockian:
                        IF ( I .GT. 0 ) CALL gethole(Fdown,S,NB,Nedown,ABS(IORBNR(1)),PHOLE,Pgdown,Pdown)
                        ! The projected fockian:
                        Fsaved = Fdown
                        QQ = EN - MATMUL(S,PHOLE)
                        Fdown = MATMUL(MATMUL(QQ,Fdown),TRANSPOSE(QQ))
                ENDIF
                ! Here we calculate the error vector/matrix according to 
                ! P. Pulay, J. Comp. Chem. 3, 556 (1982) (Eqn 4)
                ! to be used in the DIIS algorithm:
                IF ( DIISORD .NE. 0 ) THEN
                        ERRU = MATMUL(Fup,MATMUL(Pup,S)) - MATMUL(S,MATMUL(Pup,Fup))
                        ERRD = MATMUL(Fdown,MATMUL(Pdown,S)) - MATMUL(S,MATMUL(Pdown,Fdown))
                        ERRh = MATMUL(Fh,MATMUL(PHOLE,S)) - MATMUL(S,MATMUL(PHOLE,Fh))
                        ERRU = MATMUL(TRANSPOSE(SH),MATMUL(ERRU,SH))
                        ERRD = MATMUL(TRANSPOSE(SH),MATMUL(ERRD,SH))
                        ERRh = MATMUL(TRANSPOSE(SH),MATMUL(ERRh,SH))
                ENDIF
                
                CALL diaghHF( Fup,S,NB,EHFeigenup,C1,INFO1)
                CALL diaghHF( Fdown,S,NB,EHFeigendown,C2,INFO2)
                
               IF ( (INFO1 .NE. 0 .OR. INFO2 .NE. 0 ) .AND. RESET .EQ. 0 ) THEN
                        Pup   =  Pupold
                        Pdown = Pdownold
                        PHOLE = PHOLEold
                        RESET = 1
                        GOTO 10
                ENDIF
                RESET = 0               
 
                IF ( I .LT. 1 .AND. SCRATCH ) THEN
                        ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + nucE
                ELSE
                        IF ( IORBNR(1) .GT. 0 ) THEN
                                Fup = Fsaved
                        ELSE
                                Fdown = Fsaved
                        ENDIF
                        Fup   = Fup   - Vxc(1,:,:)
                        Fdown = Fdown - Vxc(2,:,:)
                        IF ( FIXNSCF .LT. 0 ) THEN
                           ETOT = 0.50d0*(SUM(H0*PT) + SUM(Fup*Pup) + SUM(Fdown*Pdown) ) + exc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ) + nucE
                           IF ( ETEMP .GT. 0.0d0  ) FTOT = ETOT - ETEMP*ENTROPY
                        ENDIF
                        
                ENDIF
                
                IF ( I .EQ. 0 .AND. POUT ) THEN
                        print*,'   =========================================================='
                        print*,'                  Entering the scf loop                      '
                        print*,'   =========================================================='
                        print*,' '
                        IF ( ETEMP .LT. 0.0d0  ) WRITE(*,'(A4,A22,A27,A27,A33)')'N','E [au]','DE [au]',' DP','DIIS'
                        IF ( ETEMP .GT. 0.0d0  ) WRITE(*,'(A4,A20,A29,A31,A27,A31)')'N','E [au]','F [au]','DF [au]',' DP','DIIS'
                ENDIF
        
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .LT. 0.0d0 ) DE = ETOT-EOLD
                IF ( I .GT. 0 .AND. FIXNSCF .LE. 0 .AND. ETEMP .GT. 0.0d0 ) DE = FTOT-FOLD
                
                IF ( STARTPRINTDIISIFO ) THEN
                        IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP,LAMDAd+LAMDAu
                        IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP,LAMDAd+LAMDAu
                ELSE
                        IF ( POUT .AND. ETEMP .GT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)'),I,ETOT,FTOT,DE,DELTAP
                        IF ( POUT .AND. ETEMP .LT. 0.0d0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20)'),I,ETOT,DE,DELTAP
                ENDIF
        
                EOLD = ETOT
                FOLD = FTOT
                I = I+1
     ENDDO

     Cup = C1
     Cdown = C2
    
     IF ( ETEMP .LT. 0.0d0 ) THEN
                Cup(:,:) = 0.0d0
                JJJ = 0
                IF ( IORBNR(1) .GT. 0 ) THEN
                        DO M=1,NB
                                VEC = C1(:,M)
                                IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Neup ) THEN
                                        JJJ = JJJ + 1
                                        Cup(:,JJJ) = VEC
                                ENDIF
                        ENDDO
                ELSE
                        DO M=1,Neup
                                Cup(:,M) = C1(:,M)
                        ENDDO
                ENDIF

                Cdown(:,:) = 0.0d0
                JJJ = 0
                IF ( IORBNR(1) .LT. 0 ) THEN
                        DO M=1,NB
                                VEC = C2(:,M)
                                IF ( DOT_PRODUCT(MATMUL(MATMUL(PHOLE,S),VEC),MATMUL(MATMUL(PHOLE,S),VEC)) .LT. 1.0E-12 .AND. JJJ .LT. Nedown ) THEN
                                        JJJ = JJJ + 1
                                        Cdown(:,JJJ) = VEC
                                ENDIF
                        ENDDO
                ELSE
                        DO M=1,Nedown
                                Cdown(:,M) = C2(:,M)
                        ENDDO
                ENDIF       
                CALL makedens(Cup,NB,Pup)
                CALL makedens(Cdown,NB,Pdown)
                Cup = C1
                Cdown = C2
     ELSE
        CALL makedensT(TOLDNe,Cup,Cdown,EHFeigenup,EHFeigendown,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
     ENDIF
     PT = Pup+Pdown

     NSCF = I

     IF ( I .LE. MAXITER .OR. FIXNSCF .GT. 0 ) THEN
             IF ( POUT ) THEN
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
                WRITE(*,'(A75,F9.6)' )'Number of electrons calculated from integrating the charge-density, Ne =  ',quadcheck(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,&
                                                                                                  & Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ)
             ENDIF
     ELSE
             print*,'---------------------------------------------------------'
             print*,'CALCULATION FAILED TO CONVERGE WITHIN ,',Tol,'au'
             print*,'           UQUANTCHEM ABORTS!                    '
             print*,'---------------------------------------------------------'
             STOP
     ENDIF
     IF ( POUT ) THEN
        OPEN(21,FILE='UHFEIGENVALUES.dat',ACTION='WRITE')
        DO I=1,NB
                WRITE(21,*)I,EHFeigenup(I),EHFeigendown(I)
        ENDDO
        CLOSE(21)
     ENDIF
    ! Here we make sure that the output equals the 
    ! free-energy.
    IF ( ETEMP .GT. 0.0d0 ) ETOT = FTOT
END SUBROUTINE DFTCHI

