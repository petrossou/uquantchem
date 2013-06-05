SUBROUTINE forces(CORRLEVEL,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsv,S,H0,Intsv,PULAY,NTOTALQUAD,Q1,Q2,Q3,LORDER,CGORDER,LQ,CGQ,force)
      USE datatypemodule
      IMPLICIT NONE
      CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
      INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
      INTEGER, INTENT(IN) :: NATOMS,Ne,NB,NRED,PULAY
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      DOUBLE PRECISION, INTENT(IN) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB),EHFeigenup(NB),EHFeigendown(NB),H0(NB,NB),S(NB,NB)
      DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
      DOUBLE PRECISION, INTENT(IN) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsv(NATOMS,3,NRED),Intsv(NRED)
      DOUBLE PRECISION :: PT(NB,NB),Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),Cu(NB,NB),Cd(NB,NB)
      DOUBLE PRECISION :: Jupgrad(NB,NB),Jdowngrad(NB,NB),Kupgrad(NB,NB),Kdowngrad(NB,NB),xcforce(NATOMS,3),dPdR(NB,NB),sumforce(3)
      DOUBLE PRECISION, INTENT(OUT) :: force(NATOMS,3)
      DOUBLE PRECISION :: Fup(NB,NB),Fdown(NB,NB),SI(NB,NB)
      INTEGER, EXTERNAL :: ijkl
      DOUBLE PRECISION, EXTERNAL :: TRACE
      INTEGER :: Neup, Nedown,I,J,K,N

      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2

        
      Cu(:,:) = 0.0d0
      Cd(:,:) = 0.0d0
      force = 0.0d0

      DO I=1,Neup
                Cu(:,I) = Cup(:,I)
      ENDDO

      DO I=1,Nedown
                Cd(:,I) = Cdown(:,I)
      ENDDO
      
      PT = Pup + Pdown
      
      CALL getJv(Pup,NB,NRED,Intsv,Jup)
      CALL getJv(Pdown,NB,NRED,Intsv,Jdown)
      CALL getKv(Pup,NB,NRED,Intsv,Kup)
      CALL getKv(Pdown,NB,NRED,Intsv,Kdown)

      DO N=1,NATOMS
                DO J=1,3

                        CALL getJv(Pup,NB,NRED,gradIntsv(N,J,:),Jupgrad)
                        CALL getJv(Pdown,NB,NRED,gradIntsv(N,J,:),Jdowngrad)
                        CALL getKv(Pup,NB,NRED,gradIntsv(N,J,:),Kupgrad)
                        CALL getKv(Pdown,NB,NRED,gradIntsv(N,J,:),Kdowngrad)
                        IF ( CORRLEVEL .EQ. 'RHF' .OR. CORRLEVEL .EQ. 'URHF' ) THEN
                                !==========================================================================
                                ! Here different ways of calculating the Pulay-force is listed and will be 
                                ! choosen by the value of the input integer PULAY
                                !==========================================================================
                                SELECT CASE (PULAY)
                                        CASE (1)
                                                ! The way it is done in "The Cook Book", page 731, or in Int. J. Quantum Chemistry: Quantum Chemistry 
                                                ! symposium 13, 225-241, see eqn (21) p. 229, or see my notes (the black note book) p. 110, eqn (14)
                                                DO I=1,Neup
                                                        force(N,J) = force(N,J) + EHFeigenup(I)*DOT_PRODUCT(Cu(:,I),MATMUL(gradS(N,J,:,:),Cu(:,I)))
                                                ENDDO
                                                DO I=1,Nedown
                                                        force(N,J) = force(N,J) + EHFeigendown(I)*DOT_PRODUCT(Cd(:,I),MATMUL(gradS(N,J,:,:),Cd(:,I)))
                                                ENDDO
                                        CASE (2)
                                                ! The way it is derived in my notes (the black note book) p. 118-119, eqn (6).
                                                IF ( N .EQ.1 .AND. J .EQ. 1 ) THEN
                                                        Fup   = H0 + Jdown - Kup   + Jup
                                                        Fdown = H0 + Jdown - Kdown + Jup
                                                        CALL invert(S,SI,NB)
                                                ENDIF
                                                force(N,J) = force(N,J) + TRACE(MATMUL(Pup,MATMUL(Fup,MATMUL(SI,gradS(N,J,:,:)))),NB)
                                                force(N,J) = force(N,J) + TRACE(MATMUL(Pdown,MATMUL(Fdown,MATMUL(SI,gradS(N,J,:,:)))),NB)
                                        CASE (3)
                                                ! The way it is derived in my notes (the black note book) p. 118-120, eqn (7), also the way 
                                                ! A. M Niklasson calculates the Pulay  force in Phys. Rev. B. 86, 174308 (2012), see Eqn (A.13)
                                                IF ( N .EQ.1 .AND. J .EQ. 1 ) THEN
                                                        Fup   = H0 + Jdown - Kup   + Jup
                                                        Fdown = H0 + Jdown - Kdown + Jup
                                                        CALL invert(S,SI,NB)
                                                ENDIF
                                                force(N,J) = force(N,J) + 0.50d0*TRACE( MATMUL(MATMUL(SI,MATMUL(Fup,Pup))     + MATMUL(Pup,MATMUL(Fup,SI)),gradS(N,J,:,:)),NB)
                                                force(N,J) = force(N,J) + 0.50d0*TRACE( MATMUL(MATMUL(SI,MATMUL(Fdown,Pdown)) + MATMUL(Pdown,MATMUL(Fdown,SI)),gradS(N,J,:,:)),NB)
                                        CASE (4)
                                                ! Here using the idempotency of the density matrix PSP = P.
                                                ! The way it is derived in my notes (the  black note book) p. 121-122, eqn (13), p. 122. Or equivalently in
                                                ! Theor. Chem. Acc. 103, 294-296, here see eqn 3-7.
                                                IF ( N .EQ.1 .AND. J .EQ. 1 ) THEN
                                                        Fup   = H0 + Jdown - Kup   + Jup
                                                        Fdown = H0 + Jdown - Kdown + Jup
                                                ENDIF
                                                force(N,J) = force(N,J) + TRACE( MATMUL(Fup,MATMUL(Pup,MATMUL(gradS(N,J,:,:),Pup))),NB)
                                                force(N,J) = force(N,J) + TRACE( MATMUL(Fdown,MATMUL(Pdown,MATMUL(gradS(N,J,:,:),Pdown))),NB)
                                        CASE DEFAULT
                                                IF ( N .EQ.1 .AND. J .EQ. 1 ) THEN
                                                        Fup   = H0 + Jdown - Kup   + Jup
                                                        Fdown = H0 + Jdown - Kdown + Jup
                                                ENDIF
                                                force(N,J) = force(N,J) + TRACE( MATMUL(Fup,MATMUL(Pup,MATMUL(gradS(N,J,:,:),Pup))),NB)
                                                force(N,J) = force(N,J) + TRACE( MATMUL(Fdown,MATMUL(Pdown,MATMUL(gradS(N,J,:,:),Pdown))),NB)
                                END SELECT
                                !===========================================================================================================

                                force(N,J) = force(N,J) - SUM(gradT(N,J,:,:)*PT) - SUM(gradV(N,J,:,:)*PT)
                                force(N,J) = force(N,J) - 0.50d0*SUM( (Jupgrad+Jdowngrad-Kupgrad)*Pup ) - 0.50d0*SUM( (Jupgrad+Jdowngrad-Kdowngrad)*Pdown )

                        ELSE IF (  CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN
                                
                                ! Here we calculate the exchange correlation force
                                IF ( N .EQ.1 .AND. J .EQ. 1  ) THEN
                                        CALL getxcforce(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,xcforce)
                                        Fup   = H0 + Jdown + Jup
                                        Fdown = H0 + Jdown + Jup
                                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                                Fup   = Fup   - 0.20d0*Kup
                                                Fdown = Fdown - 0.20d0*Kdown
                                        ENDIF
                                ENDIF

                                force(N,J) = force(N,J) + TRACE( MATMUL(Fup,MATMUL(Pup,MATMUL(gradS(N,J,:,:),Pup))),NB)
                                force(N,J) = force(N,J) + TRACE( MATMUL(Fdown,MATMUL(Pdown,MATMUL(gradS(N,J,:,:),Pdown))),NB)

                                force(N,J) = force(N,J) - SUM(gradT(N,J,:,:)*PT) - SUM(gradV(N,J,:,:)*PT)

                                IF ( CORRLEVEL .EQ. 'B3LYP'  ) THEN
                                        force(N,J) = force(N,J) - 0.50d0*SUM( (Jupgrad+Jdowngrad-0.20*Kupgrad)*Pup ) - 0.50d0*SUM( (Jupgrad+Jdowngrad-0.20*Kdowngrad)*Pdown )
                                ELSE
                                        force(N,J) = force(N,J) - 0.50d0*SUM( (Jupgrad+Jdowngrad)*Pup ) - 0.50d0*SUM( (Jupgrad+Jdowngrad)*Pdown )
                                ENDIF
                                
                                ! Adding the exchange correlation force to the total force
                                force(N,J) = force(N,J) + xcforce(N,J)
                        ENDIF

                ENDDO
                
                DO K=1,NATOMS
                        IF ( K .NE. N ) THEN
                                force(N,:) = force(N,:) + (ATOMS(N)%R - ATOMS(K)%R)*ATOMS(N)%Z*ATOMS(K)%Z/( sqrt( DOT_PRODUCT(ATOMS(N)%R - ATOMS(K)%R,ATOMS(N)%R - ATOMS(K)%R) )**3 )
                        ENDIF
                ENDDO
      ENDDO

      ! Using translational invariance, i.e SUM{FORCES} = 0
      ! to supress translational drift due to nummeerical error
      sumforce = 0.0d0
      DO I=1,NATOMS
                sumforce = sumforce + force(I,:)
      ENDDO
      ! Enforcing translational invariance
      sumforce = sumforce/NATOMS
      DO I=1,NATOMS
                force(I,:) = force(I,:) - sumforce
      ENDDO

      END SUBROUTINE forces
