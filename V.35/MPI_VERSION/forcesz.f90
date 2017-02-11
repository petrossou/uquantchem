SUBROUTINE forcesz(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,PNup,PNdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR &
& ,S,H0,Intsv,PULAY,force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL,Jup,Jdown,Kup,Kdown,Vxc,EDIR,EPROFILE,EFIELD,OMEGA,TIME,ADEF)
      USE datatypemodule
      IMPLICIT NONE
      INCLUDE "mpif.h"
      CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL,EPROFILE
      LOGICAL, INTENT(IN) :: ADEF
      INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Qstart,Qend,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend)
      INTEGER, INTENT(IN) :: NATOMS,Ne,NB,numprocessors,id,PULAY,EDIR
      INTEGER*8, INTENT(IN) :: Istart,Iend,NRED
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      TYPE(BASIS), INTENT(IN)  :: BAS
      DOUBLE PRECISION, INTENT(IN) :: Cup(NB,NB),Cdown(NB,NB),Pup(NB,NB),Pdown(NB,NB),EHFeigenup(NB),EHFeigendown(NB),H0(NB,NB),S(NB,NB)
      DOUBLE PRECISION, INTENT(IN) :: PNup(NB,NB),PNdown(NB,NB),EFIELD,OMEGA,TIME
      DOUBLE PRECISION, INTENT(IN) :: Jup(NB,NB),Jdown(NB,NB),Kup(NB,NB),Kdown(NB,NB),Vxc(2,NB,NB)
      DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
      DOUBLE PRECISION, INTENT(IN) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),gradIntsvR(NATOMS,3,Istart:Iend),Intsv(Istart:Iend)
      INTEGER, INTENT(IN) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION :: PT(NB,NB),Jupgrad(NB,NB),Jdowngrad(NB,NB),Kupgrad(NB,NB),Kdowngrad(NB,NB),Cu(NB,NB),Cd(NB,NB),gDPTENSOR(NATOMS,3,NB,NB),DPTENSOR(NB,NB),DPTENSORT(3,NB,NB)
      DOUBLE PRECISION :: Fup(NB,NB),Fdown(NB,NB),SI(Nb,NB),xcforce(NATOMS,3),sumforce(3),cl
      DOUBLE PRECISION, INTENT(OUT) :: force(NATOMS,3)
      INTEGER, EXTERNAL :: ijkl
      DOUBLE PRECISION, EXTERNAL :: TRACE
      INTEGER :: Neup, Nedown,I,J,K,N,ierr

      cl = 137.0359990740d0

      Neup   = ( Ne - MOD(Ne,2) )/2
      Nedown = ( Ne + MOD(Ne,2) )/2

      ! If a electric field is present we here calculate the gradient of the dipole tensor:
      IF ( ADEF ) THEN
              IF ( EPROFILE .EQ. 'AC' ) THEN
                      CALL graddipoletensorinhom(BAS,NATOMS,EDIR,omega,TIME,gDPTENSOR)
                      CALL dipoletensorinhom(BAS,EDIR,omega,TIME,DPTENSOR)
              ELSE
                      CALL dipoletensor(BAS,DPTENSORT)
                      IF ( EDIR .EQ. 1 ) THEN
                              CALL graddipoletensor(NATOMS,BAS,2,gDPTENSOR)
                              DPTENSOR = DPTENSORT(2,:,:)
                      ENDIF
                      IF ( EDIR .EQ. 2 ) THEN
                              CALL graddipoletensor(NATOMS,BAS,3,gDPTENSOR)
                              DPTENSOR = DPTENSORT(3,:,:)
                      ENDIF
                      IF ( EDIR .EQ. 3 ) THEN
                              CALL graddipoletensor(NATOMS,BAS,1,gDPTENSOR)
                              DPTENSOR = DPTENSORT(1,:,:)
                      ENDIF

              ENDIF
      ENDIF
        
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
      
      !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !CALL getJv(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jup)
      !CALL getJv(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Jdown)
      !CALL getKv(Pup,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kup)
      !CALL getKv(Pdown,NB,NRED,Istart,Iend,Intsv,IND1,IND2,IND3,IND4,numprocessors,id,Kdown)
      !CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

      DO N=1,NATOMS
                DO J=1,3
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                        CALL getJv(PNup,NB,NRED,Istart,Iend,gradIntsvR(N,J,:),IND1,IND2,IND3,IND4,numprocessors,id,Jupgrad)
                        CALL getJv(PNdown,NB,NRED,Istart,Iend,gradIntsvR(N,J,:),IND1,IND2,IND3,IND4,numprocessors,id,Jdowngrad)
                        CALL getKv(PNup,NB,NRED,Istart,Iend,gradIntsvR(N,J,:),IND1,IND2,IND3,IND4,numprocessors,id,Kupgrad)
                        CALL getKv(PNdown,NB,NRED,Istart,Iend,gradIntsvR(N,J,:),IND1,IND2,IND3,IND4,numprocessors,id,Kdowngrad)
                        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
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
                                                !force(N,J) = force(N,J) + TRACE(MATMUL(Pup,MATMUL(Fup,MATMUL(SI,gradS(N,J,:,:)))),NB)
                                                !force(N,J) = force(N,J) + TRACE(MATMUL(Pdown,MATMUL(Fdown,MATMUL(SI,gradS(N,J,:,:)))),NB)
                                                force(N,J) = force(N,J) + TRACE( MATMUL(MATMUL(SI,MATMUL(Fup,Pup)),gradS(N,J,:,:)),NB)
                                                force(N,J) = force(N,J) + TRACE( MATMUL(MATMUL(SI,MATMUL(Fdown,Pdown)),gradS(N,J,:,:)),NB)
                                        CASE (3)
                                                ! The way it is derived in my notes (the black note book) p. 118-120, eqn (7), also the way
                                                ! A. M Niklasson calculates the Pulay force in Phys. Rev. B. 86, 174308 (2012), see Eqn (A.13)
                                                IF ( N .EQ.1 .AND. J .EQ. 1 ) THEN
                                                        Fup   = H0 + Jdown - Kup   + Jup
                                                        Fdown = H0 + Jdown - Kdown + Jup
                                                        CALL invert(S,SI,NB)
                                                ENDIF
                                                force(N,J) = force(N,J) + 0.50d0*TRACE( MATMUL(MATMUL(SI,MATMUL(Fup,Pup))     + MATMUL(Pup,MATMUL(Fup,SI)),gradS(N,J,:,:)),NB)
                                                force(N,J) = force(N,J) + 0.50d0*TRACE( MATMUL(MATMUL(SI,MATMUL(Fdown,Pdown)) + MATMUL(Pdown,MATMUL(Fdown,SI)),gradS(N,J,:,:)),NB)
                                        CASE (4)
                                                ! Here using the idempotency of the density matrix PSP = P.
                                                ! The way it is derived in my notes (the black note book) p. 121-122, eqn (13), p. 122. Or equivalently in
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
                                force(N,J) = force(N,J) - 0.50d0*SUM((Jupgrad+Jdowngrad-Kupgrad)*(2.d0*Pup-PNup)) - 0.50d0*SUM((Jupgrad+Jdowngrad-Kdowngrad)*(2.d0*Pdown-PNdown))

                                IF ( ADEF ) THEN
                                        force(N,J) = force(N,J) - EFIELD*SUM(gDPTENSOR(N,J,:,:)*PT)
                                ENDIF

                        ELSE IF (  CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) THEN

                                ! Here we calculate the exchange correlation force
                                IF ( N .EQ.1 .AND. J .EQ. 1  ) THEN
                                        IF ( PULAY .EQ. 4 ) THEN
                                                CALL getxcforce(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,xcforce)
                                                Fup   = H0 + Jdown + Jup
                                                Fdown = H0 + Jdown + Jup
                                        ELSE
                                                ! Exc,R[2D-P]: ( Best Choice )
                                                !CALL getxcforcenipz(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,PNup,Pdown,PNdown,&
                                                !& LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,xcforce)
                                                ! Exc,R[D]:  
                                                CALL getxcforcenipz(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pup,Pdown,Pdown,&
                                                & LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,xcforce)
                                                ! Exc,R[P]: (Worst Choice )
                                                !CALL getxcforcenipz(CORRLEVEL,NATOMS,ATOMS,BAS,PNup,PNup,PNdown,PNdown,&
                                                !& LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,xcforce)

                                                !CALL getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)      
                       
                                                !CALL getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)
                                                Fup   = H0 + Jdown + Jup + Vxc(1,:,:)
                                                Fdown = H0 + Jdown + Jup + Vxc(2,:,:)
                                        ENDIF
                                        IF ( CORRLEVEL .EQ. 'B3LYP' ) THEN
                                                Fup   = Fup   - 0.20d0*Kup
                                                Fdown = Fdown - 0.20d0*Kdown
                                        ENDIF
                                ENDIF


                                !IF ( ADEF ) THEN
                                !        Fup   = Fup   + EFIELD*DPTENSOR
                                !        Fdown = Fdown + EFIELD*DPTENSOR
                                !ENDIF
                                
                                IF ( PULAY .EQ. 4 ) THEN
                                        ! Pulay force assuming idempotency
                                        force(N,J) = force(N,J) + TRACE( MATMUL(Fup,MATMUL(Pup,MATMUL(gradS(N,J,:,:),Pup))),NB)
                                        force(N,J) = force(N,J) + TRACE( MATMUL(Fdown,MATMUL(Pdown,MATMUL(gradS(N,J,:,:),Pdown))),NB)
                                ELSE
                                        ! Here we don't relay on idempotency. To be used when non-integer occupation is present, for
                                        ! instance when using Fermi-smearing of states.
                                        CALL invert(S,SI,NB)
                                        force(N,J) = force(N,J) + TRACE( MATMUL(MATMUL(SI,MATMUL(Fup,Pup)),gradS(N,J,:,:)),NB)
                                        force(N,J) = force(N,J) + TRACE( MATMUL(MATMUL(SI,MATMUL(Fdown,Pdown)),gradS(N,J,:,:)),NB)
                                ENDIF

                                force(N,J) = force(N,J) - SUM(gradT(N,J,:,:)*PT) - SUM(gradV(N,J,:,:)*PT)

                                IF ( ADEF ) THEN
                                        force(N,J) = force(N,J) - EFIELD*SUM(gDPTENSOR(N,J,:,:)*PT)
                                ENDIF

                                IF ( CORRLEVEL .EQ. 'B3LYP'  ) THEN
                                   force(N,J) = force(N,J) - 0.50d0*SUM( (Jupgrad+Jdowngrad-0.20*Kupgrad)*(2.d0*Pup-PNup) ) - 0.50d0*SUM( (Jupgrad+Jdowngrad-0.20*Kdowngrad)*(2.d0*Pdown-PNdown))
                                ELSE
                                   force(N,J) = force(N,J) - 0.50d0*SUM( (Jupgrad+Jdowngrad)*(2.d0*Pup-PNup) ) - 0.50d0*SUM( (Jupgrad+Jdowngrad)*(2.d0*Pdown-PNdown) )
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
                ! The force of the electric field on the nuclea
                IF ( ADEF ) THEN
                        IF ( EPROFILE .NE. 'AC' ) THEN
                                IF ( EDIR .EQ. 1 ) force(N,2) = force(N,2) + EFIELD*ATOMS(N)%Z
                                IF ( EDIR .EQ. 2 ) force(N,3) = force(N,3) + EFIELD*ATOMS(N)%Z
                                IF ( EDIR .EQ. 3 ) force(N,1) = force(N,1) + EFIELD*ATOMS(N)%Z
                        ELSE
                                IF ( EDIR .EQ. 1 ) force(N,2) = force(N,2) + EFIELD*sin(OMEGA*TIME + (OMEGA/cl)*ATOMS(N)%R(1))*ATOMS(N)%Z
                                IF ( EDIR .EQ. 2 ) force(N,3) = force(N,3) + EFIELD*sin(OMEGA*TIME + (OMEGA/cl)*ATOMS(N)%R(2))*ATOMS(N)%Z
                                IF ( EDIR .EQ. 3 ) force(N,1) = force(N,1) + EFIELD*sin(OMEGA*TIME + (OMEGA/cl)*ATOMS(N)%R(3))*ATOMS(N)%Z
                        ENDIF
                ENDIF
      ENDDO

      ! Using translational invariance, i.e SUM{FORCES} = 0
      ! to supress translational drift due to nummeerical error
      IF ( .not. ADEF ) THEN
                sumforce = 0.0d0
                DO I=1,NATOMS
                        sumforce = sumforce + force(I,:)
                ENDDO
                ! Enforcing translational invariance
                sumforce = sumforce/NATOMS
                DO I=1,NATOMS
                        force(I,:) = force(I,:) - sumforce
                ENDDO
      ENDIF
         

      END SUBROUTINE forcesz
