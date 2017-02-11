SUBROUTINE relax(gradS,gradT,gradV,S,H0,NB,NRED,Ne,nucE,Tol,FTol,MIX,DIISORD,DIISSTART,NATOMS,NSTEPS,NLSPOINTS,PORDER,DR,BAS,ATOMS,APPROXEE,CORRLEVEL,PRYSR,PRYSW, &
& Istarts,Iends,numprocessors,id,PULAY,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,EETOL,ETEMP,EDIR,&
& EPROFILE,EFIELDMAX,ADEF,BASAUX,ATOMSAUX,NBAUX,VRI,WRI,gradVRI,gradWRI,RIAPPROX)
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        INTEGER, INTENT(IN) :: NTOTALQUAD,LORDER,CGORDER,Qstart,Qend,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend),EDIR
        INTEGER, INTENT(IN) :: NB,Ne,NATOMS,NSTEPS,DIISORD,DIISSTART,numprocessors,id,NLSPOINTS,PORDER,PULAY,NBAUX
        INTEGER*8, INTENT(IN) :: Istarts,Iends,NRED
        INTEGER, ALLOCATABLE :: IND1(:),IND2(:),IND3(:),IND4(:)
        LOGICAL, INTENT(IN) :: APPROXEE,ADEF,RIAPPROX
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL,EPROFILE
        DOUBLE PRECISION, INTENT(INOUT) :: S(NB,NB), H0(NB,NB)
        TYPE(ATOM), INTENT(INOUT) :: ATOMS(NATOMS),ATOMSAUX(NATOMS)
        TYPE(BASIS), INTENT(INOUT)  :: BAS,BASAUX
        DOUBLE PRECISION, INTENT(IN) :: MIX,Tol,PRYSR(25,25),PRYSW(25,25),FTol,LQ(LORDER,3),CGQ(CGORDER,2),EETOL,ETEMP,EFIELDMAX
        DOUBLE PRECISION, INTENT(INOUT)  :: DR,VRI(NBAUX,NBAUX),WRI(NB,NB,NBAUX),gradVRI(NATOMS,3,NBAUX,NBAUX),gradWRI(NATOMS,3,NB,NB,NBAUX)
        DOUBLE PRECISION :: ETOT, EHFeigen(NB),EIGENVECT(NB,NB),EHFeigenup(NB),EHFeigendown(NB),Cup(NB,NB),Cdown(NB,NB),DE,EOLD
        DOUBLE PRECISION :: Pup(NB,NB),Pdown(NB,NB),P(NB,NB),DMAT(NB,NB),OMEGA,TIME
        DOUBLE PRECISION, INTENT(INOUT) :: gradS(NATOMS,3,NB,NB),gradT(NATOMS,3,NB,NB),gradV(NATOMS,3,NB,NB),NucE
        DOUBLE PRECISION :: leng, f(3),g(3),Rn(3),force(NATOMS,3),T(NB,NB),V(NB,NB),DR0,DR1,DR2,R0(NATOMS,3),R1(NATOMS,3),R2(NATOMS,3),DRMOLD,NucFieldE
        DOUBLE PRECISION :: forceold(NATOMS,3),gradold(NATOMS,3),grad(NATOMS,3),E0,E1,E2,E3,E4,const1,const2,DRM,nom,normgrad,mu,ENTROPY,DPTENSORT(3,NB,NB),DPTENSOR(NB,NB)
        DOUBLE PRECISION, ALLOCATABLE :: dInts(:),IntsvR(:),gradIntsvR(:,:,:)
        LOGICAL :: CFORCE,SCRATCH,LSEARCH,ST,SLUTA
        INTEGER :: I,J,II,JJ,KK,NLITER,JSTART,JEND,NPOINTS,POLYORDER,NSCF,ierr
        INTEGER*8 :: NONZERO,TOTALNONZERO,NONZEROO,Istart,Iend,NDIAG
        DOUBLE PRECISION, ALLOCATABLE :: ENERGY(:),DRF(:)

        OMEGA = 0.0d0
        TIME = 0.0d0

        CFORCE = .TRUE.
        I = 0
        DE = 2*Tol
        leng = 20*FTol
        EOLD = 0.0d0
        SCRATCH = .TRUE.
        NLITER = 10
        ST = .TRUE.
        NPOINTS = NLSPOINTS
        POLYORDER = PORDER
        Pup = 0.0d0
        Pdown = 0.0d0

        NDIAG = ( NB*(NB+1) )/2
        ALLOCATE(dInts(NDIAG))

        !======================================================================
        ! Calculating the ee-tensor (ij|kl)
        !=====================================================================================================================================================================
        DMAT = 1.0d0
        NONZERO = 1
        Istart = Istarts
        Iend   = Iends
        ALLOCATE(IND1(1),IND2(1),IND3(1),IND4(1))
        ! Counting the number of "non-zero" elements og the ee-tensor (ij|kl), elements satisfying |(ij|kl)| < EETOL 
        IF ( RIAPPROX ) THEN
                CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.TRUE.,id,numprocessors)
                CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.TRUE.,id,numprocessors)

                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
        ELSE
                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts)
        ENDIF
        DEALLOCATE(IND1,IND2,IND3,IND4)
        IF ( NONZERO .EQ. 0 ) THEN
            NONZEROO = NONZERO +1
        ELSE
            NONZEROO = NONZERO
        ENDIF
        ALLOCATE(IND1(NONZEROO),IND2(NONZEROO),IND3(NONZEROO),IND4(NONZEROO))
        ! creating a mapping from the non-zero index running from 1 to NONZERO on each thread to the contracted index ( see the routine ijkl.f90 for definition )
        IF ( RIAPPROX ) THEN
                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
        ELSE
                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,id,DMAT,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts)
        ENDIF
        Istart = 1
        Iend = NONZEROO
        ALLOCATE(IntsvR(Istart:Iend),gradIntsvR(NATOMS,3,Istart:Iend))
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        IF ( RIAPPROX ) THEN
                CALL eeintsRI(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,&
                & PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,.TRUE.,id,NONZERO,DMAT,dInts,NBAUX,VRI,WRI,gradVRI,gradWRI)
        ELSE
                CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,&
                & PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,.TRUE.,id,NONZERO,DMAT,dInts)
        ENDIF
        !=======================================================================================================================================================================


        ALLOCATE(ENERGY(NPOINTS),DRF(NPOINTS))

        DO WHILE ( DABS(leng) .GT. FTol .AND. I .LT. NSTEPS ) 

                !Switching to conjugate gradient:
                IF ( DABS(leng) .LT. 10.0d0*FTOL ) ST = .FALSE.
     
        !=======================================================================================================================
                
                DO J=1,NATOMS
                        R0(J,:) = ATOMS(J)%R
                ENDDO
  
                IF ( CORRLEVEL .EQ. 'RHF' ) THEN
                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                        EHFeigenup = EHFeigen
                        EHFeigendown = EHFeigen
                        Cup = EIGENVECT
                        Cdown = EIGENVECT
                        Pup = P/2.0d0
                        Pdown = P/2.0d0
                ENDIF
                
                IF ( CORRLEVEL .EQ. 'URHF' ) THEN
                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX, &
                                  & DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY)
                ENDIF
                IF ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ. 'LDA' .OR.  CORRLEVEL .EQ. 'B3LYP' ) THEN
                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,Tol,EHFeigenup,EHFeigendown,&
                        & ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY)
                ENDIF

                ! Calculating forces on atoms:
                CALL forces(IND1,IND2,IND3,IND4,Istart,Iend,NATOMS,Ne,NB,NRED,Cup,Cdown,Pup,Pdown,EHFeigenup,EHFeigendown,ATOMS,BAS,gradS,gradT,gradV,gradIntsvR,S,H0,IntsvR,PULAY,&
                & force,numprocessors,id,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,LORDER,CGORDER,LQ,CGQ,CORRLEVEL,EDIR,EPROFILE,EFIELDMAX,OMEGA,TIME,ADEF)

                leng = 0.0d0
                DO II=1,NATOMS
                        leng = leng + sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                ENDDO
                
                leng = leng/NATOMS
        !=======================================================================================================================
                EOLD = E0
                E0 = ETOT
                DE = ETOT-EOLD

                IF ( I .EQ. 0 ) THEN
                        grad = -force
                ELSE
                        const1 = 0.0d0
                        const2 = 0.0d0
                        DO J=1,NATOMS
                                const1 = const1 + DOT_PRODUCT(force(J,:)-forceold(J,:),force(J,:))
                                const2 = const2 + DOT_PRODUCT(forceold(J,:),forceold(J,:))
                        ENDDO
                        ! Uppdating the gradient according to Polack and Ribiere
                        DO J=1,NATOMS
                                IF ( ST ) THEN
                                        ! steepest deecent
                                        grad(J,:) = force(J,:) 
                                ELSE
                                        ! Conjugate gradient
                                        grad(J,:) = force(J,:) + (const1/const2)*gradold(J,:)
                                ENDIF
                        ENDDO
                ENDIF
                ! Calculating the norm of the gradient:
                normgrad = 0.0d0
                DO J=1,NATOMS
                        normgrad = normgrad + DOT_PRODUCT(grad(J,:),grad(J,:))
                ENDDO
                normgrad = (1.0d0/(1.0d0*NATOMS))*sqrt(normgrad)

                forceold = force
                gradold  = grad

        !=========================================================================================================================================
        ! Here we search for the energy minimum along the gradient grad.
        !=========================================================================================================================================
                
                
                LSEARCH = .TRUE.
                SLUTA = .FALSE.
                ENERGY(1) = E0
                DRF(1) = 0.0d0

                DO WHILE ( LSEARCH .AND. leng .GT. FTOL )
                        JSTART = 2
                        JEND   = NPOINTS

                        IF ( SLUTA ) LSEARCH = .FALSE.

                        IF ( .not. LSEARCH ) THEN
                                JSTART = 1
                                JEND   = 1
                                DRM = DABS(DRM) ! After all the step should be in the direction of the gradient
                                IF (  I .GT. 1  ) DR = 2.0d0*DRM
                                IF (  DABS(DR) .LT. 1.0E-4 ) DR = 1.0E-4
                        ENDIF

                        DO J=JSTART,JEND
                                ! Calculating new positions on atoms:
                                DO JJ=1,NATOMS
                                        g = grad(JJ,:)/normgrad
                                        IF ( LSEARCH ) THEN
                                                DRF(J) = (J-1)*DR/(NPOINTS-1)
                                                ATOMS(JJ)%R = R0(JJ,:) + g*(J-1)*DRF(J)
                                                IF ( RIAPPROX ) ATOMSAUX(JJ)%R = ATOMS(JJ)%R
                                        ELSE
                                                ATOMS(JJ)%R = R0(JJ,:) + g*DRM
                                                IF ( RIAPPROX ) ATOMSAUX(JJ)%R = ATOMS(JJ)%R
                                        ENDIF
                                ENDDO

                                ! Calculating the nuclear-nuclear repulsion energy of the
                                ! updated nuclear configuration
                                NucE = 0.0d0
                                DO II=1,NATOMS
                                        DO JJ=II+1,NATOMS
                                                Rn = ATOMS(II)%R - ATOMS(JJ)%R
                                                NucE = NucE + ATOMS(II)%Z*ATOMS(JJ)%Z/sqrt(DOT_PRODUCT(Rn,Rn))
                                        ENDDO
                                ENDDO
                                
                                ! Adding the potential energy emanating from external E-field and the atomic nuclea
                                        
                                IF ( ADEF ) THEN
                                        NucFieldE = 0.0d0
                                        DO KK=1,NATOMS
                                                IF ( EDIR .EQ. 1 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(2)
                                                ELSE IF ( EDIR .EQ. 2 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(3)
                                                ELSE IF ( EDIR .EQ. 3 ) THEN
                                                        NucFieldE = NucFieldE - EFIELDMAX*ATOMS(KK)%Z*ATOMS(KK)%R(1)
                                                ENDIF
                                        ENDDO
                                        NucE = NucE + NucFieldE
                                ENDIF

                                ! Updating the basis set:
                                CALL gettotalbasis(NATOMS,ATOMS,NB,BAS,.FALSE.)

                                ! Here the normalization of the basis functions is performed
                                CALL normalize(BAS)

                                IF ( RIAPPROX) THEN
                                        ! Updating the aux-basis set:
                                        CALL gettotalbasis(NATOMS,ATOMSAUX,NBAUX,BASAUX,.FALSE.)
                                        ! Here the normalization of the  aux-basis functions is
                                        ! performed
                                        CALL normalize(BASAUX)
                                ENDIF

                
                                ! Calculating the matrices in the basis of the updated positions:
                                CALL overlap(NATOMS,BAS,S,gradS,BAS%NBAS,.TRUE.)
                                CALL kinetic(NATOMS,BAS,T,gradT,BAS%NBAS,.TRUE.)
                                CALL potential(BAS,NATOMS,ATOMS,V,gradV,BAS%NBAS,.TRUE.,id,numprocessors)
                
                                H0 = T + V

                                IF ( ADEF ) THEN
                                        CALL dipoletensor(BAS,DPTENSORT)
                                        IF ( EDIR .EQ. 1 ) DPTENSOR = DPTENSORT(2,:,:)
                                        IF ( EDIR .EQ. 2 ) DPTENSOR = DPTENSORT(3,:,:)
                                        IF ( EDIR .EQ. 3 ) DPTENSOR = DPTENSORT(1,:,:)
                                        H0 = H0 + DPTENSOR*EFIELDMAX
                                ENDIF

                                ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                IF ( LSEARCH ) THEN
                                      !=================================================================================================
                                      ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                      !============================================================================================================================================
                                      NONZERO = 1
                                      Istart = Istarts
                                      Iend   = Iends
                                      DEALLOCATE(IND1,IND2,IND3,IND4)
                                      ALLOCATE(IND1(1),IND2(1),IND3(1),IND4(1))
                                      ! Counting the number of "non-zero" elements og the ee-tensor (ij|kl), elements satisfying |(ij|kl)| < EETOL 
                                      IF ( RIAPPROX ) THEN
                                                CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.FALSE.,id,numprocessors)
                                                CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.FALSE.,id,numprocessors)
                                                
                                                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
                                      ELSE
                                                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts)
                                      ENDIF
                                      DEALLOCATE(IND1,IND2,IND3,IND4)
                                      IF ( NONZERO .EQ. 0 ) THEN
                                          NONZEROO = NONZERO +1
                                      ELSE
                                          NONZEROO = NONZERO
                                      ENDIF
                                      ALLOCATE(IND1(NONZEROO),IND2(NONZEROO),IND3(NONZEROO),IND4(NONZEROO))
                                      ! creating a mapping from the non-zero index running from 1 to NONZERO on each thread to the contracted index
                                      ! ( see the routine ijkl.f90 for definition )
                                      IF ( RIAPPROX ) THEN
                                                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
                                      ELSE
                                                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,& 
                                                & id,Pup+Pdown,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts)
                                      ENDIF
                                      Istart = 1
                                      Iend = NONZEROO
                                      DEALLOCATE(IntsvR,gradIntsvR)
                                      ALLOCATE(IntsvR(Istart:Iend),gradIntsvR(NATOMS,3,Istart:Iend))
                                      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                                      IF ( RIAPPROX ) THEN
                                                CALL eeintsRI(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,& 
                                                  & APPROXEE,EETOL,.FALSE.,id,NONZERO,Pup+Pdown,dInts,NBAUX,VRI,WRI,gradVRI,gradWRI)
                                      ELSE
                                                CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,& 
                                                  & APPROXEE,EETOL,.FALSE.,id,NONZERO,Pup+Pdown,dInts)
                                      ENDIF
                                      !============================================================================================================================================== 
                                ELSE
                                      !=================================================================================================
                                      ! Calculation of electron repulsion tensor (ij|kl) in the basis of the updated positions:
                                      !============================================================================================================================================       
                                      NONZERO = 1
                                      Istart = Istarts
                                      Iend   = Iends
                                      DEALLOCATE(IND1,IND2,IND3,IND4)
                                      ALLOCATE(IND1(1),IND2(1),IND3(1),IND4(1))
                                      ! Counting the number of "non-zero" elements og the ee-tensor (ij|kl), elements satisfying |(ij|kl)| < EETOL 
                                      IF ( RIAPPROX ) THEN
                                                CALL calcWRI(NATOMS,NATOMS,BAS,BASAUX,WRI,gradWRI,PRYSR,PRYSW,.TRUE.,id,numprocessors)
                                                CALL calcVRI(NATOMS,NATOMS,BASAUX,VRI,gradVRI,PRYSR,PRYSW,.TRUE.,id,numprocessors)
                                                
                                                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
                                      ELSE
                                                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.FALSE.,NDIAG,NONZERO,TOTALNONZERO,dInts)
                                      ENDIF
                                      DEALLOCATE(IND1,IND2,IND3,IND4)
                                      IF ( NONZERO .EQ. 0 ) THEN
                                          NONZEROO = NONZERO +1
                                      ELSE
                                          NONZEROO = NONZERO
                                      ENDIF
                                      ALLOCATE(IND1(NONZEROO),IND2(NONZEROO),IND3(NONZEROO),IND4(NONZEROO))
                                      ! creating a mapping from the non-zero index running from 1 to NONZERO on each thread 
                                      ! to the contracted index ( see the routine ijkl.f90 for definition )
                                      IF ( RIAPPROX ) THEN
                                                CALL counteeRI(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts,NBAUX,VRI,WRI)
                                      ELSE
                                                CALL countee(BAS,IND1,IND2,IND3,IND4,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,EETOL,&
                                                & id,Pup+Pdown,.TRUE.,NDIAG,NONZEROO,TOTALNONZERO,dInts)
                                      ENDIF
                                      Istart = 1
                                      Iend = NONZEROO
                                      DEALLOCATE(IntsvR,gradIntsvR)
                                      ALLOCATE(IntsvR(Istart:Iend),gradIntsvR(NATOMS,3,Istart:Iend))
                                      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                                      IF ( RIAPPROX ) THEN
                                                CALL eeintsRI(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,&  
                                                  & APPROXEE,EETOL,.TRUE.,id,NONZERO,Pup+Pdown,dInts,NBAUX,VRI,WRI,gradVRI,gradWRI)
                                      ELSE
                                                CALL eeints(NATOMS,NATOMS,BAS,IntsvR,gradIntsvR,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istart,Iend,PRYSR,PRYSW,numprocessors,&  
                                                  & APPROXEE,EETOL,.TRUE.,id,NONZERO,Pup+Pdown,dInts)
                                      ENDIF
                                      !============================================================================================================================================== 
                               ENDIF

                                IF ( CORRLEVEL .EQ. 'RHF' .AND. LSEARCH ) THEN
                                        CALL RHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigen,ETOT,EIGENVECT,P,MIX,DIISORD,DIISSTART,NSCF,-1, &
                                                & numprocessors,id,.FALSE.,SCRATCH,.FALSE.)
                                        EHFeigenup = EHFeigen
                                        EHFeigendown = EHFeigen
                                        Cup = EIGENVECT
                                        Cdown = EIGENVECT
                                        Pup = P/2.0d0
                                        Pdown = P/2.0d0
                                ENDIF
                
                                IF ( CORRLEVEL .EQ. 'URHF' .AND. LSEARCH ) THEN
                                        CALL URHF(S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,NB,NRED,Ne,nucE,Tol/100,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1, & 
                                        & numprocessors,id,.FALSE.,SCRATCH,.FALSE.,ETEMP,ENTROPY)
                                        !IF ( id .EQ. 0 ) print*,'---->',ETOT
                                ENDIF
                                
                                IF ( ( CORRLEVEL .EQ. 'PBE' .OR. CORRLEVEL .EQ.  'LDA' .OR. CORRLEVEL .EQ. 'B3LYP' ) .AND. LSEARCH ) THEN
                                        CALL DFT(CORRLEVEL,NATOMS,ATOMS,BAS,S,H0,IntsvR,IND1,IND2,IND3,IND4,Istart,Iend,Q1,Q2,Q3,Qstart,Qend,NTOTALQUAD,NB,NRED,Ne,LORDER,CGORDER,LQ,CGQ,nucE,&
                                        & Tol,EHFeigenup,EHFeigendown,ETOT,Cup,Cdown,Pup,Pdown,MIX,DIISORD,DIISSTART,NSCF,-1,numprocessors,id,.FALSE.,SCRATCH,.FALSE.,ETEMP,mu,ENTROPY)
                                ENDIF

                                IF ( LSEARCH ) ENERGY(J) = ETOT
                       
                                ! Calculating the new position in the line-searc for the minimum
                                ! by fitting the NPOINTS energies to a  POLYORDER:th  order polynomial,
                                ! and by minimizing this polynomial in the  intervall [0,DR]
                                IF ( J .EQ. NPOINTS .AND. LSEARCH ) THEN
                                        DRMOLD = DRM
                                        CALL linesearchmin(POLYORDER,ENERGY,DRF,NPOINTS,DR,DRM)
                                        ! Safety cutoff:
                                        IF ( DRM .LT. 1.0E-10 ) DRM = DRMOLD
                                        IF ( DABS(DRM) .GT. DR ) DRM = DR
                                        SLUTA = .TRUE.
                                ENDIF

                          ENDDO

                          IF ( id .EQ. 0 ) THEN
                                  !-------------------------------------------------
                                  ! Saving the computed new atomic positions to file
                                  !-------------------------------------------------
                                  OPEN(16,FILE='ATOMPOSITIONS.dat',ACTION='WRITE')
                                  DO JJ=1,NATOMS
                                        WRITE(16,'(A4,I4,3(F15.10))')'ATOM',ATOMS(JJ)%Z,ATOMS(JJ)%R
                                  ENDDO
                                  CLOSE(16)
                                  !---------------------------------------------------------------
                          ENDIF
                          !================================================================================
                          ! Here we make sure the the mean length of the force
                          ! is consistent over all threads to avoid deadlock
                          !================================================================================
                          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                          CALL MPI_BCAST(leng,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
                          CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
                          !================================================================================
                  ENDDO
          !=========================================================================================================================================

               IF ( I .EQ. 0 .AND. id .EQ. 0 ) THEN 
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,'(A79)')'                             Relaxing the nuclear positions:                   '
                        WRITE(*,'(A77)')'            ================================================================='
                        WRITE(*,*)' '
                        WRITE(*,'(A114)')'  ITER           E [au]                        DE [au]                     <|F|> [au]                      DR [au]'
                ENDIF
                
                !leng = 0.0d0
                !DO II=1,NATOMS
                !        leng = sqrt(DOT_PRODUCT(force(II,:),force(II,:)))
                !ENDDO
                
                !leng = leng/NATOMS
                
                IF ( id .EQ. 0 ) WRITE(*,'(I4,E30.20,E30.20,E30.20,E30.20)')I,E0,DE,leng,DRM

                I = I + 1
                SCRATCH = .FALSE.

      ENDDO
        
      IF ( leng .LT. FTol .AND. id .EQ. 0 ) THEN
              WRITE(*,*)' '
              WRITE(*,'(A69)')'                 The nuclear positions have been relaxed with respect'
              WRITE(*,'(A68)')'                 to the interatomic forces. An force convergence of '
              WRITE(*,'(A27,E12.5,A24)')' F  < ',FTol,' [au] has been obtained.'
              WRITE(*,*)' ' 
              WRITE(*,'(A70)')'                The relaxed positions of the nuclea are the following:'
              WRITE(*,*)' '
              WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
              DO J=1,NATOMS
                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
              ENDDO
      ELSE
              IF ( id .EQ. 0 ) THEN
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'          Failure to relax the nuclear positions within prescribed force tolerance '
                        WRITE(*,*)' ' 
                        WRITE(*,'(A71)')'                     The last positions of the nuclea are the following:'
                        WRITE(*,*)' '
                        WRITE(*,'(A84)')'  ATOM           Rx [au]                       Ry [au]                       Rz [au]'
                        DO J=1,NATOMS
                                WRITE(*,'(I4,E30.20,E30.20,E30.20)')J,ATOMS(J)%R
                        ENDDO
                ENDIF
      ENDIF

      END SUBROUTINE relax
