SUBROUTINE eeints(NATOMS,FNATOMS,BAS,Intsv,gradIntsv,IND1,IND2,IND3,IND4,NDIAG,Istart,Iend,Istartg,Iendg,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,CFORCE,id,NONZERO,DMAT,dInts)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: id,numprocessors,NATOMS,FNATOMS
      INTEGER*8, INTENT(IN) :: Istart,Iend,Istartg,Iendg,NONZERO,NDIAG
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol,DMAT(BAS%NBAS,BAS%NBAS),dInts(NDIAG)
      LOGICAL, INTENT(IN) :: APPROXEE,CFORCE
      INTEGER, INTENT(INOUT) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(INOUT) :: Intsv(Istart:Iend),gradIntsv(FNATOMS,3,Istartg:Iendg)
      INTEGER*8 :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,GG,STRL,KK,PPP,II
      INTEGER*8 :: rcounts(numprocessors),displs(numprocessors),ierr,scount,SI,SIC(numprocessors)
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL,DMAX
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4,gradient(4,3),DMAXIMUM(6)
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
     
      IF ( NONZERO .EQ. 0 ) THEN
        ! These are dummy entries:
           IND1(1) = 1
           IND2(1) = 1
           IND3(1) = 1
           IND4(1) = 1
      ENDIF      

     Intsv(:) = 0.0d0
     gradIntsv(:,:,:) = 0.0d0

DO GG=1,1
        IF ( NONZERO .NE. 0 ) THEN
        DO NN=Istart,Iend
                        G = NN
                        I = IND1(G) 
                        J = IND2(G) 
                        K = IND3(G) 
                        L = IND4(G) 
                        IJ = (J*(J-1))/2 + I
                        KL = (L*(L-1))/2 + K
                
                        L1= BAS%PSI(I)%L(1)
                        M1= BAS%PSI(I)%L(2)
                        N1= BAS%PSI(I)%L(3)
                        A = BAS%PSI(I)%R
                        NO1 = BAS%PSI(I)%NORM 
                
                        L2= BAS%PSI(J)%L(1)
                        M2= BAS%PSI(J)%L(2)
                        N2= BAS%PSI(J)%L(3)
                        B = BAS%PSI(J)%R
                        NO2 = BAS%PSI(J)%NORM

                        L3= BAS%PSI(K)%L(1)
                        M3= BAS%PSI(K)%L(2)
                        N3= BAS%PSI(K)%L(3)
                        C = BAS%PSI(K)%R
                        NO3 = BAS%PSI(K)%NORM 
                
                        L4= BAS%PSI(L)%L(1)
                        M4= BAS%PSI(L)%L(2)
                        N4= BAS%PSI(L)%L(3)
                        D = BAS%PSI(L)%R
                        NO4 = BAS%PSI(L)%NORM
                     
                        !------------------------
                        ! Swartz inequality part
                        !------------------------
                        ! Here we use Schwartz inequality |(i,j,k,l)| <= ! [ |(i,j,i,j)|*(k,l,k,l) ]**(1/2). 
                        ! See Eqn 9.12.25 p.404 in T. Helgaker et al's Molecular Electronic-Structure Theory
                        ! We use this together with Eqn. (56) and (57) page. 86,
                        ! in Haettig, Multiscale Simulation Methods in Molecular
                        ! Sciences, J. Grotendorst, N. Attig, S. Blugel, D. Marx (Eds.),
                        ! Institute for Advanced Simulation, Forschungszentrum Julich, NIC Series, Vol. 42, ISBN 978-3-9810843-8-2,pp.  77-120, 2009.
                        
                        DMAXIMUM = (/ 2*DABS(DMAT(I,J)),2*DABS(DMAT(K,L)),DABS(DMAT(J,K)),DABS(DMAT(I,K)),DABS(DMAT(J,L)),DABS(DMAT(I,L)) /)

                        DMAX = MAXVAL(DMAXIMUM)

                        DOSUMOVERPRIMS = .FALSE.

                        IF ( ( IJ .EQ. KL .AND. sqrt(DABS(dInts(IJ)*dInts(KL)))*DMAX .GE. Tol ) .OR. ( IJ .EQ. KL .AND. .not. APPROXEE ) )  THEN
                                  Intsv(G) = dInts(IJ)
                                  IF ( CFORCE ) DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF (  APPROXEE .AND. IJ .NE. KL .AND. sqrt(DABS(dInts(IJ)*dInts(KL)))*DMAX .GE. Tol  ) THEN
                                DOSUMOVERPRIMS = .TRUE.
                        ENDIF
                        IF ( .not. APPROXEE .AND. IJ .NE. KL ) DOSUMOVERPRIMS = .TRUE.
                        !-----------------------
                        ! End of inequality part
                        !-----------------------

                        IF ( DOSUMOVERPRIMS ) THEN
                                ! SUM OVER PRIMITIVE GAUSSIANS:
                                DO M=1,BAS%PSI(I)%NPRIM
                                        al1 = BAS%PSI(I)%EXPON(M)
                                        co1 = BAS%PSI(I)%CONTRCOEFF(M)
                                        NP1 = BAS%PSI(I)%PRIMNORM(M)
                                        DO N=1,BAS%PSI(J)%NPRIM
                                                al2 = BAS%PSI(J)%EXPON(N)
                                                co2 = BAS%PSI(J)%CONTRCOEFF(N)
                                                NP2 = BAS%PSI(J)%PRIMNORM(N)
                                                DO P=1,BAS%PSI(K)%NPRIM
                                                        al3 = BAS%PSI(K)%EXPON(P)
                                                        co3 = BAS%PSI(K)%CONTRCOEFF(P)
                                                        NP3 = BAS%PSI(K)%PRIMNORM(P)
                                                        DO Q=1,BAS%PSI(L)%NPRIM
                                                                al4 = BAS%PSI(L)%EXPON(Q)
                                                                co4 = BAS%PSI(L)%CONTRCOEFF(Q)
                                                                NP4 = BAS%PSI(L)%PRIMNORM(Q)
                                                                TERM = NO1*NO2*NO3*NO4*NP1*NP2*NP3*NP4*co1*co2*co3*co4
                                                                ! Here we check if all the orbitals i,j,k,l are centered on the
                                                                ! same atom. If so we can use the pre calculated weights and
                                                                ! roots of the Rys polynomals in order to calculate (ij|kl)
                                                                IF ( DOT_PRODUCT(A-B,A-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-B,C-B) .EQ. 0.0d0 .AND.  DOT_PRODUCT(C-D,C-D) .EQ. 0.0d0 ) THEN
                                                                        SAMECITE = .TRUE.
                                                                ELSE
                                                                        SAMECITE = .FALSE.
                                                                ENDIF
                                                                
                                                                !-----------------------------------------------
                                                                ! Here we calculate the two electron integrals
                                                                !-----------------------------------------------
                                                                IF ( KL .NE. IJ ) THEN
                                                                    Intsv(G) = Intsv(G) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)
                                                                ENDIF
                                                                !-------------------------------------------------------------
                                                                ! Here we calculate the gradient of the two electron integrals
                                                                !-------------------------------------------------------------

                                                                IF ( CFORCE ) THEN
                                                                        CALL gradprimeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW,gradient)
                                                                        gradIntsv(BAS%PSI(I)%ATYPE,:,G) = gradIntsv(BAS%PSI(I)%ATYPE,:,G) + TERM*gradient(1,:) 
                                                                        gradIntsv(BAS%PSI(J)%ATYPE,:,G) = gradIntsv(BAS%PSI(J)%ATYPE,:,G) + TERM*gradient(2,:)
                                                                        gradIntsv(BAS%PSI(K)%ATYPE,:,G) = gradIntsv(BAS%PSI(K)%ATYPE,:,G) + TERM*gradient(3,:)
                                                                        gradIntsv(BAS%PSI(L)%ATYPE,:,G) = gradIntsv(BAS%PSI(L)%ATYPE,:,G) + TERM*gradient(4,:)
                                                                ENDIF
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDIF
       ENDDO
       ENDIF
       CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
              
       
ENDDO
END SUBROUTINE eeints
