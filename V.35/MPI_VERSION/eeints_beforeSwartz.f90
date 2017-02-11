SUBROUTINE eeints(BAS,Intsv,IND1,IND2,IND3,IND4,NRED,Istart,Iend,PRYSR,PRYSW,numprocessors,APPROXEE,Tol,id)
      ! This subroutine calculates the 
      ! electron-electron repulsion tensor.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: primeeintegral
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: NRED,id,numprocessors,Istart,Iend
      DOUBLE PRECISION, INTENT(IN) :: PRYSR(25,25),PRYSW(25,25),Tol
      LOGICAL, INTENT(IN) :: APPROXEE
      INTEGER, INTENT(INOUT) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(INOUT) :: Intsv(Istart:Iend)
      INTEGER :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART
      INTEGER :: rcounts(numprocessors),displs(numprocessors),ierr,scount
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS,SAMECITE
      ! First we calculate the mapping from one index, g, onto four
      ! ndexes (I,J,K,L), i.e IMAP:g ----->(I,J,K,L) 
      
      P = 0
      DO J=1,BAS%NBAS
        DO I=1,J
            DO L=1,BAS%NBAS
              DO K=1,L
                M = (J*(J-1))/2 + I
                N = (L*(L-1))/2 + K
                P = P + 1
                IF ( P .LE. Iend .AND. P .GE. Istart ) THEN
                             IMAP(P,1) = I
                             IND1(P) = I
                             IMAP(P,2) = J
                             IND2(P) = J
                             IMAP(P,3) = K
                             IND3(P) = K
                             IMAP(P,4) = L
                             IND4(P) = L
                             IMAP(P,5) = M
                             IMAP(P,6) = N
                ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
      
     Intsv(:) = 0.0d0
!==================================================
!==================================================
        DO NN=Istart,Iend
                        G = NN
                        I = IMAP(G,1) 
                        J = IMAP(G,2) 
                        K = IMAP(G,3) 
                        L = IMAP(G,4) 
                        IJ = IMAP(G,5)
                        KL = IMAP(G,6) 
                
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
                        
                        DOSUMOVERPRIMS = .TRUE.
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
                                                                Intsv(G) = Intsv(G) + TERM*primeeintegral(L1,M1,N1,A,al1,L2,M2,N2,B,al2,L3,M3,N3,C,al3,L4,M4,N4,D,al4,SAMECITE,PRYSR,PRYSW)
                                                        ENDDO
                                                ENDDO
                                        ENDDO
                                ENDDO
                        ENDIF
        ENDDO
END SUBROUTINE eeints
