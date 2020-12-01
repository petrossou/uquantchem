SUBROUTINE doubleoverlap(NATOMS,BAS,Intsv,IND1,IND2,IND3,IND4,Istart,Iend,numprocessors,id,NONZERO)
      ! This subroutine calculates the 
      ! double overlap tesor to be used to calculate the relativistic Dirac term/correction emanating from the Hartree electron
      ! repulsion potential.
      USE datatypemodule
      USE OMP_LIB
      IMPLICIT NONE
      INCLUDE "mpif.h"
      DOUBLE PRECISION, EXTERNAL :: doubleprimoverlap
      TYPE(BASIS), INTENT(IN) :: BAS
      INTEGER, INTENT(IN) :: id,numprocessors,NATOMS
      INTEGER*8, INTENT(IN) :: Istart,Iend,NONZERO
      INTEGER, INTENT(INOUT) :: IND1(Istart:Iend),IND2(Istart:Iend),IND3(Istart:Iend),IND4(Istart:Iend)
      DOUBLE PRECISION, INTENT(INOUT) :: Intsv(Istart:Iend)
      INTEGER*8 :: IMAP(Istart:Iend,6),N0p(numprocessors)
      INTEGER*8 :: I,J,K,L,M,N,NN,MM,P,Q,G,IJ,KL,L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4,NSTART,GG,STRL,KK,PPP,II
      INTEGER*8 :: rcounts(numprocessors),displs(numprocessors),ierr,scount,SI,SIC(numprocessors)
      DOUBLE PRECISION :: NO1,NO2,NO3,NO4,NP1,NP2,NP3,NP4,TERM,N0,INTEGRAAL
      DOUBLE PRECISION :: A(3),B(3),C(3),D(3),al1,al2,al3,al4,co1,co2,co3,co4
      INTEGER, EXTERNAL :: ijkl
      LOGICAL :: DOSUMOVERPRIMS
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION, PARAMETER :: lightspeed = 137.035999084
      
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
                                                                
                                                                !-----------------------------------------------------------
                                                                ! Here we calculate the double overlap integrals < a b c d >
                                                                !------------------------------------------------------------
                                                                IF ( KL .NE. IJ ) THEN
                                                                    Intsv(G) = Intsv(G) + TERM*doubleprimoverlap(L1,M1,N1,A,al1,L2,M2,N2,B,al2, &
                                                                               & L3,M3,N3,C,al3,L4,M4,N4,D,al4)
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
END SUBROUTINE doubleoverlap
