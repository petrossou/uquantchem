FUNCTION doubleprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,L3,M3,N3,C,alpha3,L4,M4,N4,D,alpha4)
        ! This function calculates the overlap integral
        ! S = S_x*S_y*S_z as given by the expressions 
        ! on the pages 234-238 in the Cook book
        IMPLICIT NONE
        INTEGER :: L1,M1,N1,L2,M2,N2,L3,M3,N3,L4,M4,N4
        DOUBLE PRECISION :: doubleprimoverlap
        DOUBLE PRECISION :: alpha1,alpha2,alpha3,alpha4,A(3),B(3),C(3),D(3),PA(3),PB(3),P(3),QC(3),QD(3),Q(3),RP(3),RQ(3),CC
        DOUBLE PRECISION :: Sx,Sy,Sz,gama,delta,xi
        INTEGER :: I,J,K
        DOUBLE PRECISION, EXTERNAL :: fj,dfac
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0

        Sx = 0.0d0
        Sy = 0.0d0
        Sz = 0.0d0
        
        gama = alpha1 + alpha2
        delta = alpha3 + alpha4
        xi = gama + delta

        P = (A*alpha1 + B*alpha2)/gama
        Q = (C*alpha3 + D*alpha4)/delta
        R = (P*gama + Q*delta)/xi

        CC = (alpha1*alpha2/gama)*DOT_PRODUCT(A-B,A-B) + (alpha3*alpha4/delta)*DOT_PRODUCT(C-D,C-D) + (gama*delta/xi)*DOT_PRODUCT(P-Q,P-Q)

        PA = P-A
        PB = P-B
        QC = Q-C
        QD = Q-D
        RP = R-P
        RQ = R-Q

        DO I=0,L1+L2
                DO J=0,L3+L4
                        DO K=0,(I+J-MOD(I+J,2))/2
                                Sx = Sx+fj(I,L1,L2,PA(1),PB(1))*fj(J,L3,L4,QC(1),QD(1))*fj(2*K,I,J,RP(1),RQ(1))*dfac(2*K-1)/((2.0d0*xi)**K)
                        ENDDO
                ENDDO
        ENDDO
        Sx = sqrt(pi/xi)*Sx
                        
        DO I=0,M1+M2
                DO J=0,M3+M4
                        DO K=0,(I+J-MOD(I+J,2))/2
                                Sy = Sy+fj(I,M1,M2,PA(2),PB(2))*fj(J,M3,M4,QC(2),QD(2))*fj(2*K,I,J,RP(2),RQ(2))*dfac(2*K-1)/((2.0d0*xi)**K)
                        ENDDO
                ENDDO
        ENDDO
        Sy = sqrt(pi/xi)*Sy
        
        DO I=0,N1+N2
                DO J=0,N3+N4
                        DO K=0,(I+J-MOD(I+J,2))/2
                                Sz = Sz+fj(I,N1,N2,PA(3),PB(3))*fj(J,N3,N4,QC(3),QD(3))*fj(2*K,I,J,RP(3),RQ(3))*dfac(2*K-1)/((2.0d0*xi)**K)
                        ENDDO
                ENDDO
        ENDDO
        Sz = sqrt(pi/xi)*Sz
       
        doubleprimoverlap = exp(-CC)*Sx*Sy*Sz

 END FUNCTION doubleprimoverlap
