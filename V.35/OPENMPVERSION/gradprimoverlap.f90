SUBROUTINE gradprimoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,gradient)
        ! This function calculates the gradien overlap integral
        ! grad(S), with respect to the nuclear coordinate of 
        ! orbital 1 and 2. The output 
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: L1,M1,N1,L2,M2,N2
        DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,A(3),B(3)
        DOUBLE PRECISION, EXTERNAL :: primoverlap
        DOUBLE PRECISION, INTENT(OUT) :: gradient(2,3)
        
        ! x-coordinat with respect to atom 1:
        gradient(1,1)  = 0.0d0
        IF ( L1 .GT. 0 ) gradient(1,1) = -L1*primoverlap(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)
        gradient(1,1)  = gradient(1,1) + 2*alpha1*primoverlap(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2)

        ! x-coordinat with respect to atom 2:
        gradient(2,1)  = 0.0d0
        IF ( L2 .GT. 0 ) gradient(2,1) = -L2*primoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2)
        gradient(2,1)  = gradient(2,1) + 2*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2)
        
        ! y-coordinat with respect to atom 1:
        gradient(1,2)  = 0.0d0
        IF ( M1 .GT. 0 ) gradient(1,2) = -M1*primoverlap(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2)
        gradient(1,2)  = gradient(1,2) + 2*alpha1*primoverlap(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2)
        
        ! y-coordinat with respect to atom 2:
        gradient(2,2)  = 0.0d0
        IF ( M2 .GT. 0 ) gradient(2,2) = -M2*primoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2)
        gradient(2,2)  = gradient(2,2) + 2*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2)
        
        ! z-coordinat with respect to atom 1:
        gradient(1,3)  = 0.0d0
        IF ( N1 .GT. 0 ) gradient(1,3) = -N1*primoverlap(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2)
        gradient(1,3)  = gradient(1,3) + 2*alpha1*primoverlap(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2)
        
        ! z-coordinat with respect to atom 2:
        gradient(2,3)  = 0.0d0
        IF ( N2 .GT. 0 ) gradient(2,3) = -N2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2)
        gradient(2,3)  = gradient(2,3) + 2*alpha2*primoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2)
        
        END SUBROUTINE gradprimoverlap
