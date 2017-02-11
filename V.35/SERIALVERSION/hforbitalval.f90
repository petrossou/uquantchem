FUNCTION hforbitalval(BAS,C,N,R)
        ! Calculates the value of the N:th Hartree-Fock orbital PHI_N(R)
        ! at the spatial point R.  C(:,N) = the expansion coefficients of the 
        ! N:th HF-orbital. See QDMC-notes p.9 Eqn.24
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: hforbitalval
        DOUBLE PRECISION :: R(3)
        INTEGER :: N
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: C(BAS%NBAS,BAS%NBAS),RR(3),RB,stoval,P(6)
        INTEGER :: I
        DOUBLE PRECISION, EXTERNAL :: basfunkval
        
        hforbitalval = 0.0d0
        
        DO I=1,BAS%NBAS
                RR = R-BAS%PSI(I)%R
                RB = sqrt(DOT_PRODUCT(RR,RR))
                IF ( RB .GT. BAS%PSI(I)%RC + BAS%PSI(I)%DR .OR. .not. BAS%PSI(I)%HYBRID ) THEN
                        hforbitalval  = hforbitalval + C(I,N)*basfunkval(BAS%PSI(I),R)
                ENDIF
                IF ( RB .LE. BAS%PSI(I)%RC .AND. BAS%PSI(I)%HYBRID ) THEN
                        stoval =  C(I,N)*BAS%PSI(I)%NORM*BAS%PSI(I)%A*( RR(1)**BAS%PSI(I)%L(1) )*( RR(2)**BAS%PSI(I)%L(2) )*( RR(3)**BAS%PSI(I)%L(3) )*exp(-BAS%PSI(I)%STOEXPON*RB)
                        hforbitalval  = hforbitalval + stoval
                        stoval = C(I,N)*BAS%PSI(I)%NORM*BAS%PSI(I)%BETA*RB*( RR(1)**BAS%PSI(I)%L(1) )*( RR(2)**BAS%PSI(I)%L(2) )*( RR(3)**BAS%PSI(I)%L(3) )
                        hforbitalval  = hforbitalval + stoval
                        stoval = C(I,N)*BAS%PSI(I)%NORM*BAS%PSI(I)%GAMA*( RR(1)**BAS%PSI(I)%L(1) )*( RR(2)**BAS%PSI(I)%L(2) )*( RR(3)**BAS%PSI(I)%L(3) )
                        hforbitalval  = hforbitalval + stoval
                ENDIF
                IF ( RB .GT. BAS%PSI(I)%RC .AND. RB .LE. BAS%PSI(I)%RC + BAS%PSI(I)%DR .AND. BAS%PSI(I)%HYBRID ) THEN
                        P = BAS%PSI(I)%P
                        stoval =  C(I,N)*BAS%PSI(I)%NORM*( P(1) + P(2)*RB + (1.0d0/2.0d0)*P(3)*RB**2 + (1.0d0/6.0d0)*P(4)*RB**3 + (1.0d0/12.0d0)*P(5)*RB**4 + (1.0d0/20.0d0)*P(6)*RB**5 )
                        hforbitalval  = hforbitalval + stoval
                ENDIF
        ENDDO
END FUNCTION hforbitalval

