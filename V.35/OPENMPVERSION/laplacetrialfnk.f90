FUNCTION laplacetrialfnk(N,N3,Cup,Cdown,BAS,b,c,r)
        ! This function calculates the laplacian of the trial function used in the
        ! DQMC importance sampling. 
        ! N = Total number of electrons, Cup = Expansion coefficients to 
        ! the Hartree fock spin-up orbitals, Cdown = Expansion coefficients 
        ! to Hartree fock spin-down orbitals.
        ! BAS = The basis used to solve the HF-equations
        ! b,c = Jastrow coefficients
        ! r = array containing the coordinates of all the electrons, i.e a Nx3
        ! See Eqn.(29) p.10 in my DQMC notes.
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: laplacetrialfnk
        INTEGER  :: N,N3
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: b,c,r(N3)
        DOUBLE PRECISION :: LSD,LJS,TEST
        INTEGER :: I,M,NS
        DOUBLE PRECISION :: gradD(3),gradJ(3),SU,SD,VECT1(3),VECT2(3)
        DOUBLE PRECISION, EXTERNAL :: slaterdet,laplacejastrow,laplaceslaterdet,jastrow
        DOUBLE PRECISION, EXTERNAL :: laplacehforbitalval,hforbitalval
        EXTERNAL :: gradslaterdet, gradjastrow
        
        LJS = 0.0d0
        LSD = 0.0d0

        M = ( N - MOD(N,2) )/2 

        NS = ( N - MOD(N,2) )/2
        SU = slaterdet(N,N3,NS,BAS,Cup,r,.TRUE.)
        
        NS = ( N + MOD(N,2) )/2
        SD = slaterdet(N,N3,NS,BAS,Cdown,r,.FALSE.)
        
        laplacetrialfnk = 0.0d0

        DO I=1,N
                call gradjastrow(I,N,N3,r,b,c,gradJ)
                IF ( I .LE. M ) THEN
                        NS = ( N - MOD(N,2) )/2
                        call gradslaterdet(I,N,N3,NS,BAS,Cup,r,.TRUE.,gradD)
                        LSD = laplaceslaterdet(I,N,N3,NS,BAS,Cup,r,.TRUE.)
                        LJS = laplacejastrow(I,N,N3,r,b,c)
                        !laplacetrialfnk = laplacetrialfnk + LSD*SD
                        laplacetrialfnk = laplacetrialfnk + (LSD + 2*DOT_PRODUCT(gradD,gradJ)+ SU*LJS + SU*DOT_PRODUCT(gradJ,gradJ) )*SD
                ELSE
                        NS = ( N + MOD(N,2) )/2
                        call gradslaterdet(I,N,N3,NS,BAS,Cdown,r,.FALSE.,gradD)
                        LSD = laplaceslaterdet(I,N,N3,NS,BAS,Cdown,r,.FALSE.)
                        LJS = laplacejastrow(I,N,N3,r,b,c)
                        !laplacetrialfnk = laplacetrialfnk + LSD*SU
                        laplacetrialfnk = laplacetrialfnk + (LSD + 2*DOT_PRODUCT(gradD,gradJ)+ SD*LJS + SD*DOT_PRODUCT(gradJ,gradJ) )*SU
                ENDIF
        ENDDO
        laplacetrialfnk = laplacetrialfnk*exp(jastrow(N,N3,r,b,c))
END FUNCTION laplacetrialfnk
