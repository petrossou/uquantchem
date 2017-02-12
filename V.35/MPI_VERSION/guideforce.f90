SUBROUTINE guideforce(N,N3,Cup,Cdown,BAS,b,c,r,force)
        ! This subroutine calculates the guiding force used in the
        ! DQMC importance sampling. The force is a 3N dimensional array
        ! N = Total number of electrons, Cup = Expansion coefficients to 
        ! the Hartree fock spin-up orbitals, Cdown = Expansion coefficients 
        ! to Hartree fock spin-down orbitals.
        ! BAS = The basis used to solve the HF-equations
        ! b,c = Jastrow coefficients
        ! r = array containing the coordinates of all the electrons, i.e a
        ! 3*N-long array
        ! array. force = the 3N diminsional force array.
        ! See Eqn.(18) p.7 in my DQMC notes.
        USE datatypemodule
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N,N3
        TYPE(BASIS), INTENT(IN) :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: b,c,r(N3)
        DOUBLE PRECISION, INTENT(OUT) :: force(N3)
        INTEGER :: I,M,NS
        DOUBLE PRECISION :: gradD(3),gradJ(3),SDET
        DOUBLE PRECISION, EXTERNAL :: slaterdet
        EXTERNAL :: gradslaterdet, gradjastrow
        
        M = ( N - MOD(N,2) )/2 

        force = 0.0d0

        DO I=1,N
                IF ( I .LE. M ) THEN
                        NS = ( N - MOD(N,2) )/2
                        call gradslaterdet(I,N,N3,NS,BAS,Cup,r,.TRUE.,gradD)
                        IF ( I .EQ. 1 ) SDET = slaterdet(N,N3,NS,BAS,Cup,r,.TRUE.)
                ELSE
                        NS = ( N + MOD(N,2) )/2
                        call gradslaterdet(I,N,N3,NS,BAS,Cdown,r,.FALSE.,gradD)
                        IF ( I .EQ. M+1 ) SDET = slaterdet(N,N3,NS,BAS,Cdown,r,.FALSE.)
                ENDIF
                call gradjastrow(I,N,N3,r,b,c,gradJ)
                IF ( SDET .NE. 0.0d0 ) THEN
                     force(1+3*(I-1)) = ( gradD(1)/SDET  + gradJ(1) )
                     force(2+3*(I-1)) = ( gradD(2)/SDET  + gradJ(2) )
                     force(3+3*(I-1)) = ( gradD(3)/SDET  + gradJ(3) )
                ELSE
                     force(1+3*(I-1)) = 0.0d0
                     force(2+3*(I-1)) = 0.0d0
                     force(3+3*(I-1)) = 0.0d0
                ENDIF
        ENDDO
END SUBROUTINE guideforce
