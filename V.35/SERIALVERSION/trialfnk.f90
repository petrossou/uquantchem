FUNCTION trialfnk(N,N3,Cup,Cdown,BAS,b,c,r)
        ! This function calculates the trial function used in the
        ! DQMC importance sampling. 
        ! N = Total number of electrons, Cup = Expansion coefficients to 
        ! the Hartree fock spin-up orbitals, Cdown = Expansion coefficients 
        ! to Hartree fock spin-down orbitals.
        ! BAS = The basis used to solve the HF-equations
        ! b,c = Jastrow coefficients
        ! r = array containing the coordinates of all the electrons, i.e a 3*N-long array
        ! See Eqn.(16) p.7 in my DQMC notes.
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: trialfnk
        INTEGER  :: N,N3
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: b,c,r(N3),VECT(3)
        INTEGER :: NS
        DOUBLE PRECISION :: SU,SD
        DOUBLE PRECISION, EXTERNAL :: slaterdet,jastrow

        NS = ( N - MOD(N,2) )/2
        SU = slaterdet(N,N3,NS,BAS,Cup,r,.TRUE.)
        
        NS = ( N + MOD(N,2) )/2
        SD = slaterdet(N,N3,NS,BAS,Cdown,r,.FALSE.)

        trialfnk = SU*SD*exp(jastrow(N,N3,r,b,c))
END FUNCTION trialfnk
