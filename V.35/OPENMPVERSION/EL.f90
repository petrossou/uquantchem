FUNCTION EL(N,N3,NATOMS,Cup,Cdown,BAS,ATOMS,b,c,r)
        ! This function calculates the Local energy  of the trial function used in the
        ! DQMC importance sampling. 
        ! N = Total number of electrons, Cup = Expansion coefficients to 
        ! the Hartree fock spin-up orbitals, Cdown = Expansion coefficients 
        ! to Hartree fock spin-down orbitals.
        ! BAS = The basis used to solve the HF-equations
        ! b,c = Jastrow coefficients
        ! r = array containing the coordinates of all the electrons, i.e a Nx3
        ! See Eqn.(7) p.6 in my DQMC notes.
        ! NATOMS = Number of atoms in molecule, ATOMS = Atoms (Position and Z-number)
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: EL
        INTEGER  :: N,N3,NATOMS
        TYPE(BASIS) :: BAS
        TYPE(ATOM) :: ATOMS(NATOMS)
        DOUBLE PRECISION :: Cup(BAS%NBAS,BAS%NBAS),Cdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: b,c,r(N3)
        DOUBLE PRECISION, EXTERNAL :: laplacetrialfnk, trialfnk
        INTEGER :: I,J,K
        DOUBLE PRECISION :: rij, V(3), U(3), W(3),POTE
        
        POTE = 0.0d0

        ! Calculation of the electron-electron contribution 
        ! to the local energy and the electron-nucleon contribution
        
        DO I=1,N
             V(1) = r(3*(I-1) + 1 )
             V(2) = r(3*(I-1) + 2 )
             V(3) = r(3*(I-1) + 3 )
             DO K=1,NATOMS
                        W = V - ATOMS(K)%R
                        POTE = POTE - (1.0d0*ATOMS(K)%Z)/( sqrt(DOT_PRODUCT(W,W)) )
             ENDDO
             DO J =I+1,N
                        U(1) = r(3*(J-1) + 1 )
                        U(2) = r(3*(J-1) + 2 )
                        U(3) = r(3*(J-1) + 3 )
                        rij = sqrt(DOT_PRODUCT(V-U,V-U))
                        POTE = POTE + 1.0d0/rij
             ENDDO
      ENDDO
                        
      ! Adding the local Kinetic energy

      EL = POTE  - 0.5d0*laplacetrialfnk(N,N3,Cup,Cdown,BAS,b,c,r)/trialfnk(N,N3,Cup,Cdown,BAS,b,c,r)
END FUNCTION EL
