FUNCTION jastrow(N,N3,r,b,c)
        ! This function calculates the (Total) Jastrow factor.
        ! N = total number of electrons
        ! r = coordinates of all electrons, r is a 3*N-long array
        ! b,c = jastrow parameters. See Eqn.(15), p.7 in my QDMC notes
        IMPLICIT NONE
        DOUBLE PRECISION :: jastrow
        INTEGER :: N,N3
        DOUBLE PRECISION :: r(N3),b,c
        DOUBLE PRECISION, EXTERNAL :: jastrowup,jastrowdown,jastrowud

        jastrow = jastrowup(N,r,b,c) + jastrowdown(N,r,b,c) + jastrowud(N,r,b,c)

END FUNCTION jastrow
