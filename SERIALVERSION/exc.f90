FUNCTION exc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ)
        ! This function calculates the exchange correlation energy
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: exc
        CHARACTER(LEN=20) :: CORRLEVEL
        INTEGER :: NATOMS,LORDER,CGORDER
        TYPE(ATOM) :: ATOMS(NATOMS)
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION :: r(3),rm,rad,p,ptot,x,faktor,norm,rcutt
        DOUBLE PRECISION :: Vxcr(2,BAS%NBAS,BAS%NBAS),gVxcr(NATOMS,2,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi,excdr,rho
        INTEGER :: n,m,I,J
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
       
        rcutt = 40.0d0 ! cuttoff radius in au for numerical quadrature
        
        exc = 0.0d0
        DO n=1,CGORDER                  ! Loop over chebyshev-gauss radial quadrature of second order.
                DO m=1,LORDER           ! Loop over Lebedev quadrature over the surface of the sphere.
                        DO I=1,NATOMS   ! Loop over fuzzy voronoi polyheadra centered on each of the atoms.
                                x = CGQ(n,1)
                                rm = atomicradii(ATOMS(I)%Z)
                                IF ( ATOMS(I)%Z .NE. 1 ) rm = 0.50d0*rm
                                !rad = rm*( 1.0d0 + x )/( 1.0d0 - x )  ! (25)
                                rad = -rm*LOG(1.0d0 - ((x+1.0d0)/2.0d0)**4)
                                IF ( rad .LE. rcutt ) THEN 
                                        r(1) =  rad*sin(LQ(m,2))*cos(LQ(m,1))
                                        r(2) =  rad*sin(LQ(m,2))*sin(LQ(m,1))
                                        r(3) =  rad*cos(LQ(m,2))
                                        r = r + ATOMS(I)%R
                                        p = pvoronoi(I,NATOMS,ATOMS,r)
                                        ! This factor comes from the variable change r ! --> x, 
                                        ! and from using chebyshev-gauss radial quadrature of second order
                                        !faktor = 2.0d0*rm/(sqrt(1.0d0 - x**2)*(1.0d0 - x )**2 )
                                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                                        exc = exc  + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*excdr(CORRLEVEL,BAS,Pup,Pdown,r)*rad**2
                                ENDIF
                         ENDDO
                ENDDO
        ENDDO
                        
END FUNCTION exc

