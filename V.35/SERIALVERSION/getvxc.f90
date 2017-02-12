SUBROUTINE getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,Vxc)
        ! This subroutine calculates the exchange correlation potential 
        ! matrix elements, <i|Vxc|j>, and the nuclear gradients of these matrix
        ! elements, g<i|Vxc|j>, using the resulotion of the identity into fuzzy
        ! polyheadra prescribed by A. D. Becke, J. Chem. Phys. 88, 2547 (1988),
        ! together with Gauss-Chebysshev and Lebedev quadrature.
        ! All numbers in parenthsesis appearing as comments are refering to
        ! equations in  J. Chem. Phys. 88, 2547 (1988).
        USE datatypemodule
        IMPLICIT NONE
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NATOMS,LORDER,CGORDER
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION, INTENT(OUT) :: Vxc(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: r(3),rm,rad,p,ptot,x,faktor,rcutt
        DOUBLE PRECISION :: Vxcr(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi
        INTEGER :: n,m,I,J
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        
        rcutt = 100.0d0 ! cuttoff radius in au for numerical quadrature
        
        Vxc(:,:,:) = 0.0d0
        
        DO n=1,CGORDER                  ! Loop over chebyshev-gauss radial quadrature of second order.
                DO m=1,LORDER           ! Loop over Lebedev quadrature over the surface of the sphere.
                        DO I=1,NATOMS   ! Loop over fuzzy voronoi polyheadra centered on each of the atoms.
                                x = CGQ(n,1)
                                rm = atomicradii(ATOMS(I)%Z)
                                IF ( ATOMS(I)%Z .NE. 1 ) rm = 0.50d0*rm
                                !rad = rm*( 1.0d0 + x )/( 1.0d0 - x )    ! (25)
                                rad = -rm*LOG(1.0d0 - ((x+1.0d0)/2.0d0)**4)
                                IF ( rad .LE. rcutt ) THEN  ! Cut off radius of 10*rm  is employed
                                        r(1) =  rad*sin(LQ(m,2))*cos(LQ(m,1))
                                        r(2) =  rad*sin(LQ(m,2))*sin(LQ(m,1))
                                        r(3) =  rad*cos(LQ(m,2))
                                        r = r + ATOMS(I)%R
                                        p = pvoronoi(I,NATOMS,ATOMS,r)
                                        CALL getvxcr(CORRLEVEL,NATOMS,BAS,Pup,Pdown,gradS,r,Vxcr)
                                        ! This factor comes from the variable change r ! --> x, 
                                        ! and from using chebyshev-gauss radial quadrature of second order
                                        !faktor = 2.0d0*rm/(sqrt(1.0d0 - x**2)*(1.0d0 - x )**2 )
                                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                                        Vxc  = Vxc    + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*Vxcr*rad**2
                                ENDIF
                         ENDDO
                ENDDO
        ENDDO
                        
END SUBROUTINE getvxc

