SUBROUTINE getvxc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,Vxc)
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
        INTEGER, INTENT(IN) :: NATOMS,LORDER,CGORDER,NTOTALQUAD,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION, INTENT(OUT) :: Vxc(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: rm,rad,p,ptot,x,faktor,rcutt,Vxctemp(2,BAS%NBAS,BAS%NBAS,NTOTALQUAD)
        DOUBLE PRECISION, ALLOCATABLE :: r(:),Vxcr(:,:,:)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi
        INTEGER :: n,m,I,J,id,nthreads
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        INTEGER :: omp_get_thread_num,omp_get_num_threads


        rcutt = 100.0d0 ! cuttoff radius in au for numerical quadrature
        
        
        Vxc = 0.0d0
        !$OMP PARALLEL &
        !$OMP PRIVATE(J,n,m,I,x,rm,rad,r,faktor,Vxcr,p) 
        ALLOCATE(r(3),Vxcr(2,BAS%NBAS,BAS%NBAS))
        !$OMP DO 
        DO J=1,NTOTALQUAD
                n = Q1(J)
                m = Q2(J)
                I = Q3(J)
                x = CGQ(n,1)
                rm = atomicradii(ATOMS(I)%Z)
                IF ( ATOMS(I)%Z .NE. 1 ) rm = 0.50d0*rm
                rad = -rm*LOG(1.0d0 - ((x+1.0d0)/2.0d0)**4)
                IF ( rad .LE. rcutt ) THEN  ! Cut off radius of 10*rm  is employed
                        r(1) =  rad*sin(LQ(m,2))*cos(LQ(m,1))
                        r(2) =  rad*sin(LQ(m,2))*sin(LQ(m,1))
                        r(3) =  rad*cos(LQ(m,2))
                        r = r + ATOMS(I)%R
                        p = pvoronoi(I,NATOMS,ATOMS,r)
                        ! This factor comes from the variable change r ! --> x, 
                        ! and from using chebyshev-gauss radial quadrature of second order
                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                        CALL getvxcr(CORRLEVEL,NATOMS,BAS,Pup,Pdown,gradS,r,Vxcr)
                        !Vxc = Vxc +  4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*Vxcr*rad**2
                        Vxctemp(:,:,:,J) =  4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*Vxcr*rad**2
                ENDIF
        ENDDO
        !$OMP END DO
        DEALLOCATE(r,Vxcr)
        !$OMP END PARALLEL
        
        !========================================================================
        ! Since the results depend on the order of summation we have to do this
        ! to obtain consistent results
        !========================================================================
        DO J=1,NTOTALQUAD
                Vxc = Vxc + Vxctemp(:,:,:,J)
        ENDDO
END SUBROUTINE getvxc

