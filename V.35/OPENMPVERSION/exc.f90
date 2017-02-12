FUNCTION exc(CORRLEVEL,NATOMS,ATOMS,NTOTALQUAD,Q1,Q2,Q3,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ)
        ! This function calculates the exchange correlation energy
        USE datatypemodule
        IMPLICIT NONE
        DOUBLE PRECISION :: exc
        CHARACTER(LEN=20) :: CORRLEVEL
        INTEGER :: NATOMS,LORDER,CGORDER,NTOTALQUAD,Q1(NTOTALQUAD),Q2(NTOTALQUAD),Q3(NTOTALQUAD)
        TYPE(ATOM) :: ATOMS(NATOMS)
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION :: rm,rad,p,ptot,x,faktor,norm,rcutt
        DOUBLE PRECISION, ALLOCATABLE :: r(:)
        DOUBLE PRECISION :: Vxcr(2,BAS%NBAS,BAS%NBAS),tempexc(NTOTALQUAD)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi,excdr,rho
        INTEGER :: n,m,I,J,id,nthreads
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        INTEGER :: omp_get_thread_num,omp_get_num_threads


        rcutt = 100.0d0 
        
        exc = 0.0d0
        !$OMP PARALLEL &
        !$OMP PRIVATE(J,n,m,I,x,rm,rad,r,faktor,p)
        ALLOCATE(r(3))
        !$OMP DO 
        DO J=1,NTOTALQUAD
                n = Q1(J)
                m = Q2(J)
                I = Q3(J)
                x = CGQ(n,1)
                rm = atomicradii(ATOMS(I)%Z)
                IF ( ATOMS(I)%Z .NE. 1 ) rm = 0.50d0*rm
                rad = -rm*LOG(1.0d0 - ((x+1.0d0)/2.0d0)**4)
                IF ( rad .LE. rcutt ) THEN 
                        r(1) =  rad*sin(LQ(m,2))*cos(LQ(m,1))
                        r(2) =  rad*sin(LQ(m,2))*sin(LQ(m,1))
                        r(3) =  rad*cos(LQ(m,2))
                        r = r + ATOMS(I)%R
                        p = pvoronoi(I,NATOMS,ATOMS,r)
                        ! This factor comes from the variable change r ! --> x, 
                        ! and from using chebyshev-gauss radial quadrature of second order
                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                        !exc = exc + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*excdr(CORRLEVEL,BAS,Pup,Pdown,r)*rad**2
                        tempexc(J) = 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*excdr(CORRLEVEL,BAS,Pup,Pdown,r)*rad**2
                ENDIF
        ENDDO
        !$OMP END DO
        DEALLOCATE(r)
        !$OMP END PARALLEL

        !========================================================================
        ! Since the results depend on the order of summation we have to do this
        ! to obtain consistent results
        !========================================================================
        DO J=1,NTOTALQUAD
                exc = exc + tempexc(J)
        ENDDO

        END FUNCTION exc
