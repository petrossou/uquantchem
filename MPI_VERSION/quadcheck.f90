FUNCTION quadcheck(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id)
        ! This function calculates the exchange correlation energy
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        DOUBLE PRECISION :: quadcheck
        CHARACTER(LEN=20) :: CORRLEVEL
        INTEGER :: NTOTALQUAD,Qstart,Qend
        INTEGER :: NATOMS,LORDER,CGORDER,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend)
        INTEGER :: numprocessors,id
        TYPE(ATOM) :: ATOMS(NATOMS)
        TYPE(BASIS) :: BAS
        DOUBLE PRECISION :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION :: r(3),rm,rad,p,ptot,x,faktor,norm,rcutt
        DOUBLE PRECISION :: quadcheckthread
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi,excdr,rho
        INTEGER :: n,m,I,J,ierr
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
       
        rcutt = 20.0d0 ! cuttoff radius in au for numerical quadrature
        
        quadcheck = 0.0d0
        quadcheckthread = 0.0d0
        DO J=Qstart,Qend
                                n = Q1(J)
                                m = Q2(J)
                                I = Q3(J)
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
                                        quadcheckthread = quadcheckthread + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*rho(BAS,Pup+Pdown,r)*rad**2
                                ENDIF
        ENDDO
        CALL MPI_REDUCE(quadcheckthread,quadcheck,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(quadcheck,1,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END FUNCTION quadcheck

