SUBROUTINE getvxc(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,Vxc)
        ! This subroutine calculates the exchange correlation potential 
        ! matrix elements, <i|Vxc|j>, and the nuclear gradients of these matrix
        ! elements, g<i|Vxc|j>, using the resulotion of the identity into fuzzy
        ! polyheadra prescribed by A. D. Becke, J. Chem. Phys. 88, 2547 (1988),
        ! together with Gauss-Chebysshev and Lebedev quadrature.
        ! All numbers in parenthsesis appearing as comments are refering to
        ! equations in  J. Chem. Phys. 88, 2547 (1988).
        USE datatypemodule
        IMPLICIT NONE
        INCLUDE "mpif.h"
        CHARACTER(LEN=20), INTENT(IN) :: CORRLEVEL
        INTEGER, INTENT(IN) :: NTOTALQUAD,Qstart,Qend
        INTEGER, INTENT(IN) :: NATOMS,LORDER,CGORDER,Q1(Qstart:Qend),Q2(Qstart:Qend),Q3(Qstart:Qend),numprocessors,id
        TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
        TYPE(BASIS), INTENT(IN)  :: BAS
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION, INTENT(OUT) :: Vxc(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: r(3),rm,rad,p,ptot,x,faktor,rcutt
        DOUBLE PRECISION :: Vxcr(2,BAS%NBAS,BAS%NBAS),Vxcsend(2,BAS%NBAS,BAS%NBAS),Vxcrecieve(2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: Vxccollect(numprocessors,2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi
        INTEGER :: n,m,I,J,ierr,status(MPI_STATUS_SIZE),tag,NN,NB
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0

        NB = BAS%NBAS
        rcutt = 20.0d0 ! cuttoff radius in au for numerical quadrature
        tag = 1
 
        Vxc = 0.0d0
        Vxcsend = 0.0d0
        DO J=Qstart,Qend
        !DO J=1,NTOTALQUAD
               n = Q1(J)
               m = Q2(J)
               I = Q3(J)
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
                                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                                        CALL getvxcr(CORRLEVEL,NATOMS,BAS,Pup,Pdown,r,Vxcr)
                                        ! This factor comes from the variable change r ! --> x, 
                                        ! and from using chebyshev-gauss radial quadrature of second order
                                        !faktor = 2.0d0*rm/(sqrt(1.0d0 - x**2)*(1.0d0 - x )**2 )
                                        Vxcsend = Vxcsend  + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*Vxcr*rad**2
                                        !Vxc = Vxc  + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*Vxcr*rad**2
               ENDIF
        ENDDO

        CALL MPI_REDUCE(Vxcsend,Vxc,2*BAS%NBAS*BAS%NBAS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(Vxc,2*BAS%NBAS*BAS%NBAS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
END SUBROUTINE getvxc

