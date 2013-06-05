SUBROUTINE getxcforce(CORRLEVEL,NATOMS,ATOMS,BAS,Pup,Pdown,gradS,LORDER,CGORDER,LQ,CGQ,NTOTALQUAD,Q1,Q2,Q3,Qstart,Qend,numprocessors,id,xcforce)
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
        DOUBLE PRECISION, INTENT(IN) :: Pup(BAS%NBAS,BAS%NBAS), Pdown(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, INTENT(IN) :: LQ(LORDER,3),CGQ(CGORDER,2)
        DOUBLE PRECISION, INTENT(OUT) :: xcforce(NATOMS,3)
        DOUBLE PRECISION :: r(3),rm,rad,p,ptot,x,faktor,rcutt,sumforce(3),gradSS(NATOMS,3,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION :: excforcedens(NATOMS,3),xcforcesend(NATOMS,3),gradp(NATOMS,3)
        DOUBLE PRECISION :: Vxccollect(numprocessors,2,BAS%NBAS,BAS%NBAS)
        DOUBLE PRECISION, EXTERNAL :: atomicradii,pvoronoi,excdr
        INTEGER :: JJ,n,m,I,J,ierr,status(MPI_STATUS_SIZE),tag,NN,NB
        DOUBLE PRECISION, PARAMETER :: PI = 3.1415926535897932384626433832795028841970d0
        
        NB = BAS%NBAS
        rcutt = 40.0d0 ! cuttoff radius in au for numerical quadrature
        tag = 1
 
        xcforcesend = 0.0d0
        xcforce = 0.0d0
        DO J=Qstart,Qend
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
                                        !CALL gradoverlap(NATOMS,ATOMS,BAS,SHIFT,gradSS)
                                        CALL dftforcedens(CORRLEVEL,NATOMS,BAS,Pup,Pdown,gradS,r,I,excforcedens)
                                        CALL gradpvoronoi(I,NATOMS,ATOMS,r,gradp)
                                        ! This factor comes from the variable change r ! --> x, 
                                        ! and from using chebyshev-gauss radial quadrature of second order
                                        excforcedens(I,:) = 0.0d0
                                        gradp(I,:) = 0.0d0
                                        faktor = 2.0d0*rm*((x+1.0d0)/2.0d0)**3/(sqrt(1.0d0 - x**2)*(1.0d0 - ((x+1.0d0)/2.0d0)**4 ) )
                                        xcforcesend  = xcforcesend + 4*PI*p*CGQ(n,2)*LQ(m,3)*faktor*excforcedens*rad**2
                                        xcforcesend  = xcforcesend - 4*PI*CGQ(n,2)*LQ(m,3)*faktor*excdr(CORRLEVEL,BAS,Pup,Pdown,r)*gradp*rad**2
                                        ! Using translational invariance to
                                        ! calculate the oncite contribution to the force
                                        sumforce = 0.0d0
                                        DO JJ=1,NATOMS
                                                IF ( JJ .NE. I ) sumforce = sumforce - xcforcesend(JJ,:)
                                        ENDDO 
                                        xcforcesend(I,:) = sumforce
               ENDIF
        ENDDO
        CALL MPI_REDUCE(xcforcesend,xcforce,3*NATOMS,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(xcforce,3*NATOMS,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        
END SUBROUTINE getxcforce

