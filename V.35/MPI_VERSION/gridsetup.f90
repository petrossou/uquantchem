     subroutine gridsetup(nproc,nprow,npcol)
!
! This subroutine factorizes the number of processors (nproc)
! into nprow and npcol, that are the sizes of the 2d processors mesh.
!
      implicit none
      integer nproc,nprow,npcol
      integer sqrtnp,i 

      sqrtnp = int( sqrt( dble(nproc) ) + 1 )
      do i=1,sqrtnp
        if(mod(nproc,i).eq.0) nprow = i
      end do
      npcol = nproc/nprow

      return
      end
