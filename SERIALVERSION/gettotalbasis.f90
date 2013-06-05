SUBROUTINE gettotalbasis(NATOMS,ATOMS,NB,BAS,COUNTER)
      ! If COUNTER=.TRUE. this subroutine calculates the 
      ! size of the total basis set from the local basis sets centered 
      ! at the different atoms
      ! If COUNTER=.FALSE. this subroutine returns the total basis set
      USE datatypemodule
      IMPLICIT NONE
      INTEGER, INTENT(INOUT) :: NB
      TYPE(BASIS), INTENT(INOUT) :: BAS
      INTEGER, INTENT(IN) :: NATOMS
      TYPE(ATOM), INTENT(IN) :: ATOMS(NATOMS)
      LOGICAL, INTENT(IN) :: COUNTER
      INTEGER :: I,J,K,M,N,L(3)
      INTEGER :: PMAP(3,3),DMAP(6,3),FMAP(10,3),GMAP(15,3),HMAP(21,3),IMAP(28,3),JMAP(36,3),KMAP(42,3)
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION :: al
      DOUBLE PRECISION, EXTERNAL :: dfac
      N = 0
      IF (COUNTER) THEN
              NB = 0
              DO I=1,NATOMS
                        DO J=1,ATOMS(I)%NBAS
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 0 ) NB = NB + 1
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 1 ) NB = NB + 3
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 2 ) NB = NB + 6
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 3 ) NB = NB + 10
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 4 ) NB = NB + 15
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 5 ) NB = NB + 21
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 6 ) NB = NB +28
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 7 ) NB = NB +36
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 8 ) NB = NB +42
                        ENDDO
                ENDDO
       ELSE
               !--Permutations of p-orbitals----
                PMAP(1,:) = (/1,0,0/)
                PMAP(2,:) = (/0,1,0/)
                PMAP(3,:) = (/0,0,1/)
                !--Permutations of d-orbitals----
                DMAP(1,:) = (/2,0,0/)
                DMAP(2,:) = (/1,1,0/)
                DMAP(3,:) = (/1,0,1/)
                DMAP(4,:) = (/0,2,0/)
                DMAP(5,:) = (/0,1,1/)
                DMAP(6,:) = (/0,0,2/)
                !--Permutations of f-orbitals----
                FMAP(1,:) = (/3,0,0/)
                FMAP(2,:) = (/0,3,0/)
                FMAP(3,:) = (/0,0,3/)
                FMAP(4,:) = (/2,1,0/)
                FMAP(5,:) = (/1,2,0/)
                FMAP(6,:) = (/2,0,1/)
                FMAP(7,:) = (/1,0,2/)
                FMAP(8,:) = (/0,2,1/)
                FMAP(9,:) = (/0,1,2/)
                FMAP(10,:) = (/1,1,1/)
                !--Permutations of g-orbitals----
                GMAP(1,:) = (/4,0,0/)
                GMAP(2,:) = (/0,4,0/)
                GMAP(3,:) = (/0,0,4/)
                GMAP(4,:) = (/3,1,0/)
                GMAP(5,:) = (/1,3,0/)
                GMAP(6,:) = (/3,0,1/)
                GMAP(7,:) = (/1,0,3/)
                GMAP(8,:) = (/0,3,1/)
                GMAP(9,:) = (/0,1,3/)
                GMAP(10,:) = (/2,2,0/)
                GMAP(11,:) = (/2,0,2/)
                GMAP(12,:) = (/0,2,2/)
                GMAP(13,:) = (/2,1,1/)
                GMAP(14,:) = (/1,2,1/)
                GMAP(15,:) = (/1,1,2/)
                !--Permutations of h-orbitals----
                HMAP(1,:) = (/5,0,0/)
                HMAP(2,:) = (/0,5,0/)
                HMAP(3,:) = (/0,0,5/)
                HMAP(4,:) = (/4,1,0/)
                HMAP(5,:) = (/1,4,0/)
                HMAP(6,:) = (/4,0,1/)
                HMAP(7,:) = (/1,0,4/)
                HMAP(8,:) = (/0,4,1/)
                HMAP(9,:) = (/0,1,4/)
                HMAP(10,:) = (/3,2,0/)
                HMAP(11,:) = (/2,3,0/)
                HMAP(12,:) = (/3,0,2/)
                HMAP(13,:) = (/2,0,3/)
                HMAP(14,:) = (/0,3,2/)
                HMAP(15,:) = (/0,2,3/)
                HMAP(16,:) = (/2,2,1/)
                HMAP(17,:) = (/2,1,2/)
                HMAP(18,:) = (/1,2,2/)
                HMAP(19,:) = (/1,1,3/)
                HMAP(20,:) = (/1,3,1/)
                HMAP(21,:) = (/3,1,1/)
                !--Permutations of i-orbitals----
                IMAP(1,:) = (/6,0,0/)
                IMAP(2,:) = (/0,6,0/)
                IMAP(3,:) = (/0,0,6/)
                IMAP(4,:) = (/5,1,0/)
                IMAP(5,:) = (/1,5,0/)
                IMAP(6,:) = (/5,0,1/)
                IMAP(7,:) = (/1,0,5/)
                IMAP(8,:) = (/0,5,1/)
                IMAP(9,:) = (/0,1,5/)
                IMAP(10,:) = (/4,2,0/)
                IMAP(11,:) = (/2,4,0/)
                IMAP(12,:) = (/4,0,2/)
                IMAP(13,:) = (/2,0,4/)
                IMAP(14,:) = (/0,4,2/)
                IMAP(15,:) = (/0,2,4/)
                IMAP(16,:) = (/3,3,0/)
                IMAP(17,:) = (/3,0,3/)
                IMAP(18,:) = (/0,3,3/)
                IMAP(19,:) = (/4,1,1/)
                IMAP(20,:) = (/1,4,1/)
                IMAP(21,:) = (/1,1,4/)
                IMAP(22,:) = (/3,2,1/)
                IMAP(23,:) = (/2,3,1/)
                IMAP(24,:) = (/3,1,2/)
                IMAP(25,:) = (/2,1,3/)
                IMAP(26,:) = (/1,3,2/)
                IMAP(27,:) = (/1,2,3/)
                IMAP(28,:) = (/2,2,2/)
                !--Permutations of j-orbitals----
                JMAP(1,:) = (/7,0,0/)
                JMAP(2,:) = (/0,7,0/)
                JMAP(3,:) = (/0,0,7/)
                JMAP(4,:) = (/6,1,0/)
                JMAP(5,:) = (/1,6,0/)
                JMAP(6,:) = (/6,0,1/)
                JMAP(7,:) = (/1,0,6/)
                JMAP(8,:) = (/0,6,1/)
                JMAP(9,:) = (/0,1,6/)
                JMAP(10,:) = (/5,2,0/)
                JMAP(11,:) = (/2,5,0/)
                JMAP(12,:) = (/5,0,2/)
                JMAP(13,:) = (/2,0,5/)
                JMAP(14,:) = (/0,5,2/)
                JMAP(15,:) = (/0,2,5/)
                JMAP(16,:) = (/4,3,0/)
                JMAP(17,:) = (/3,4,0/)
                JMAP(18,:) = (/4,0,3/)
                JMAP(19,:) = (/3,0,4/)
                JMAP(20,:) = (/0,4,3/)
                JMAP(21,:) = (/0,3,4/)
                JMAP(22,:) = (/5,1,1/)
                JMAP(23,:) = (/1,5,1/)
                JMAP(24,:) = (/1,1,5/)
                JMAP(25,:) = (/4,2,1/)
                JMAP(26,:) = (/2,4,1/)
                JMAP(27,:) = (/4,1,2/)
                JMAP(28,:) = (/2,1,4/)
                JMAP(29,:) = (/1,4,2/)
                JMAP(30,:) = (/1,2,4/)
                JMAP(31,:) = (/3,2,2/)
                JMAP(32,:) = (/2,3,2/)
                JMAP(33,:) = (/2,2,3/)
                JMAP(34,:) = (/3,3,1/)
                JMAP(35,:) = (/3,1,3/)
                JMAP(36,:) = (/1,3,3/)
                !--Permutations of h-orbitals----
                KMAP(1,:) = (/8,0,0/)
                KMAP(2,:) = (/0,8,0/)
                KMAP(3,:) = (/0,0,8/)
                KMAP(4,:) = (/7,1,0/)
                KMAP(5,:) = (/1,7,0/)
                KMAP(6,:) = (/7,0,1/)
                KMAP(7,:) = (/1,0,7/)
                KMAP(8,:) = (/0,7,1/)
                KMAP(9,:) = (/0,1,7/)
                KMAP(10,:) = (/6,2,0/)
                KMAP(11,:) = (/2,6,0/)
                KMAP(12,:) = (/6,0,2/)
                KMAP(13,:) = (/2,0,6/)
                KMAP(14,:) = (/0,6,2/)
                KMAP(15,:) = (/0,2,6/)
                KMAP(16,:) = (/5,3,0/)
                KMAP(17,:) = (/3,5,0/)
                KMAP(18,:) = (/5,0,3/)
                KMAP(19,:) = (/3,0,5/)
                KMAP(20,:) = (/0,5,3/)
                KMAP(21,:) = (/0,3,5/)
                KMAP(22,:) = (/4,4,0/)
                KMAP(23,:) = (/4,0,4/)
                KMAP(24,:) = (/0,4,4/)
                KMAP(25,:) = (/6,1,1/)
                KMAP(26,:) = (/1,6,1/)
                KMAP(27,:) = (/1,1,6/)
                KMAP(28,:) = (/5,2,1/)
                KMAP(29,:) = (/2,5,1/)
                KMAP(30,:) = (/5,1,2/)
                KMAP(31,:) = (/2,1,5/)
                KMAP(32,:) = (/1,5,2/)
                KMAP(33,:) = (/1,2,5/)
                KMAP(34,:) = (/4,2,2/)
                KMAP(35,:) = (/2,4,2/)
                KMAP(36,:) = (/2,2,4/)
                KMAP(37,:) = (/4,3,1/)
                KMAP(38,:) = (/3,4,1/)
                KMAP(39,:) = (/4,1,3/)
                KMAP(40,:) = (/3,1,4/)
                KMAP(41,:) = (/1,4,3/)
                KMAP(42,:) = (/1,3,4/)
                
                DO I=1,NATOMS
                        DO J=1,ATOMS(I)%NBAS
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 0 ) THEN
                                        N = N+1
                                        BAS%PSI(N) = ATOMS(I)%PSI(J)
                                        BAS%PSI(N)%ATYPE = I
                                        BAS%PSI(N)%R = ATOMS(I)%R
                                        DO M=1,BAS%PSI(N)%NPRIM
                                                al = BAS%PSI(N)%EXPON(M)
                                                L = BAS%PSI(N)%L
                                                BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                        ENDDO
                                        BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 1 ) THEN
                                        DO K=1,3
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = PMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 2 ) THEN
                                        DO K=1,6
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = DMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 3 ) THEN
                                        DO K=1,10
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = FMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 4 ) THEN
                                        DO K=1,15
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = GMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 5 ) THEN
                                        DO K=1,21
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = HMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 6 ) THEN
                                        DO K=1,28
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = IMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 7 ) THEN
                                        DO K=1,36
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = JMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                                IF ( ATOMS(I)%PSI(J)%L(1) .EQ. 8 ) THEN
                                        DO K=1,42
                                                N=N+1
                                                BAS%PSI(N) = ATOMS(I)%PSI(J)
                                                BAS%PSI(N)%ATYPE = I
                                                BAS%PSI(N)%R = ATOMS(I)%R
                                                BAS%PSI(N)%L = KMAP(K,:)
                                                DO M=1,BAS%PSI(N)%NPRIM
                                                        al = BAS%PSI(N)%EXPON(M)
                                                        L = BAS%PSI(N)%L
                                                        BAS%PSI(N)%PRIMNORM(M)=(2.0d0*al/pi)**(0.750d0)*al**(0.50d0*SUM(L))*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                                ENDDO
                                                BAS%PSI(N)%A = BAS%PSI(N)%A*2.0d0**(SUM(L))/sqrt(dfac(2*L(1)-1)*dfac(2*L(2)-1)*dfac(2*L(3)-1))
                                         ENDDO
                                ENDIF
                         ENDDO
               ENDDO
       ENDIF
END SUBROUTINE gettotalbasis
                                




