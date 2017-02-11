SUBROUTINE hessianrhohole(NATOMS,BAS,P,PHOLE,gradS,r,BECKECENTER,hr)
      ! This subroutine calculates the hessian matrix 
      ! of the charge density, hr(i,j,k) = d^2rho/dR(k)_i*dr_j
      ! where R(k)_i = i:th coordinate of the nucleus positioned 
      ! at R(k)_i, r_j = the j:th coordinate of the spatial electron 
      ! coordinate r.
      USE datatypemodule
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS,BECKECENTER
      TYPE(BASIS), INTENT(IN) :: BAS
      DOUBLE PRECISION, INTENT(IN) :: P(BAS%NBAS,BAS%NBAS),gradS(NATOMS,3,BAS%NBAS,BAS%NBAS),r(3),PHOLE(BAS%NBAS,BAS%NBAS)
      DOUBLE PRECISION, INTENT(OUT) :: hr(3,3,NATOMS)
      INTEGER :: I,J,K
      DOUBLE PRECISION :: grad(3),dPHImx(BAS%NBAS),dPHImy(BAS%NBAS),dPHImz(BAS%NBAS),hb(3,3)
      DOUBLE PRECISION :: dPHInx(BAS%NBAS),dPHIny(BAS%NBAS),dPHInz(BAS%NBAS)
      DOUBLE PRECISION :: dPHImxx(BAS%NBAS),dPHImyy(BAS%NBAS),dPHImzz(BAS%NBAS)
      DOUBLE PRECISION :: dPHImyx(BAS%NBAS),dPHImzx(BAS%NBAS),dPHImzy(BAS%NBAS)
      DOUBLE PRECISION :: dP(BAS%NBAS,BAS%NBAS,3),PHIn(BAS%NBAS)
      DOUBLE PRECISION, EXTERNAL :: basfunkval

      DO I=1,NATOMS
         DO J=1,BAS%NBAS

                IF ( I .EQ. 1 ) THEN
                        CALL gradbasfunkval(BAS%PSI(J),r,grad)
                        dPHImx(J) = grad(1)
                        dPHImy(J) = grad(2)
                        dPHImz(J) = grad(3)
                        PHIn(J) = basfunkval(BAS%PSI(J),r)
                ENDIF
                
                dPHImxx(J) = 0.0d0
                dPHImyy(J) = 0.0d0
                dPHImzz(J) = 0.0d0

                dPHImyx(J) = 0.0d0
                dPHImzx(J) = 0.0d0
                dPHImzy(J) = 0.0d0

                IF ( BAS%PSI(J)%ATYPE .EQ. BECKECENTER ) THEN
                        dPHImxx(J) = 0.0d0
                        dPHImyy(J) = 0.0d0
                        dPHImzz(J) = 0.0d0

                        dPHImyx(J) = 0.0d0
                        dPHImzx(J) = 0.0d0
                        dPHImzy(J) = 0.0d0
                ELSE
                        IF ( I .EQ. BECKECENTER ) THEN
                                CALL hessianbasfunk(BAS%PSI(J),r,hb)
                                dPHImxx(J) = -hb(1,1)
                                dPHImyy(J) = -hb(2,2)
                                dPHImzz(J) = -hb(3,3)

                                dPHImyx(J) = -hb(2,1)
                                dPHImzx(J) = -hb(3,1)
                                dPHImzy(J) = -hb(3,2)
                        ENDIF
                        IF ( BAS%PSI(J)%ATYPE .EQ. I ) THEN
                                CALL hessianbasfunk(BAS%PSI(J),r,hb)
                                dPHImxx(J) = hb(1,1)
                                dPHImyy(J) = hb(2,2)
                                dPHImzz(J) = hb(3,3)

                                dPHImyx(J) = hb(2,1)
                                dPHImzx(J) = hb(3,1)
                                dPHImzy(J) = hb(3,2)
                        ENDIF
                ENDIF
                
                dPHInx(J) = 0.0d0
                dPHIny(J) = 0.0d0
                dPHInz(J) = 0.0d0
                
                IF ( BAS%PSI(J)%ATYPE .EQ. BECKECENTER ) THEN
                        dPHInx(J) = 0.0d0
                        dPHIny(J) = 0.0d0
                        dPHInz(J) = 0.0d0
                ELSE
                        IF ( I .EQ. BECKECENTER ) THEN
                                dPHInx(J) = -dPHImx(J)
                                dPHIny(J) = -dPHImy(J)
                                dPHInz(J) = -dPHImz(J)
                        ENDIF
                        IF ( BAS%PSI(J)%ATYPE .EQ. I ) THEN
                                dPHInx(J) = dPHImx(J)
                                dPHIny(J) = dPHImy(J)
                                dPHInz(J) = dPHImz(J)
                        ENDIF
                ENDIF

         ENDDO

         dP(:,:,1) = -MATMUL(P,MATMUL(gradS(I,1,:,:),P))+MATMUL(PHOLE,MATMUL(gradS(I,1,:,:),PHOLE))
         dP(:,:,2) = -MATMUL(P,MATMUL(gradS(I,2,:,:),P))+MATMUL(PHOLE,MATMUL(gradS(I,2,:,:),PHOLE))
         dP(:,:,3) = -MATMUL(P,MATMUL(gradS(I,3,:,:),P))+MATMUL(PHOLE,MATMUL(gradS(I,3,:,:),PHOLE))

         hr(1,1,I) = 2*( DOT_PRODUCT(dPHImx,MATMUL(dP(:,:,1),PHIn)) - DOT_PRODUCT(dPHImxx,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImx,MATMUL(P-PHOLE,dPHInx)) )
         hr(2,1,I) = 2*( DOT_PRODUCT(dPHImx,MATMUL(dP(:,:,2),PHIn)) - DOT_PRODUCT(dPHImyx,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImx,MATMUL(P-PHOLE,dPHIny)) )
         hr(3,1,I) = 2*( DOT_PRODUCT(dPHImx,MATMUL(dP(:,:,3),PHIn)) - DOT_PRODUCT(dPHImzx,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImx,MATMUL(P-PHOLE,dPHInz)) )

         hr(1,2,I) = 2*( DOT_PRODUCT(dPHImy,MATMUL(dP(:,:,1),PHIn)) - DOT_PRODUCT(dPHImyx,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImy,MATMUL(P-PHOLE,dPHInx)) )
         hr(2,2,I) = 2*( DOT_PRODUCT(dPHImy,MATMUL(dP(:,:,2),PHIn)) - DOT_PRODUCT(dPHImyy,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImy,MATMUL(P-PHOLE,dPHIny)) )
         hr(3,2,I) = 2*( DOT_PRODUCT(dPHImy,MATMUL(dP(:,:,3),PHIn)) - DOT_PRODUCT(dPHImzy,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImy,MATMUL(P-PHOLE,dPHInz)) )

         hr(1,3,I) = 2*( DOT_PRODUCT(dPHImz,MATMUL(dP(:,:,1),PHIn)) - DOT_PRODUCT(dPHImzx,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImz,MATMUL(P-PHOLE,dPHInx)) )
         hr(2,3,I) = 2*( DOT_PRODUCT(dPHImz,MATMUL(dP(:,:,2),PHIn)) - DOT_PRODUCT(dPHImzy,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImz,MATMUL(P-PHOLE,dPHIny)) )
         hr(3,3,I) = 2*( DOT_PRODUCT(dPHImz,MATMUL(dP(:,:,3),PHIn)) - DOT_PRODUCT(dPHImzz,MATMUL(P-PHOLE,PHIn)) - DOT_PRODUCT(dPHImz,MATMUL(P-PHOLE,dPHInz)) )

       ENDDO
       !WRITE(*,'(3(F15.10,2X))'),hr(1,:,2)
       !WRITE(*,'(3(F15.10,2X))'),hr(2,:,2)
       !WRITE(*,'(3(F15.10,2X))'),hr(3,:,2)
       !print*,'======================'
       !WRITE(*,'(3(F15.10,2X))'),hb(1,:)
       !WRITE(*,'(3(F15.10,2X))'),hb(2,:)
       !WRITE(*,'(3(F15.10,2X))'),hb(3,:)
END SUBROUTINE hessianrhohole

