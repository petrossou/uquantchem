SUBROUTINE gradpriminhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,d,omega,t,gDpole)
        ! This subroutine calculates the nuclear gradients of the three dipole tensor 
        ! elements between two primitive gaussian basis-functions 
        ! in the case of an inhomogeneous electric-field, i.e 
        ! E(r,t) = E0*sin( omega*t - (omega/c)*r ), d =x or y or z
        ! the out-put is given by ( if POVERLAP = .TRUE. ):
        ! Dpole = <L1,M1,N1,A,alpha1|E(x,t)*y|L2,M2,N2,B,alpha2>, r = x,  if d = 1
        ! Dpole = <L1,M1,N1,A,alpha1|E(y,t)*z|L2,M2,N2,B,alpha2>, r = y,  if d = 2
        ! Dpole = <L1,M1,N1,A,alpha1|E(z,t)*x|L2,M2,N2,B,alpha2>, r = z,  if d = 3
        ! the out-put is given by ( if POVERLAP = .FALSE. ):
        ! Dpole = <L1,M1,N1,A,alpha1|E(x,t)*(y-B(2))|L2,M2,N2,B,alpha2>, r = x,  if d = 1
        ! Dpole = <L1,M1,N1,A,alpha1|E(y,t)*(z-B(3))|L2,M2,N2,B,alpha2>, r = y,  if d = 2
        ! Dpole = <L1,M1,N1,A,alpha1|E(z,t)*(x-B(1))|L2,M2,N2,B,alpha2>, r = z,  if d = 3
        ! See my notes in the blue note-book pages 13-22. Especially 
        ! Equation (80) on page 22.
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: L1,M1,N1,L2,M2,N2,d
        DOUBLE PRECISION, INTENT(OUT) :: gDpole(2,3)
        DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,A(3),B(3),omega,t
        LOGICAL ::  POVERLAP
        DOUBLE PRECISION :: Sx,Sy,Sz,Sc,gama,phi,kv,cl,PA(3),PB(3),P(3),C,Dpole
        INTEGER :: I,J,K,n,m
        DOUBLE PRECISION, EXTERNAL :: fj,dfac,binomfac
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
        
        ! Speed of light in atomic units: c = 1/fine-structure-constant

        cl = 137.0359990740d0

        gDpole = 0.0d0

        IF ( d .EQ. 1 ) THEN
                !x-component:
                IF ( L1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,1) = gDpole(1,1) -L1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,1) = gDpole(1,1) + 2*alpha1*Dpole

                IF ( L2 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(2,1) = gDpole(2,1) -L2*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,1) = gDpole(2,1) + 2*alpha2*Dpole
                
                !y-component:
                IF ( M1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,2) = gDpole(1,2) -M1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,2) = gDpole(1,2) + 2*alpha1*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) -(M2+1)*Dpole
                
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) + 2*alpha2*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) + Dpole
                
                !z-component:
                IF ( N1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,3) = gDpole(1,3) -N1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,3) = gDpole(1,3) + 2*alpha1*Dpole

                IF ( L2 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(2,3) = gDpole(2,3) -N2*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,3) = gDpole(2,3) + 2*alpha2*Dpole
        ENDIF
        IF ( d .EQ. 2 ) THEN
                !x-component:
                IF ( L1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,1) = gDpole(1,1) -L1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,1) = gDpole(1,1) + 2*alpha1*Dpole

                IF ( L2 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(2,1) = gDpole(2,1) -L2*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,1) = gDpole(2,1) + 2*alpha2*Dpole
                
                !y-component:
                IF ( M1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,2) = gDpole(1,2) -M1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,2) = gDpole(1,2) + 2*alpha1*Dpole
        
                IF ( M2 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(2,2) = gDpole(2,2) -M2*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) + 2*alpha2*Dpole

                !z-component:
                IF ( N1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,3) = gDpole(1,3) -N1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,3) = gDpole(1,3) + 2*alpha1*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,3) = gDpole(2,3) -(N2+1)*Dpole
                
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,3) = gDpole(2,3) + 2*alpha2*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,3) = gDpole(2,3) + Dpole
        ENDIF
        IF ( d .EQ. 3 ) THEN
                !x-component:
                IF ( L1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1-1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,1) = gDpole(1,1) -L1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1+1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,1) = gDpole(1,1) + 2*alpha1*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,1) = gDpole(2,1) -(L2+1)*Dpole
                
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2+1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,1) = gDpole(2,1) + 2*alpha2*Dpole
                
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2-1,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,1) = gDpole(2,1) + Dpole
                
                !y-component:
                IF ( M1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1-1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,2) = gDpole(1,2) -M1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1+1,N1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,2) = gDpole(1,2) + 2*alpha1*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) -(M2+1)*Dpole
                
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2+1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) + 2*alpha2*Dpole

                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2-1,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,2) = gDpole(2,2) + Dpole
                
                !z-component:
                IF ( N1 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1-1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(1,3) = gDpole(1,3) -N1*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1+1,A,alpha1,L2,M2,N2,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(1,3) = gDpole(1,3) + 2*alpha1*Dpole

                IF ( L2 .GE. 1 ) THEN
                        CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2-1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                        gDpole(2,3) = gDpole(2,3) -N2*Dpole
                ENDIF
                CALL inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2+1,B,alpha2,.FALSE.,d,omega,t,Dpole)
                gDpole(2,3) = gDpole(2,3) + 2*alpha2*Dpole
        ENDIF


 END SUBROUTINE gradpriminhomefieldoverlap
