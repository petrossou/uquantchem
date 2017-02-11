SUBROUTINE inhomefieldoverlap(L1,M1,N1,A,alpha1,L2,M2,N2,B,alpha2,POVERLAP,d,omega,t,Dpole)
        ! This subroutine calculates the three dipole tensor 
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
        DOUBLE PRECISION, INTENT(OUT) :: Dpole
        DOUBLE PRECISION, INTENT(IN) :: alpha1,alpha2,A(3),B(3),omega,t
        LOGICAL, INTENT(IN) ::  POVERLAP
        DOUBLE PRECISION :: Sx,Sy,Sz,Sc,gama,phi,kv,cl,PA(3),PB(3),P(3),C
        INTEGER :: I,J,K,n,m
        DOUBLE PRECISION, EXTERNAL :: fj,dfac,binomfac
        DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
        
        ! Speed of light in atomic units: c = 1/fine-structure-constant

        cl = 137.0359990740d0

        Sx = 0.0d0
        Sy = 0.0d0
        Sz = 0.0d0
        Sc = 0.0d0

        phi = omega*t 
        kv = omega/cl
        gama = alpha1 + alpha2

        P = (A*alpha1 + B*alpha2)/gama

        C = (alpha1*alpha2/gama)*DOT_PRODUCT(A-B,A-B)

        PA = P-A
        PB = P-B

        IF ( d .EQ. 1 ) THEN

                phi = phi +kv*P(1)

                DO n=0,L1+L2
                        IF ( MOD(n,2) .EQ. 0 ) THEN
                                DO m=0,(n-MOD(n,2))/2
                                        Sx = Sx+fj(n,L1,L2,PA(1),PB(1))*( binomfac(n,2*m)*( (kv/(2*gama))**(n-2*m) )*cos( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*sin(phi)
                                ENDDO
                        ELSE
                                DO m=0,(n-MOD(n,2))/2
                                        Sx = Sx+fj(n,L1,L2,PA(1),PB(1))*( binomfac(n,2*m)*( (kv/(2*gama))**(n-2*m) )*sin( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*cos(phi)
                                ENDDO
                        ENDIF
                ENDDO
                Sx = sqrt(pi/gama)*Sx*EXP(-kv**2/(4*gama))

                DO J=0,(M1+M2+1-MOD(M1+M2+1,2))/2
                        Sy = Sy+fj(2*J,M1,M2+1,PA(2),PB(2))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sy = sqrt(pi/gama)*Sy
        
                DO J=0,(N1+N2-MOD(N1+N2,2))/2
                        Sz = Sz+fj(2*J,N1,N2,PA(3),PB(3))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sz = sqrt(pi/gama)*Sz
                
                DO J=0,(M1+M2-MOD(M1+M2,2))/2
                        Sc = Sc+fj(2*J,M1,M2,PA(2),PB(2))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sc = sqrt(pi/gama)*Sc

                IF (POVERLAP ) Sy = Sy + B(2)*Sc
        ENDIF
        
        IF ( d .EQ. 2 ) THEN

                phi = phi + kv*P(2)

                DO J=0,(L1+L2-MOD(L1+L2,2))/2
                        Sx = Sx+fj(2*J,L1,L2,PA(1),PB(1))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sx = sqrt(pi/gama)*Sx

                DO n=0,M1+M2
                        IF ( MOD(n,2) .EQ. 0 ) THEN
                                DO m=0,(n-MOD(n,2))/2
                                        Sy = Sy+fj(n,M1,M2,PA(2),PB(2))*binomfac(n,2*m)*( ( (kv/(2*gama))**(n-2*m) )*cos( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*sin(phi)
                                ENDDO
                        ELSE
                                DO m=0,(n-MOD(n,2))/2
                                        Sy = Sy+fj(n,M1,M2,PA(2),PB(2))*binomfac(n,2*m)*( ( (kv/(2*gama))**(n-2*m) )*sin( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*cos(phi)
                                ENDDO
                        ENDIF
                ENDDO
                Sy = sqrt(pi/gama)*Sy*EXP(-kv**2/(4*gama))
        
                DO J=0,(N1+N2+1-MOD(N1+N2+1,2))/2
                        Sz = Sz+fj(2*J,N1,N2+1,PA(3),PB(3))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sz = sqrt(pi/gama)*Sz
                
                DO J=0,(N1+N2-MOD(N1+N2,2))/2
                        Sc = Sc+fj(2*J,N1,N2,PA(3),PB(3))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sc = sqrt(pi/gama)*Sc

                IF ( POVERLAP ) Sz = Sz + B(3)*Sc
        ENDIF

        IF ( d .EQ. 3 ) THEN
                
                phi = phi + kv*P(3)
                
                DO J=0,(L1+L2+1-MOD(L1+L2+1,2))/2
                        Sx = Sx+fj(2*J,L1,L2+1,PA(1),PB(1))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sx = sqrt(pi/gama)*Sx

                DO J=0,(M1+M2-MOD(M1+M2,2))/2
                        Sy = Sy+fj(2*J,M1,M2,PA(2),PB(2))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sy = sqrt(pi/gama)*Sy
        
                DO n=0,N1+N2
                        IF ( MOD(n,2) .EQ. 0 ) THEN
                                DO m=0,(n-MOD(n,2))/2
                                        Sz = Sz+fj(n,N1,N2,PA(3),PB(3))*binomfac(n,2*m)*( ( (kv/(2*gama))**(n-2*m) )*cos( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*sin(phi)
                                ENDDO
                        ELSE
                                DO m=0,(n-MOD(n,2))/2
                                        Sz = Sz+fj(n,N1,N2,PA(3),PB(3))*binomfac(n,2*m)*( ( (kv/(2*gama))**(n-2*m) )*sin( (pi/2.0d0)*(n-2*m) )*dfac(2*m-1)/((2.0d0*gama)**m) )*cos(phi)
                                ENDDO
                        ENDIF
                ENDDO
                Sz = sqrt(pi/gama)*Sz*EXP(-kv**2/(4*gama))
                
                DO J=0,(L1+L2-MOD(L1+L2,2))/2
                        Sc = Sc+fj(2*J,L1,L2,PA(1),PB(1))*dfac(2*J-1)/((2.0d0*gama)**J)
                ENDDO
                Sc = sqrt(pi/gama)*Sc

                IF ( POVERLAP ) Sx = Sx + B(1)*Sc
        ENDIF
        
        Dpole = exp(-C)*Sx*Sy*Sz

 END SUBROUTINE inhomefieldoverlap
