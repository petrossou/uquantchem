SUBROUTINE stoexponent(PSI,deltar,LAMDA,A,GAMA,BETA,POLY,rc,USEGTO,CORRALCUSP)
        ! Finds the exponent, LAMDA, and coeficient, A,  of a slater type orbital, A*exp(-LAMDA*r) 
        ! by fiting it to the gaussian basis function, gto(r), at the inflexion
        ! point rc, i.e at the point where gto'(rc) = minimum. The
        ! exponent, LAMDA, and coeficient, A are calculated from the
        ! continuity conditions:
        ! (1)  gto(rc) = A*exp(-LAMDA*rc)
        ! (2)  gto'(rc) = -A*LAMDA*exp(-LAMDA*rc)
        ! leading to LAMDA = - gto'(rc)/gto(rc), and A = gto(rc)*exp(-gto'(rc)*rc/gto(rc))
        USE datatypemodule
        IMPLICIT NONE
        TYPE(BASFUNCT), INTENT(IN) :: PSI
        DOUBLE PRECISION, INTENT(IN) :: deltar
        DOUBLE PRECISION, INTENT(INOUT) :: rc
        DOUBLE PRECISION, INTENT(OUT) :: LAMDA,A,GAMA,BETA,POLY(6)
        LOGICAL, INTENT(IN) :: CORRALCUSP
        LOGICAL, INTENT(OUT) :: USEGTO
        INTEGER :: I,J,SEARCHMESH,NPC
        DOUBLE PRECISION :: DR, r1,r2,signn,NODE,P(4),rcposible(3),C,K1,K2
        DOUBLE PRECISION, EXTERNAL :: gto, gtop, gtopp,gtoppp
        LAMDA = 0.0d0
        GAMA = 0.0d0
        BETA = 0.0d0
        POLY = 0.0d0
        A = 0.0d0
        NPC = 1
        IF ( CORRALCUSP ) NPC = 0
        IF ( PSI%NPRIM .GT. NPC ) THEN
                !=============================================
                ! SEARCHING FOR THE LOCATION OF THE FIRST NODE
                !=============================================
                USEGTO = .FALSE.
                NODE = -1000.0d0
                SEARCHMESH = 50000
                DR = 5.0d0/SEARCHMESH
                r1 = DR
                r2 = r1 + DR
                signn = 1.0d0
                DO WHILE ( signn .GT. 0.0d0 .AND. I .LT. SEARCHMESH ) 
                         signn = gto(PSI,r1)*gto(PSI,r2)
                         IF ( signn .LE. 0.0d0 ) NODE = 0.50d0*(r1+r2)
                         r1 = r1+DR
                         r2 = r2 +DR
                ENDDO
                IF ( signn .LT. 0.0d0 ) THEN
                        rcposible(1) = NODE
                ELSE
                        rcposible(1) = -1.0d0
                ENDIF
                !======================================================
                ! SEARCHING FOR THE INFLEXION POINT WHERE gto'(r) = min
                !======================================================
                SEARCHMESH = 50000
                DR = 5.0d0/SEARCHMESH
                r1 = DR
                r2 = r1 + DR
                signn = 1.0d0
                DO WHILE ( signn .GT. 0.0d0 .AND. r2 .LT. 5.0d0 ) 
                        signn = gtopp(PSI,r1)*gtopp(PSI,r2)
                        r1 = r1 + DR
                        r2 = r2 + DR
                ENDDO
                IF ( signn .LT. 0.0d0 ) THEN
                        rcposible(2) = 0.5d0*(r1 + r2 - 2*DR )
                ELSE
                        rcposible(2) = -1.0d0
                ENDIF
                !=======================================================================
                ! SEARCHING FOR THE POINT FURTHEST AWAY FROM r=0 FOR WHICH gto'''(r) = 0 
                !=======================================================================
                SEARCHMESH = 10000
                DR = 1.0d0/SEARCHMESH
                r1 = 1.0d0
                r2 = r1 - DR
                signn = 1.0d0
                DO WHILE ( signn .GT. 0.0d0 .AND. r2 .GE. 0.0d0 )
                        signn = gtoppp(PSI,r1)*gtoppp(PSI,r2)
                        r1 = r1 - DR
                        r2 = r2 - DR
                ENDDO
                
                IF ( signn .LT. 0.0d0 ) THEN
                        rcposible(3) = 0.5d0*(r1 + r2 + 2*DR )
                ELSE
                        rcposible(3) = -1.0d0
                ENDIF
                
                ! WE USE A FIXED VALUE FOR rc!
                ! Default is rc =0.10d0
                
                
                IF ( rcposible(1) .LT. rc  .AND. rcposible(1) .GT.  0.0d0 ) THEN 
                        WRITE(*,*)'CUSP FITTING RADIUS CHANGED TO (1/2)*(distance to first node)'
                        WRITE(*,*)'SINCE LOCATION OF NODE =',rcposible(1),'<rc=',rc
                        rc = rcposible(1)/2.0d0
                ELSE
                        call leastsq(rc,1000,PSI,A,LAMDA,C)
                        call interpol(rc,deltar,PSI,A,LAMDA,P)

                        K1 = gtop(PSI,rc+deltar) - P(1)*(rc+deltar) - (P(2)/2.0d0)*(rc+deltar)**2 - (P(3)/3.0d0)*(rc+deltar)**3 - (P(4)/4.0d0)*(rc+deltar)**4
                        K2 = gto(PSI,rc+deltar) - K1*(rc+deltar) - (P(1)/2.0d0)*(rc+deltar)**2 - (P(2)/6.0d0)*(rc+deltar)**3 - (P(3)/12.0d0)*(rc+deltar)**4 - (P(4)/20.0d0)*(rc+deltar)**5

                        BETA = K1       +      P(1)*rc      +      (P(2)/2.0d0)*rc**2 + (P(3)/3.0d0)*rc**3 +  (P(4)/4.0d0)*rc**4 + A*LAMDA*EXP(-LAMDA*rc)
                        GAMA = K2 + K1*rc   + (P(1)/2.0d0)*rc**2 + (P(2)/6.0d0)*rc**3 + (P(3)/12.0d0)*rc**4 + (P(4)/20.0d0)*rc**5 - BETA*rc - A*EXP(-LAMDA*rc)

                        POLY(1) = K2
                        POLY(2) = K1
                        POLY(3) = P(1)
                        POLY(4) = P(2)
                        POLY(5) = P(3)
                        POLY(6) = P(4)
                ENDIF

        ELSE
                USEGTO = .TRUE.
        ENDIF
END SUBROUTINE stoexponent
