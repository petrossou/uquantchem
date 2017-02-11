SUBROUTINE gethole(F,S,NB,Ne,ORBINDEX,PHOLE,P,Pp)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Ne,ORBINDEX
      DOUBLE PRECISION, INTENT(IN) :: F(NB,NB), S(NB,NB),P(NB,NB),Pp(NB,NB)
      DOUBLE PRECISION, INTENT(INOUT) :: PHOLE(NB,NB)
      DOUBLE PRECISION :: Qs(NB,NB),Eholeigen(NB),DP,PTOL,PHOLEold(NB,NB),FTRANSF(NB,NB),EN(NB,NB)
      DOUBLE PRECISION :: U(NB,NB),C(NB,NB),MIX,Q(NB,NB),Ch(NB,NB),C2(NB,NB),VEC(NB),DELTAP,Psave(NB,NB)
      INTEGER :: I,J,K,INFO,MAXITER
      
      Psave = PHOLE
      PTOL = 1.0E-8
      MAXITER = 500
      MIX = 0.50d0

      EN = 0.0d0
      DO I=1,NB
        EN(I,I) = 1.0d0
      ENDDO


      I = 0
      PHOLEold = PHOLE
      DELTAP = 2.0d0*PTOL
      DO WHILE ( I .LT. MAXITER .AND. DELTAP .GT. PTOL ) 
                Qs = EN - Pp - PHOLE
                U = MATMUL(MATMUL(S,P),(EN - MATMUL(S,Qs)))
                FTRANSF = MATMUL(U,MATMUL(F,TRANSPOSE(U)))
                CALL diaghHF( FTRANSF,S,NB,Eholeigen,Ch,INFO)
                IF ( INFO .EQ. 0 ) THEN
                        C = 0.0d0
                        IF ( NB+(ORBINDEX-Ne) .LE. NB .AND. NB+(ORBINDEX-Ne) .GE. 1 ) THEN
                                C(:,NB+(ORBINDEX-Ne)) = Ch(:,NB+(ORBINDEX-Ne))
                                CALL makedens(C,NB,PHOLE)
                        ENDIF
                ENDIF
                DELTAP = sqrt(DOT_PRODUCT(reshape(PHOLE-PHOLEold,(/NB**2/)),reshape(PHOLE-PHOLEold,(/NB**2/))))
                PHOLEold = PHOLE
                PHOLE = PHOLE*MIX + PHOLEold*(1.0d0-MIX)
                I = I+1
      ENDDO

      IF ( I .EQ. MAXITER ) PHOLE = Psave
END SUBROUTINE gethole


