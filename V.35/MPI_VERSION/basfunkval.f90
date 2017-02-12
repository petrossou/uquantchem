FUNCTION basfunkval(PSI,R)
      ! Calculates the value of a real basis-function at the spatial
      ! point R. I.e, The function = Psi(R)
      USE datatypemodule
      IMPLICIT NONE
      DOUBLE PRECISION :: basfunkval
      TYPE(BASFUNCT) :: PSI
      DOUBLE PRECISION :: R(3),RP(3),RP2,faktor
      INTEGER :: I
     
      basfunkval = 0.0d0
      RP = R - PSI%R
      RP2 = DOT_PRODUCT(RP,RP)
      faktor = ( RP(1)**PSI%L(1) )*( RP(2)**PSI%L(2) )*( RP(3)**PSI%L(3) )
      
      IF ( DABS(faktor)*EXP(-MINVAL(PSI%EXPON)*RP2) .GT. 1.0E-10 ) THEN
          DO  I=1,PSI%NPRIM
             basfunkval = basfunkval + PSI%CONTRCOEFF(I)*PSI%PRIMNORM(I)*PSI%NORM*faktor*EXP(-PSI%EXPON(I)*RP2)
          ENDDO
      ENDIF
     
       END FUNCTION basfunkval

