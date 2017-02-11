SUBROUTINE makedensT(TOLDNe,Cup,Cdown,Eu,Ed,ETEMP,NB,Ne,Pup,Pdown,mu,ENTROPY)
      ! This subroutine calculates the density matrices, Pup, Pdown, at a given
      ! temperature employing the fermi-dirac occupation of the electronic
      ! states at the temperature, ETEMP. Appart from the density matrices the 
      ! chemical potential, mu, and the entropy, ENTROPY, are also given as output.
      ! Shalf, is the Lowdin orhogonalization matrix that enables us to work 
      ! in the orthogonal representation. Ne = number of electrons.
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NB,Ne
      DOUBLE PRECISION, INTENT(INOUT) :: Cup(NB,NB),Cdown(NB,NB),mu
      DOUBLE PRECISION, INTENT(IN) ::  Eu(NB),Ed(NB),ETEMP
      DOUBLE PRECISION, INTENT(OUT) :: Pup(NB,NB),Pdown(NB,NB),ENTROPY
      DOUBLE PRECISION :: Pups(NB,NB),Pdowns(NB,NB),COEFFus(NB,NB),COEFFds(NB,NB),OCC,DOCC,OCCsave,DDOCC
      DOUBLE PRECISION :: SHI(NB,NB),DNel,DNels,COEFFu(NB,NB),COEFFd(NB,NB),ID(NB,NB),PT(NB,NB),muo,fu,fd,musave
      INTEGER :: I,J,K,Neu,Ned,ITER
      DOUBLE PRECISION, PARAMETER :: KB = 0.000003166811524300030d0 !  Boltzmanns constant in au/K
      DOUBLE PRECISION, INTENT(IN) :: TOLDNe                        !  Tolerance for the Newon-Raphsson search of the chemical constant
      DOUBLE PRECISION, EXTERNAL :: TRACE


      ! Calculating the index's of the highest occupied spin-orbitals
      Neu  = ( Ne - MOD(Ne,2) )/2
      Ned  = ( Ne + MOD(Ne,2) )/2
      
      DNel = 1.0d0
      DNels = 100000000.0

      ! Calculating the starting guess for the chemical potential:
      IF ( MOD(Ne,2) .EQ. 0 ) THEN
           mu = 0.250d0*(Eu(Neu) + Eu(Neu+1) + Ed(Ned) + Ed(Ned+1) )
      ELSE
           mu = 0.50d0*(Eu(Ned) + Ed(Ned))
      ENDIF

      muo = mu
      
      ITER = 0

      ! Here starts the Newton-Raphsson routine to find the chemical constant 
      ! and eventually the density matrices.

      DO WHILE ( DABS(DNel) .GT.  TOLDNe .AND. ITER .LE. 100000 )
                ! CAlculating the density matrices and the entropy
                !ENTROPY = 0.0d0
                OCC = 0.0d0
                DOCC = 0.0d0
                DDOCC = 0.0d0
                DO I=0,NB-1
                        fu = 1.0d0/( 1.0d0 + EXP( (Eu(NB-I)-mu)/(KB*ETEMP) ) )
                        fd = 1.0d0/( 1.0d0 + EXP( (Ed(NB-I)-mu)/(KB*ETEMP) ) )
                        OCC = OCC + fu + fd
                        DOCC = DOCC + (1.0d0/(KB*ETEMP))*( fu*(1.0d0-fu) + fd*(1.0d0-fd) )
                ENDDO
                
                ! Calculating the new chemical potential through Newton-Raphsson:
                ! ( See my notes in the blue note book, p.37 equation 127)
                DNel =  OCC - Ne 
                
                ! Saving the density matrix if the approximation of the 
                ! chemical potential has been improved
                IF ( DABS(DNel) .LE. DABS(DNels) ) THEN
                        musave = mu
                        DNels = DNel
                        OCCsave = OCC
                ENDIF
                
                IF ( DABS(DOCC) .GT. 0.0d0 ) THEN 
                        mu = muo - DNel/DOCC
                        muo = mu
                ENDIF

                ITER = ITER + 1
                
      ENDDO

      ! Here we transform the density-matrices to there non-ortogonal representation:
      ! (See Eqn 113, p.32 in my blue notebook.)
      mu = musave
      ENTROPY = 0.0d0
      
      DO I=1,NB
            fu = 1.0d0/( 1.0d0 + EXP( (Eu(I)-mu)/(KB*ETEMP) ) )
            fd = 1.0d0/( 1.0d0 + EXP( (Ed(I)-mu)/(KB*ETEMP) ) )
            Cup(:,I) = Cup(:,I)*fu**(1.0d0/2.0d0)
            Cdown(:,I) = Cdown(:,I)*fd**(1.0d0/2.0d0)
            if ( fu .NE. 1.0d0 ) ENTROPY = ENTROPY -KB*(  (1.0d0-fu)*log(DABS(1.0d0-fu)) )
            if ( fd .NE. 1.0d0 ) ENTROPY = ENTROPY -KB*(  (1.0d0-fd)*log(DABS(1.0d0-fd)) )
      ENDDO

      DO I=0,NB-1
            fu = 1.0d0/( 1.0d0 + EXP( (Eu(NB-I)-mu)/(KB*ETEMP) ) )
            fd = 1.0d0/( 1.0d0 + EXP( (Ed(NB-I)-mu)/(KB*ETEMP) ) )
            if ( fu .NE. 0.0d0 ) ENTROPY = ENTROPY -KB*( fu*log(DABS(fu))  )
            if ( fd .NE. 0.0d0 ) ENTROPY = ENTROPY -KB*( fd*log(DABS(fd))  )
      ENDDO
      
      ! Finally calculating the temperatur smeared density matrices
      Pup = MATMUL(Cup,TRANSPOSE(Cup))
      Pdown = MATMUL(Cdown,TRANSPOSE(Cdown))
      PT = Pup + Pdown
END SUBROUTINE makedensT
           


