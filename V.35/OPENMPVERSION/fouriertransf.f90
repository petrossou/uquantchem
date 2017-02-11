SUBROUTINE fouriertransf(NPOINTS,TIMESTEP,MOMENT,EMAX)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOINTS
      DOUBLE PRECISION, INTENT(IN) :: TIMESTEP, MOMENT(NPOINTS),EMAX
      DOUBLE PRECISION :: DELTAOMEGA,OMEGA(NPOINTS),T,LAMDA
      INTEGER :: I,J
      DOUBLE PRECISION :: F(2,NPOINTS),MEAN,SS(NPOINTS)
      DOUBLE PRECISION, PARAMETER :: pi = 3.1415926535897932384626433832795028841970d0
      DOUBLE PRECISION :: SMEARMATRIX(NPOINTS,NPOINTS)

      DELTAOMEGA = 2*pi/(TIMESTEP*NPOINTS)

      MEAN = SUM(MOMENT)/(1.0d0*NPOINTS)

      LAMDA = ((TIMESTEP*NPOINTS/(2*pi))**2)*LOG(2.0d0)

      OPEN(444,FILE='ABSSPECTRUM.dat',ACTION='WRITE')
      
      DO I =0,NPOINTS-1
                OMEGA(I+1) = DELTAOMEGA*I
                F(1,I+1) = 0.0d0
                F(2,I+1) = 0.0d0
                DO J=0,NPOINTS-1
                        SMEARMATRIX(I+1,J+1)=EXP(-LAMDA*(DELTAOMEGA*(I-J))**2)
                        T = J*TIMESTEP
                        IF ( EMAX .NE. 0.0d0 ) THEN
                                F(1,I+1) = F(1,I+1) + 2.0d0*TIMESTEP*sin(OMEGA(I+1)*T)*(MOMENT(J+1)-MEAN)/EMAX
                                F(2,I+1) = F(2,I+1) + 2.0d0*TIMESTEP*cos(OMEGA(I+1)*T)*(MOMENT(J+1)-MEAN)/EMAX
                        ELSE
                                F(1,I+1) = F(1,I+1) + 2.0d0*TIMESTEP*sin(OMEGA(I+1)*T)*(MOMENT(J+1)-MEAN)
                                F(2,I+1) = F(2,I+1) + 2.0d0*TIMESTEP*cos(OMEGA(I+1)*T)*(MOMENT(J+1)-MEAN)
                        ENDIF
                ENDDO
      ENDDO   

      SS = MATMUL(SMEARMATRIX,F(1,:))
      
      DO I =0,(NPOINTS-MOD(NPOINTS,2))/2
            WRITE(444,'(3(E30.20))')OMEGA(I+1),DABS(F(1,I+1)),DABS(SS(I+1))
      ENDDO
     
      CLOSE(444)

END SUBROUTINE fouriertransf


