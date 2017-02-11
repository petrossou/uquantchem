MODULE datatypemodule
      IMPLICIT NONE
      TYPE BASFUNCT
              INTEGER :: NPRIM
              INTEGER :: ATYPE
              LOGICAL :: HYBRID
              DOUBLE PRECISION :: CONTRCOEFF(100)
              DOUBLE PRECISION :: R(3)
              DOUBLE PRECISION  :: EXPON(100)
              DOUBLE PRECISION :: PRIMNORM(100)
              INTEGER :: L(3)
              DOUBLE PRECISION :: NORM
              DOUBLE PRECISION  :: STOEXPON
              DOUBLE PRECISION :: A
              DOUBLE PRECISION :: BETA
              DOUBLE PRECISION :: GAMA
              DOUBLE PRECISION :: RC
              DOUBLE PRECISION :: DR
              DOUBLE PRECISION :: P(6)
      END TYPE BASFUNCT

      TYPE ATOM
              INTEGER :: Z
              DOUBLE PRECISION :: R(3)
              DOUBLE PRECISION :: M
              INTEGER :: NBAS
              TYPE(BASFUNCT), ALLOCATABLE :: PSI(:)
      END TYPE ATOM

      TYPE BASIS
              INTEGER :: NBAS
              TYPE(BASFUNCT), ALLOCATABLE :: PSI(:)
      END TYPE BASIS
END MODULE datatypemodule
              
