##########################################################
#  Compilation instructions for the uquantchem V.28 code #
##########################################################

To compile the uquantchem code you need to have the 
gfortran compiler installed on your system. The lapack 
and Blas libraries needed will be compiled together 
with the code.

In the below instructions it is assumed that the 
UQUANTCHEM.tar.gz package is located at the arbitrary 
directory "HOME"

#==============================
# Compiling the openmp-version:
#==============================
To compile the openmp-version of the code open a console 
and execute the following commands on the command line:

# Go to the HOME directory:
(1)  cd HOME

# Unzip the package:
(2) gunzip UQUANTCHEM.tar.gz 

# unpack the package:
(3) tar -xvf UQUANTCHEM.tar

# Change directory to the src directory where the source 
# code is located:
(4) cd UQUANTCHEM/OPENMPVERSION/

# Compile the code together with the lapack and blas libraries.
# this might take a while
(5) make all


If successfull you should have ended up with an executable
named "uquantchem.omp" located in HOME/UQUANTCHEM/src

#==============================
# Compiling the serial version:
#==============================
In order to compile the serial version follow the exact same
steps as are listed for the compilation of the openmp version,
except for step (4) where you issue the command:

cd UQUANTCHEM/SERIALVERSION/

instead.

#===========================
# Compiling the mpi version:
#===========================

# Go to the HOME directory:
(1)  cd HOME

# Unzip the package:
(2) gunzip UQUANTCHEM.tar.gz 

# unpack the package:
(3) tar -xvf UQUANTCHEM.tar

# Change directory to the src directory where the source 
# code is located:
(4) cd UQUANTCHEM/MPI_VERSION/

# Edit the "Makefile" so that it complies with the system
# specifications such as the type of mpi compiler and the 
# paths to the blas, lapack and scalapack libraries.
(5) Edit Makefile

# Compile the code:
(6) make 

############################################################
#  The prepared test calculations comming with the package #
############################################################

There are 7 different prepared test calculations that can 
be run. The input-files and respective output files 
are located in the subdirectories
TESTS/RUN1,TESTS/RUN2, ... , TESTS/RUN7
of the UQUANTCHEM/ directory. The respective pre-calulated 
output files are named OUT.1.dat, OUT.2.dat, ... , OUT.7.dat.

In order to run these test calculations for the openmpi version
issue the following comands: 
(for example the test calculation in the dir TESTS/RUN1)

(1) cd UQUANTCHEM/TESTS/RUN1

(2) ./run
 
##############################################################

Below follow a short description of these test calculations:

RUN1: An unrestricted Hartree-Fock total energy calculation of 
      a water molecule using a 6-31G** basis-set.

RUN2: A PBE total energy calculation of 
      a water molecule using a 6-31G** basis-set.

RUN3: A B3LYP total energy and force calculation of 
      a water molecule using a 6-31G** basis-set.

RUN4: A MP2  total energy calculation of 
      a water molecule using a 6-31G** basis-set.

RUN5: A CISD  total energy calculation of 
      a Be atom using a 6-31G** basis-set.

RUN6: DQMC total energy calculation of 
      a He atom using a cc-pVTZ basis-set.

RUN7: Unrestricted Hartree-Fock Geometry optimization of a
      a water molecule using a STO-3G basis-set.

##################################################################
If you instead want to run the serial version when performing 
the test calculations, instead of the command "./run" execute 
the command ../../SERIALVERSION/uquantchem.s

