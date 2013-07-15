#################################################################################
# This makefile is for building the MPI version of uquantchem on PDS's Lindgren #
#     matter. Before running make load the appropriate libraries by running:    #
#               module swap PrgEnv-pgi PrgEnv-intel                             #
#              using the ifort mkl lapack library, BLACS and SCALAPACK          #
#              ( good debugg flagg: -check all -traceback )                     #
#               Petros Souvatzis: petros.souvatzis@fysik.uu.se                  #
#                       Wed Apr 25 09:23:31 CEST 2012                           #
#################################################################################

FC = ftn

FFLAGS = 

F77FLAGS = 

LIB =

OFLAG  = 

SOURCE = datatypemodule.o exchcorrmodule.o random.o Afunc.o RHF.o Fnurec.o Fnu.o CISD.o getKv.o getK.o getJv.o getJ.o gammaf.o fj.o fac.o  \
	 gto.o gtop.o gtopp.o eeints.o diaghHF.o diagh.o dfac.o checkinputfile.o chargedenssave.o exciteornot.o gridsetup.o diagscalapack.o \
	 hforbitalsave.o gtoppp.o binomfac.o basfunkval.o URHF.o primeeintegral.o potential.o overlap.o oneint.o normalize.o ROOTWEIGHT5.o \
	 leastsq.o makedens.o kinetic.o ijkl.o homolumosave.o gettotalbasis.o rysquad.o rho.o readin.o ROOTWEIGHT4.o blockset.o \
	 interpol.o gradbasfunkval.o laplacebasfunkval.o hforbitalval.o gradhforbitalval.o laplacehforbitalval.o det.o \
	 stoexponent.o slaterdet.o laplaceslaterdet.o gradslaterdet.o jastrowup.o jastrowdown.o jastrowud.o jastrow.o \
	 DIIS.o preparedifusion.o gradjastrow.o laplacejastrow.o guideforce.o laplacetrialfnk.o trialfnk.o EL.o dqmc.o \
	 findclosestatom.o ROOTWEIGHTMAX3.o readbasis.o primpotential.o primoverlap.o primkinetic.o MP2.o uquantchem.o \
	 gradprimoverlap.o gradprimkinetic.o gradprimpotential.o gradprimeeintegral.o forces.o dAfunc.o relax.o linesearchmin.o \
	 massa.o moleculardynamics.o moleculardynamicssoft.o invert.o vxc.o vxcalt.o dvxc.o ngbasfunkval.o nggradbasfunkval.o nglaplacebasfunkval.o \
	 getvxcr.o atomicradii.o lebedev.o chebgauss.o sofmu.o pvoronoi.o getvxc.o excdens.o excdr.o exc.o DFT.o quadcheck.o TRACE.o \
	 countee.o hessianbasfunk.o hessianrho.o dftforcedens.o getxcforce.o pnonorm.o tmu.o gradpvoronoi.o relaxf.o

uquantchem:  $(SOURCE) 
	$(FC) $(FFLAGS) $(SOURCE) $(LIB) -o uquantchem.mpi

clean: 
	-/bin/rm *.o *~ *.mod uquantchem.mpi

%.o: %.f
	$(FC) $(LIB) $(F77FLAGS) $(OFLAG) -c $<

%.o: %.f90
	$(FC) $(LIB) $(FFLAGS) $(OFLAG) -c $<

lebedev.o: lebedev.f90
	$(FC) -O0 -c lebedev.f90

