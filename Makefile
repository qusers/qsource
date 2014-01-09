# makefile for the Q package:
# Qprep5/Qdyn5/Qdyn5p/Qfep5/Qdum5/Qcalc5
# version 5.01, 2003-08-26

FC = mpiifort
FFLAGS = -O0 -g -check uninit -traceback
LD = $(FC)
LDFLAGS = -traceback

QfepSource = qfep.f90 mpiglob.f90 nrgy.f90 misc.f90 parse.f90

QprepSource = q_prep.f90 topo.f90 misc.f90 mpiglob.f90 parse.f90 prefs.f90 prep.f90 prmfile.f90 index.f90 mask.f90 trj.f90 sizes.f90 avetr.f90

QcalcSource = calc_base.f90 calc_chemscore.f90 calc_fit.f90 calc_geom.f90 calc_pmfscore.f90 calc_com_ke.f90 calc_com.f90 calc_rdf.f90 calc_rms.f90 calc_rmsf.f90 calc_entropy.f90 calc_nb.f90 calc_xscore.f90 eigen.f90 index.f90 mask.f90 maskmanip.f90 misc.f90 mpiglob.f90 nrgy.f90 parse.f90 prmfile.f90 qatom.f90 qcalc.f90 sizes.f90 topo.f90 trj.f90

QdynSource = md.f90 mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 qdyn.f90 sizes.f90 topo.f90 trj.f90 index.f90

QdumSource = md_dum.f90 mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 qdyn.f90 sizes.f90 topo.f90 trj.f90 index.f90

QdynpSource = md_mpi.f90 mask.f90 misc.f90 mpiglob.f90 nrgy.f90 prmfile.f90 qatom.f90 qdyn_mpi.f90 sizes.f90 topo.f90 trj.f90 index.f90

QfepObjects = $(QfepSource:.f90=.o)
QprepObjects = $(QprepSource:.f90=.o)
QcalcObjects = $(QcalcSource:.f90=.o)
QdynObjects = $(QdynSource:.f90=.o)
QdumObjects = $(QdumSource:.f90=.o)
QdynpObjects = $(QdynpSource:.f90=.o)

PROGRAMS = Qfep5 Qdyn5p Qdyn5 Qprep5 Qcalc5 Qdum5

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

###########################################################
#
# pseudo-targets
#
###########################################################

all:	$(PROGRAMS)

clean:
	rm -f *.o *.F90 *.mod *.M *.kmo *.il $(PROGRAMS)

###########################################################
#
# real build targets: programs
#
###########################################################

Qfep5: $(QfepObjects)
	$(LD) $(LDFLAGS) -o $@ $(QfepObjects) $(FLIBS)

Qprep5: $(QprepObjects)
	$(LD) $(LDFLAGS) -o $@ $(QprepObjects) $(FLIBS)

Qcalc5: $(QcalcObjects)
	$(LD) $(LDFLAGS) -o $@ $(QcalcObjects) $(FLIBS)

Qdyn5:	$(QdynObjects)
	$(LD) $(LDFLAGS) -o $@ $(QdynObjects) $(FLIBS)

Qdum5: $(QdumObjects)
	$(LD) $(LDFLAGS) -o $@ $(QdumObjects) $(FLIBS)

Qdyn5p: $(QdynpObjects)
	$(LD) $(LDFLAGS) -o $@ $(QdynpObjects) $(FLIBS)
