###############################################################################
#                       Q v.5 Makefile for Linux                              #
###############################################################################

### -i-static -static
# Put your favorite compilers and flags here.
FC=		ifort
CC=		gcc
MPIFC=		mpif90
FC_OPT=		-O3 -ip -ipo -unroll -i-static 
CC_OPT=		-O3
#MPIFLAGS=	-f90=${FC}
OMPI_FC=	${FC}
OMPI_CC=	${CC}	
FFLAGS=		-fstack-protector ${FC_OPT}
CCFLAGS=	${CC_OPT}
OMPI_FCFLAGS=	${FFLAGS}
OMPI_CXXFLAGS=	${CCFLAGS}
FLIBS=		-Vaxlib
FPP_FLG=	-fpp  
# 
###

###
# Display options
default: what
###

###
# Targets
all:	Qfep5 Qprep5 Qdyn5 Qdum5 Qcalc5

debug:	
	@make FFLAGS='-debug -traceback -g -Wcheck'         \
	Qfep5 Qprep5 Qdyn5 Qdum5 Qcalc5

mpi:	Qdyn5p

mpi_debug:	
	@make FFLAGS='-debug -traceback -g -Wcheck'         \
	Qdyn5p

clean:
	-rm -f *.o *.F90 *.mod *.M *.kmo *.il

nuke:
	-rm -f *.o *.F90 *.mod *.M *.kmo *.il Qfep5 Qdyn5p Qdyn5 Qprep5 Qcalc5 Qdum5

Qcalc5 Qdyn5 Qdum5 Qdyn5p Qprep5 Qfep5: misc.o mpiglob.o 

Qcalc5 Qdyn5 Qdum5 Qdyn5p Qprep5: mask.o prmfile.o sizes.o topo.o trj.o index.o 

Qcalc5 Qdyn5 Qdum5 Qdyn5p Qfep5: nrgy.o

Qcalc5 Qdyn5 Qdum5 Qdyn5p: qatom.o  

Qcalc5 Qprep5 Qfep5: parse.o 

Qcalc5 Qprep5: maskmanip.o 

Qdyn5 Qdum5: qdyn.o

Qprep5: q_prep.o prefs.o prep.o avetr.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 

Qfep5: qfep.o 
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 

Qcalc5: calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o \
        calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o \
        calc_entropy.o calc_nb.o calc_xscore.o eigen.o qcalc.o 
	@cp qcalc.f90 qcalc.F90
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 

Qdyn5 : md.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 

Qdum5 : md_dum.o 
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@ 

Qdyn5p : md_mpi.o qdyn_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@ 

###
# Object modules
avetr.o: avetr.f90 prep.o
	${FC} ${FFLAGS} -c $<

calc_base.o:calc_base.f90 topo.o
	${FC} ${FFLAGS} -c $<

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${FC} ${FFLAGS} -c $<

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_geom.o:calc_geom.f90 calc_base.o
	${FC} ${FFLAGS} -c $<

calc_nb.o:calc_nb.f90 calc_base.o maskmanip.o parse.o
	${FC} ${FFLAGS} -c $<

calc_pmfscore.o: calc_pmfscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_xscore.o: calc_xscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_rdf.o:calc_rdf.f90 calc_base.o parse.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rms.o:calc_rms.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com_ke.o:calc_com_ke.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com.o:calc_com.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rmsf.o:calc_rmsf.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

eigen.o:eigen.f90
	${FC} ${FFLAGS} -c $<

invsqrt_q.o:invsqrt_q.c
	$(CC) $(CCFLAGS) -c $<

index.o:index.f90
	${FC} ${FFLAGS} -c $<

mask.o:	mask.f90 topo.o
	${FC} ${FFLAGS} -c $<

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${FC} ${FFLAGS} -c $<

md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o 
	@cp md.f90 md.F90
	${FC} ${FFLAGS} ${FPP_FLG} -c md.F90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o topo.o
	@cp md.f90 md.F90
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.F90 -o md_dum.o 

md_mpi.o: md.f90 mpiglob.o qatom.o sizes.o topo.o trj.o
	@cp md.f90 md.F90
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.F90 -o md_mpi.o 

misc.o: misc.f90 sizes.o
	${FC} ${FFLAGS} -c $<

mpiglob.o: mpiglob.f90 sizes.o nrgy.o
	${FC} ${FFLAGS} -c $<

nrgy.o: nrgy.f90 sizes.o
	${FC} ${FFLAGS} -c $<

parse.o: parse.f90 misc.o
	${FC} ${FFLAGS} -c $<

prefs.o: prefs.f90
	${FC} ${FFLAGS} -c $<

prep.o: prep.f90 maskmanip.o sizes.o parse.o prmfile.o trj.o index.o prefs.o 
	${FC} ${FFLAGS} -c $<

prmfile.o: prmfile.f90 misc.o mpiglob.o
	${FC} ${FFLAGS} -c $<

q_prep.o: q_prep.f90 prep.o avetr.o
	${FC} ${FFLAGS} -c $<

qatom.o: qatom.f90 misc.o nrgy.o prmfile.o sizes.o index.o topo.o
	${FC} ${FFLAGS} -c $<

qcalc.o: qcalc.f90 calc_chemscore.o calc_pmfscore.o calc_xscore.o trj.o calc_base.o calc_rms.o calc_fit.o calc_geom.o
	@cp qcalc.f90 qcalc.F90 
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.F90

qdyn.o: qdyn.f90 md.o mpiglob.o
	@cp qdyn.f90 qdyn.F90
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn.F90

qdyn_dum.o: qdyn.f90 md.o mpiglob.o
	@cp qdyn.f90 qdyn.F90
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c qdyn.F90 -o qdyn_dum.o 

qdyn_mpi.o: qdyn.f90 md_mpi.o mpiglob.o
	@cp qdyn.f90 qdyn.F90
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
        -c qdyn.F90 -o qdyn_mpi.o

qfep.o: qfep.f90 nrgy.o parse.o
	@cp qfep.f90 qfep.F90
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.F90

sizes.o: sizes.f90
	${FC} ${FFLAGS} -c $<

topo.o:  topo.f90 misc.o mpiglob.o sizes.o
	${FC} ${FFLAGS} -c $<

trj.o:  trj.f90 mask.o misc.o 
	${FC} ${FFLAGS} -c $<
#
###

###
# Build instructions
what:
	@echo "usage: make <target>"
	@echo
	@echo "<target> is one of:"
	@echo "all              Everything except for the MPI parallel version of Q (Qdyn5p)"
	@echo "debug            'all' but with debug information (stacktraces etc)."
	@echo "mpi              Qdyn5p using the currently loaded mpi library"
	@echo "mpi_debug        'mpi' but with debug information (stacktraces etc)."
	@echo
#
###
