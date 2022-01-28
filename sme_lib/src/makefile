#
# Makefile for smelib src
#

SMELIB_DIR = ../..
include $(SMELIB_DIR)/makefile.def 

SME_LIB_A = ../../lib/smelib.a

SRCS = 	\
	     Common.cc \
             Grid.cc \
             Region.cc \
             Covariance.cc \
             ExpCovariance.cc \
             GaussCovariance.cc \
             Perm.cc \
             Fluid.cc \
             Control.cc \
             Well.cc \
             Eqn.cc \
             Solver.cc \
             P0Eqn.cc \
             CYPEqn.cc \
             P2Eqn.cc \
             CPPEqn.cc \
             FullLUSolver.cc \
             LapackBandSolver.cc \
             LapackFullSolver.cc \
             LapackSymBandSolver.cc \
             Velocity.cc \
             BLSolution.cc \
             FractionFlow.cc \
             Saturation.cc \
             Flow.cc \
             Production.cc \
             Domain.cc \
             Block.cc \
             Particle.cc \
             ParticleTrack.cc \
             Point.cc \
             WellArc.cc \
             MCTrvl.cc \
             gauleg.cc \
             pdf.cc \
             MC2SME.cc \
             MasterPoint.cc \
             KrigWeigh.cc \
             ran3.cc \
	     quadrature.cc \

OBJS    =  \
	     Common.o \
             Grid.o \
             Region.o \
             Covariance.o \
             ExpCovariance.o \
             GaussCovariance.o \
             Perm.o \
             Fluid.o \
             Control.o \
             Well.o \
             Eqn.o \
             Solver.o \
             P0Eqn.o \
             CYPEqn.o \
             P2Eqn.o \
             CPPEqn.o \
             FullLUSolver.o \
             LapackBandSolver.o \
             LapackFullSolver.o \
             LapackSymBandSolver.o \
             Velocity.o \
             BLSolution.o \
             FractionFlow.o \
             Saturation.o \
             Flow.o \
             Production.o \
             Domain.o \
             Block.o \
             Particle.o \
             ParticleTrack.o \
             Point.o \
             WellArc.o \
             MCTrvl.o \
             gauleg.o \
             pdf.o \
             MC2SME.o \
             MasterPoint.o \
             KrigWeigh.o \
             ran3.o \
	     quadrature.o \

all: lib

srcs: $(SRCS)
lib: $(SME_LIB_A)

$(SME_LIB_A) : $(OBJS)
	@ ar rv $(SME_LIB_A) $?
	@ case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(SME_LIB_A);; esac

clean: 
	/bin/rm -f $(OBJS) 

wipe:
	/bin/rm -f $(OBJS) $(SME_LIB_A) 
