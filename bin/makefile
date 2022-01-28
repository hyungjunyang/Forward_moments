# 
# Makefile for smelib testing
#

SMELIB_DIR = ..
include $(SMELIB_DIR)/makefile.def

SRCS		= Main.cc 
EXECS		= SME3D
OBJS		= Main.o

all: progs

progs: $(EXECS)

clean:
	/bin/rm -f $(OBJS) $(EXECS) out *.out

wipe: 
	/bin/rm -f $(OBJS) $(EXECS) out *.out

Main.o : Main.cc  

SME3D : Main.o
	$(CCC) $(CCCFLAGS)  -o SME3D Main.o $(LDFLAGS) $(CLUB) $(LAPACK_BLAS_LIBS) $(FOBJS)

