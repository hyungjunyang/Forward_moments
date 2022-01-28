
# Important: use optimization (base optimization as a minimum)
#            flags for compilation.
#
# using flux machine to build lib
FC = /vend/pgi/5.2-4/linux86-64/bin/pgf90

LM_LICENSE_FILE=/vend/pgi/license.dat-flux

FFLAG = -c -O4 -tpp7
FFLAG = -c -O4

#SMELIB_DIR = ../..
MINIM_LIB_A = ../lib/minimlib.a

OBJS = BUILDF.o \
       BUILDFD.o \
       MINIM.o \
       BUILDF_Var.o \
       simp1.o \
       simp2.o \
       simp3.o \
       simplx.o

.f.o:
	$(FC) $(FFLAG) $*.f

#	rm -f       libminim.a
#	ar crv      libminim.a $(OBJS)
#	ar tl libminim.a

lib: $(MINIM_LIB_A)
$(MINIM_LIB_A) : $(OBJS)
	@ ar rv $(MINIM_LIB_A) $?
	@ case x$(HASRANLIB) in xt ) echo ranlib; ranlib $(MINIM_LIB_A);; esac
	
#
clean:
	/bin/rm -f $(OBJS)	

wipe:
	/bin/rm -f $(OBJS) $(MINIM_LIB_A) 
