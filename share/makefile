OPTIFLAGS= -O3 -fomit-frame-pointer -ffast-math -funroll-all-loops -fstrength-reduce -fthread-jumps -frerun-cse-after-loop -frerun-loop-opt -fgcse -fexpensive-optimizations

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

CXXFLAGS= -c -Wall $(ROOTCFLAGS) -fPIC -I$(HOME)/Avec/

LDFLAGS       = 

SOFLAGS       = -shared

LD            = g++

LIBS          = $(ROOTLIBS)

HDRS          = LAS_alpar.h LASGlobalData.h LAS_globaldata_tools.h LAS_histo.h LAS_RDC_tools.h LAS_Tec_Reconstruction.h LAS_basic_tools.h LAS_control_plots.h LASGlobalDataLoop.h LAS_stabil.h LAS_vectorfloat_tools.h LAS_At_rec.h

SRCS          = LASGlobalDataLoop.cc LAS_alpar.C LAS_basic_tools.C LAS_globaldata_tools.C LAS_stabil.C LAS_control_plots.C LAS_histo.C LAS_RDC_tools.C LAS_Tec_Reconstruction.C LAS_vectorfloat_tools.C

OBJS          = \
                LASDict.o\
                LAS_vectorfloat_tools.o\
                LAS_basic_tools.o\
                LAS_globaldata_tools.o\
                LASGlobalDataLoop.o\
                LAS_alpar.o\
                LAS_stabil.o\
                LAS_control_plots.o\
                LAS_histo.o\
                LAS_Tec_Reconstruction.o\
                LAS_At_rec.o\
                LAS_RDC_tools.o

SHARED_OBJS   = LAS_analysis.so

PROGRAM       = 

all:            $(SHARED_OBJS)

clean:;		rm -f *.o *.d *.so core LASDict.*

# General rule for compiling C files
%.o: ../%.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(CXXFLAGS) -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# General rule for compiling cc files
%.o: ../%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(CXXFLAGS) -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

LASDict.C: LAS_globaldata_tools.h LAS_alpar.h LinkDef.h
	@echo "Generating dictionary ..."
#	rootcint -f LASDict.C -c -I$(HOME)/Avec/ -p LAS_globaldata_tools.h LAS_RDC_tools.h LAS_alpar.h Avec2D.h LinkDef.h
	rootcint -f LASDict.C -c -I$(HOME)/Avec/ -p $(HDRS) Avec.h Avec2D.h LinkDef.h


LAS_analysis.so:	$(OBJS)
		$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@
		@echo "$@ done"

