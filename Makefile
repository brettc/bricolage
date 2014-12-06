# This is a very naive Makefile for development purposes on my Mac.
#
# As far as I have found (yet), there is no IDEAL solution to developing
# complex cross platform cython extensions with C++ libraries and multiple
# cython files. Hence this dumb, platform specific, totally specialised
# solution that works well with my personal development environment for now. It
# goes against all advice about writing Makefiles, as it makes everything
# explicit. It is a lot easier to see what is going on that way, though it
# clearly won't scale well.
#
# CURRENT PLATFORM: MacOSX 10.10.1, using Anaconda python 2.7

PYTHONDIR=/Users/Brett/anaconda
PYSRC=./organismal
CPPSRC=./src

LIBINC=-L$(PYTHONDIR)/lib -L$(PYSRC)
CYTHON=$(PYTHONDIR)/bin/cython

CCFLAGS = \
		 -fno-strict-aliasing -arch x86_64 -fPIC -DNDEBUG -g -fwrapv -O3 \
		 -Wall -Wstrict-prototypes -Wno-unused-function \
		 -stdlib=libc++ -std=c++11 -mmacosx-version-min=10.8 \
		 -MMD

INCLUDES = \
		-I. \
		-I$(PYTHONDIR)/include \
		-I$(PYTHONDIR)/lib/python2.7/site-packages/numpy/core/include \
		-I$(PYTHONDIR)/include/python2.7 

# NOT the same thing
PYEXT_FLAGS=-bundle -undefined dynamic_lookup -arch x86_64
DYLIB_FLAGS=-dynamiclib -undefined dynamic_lookup -arch x86_64

LIBS=-lpython2.7 -lstdc++

GRN_SRCS = $(wildcard $(CPPSRC)/*.cpp)
GRN_OBJS = $(GRN_SRCS:.cpp=.o)
CY_SRCS = $(wildcard $(PYSRC)/*.pyx)
CY_EXTS = $(CY_SRCS:.pyx=.so)

GRN_DYLIB = $(PYSRC)/libgrn.dylib

all: $(CY_EXTS)

# Build the shared libary of all c++ code
$(GRN_DYLIB): $(GRN_OBJS)
	$(CC) $(DYLIB_FLAGS) $(LIBINC) $(GRN_OBJS) -o $(GRN_DYLIB)

# ---- Automatic rules
# Make objects for dynamic library
$(CPPSRC)/%.o : $(CPPSRC)/%.cpp 
	$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

# Make python extensions from cython objects
$(PYSRC)/%.so: $(PYSRC)/%.o $(GRN_DYLIB)
	$(CC) $(PYEXT_FLAGS) $(LIBINC) $< -o $@ -lgrn

# Make objects for all cython output
$(PYSRC)/%.o : $(PYSRC)/%.cpp 
	$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

# Make cpp files from cython
$(PYSRC)/%.cpp : $(PYSRC)/%.pyx $(PYSRC)/core_ext.pxd $(PYSRC)/grn.pxd $(PYSRC)/utility.pxd 
	$(CYTHON) --cplus $< -o $@

# Add the Header Dependencies generated by -MMD
-include $(wildcard $(CPPSRC)/*.d)
-include $(wildcard $(PYSRC)/*.d)

clean: 
	rm -f $(CY_EXTS) $(GRN_DYLIB)
	rm -f $(CPPSRC)/*.d $(PYSRC)/*.d
	rm -f $(CPPSRC)/*.o $(PYSRC)/*.o
	rm -f $(PYSRC)/*.cpp
	rm -f $(PYSRC)/*.pyc

cleanall:
	rm -rf **/*.o

.PHONY: all clean cleanlib
