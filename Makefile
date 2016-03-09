# This is a very naive Makefile for development purposes on my Mac.
#
# As far as I have found (yet), there is no IDEAL solution to developing
# complex cross-platform cython extensions that include C++ libraries and
# multiple cython files. Hence this dumb, platform specific, totally
# specialised solution that works well with my personal development environment
# for now. I'm sure it goes against lots of advice about writing Makefiles, as
# it makes everything explicit. But it WORKS, and it is a lot easier to see
# what is going on that way (though it clearly won't scale well).
#
# CURRENT PLATFORM: MacOSX 10.10.1, using Anaconda python 2.7

PYTHONDIR=/Users/Brett/anaconda
PYSRC=./bricolage
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
		-I /usr/local/include \
		-I$(PYTHONDIR)/include \
		-I$(PYTHONDIR)/lib/python2.7/site-packages/numpy/core/include \
		-I$(PYTHONDIR)/include/python2.7

# NOT the same thing
PYEXT_FLAGS=-bundle -undefined dynamic_lookup -arch x86_64
DYLIB_FLAGS=-dynamiclib -undefined dynamic_lookup -arch x86_64

LIBS=-lpython2.7 -lstdc++

CPP_SRCS = $(wildcard $(CPPSRC)/*.cpp)
CPP_OBJS = $(CPP_SRCS:.cpp=.o)
CY_SRCS = $(wildcard $(PYSRC)/*.pyx)
CY_EXTS = $(CY_SRCS:.pyx=.so)
CY_PXDS = $(wildcard $(PYSRC)/*.pxd) $(wildcard $(CPPSRC)/*.pxd)

# We need to give our library a name
CPP_LIBNAME = grn
GRN_DYLIB_NAME = lib$(CPP_LIBNAME).dylib
GRN_DYLIB = $(PYSRC)/$(GRN_DYLIB_NAME)

all: $(CY_EXTS)

# all: $(CPP_OBJS)

# This is useful for invoking via vim
cython: $(CY_SRCS:.pyx=.cpp)


shared: $(GRN_DYLIB)


# Build the shared libary of all c++ code
# NOTE: need to change the "install-name" of the dylib so that it loads
# relative to the binaries that will be using it (the python extensions)
$(GRN_DYLIB): $(CPP_OBJS)
	$(CC) $(DYLIB_FLAGS) $(LIBINC) $(CPP_OBJS) -o $(GRN_DYLIB)
	install_name_tool -id "@loader_path/$(GRN_DYLIB_NAME)" $(GRN_DYLIB)

# otool shit?
# $(CC) $(DYLIB_FLAGS) $(LIBINC) $(CPP_OBJS) -o libgrn.dylib

# ---- Automatic rules
# Make objects for dynamic library

$(CPPSRC)/%.o : $(CPPSRC)/%.cpp 
	$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

# Make python extensions from cython objects
$(PYSRC)/%.so: $(PYSRC)/%.o $(GRN_DYLIB)
	$(CC) $(PYEXT_FLAGS) $(LIBINC) $< -o $@ -l$(CPP_LIBNAME)

# Make objects for all cython output
$(PYSRC)/%.o : $(PYSRC)/%.cpp 
	$(CC) $(CCFLAGS) $(INCLUDES) -c $< -o $@

# Make cpp files from cython
# Overkill -- but better to make sure everything rebuilds on pxd change. 
$(PYSRC)/%.cpp : $(PYSRC)/%.pyx $(CY_PXDS)
	$(CYTHON) --include-dir $(CPPSRC) --cplus $< -o $@

# Add the Header Dependencies generated by -MMD
-include $(wildcard $(CPPSRC)/*.d)
-include $(wildcard $(PYSRC)/*.d)

# Manually add dependencies for your pyx
# $(PYSRC)/threshold3.pyx : 

clean: 
	rm -f $(CY_EXTS) $(GRN_DYLIB) $(MOVED_GRN_DYLIB)
	rm -f $(CPPSRC)/*.d $(PYSRC)/*.d
	rm -f $(CPPSRC)/*.o $(PYSRC)/*.o
	rm -f $(PYSRC)/*.cpp
	rm -f $(PYSRC)/*.pyc
	rm -f $(PYSRC)/*.dylib
	rm -f $(PYSRC)/*.so

cleanall:
	rm -f **/*.o
	rm -f **/*.d
	rm -f **/*.pyc
	rm -f $(PYSRC)/*.cpp

.PHONY: all clean cleanall 
