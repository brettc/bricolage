# This is a very naive Makefile for development purposes on my Mac.
#
# As far as I have found (yet), there is no IDEAL solution to developing
# complex cython extensions, hence this dumb, platform specific, totally
# specialised solution that works well with my personal development environment
# for now. It goes against all advice about writing Makefiles, as it makes
# everything explicit.
#
# CURRENT PLATFORM: MacOSX 10.10.1, using Anaconda python 2.7

PYTHONDIR=/Users/Brett/anaconda
PYSRC=./organismal
CPPSRC=./src

LIBINC=-L$(PYTHONDIR)/lib -L$(PYSRC)
CYTHON=$(PYTHONDIR)/bin/cython

CFLAGS = \
		 -fno-strict-aliasing -arch x86_64 -fPIC -DNDEBUG -g -fwrapv -O3 \
		 -Wall -Wstrict-prototypes -Wno-unused-function \
		 -stdlib=libc++ -std=c++11 -mmacosx-version-min=10.8

INCLUDES = \
		-I$(CPPSRC) \
		-I$(PYTHONDIR)/include \
		-I$(PYTHONDIR)/lib/python2.7/site-packages/numpy/core/include \
		-I$(PYTHONDIR)/include/python2.7 

PYEXT_FLAGS=-bundle -undefined dynamic_lookup -arch x86_64
DYLIB_FLAGS=-dynamiclib -undefined dynamic_lookup -arch x86_64

LIBS=-lpython2.7 -lstdc++

PYEXT_SOURCES = \
		  $(PYSRC)/pubsub2_ext.cpp

GRN_SOURCES = \
		  $(CPPSRC)/pubsub2_c.cpp \
		  $(CPPSRC)/scheme_cooperative.cpp

PYEXT_OBJS = $(PYEXT_SOURCES:.cpp=.o)
GRN_OBJECTS = $(GRN_SOURCES:.cpp=.o)

PUBSUB_EXT = $(PYSRC)/pubsub2_ext.so 
GRN_DYLIB = $(PYSRC)/libgrn.dylib

all: $(PUBSUB_EXT) $(GRN_DYLIB)

$(GRN_DYLIB): $(GRN_OBJECTS)
	$(CC) $(DYLIB_FLAGS) $(LIBINC) $(GRN_OBJECTS) -o $(GRN_DYLIB)

$(PUBSUB_EXT): $(PYEXT_OBJS) $(GRN_DYLIB)
	$(CC) $(PYEXT_FLAGS) $(LIBINC) $(PYEXT_OBJS) -o $(PUBSUB_EXT) -lgrn

# TODO: Add your includes here
$(PYSRC)/%.o : $(PYSRC)/%.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(CPPSRC)/%.o : $(CPPSRC)/%.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(PYSRC)/pubsub2_ext.cpp : $(PYSRC)/pubsub2_ext.pyx $(PYSRC)/*.pxd
	$(CYTHON) --cplus $(PYSRC)/pubsub2_ext.pyx -o $(PYSRC)/pubsub2_ext.cpp

# Add the Header Dependencies separately
$(PYSRC)/pubsub2_ext.o: $(CPPSRC)/pubsub2_c.h
$(PYSRC)/pubsub2_c.o: $(CPPSRC)/pubsub2_c.h

clean: 
	rm -f $(PYSRC)/*.o $(PYSRC)/*.so $(PYSRC)/*.dylib 
	rm -f $(PYEXT_SOURCES)
	rm -f $(CPPSRC)/*.o 


.PHONY: all clean cleanlib
