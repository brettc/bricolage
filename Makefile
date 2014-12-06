# This is a very naive Makefile for development purposes on my Mac.
#
# As far as I can figure, there is no IDEAL solution to developing complex
# cython extensions, hence this dumb, platform specific, totally specialised
# solution that works well with my personal development environment for for
# now.

CC=gcc

PYTHONDIR=/Users/Brett/anaconda
LIBINC=-L$(PYTHONDIR)/lib
CYTHON=$(PYTHONDIR)/bin/cython
ORG=./organismal

CFLAGS=-fno-strict-aliasing -arch x86_64 -DNDEBUG -g -fwrapv -O3 \
	   -Wall -Wstrict-prototypes -Wno-unused-function \
	   -stdlib=libc++ -std=c++11 -mmacosx-version-min=10.8

INCLUDES=-I$(PYTHONDIR)/include \
		 -I$(PYTHONDIR)/lib/python2.7/site-packages/numpy/core/include \
		 -I$(PYTHONDIR)/include/python2.7 

LDFLAGS=-bundle -undefined dynamic_lookup -arch x86_64

# OBJECTS= $(SOURCES:%.cpp=$(OBJ_DIR)/%.o)
SOURCES = $(ORG)/pubsub2_ext.cpp \
		  $(ORG)/pubsub2_c.cpp \
		  $(ORG)/scheme_cooperative.cpp

OBJECTS = $(SOURCES:.cpp=.o)

PUBSUB_SO = $(ORG)/pubsub2_ext.so 

all: $(PUBSUB_SO)

$(PUBSUB_SO): $(OBJECTS)
	$(CC) $(LDFLAGS) $(LIBINC) $(OBJECTS) -o $(PUBSUB_SO)

# TODO: Add your includes here
$(ORG)/%.o : $(ORG)/%.cpp 
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(ORG)/pubsub2_ext.cpp : $(ORG)/pubsub2_ext.pyx $(ORG)/*.pxd
	$(CYTHON) --cplus $(ORG)/pubsub2_ext.pyx -o $(ORG)/pubsub2_ext.cpp


# Add the Header Dependencies separately
$(ORG)/pubsub2_ext.o: $(ORG)/pubsub2_c.h
$(ORG)/pubsub2_c.o: $(ORG)/pubsub2_c.h

clean: 
	rm -f $(ORG)/*.o $(ORG)/*.so $(ORG)/pubsub2_ext.cpp

.PHONY: all clean
