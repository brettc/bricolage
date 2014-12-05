CC=gcc
PYTHONDIR=/Users/Brett/anaconda
CFLAGS=-fno-strict-aliasing -arch x86_64 -DNDEBUG -g -fwrapv -O3 -Wall -Wstrict-prototypes -Wno-unused-function -stdlib=libc++ -std=c++11 -mmacosx-version-min=10.8
INCLUDES=-I$(PYTHONDIR)/include -I$(PYTHONDIR)/lib/python2.7/site-packages/numpy/core/include -I$(PYTHONDIR)/include/python2.7 
LDFLAGS=-bundle -undefined dynamic_lookup -arch x86_64
LIBINC=-L$(PYTHONDIR)/lib
CYTHON=$(PYTHONDIR)/bin/cython

SOURCES = pubsub2_ext.cpp pubsub2_c.cpp scheme_cooperative.cpp
OBJECTS = $(SOURCES:.cpp=.o)

all: pubsub2_ext.so

# pubsub2_ext.so: pubsub2_ext.o pubsub2_c.o scheme_cooperative.o
pubsub2_ext.so: $(OBJECTS)
	$(CC) $(LDFLAGS) $(LIBINC) $(OBJECTS) -o pubsub2_ext.so

%.o : %.cpp pubsub2_c.h scheme_cooperative.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

pubsub2_ext.cpp : pubsub2_ext.pyx *.pxd
	$(CYTHON) --cplus pubsub2_ext.pyx -o pubsub2_ext.cpp

clean: 
	rm -f *.o *.so pubsub2_ext.cpp

.PHONY: all clean
