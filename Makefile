include config.mk

lflags=-L. `gsl-config --libs`
iflags=`gsl-config --cflags`

# Lists of files to be built
#objs=
#src=$(patsubst %.o,%.cc,$(objs))
execs=dist_model

all: $(execs)

#include Makefile.dep

#depend:
#	$(cxx) -MM $(src) >Makefile.dep

clean:
	rm $(objs) $(execs)

dist_model: dist_model.cc
	$(cxx) $(cflags) $(iflags) -o $@ $< $(lflags)

%.o: %.cc
	$(cxx) $(cflags) -c $<

.PHONY: clean depend
