
OBJS = block_data.o channelex.o eg.o ev_nuc.o fermi_d.o prim_kine.o yn3.o yn5i.o yn6.o gamsu2.o inG.o

OBJSM = genbosG.o 

FOR = gfortran -c -fPIC -fno-automatic -O2  -std=legacy -g -Wno-conversion
LIBS = -lgfortran -lgenbos -lmcern -lkernlib -lm

CDIR = $(HOME)/cern

LDIR = -L$(CDIR)/lib -L.

MCERN = mcern/genbod.o mcern/pdk.o mcern/rotes2.o mcern/norranux.o mcern/rnorthux.o

all:    tester genbosG libgenbos.so libmcern.so

%.o:	%.F *.inc
	${FOR} $<

genbosG: $(OBJS) $(OBJSM) libgenbos.so libmcern.so
	gfortran $(OBJS) $(OBJSM) $(LDIR) $(LIBS) -o genbosG

libgenbos.so: $(OBJS)
	gfortran -shared -o libgenbos.so $(OBJS)

libmcern.so: $(MCERN)
	make -C mcern
	gfortran -shared -o libmcern.so $(MCERN)

tester:	tester.c libgenbos.so libmcern.so
	gcc tester.c -g $(LDIR) $(LIBS) -o tester

install: 
	cp -a libgenbos.so libmcern.so $(CDIR)/lib
	cp -a genbos.h genbos.hh $(CDIR)/include

clean:  
	/bin/rm *.o *.so genbosG tester mcern/*.o mcern/*.so
