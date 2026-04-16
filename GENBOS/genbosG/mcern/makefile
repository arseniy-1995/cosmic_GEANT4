
FOR = gfortran -c -fPIC -fno-automatic -O2  -std=legacy -g -Wno-conversion

MCERN = genbod.o pdk.o rotes2.o norranux.o rnorthux.o

all:    libmcern.so

%.o:	%.F *.h
	${FOR} $<

libmcern.so: $(MCERN)
	gfortran -shared -o libmcern.so $(MCERN)

clean:
	/bin/rm *.so

