
OBJS = block_data.o channelex.o eg.o ev_nuc.o fermi_d.o prim_kine.o yn3.o yn5i.o yn6.o gamsu2.o inG.o

OBJSM = genbosG.o 

FOR = gfortran -fno-automatic -O2

#LIBS = -lkernlib -lphtools -lmathlib -lgfortran -lm
#LIBS = -lkernlib -lnsl -lcrypt -ldl  -lphtools -lmathlib -lgeant321 -lpacklib -lgfortran -lm -lnsl -lphtools /usr/local/lib/cernlib/2023/lib/libpacklib.a -lnsl -lcrypt -ldl

#LIBS = -lgfortran -L/usr/local/lib/cernlib/2023/lib -lkernlib -lphtools -lmathlib -lpacklib -lm -lgfortran
#LIBS = -lgfortran -lkernlib -lnsl -lcrypt -ldl  -lphtools -lmathlib -lgeant321 -lpacklib -lgfortran -lm /usr/local/lib/cernlib/2023/lib/libpacklib.a -lnsl -lphtools
#LIBS = -L/usr/local/lib/cernlib/2023/lib -Wl,-static -Wl,-dy -lnsl -lcrypt -ldl

LIBS = -lmathlib -lgfortran -L/usr/local/lib/cernlib/2023/lib -lpacklib -lphtools -lkernlib -lm

#LIBS += -lpacklib -lkernlib -lgeant321 -lphtools -lmathlib

#LIBS += /usr/local/lib/cernlib/2023/lib/libpacklib.a -lnsl -lcrypt -ldl

all:    genbosG libgenbos.a

%.o:	%.F *.inc
	${FOR} -c $<

genbosG: $(OBJS) $(OBJSM)
	${FOR} $(OBJS) $(OBJSM) $(LIBS) -o genbosG

libgenbos.a: $(OBJS)
	ar cr libgenbos.a $(OBJS)
	
tester:	tester.c libgenbos.a
	gcc tester.c libgenbos.a $(LIBS) -o tester


clean:  
	/bin/rm *.o *.so *.a genbosG
