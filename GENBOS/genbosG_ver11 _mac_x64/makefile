
OBJS = block_data.o channelex.o eg.o ev_nuc.o fermi_d.o prim_kine.o yn3.o yn5i.o yn6.o gamsu2.o inG.o

OBJSM = genbosG.o 

FOR = gfortran -fno-automatic -O2

LIBS = -lkernlib -lphtools -lmathlib -lgfortran -lm

#LIBS =   -lpacklib -lkernlib -lgeant321 -lphtools -lmathlib

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
	/bin/rm *.o genbosG
