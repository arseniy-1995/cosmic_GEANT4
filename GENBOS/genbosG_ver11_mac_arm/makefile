
OBJS = block_data.o channelex.o eg.o ev_nuc.o fermi_d.o prim_kine.o yn3.o yn5i.o yn6.o gamsu2.o inG.o

OBJSM = genbosG.o 

#FOR = gfortran

#FOR = g++ -lstdc++ -lgfortran

FOR = gfortran -fno-automatic -O2

#FOR = f2c

LIBS = -lgfortran -L/usr/local/lib/cernlib/2023/lib -lkernlib -lphtools -lmathlib -lpacklib -lm -lgfortran

#LIBS =   -lpacklib -lkernlib -lgeant321 -lphtools -lmathlib -lgfortran /opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/current/libgfortran.a

#all:    genbosG libgenbos.a libgenbos.o libgenbos.so

all:    genbosG libgenbos.a

%.o:	%.F *.inc
	${FOR} -c $<

genbosG: $(OBJS) $(OBJSM)
	${FOR} $(OBJS) $(OBJSM) $(LIBS) -o genbosG

libgenbos.a: $(OBJS)
	ar cr libgenbos.a $(OBJS)

#libgenbos.o: $(OBJS)
#	c libgenbos.o $(OBJS)

#libgenbos.so: $(OBJS)
#	-shared libgenbos.so $(OBJS)

	
tester:	tester.c libgenbos.a
	gcc tester.c libgenbos.a $(LIBS) -o tester


clean:  
	/bin/rm *.o *.so *.a genbosG
