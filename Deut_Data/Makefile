#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

#PROGRAM = GetXS_Deut_Plat
PROGRAM = GetXS_Deut
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include 
RTLIBS := -L${ROOTSYS}/lib -lCore -lgfortran # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore /usr/lib64/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore  /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
all: $(PROGRAM)

#GFORTRAN := f77
GFORTRAN := gfortran
#FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

$(PROGRAM): $(PROGRAM).o 
	g++ -o $(PROGRAM) $(PROGRAM).o $(CFLAGS) $(LIBS) $(RTLIBS)

$(PROGRAM).o: $(PROGRAM).C
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C ${INCLUDE}

clean:
	rm -f *.o ./SRC/*.o
