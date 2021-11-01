TARGET=read_fort_cnt
FFLAGS=-O3

all: $(TARGET)
read_fort_cnt : $(TARGET).f90
	f2py --fcompiler=gfortran -m $(TARGET) -c --f90flags=$(FFLAGS) $<
