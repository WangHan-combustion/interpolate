
#FLAGS  =  -O5 # -mcmodel=medium -shared-intel #-finline-limit=500 -ipo  #-g -fbacktrace -ffpe-trap=zero,overflow,underflow #-p
#FLAGS  = -O0 -g #-fbacktrace -ffpe-trap=zero,overflow,underflow #-p
FLAGSC  = $(FLAGS) $(DBG) -c
PROF =  

#DECOMP = /home/uttiya/Libraries/2decomp_fft
#DECOMP = /home/simone/LIBRARIES/2decomp_fft
DECOMP = /home00/simone/LIBRARIES/2decomp_fft

INC = -I$(DECOMP)/include

#LIB =  -L/home00/renep/Libraries -lblas -llapack -L$(DECOMP)/lib -l2decomp_fft 
#LIB =  -L$(DECOMP)/lib -l2decomp_fft /home/simone/LIBRARIES/liblapack.so /home/simone/LIBRARIES/libblas.so  -lc -lstdc++
LIB =  -L$(DECOMP)/lib -l2decomp_fft /home00/simone/LIBRARIES/liblapack.so /home00/simone/LIBRARIES/libblas.so -lc -lstdc++
 
#LIB = -lacml -L$(DECOMP)/lib -l2decomp_fft 

# ifort compiler 
COMP  = /opt/openmpi/bin/mpif90 -mcmodel=large -r8 -132
EXE   = g++ 
#COMP  = /opt/mpich/bin/mpif90 -r8 -132
#COMP = mpif90 -mcmodel=large -r8 -132 -g
#COMP = mpif90 -r8 -132 
#COMP = mpiifort -r8 -132 -g

# compiler options for gfortran 
#COMP = mpif90 -fdefault-real-8 -ffixed-line-length-none  
#COMP = gfortran -fdefault-real-8 -ffixed-line-length-none  

PROGRAM = Interpolate 
OBJS = main_comp.o main.o spline.o


all: $(PROGRAM)

$(PROGRAM): $(OBJS) 
	$(COMP) $(FLAGS) $(OBJS) $(LIB)  -o $(PROGRAM)

main_comp.o : main_comp.f90 param.txt 
	$(COMP) $(INC) $(FLAGSC) main_comp.f90
main.o      : main.cpp  
	$(EXE) $(INC) $(FLAGSC) main.cpp
spline.o    : spline.cpp  
	$(EXE) $(INC) $(FLAGSC) spline.cpp

clean:
	$(RM) a.out core* slurm* *.mod $(PROGRAM) $(OBJS) *.o *.mod new* old*


