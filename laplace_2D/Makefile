
# For the time being this has no CUDA and just compiles the "gold" code
#INC	:= -I$(CUDA_HOME)/include -I.
#LIB	:= -L$(CUDA_HOME)/lib64 -lcudart

#NVCCFLAGS	:= -lineinfo -arch=sm_35 --ptxas-options=-v --use_fast_math
# NVCCFLAGS	:= -lineinfo -arch=sm_50 --ptxas-options=-v --use_fast_math

all: 		laplace2d 

laplace2d:	laplace2d.cpp CPU_jacobi_laplace2d.cpp CPU_norm2d.cpp Makefile
		g++ laplace2d.cpp CPU_jacobi_laplace2d.cpp CPU_norm2d.cpp -o laplace2d \
#		$(INC) $(NVCCFLAGS) $(LIB)

clean:
		rm -f laplace2d 
