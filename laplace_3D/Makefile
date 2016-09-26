
# For the time being this has no CUDA and just compiles the "gold" code
#INC	:= -I$(CUDA_HOME)/include -I.
#LIB	:= -L$(CUDA_HOME)/lib64 -lcudart

#NVCCFLAGS	:= -lineinfo -arch=sm_35 --ptxas-options=-v --use_fast_math
# NVCCFLAGS	:= -lineinfo -arch=sm_50 --ptxas-options=-v --use_fast_math

all: 		laplace3d 

laplace3d:	laplace3d.cpp laplace3d_gold.cpp Makefile
		g++ laplace3d.cpp laplace3d_gold.cpp -o laplace3d \
#		$(INC) $(NVCCFLAGS) $(LIB)

clean:
		rm -f laplace3d 
