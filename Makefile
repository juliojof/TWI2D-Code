CXX      = g++
CXX2     = mpicc 
CCFLAGS  = -lm -lfftw3 -fopenmp -lmpi
#NVCC     = /usr/local/cuda-11.2/bin/nvcc
#NVCC     = /usr/local/cuda-10.2/bin/nvcc
NVCC     = /usr/local/cuda-10.1/bin/nvcc
NVFLAGS  = -gencode -arch=compute_75, code=compute_75

Play: Play.o PropagationFunctions_GPU.o
	$(CXX) -v -o Play Play.o PropagationFunctions_GPU.o -L/usr/local/cuda-10.1/lib64 -lcudart -lcufft $(CCFLAGS) 
#	$(CXX) -v -o Play Play.o PropagationFunctions_GPU.o -L/usr/local/cuda-10.2/lib64 -lcudart -lcufft $(CCFLAGS)
#	$(CXX) -v -o Play Play.o PropagationFunctions_GPU.o -L/usr/local/cuda-11.2/lib64 -lcudart -lcufft $(CCFLAGS)


Play.o: Play.c
	$(CXX2) $(CCFLAGS) -v -c Play.c -o Play.o -I/usr/local/cuda-10.1/include 
#	$(CXX2) $(CCFLAGS) -v -c Play.c -o Play.o -I/usr/local/cuda-10.2/include 
#	$(CXX2) $(CCFLAGS) -v -c Play.c -o Play.o -I/usr/local/cuda-11.2/include 



PropagationFunctions_GPU.o: PropagationFunctions_GPU.cu
	$(NVCC) -v -c PropagationFunctions_GPU.cu -o PropagationFunctions_GPU.o


clean :
	rm *.o