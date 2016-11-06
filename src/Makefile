## please change the following paths according to your environment
INCLUDE	= -I/usr/local/include/igraph -I/usr/local/include/eigen3
LIB	= -L/usr/local/lib -ligraph

main: graphKernels_test.cc graphKernels.h
	g++ -O3 $(INCLUDE) graphKernels_test.cc $(LIB) -o gkernel
