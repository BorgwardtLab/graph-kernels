# Graph Kernels
A fast C++ implementation of graph kernels including:
* simple kernels between vertex and/or edge label histograms,
* random walk kernels (popular baselines), and
* the Weisfeiler-Lehman graph kernel (state-of-the-art).

Please see the following paper for more details:
* Sugiyama, M., Borgwardt, K. M.: **Halting in Random Walk Kernels**, *Advances in Neural Information Processing Systems (NIPS 2015)*, 2015 [[PDF]](https://papers.nips.cc/paper/5688-halting-in-random-walk-kernels.pdf)

Other implementations of graph kernels and graph benchmark datasets are available [here](https://www.bsse.ethz.ch/mlcb/research/machine-learning/graph-kernels.html).

## Usage
### In your program
You can compute the full kernel matrix by calling the function `graphKernelMatrix`.
To use it, you just need to include the header file "graphKernels.h" in your program.
The [Eigen](http://eigen.tuxfamily.org) and [igraph](http://igraph.org/c/) libraries are needed.

The main function `graphKernelMatrix` is defined as:
```
void graphKernelMatrix(vector<igraph_t>& g, vector<double>& par,
                       string& kernel_name, MatrixXd& K);
```
* `g`: a vector of input graphs
* `par`: a vector of parameters
* `kernel_name`: a string to specify a kernel (see the list below)
* `K`: the full kernel matrix will be returned here

### In terminal
To try the code, we also provide a graph benchmark dataset "mutag" and a test code "graphKernels_test.cc", which includes input and output interface for graph files.

For example:
```
$ cd src/cc
$ make
$ ./gkernel -k kR -p 1,2,1 -i ../../data/mutag.list -g ../../data/mutag/ -o mutag_kR.kernel
> Reading files ... end
> Information:
  Kernel:    k-step random walk kernel
  Parameter: k = 3, lambda = (1, 2, 1)
  Number of graphs: 188
  The average number of vertices: 17.9309
  The maximum number of vertices: 28
  The average number of edges:    19.7926
  The maximum number of edges:    33
> Computing the kernel matrix ... end
  Runtime for the kernel matrix computation: 2.9501 (sec.)
> Writing the kernel matrix to "mutag_kR.kernel" ... end
$ ./gkernel -k WL -p 5 -i ../../data/mutag.list -g ../../data/mutag/ -o mutag_WL.kernel
> Reading files ... end
> Information:
  Kernel:    Weisfeiler-Lehman kernel
  Parameter: h = 5
  Number of graphs: 188
  The average number of vertices: 17.9309
  The maximum number of vertices: 28
  The average number of edges:    19.7926
  The maximum number of edges:    33
> Computing the kernel matrix ... end
  Runtime for the kernel matrix computation: 0.00567007 (sec.)
> Writing the kernel matrix to "mutag_WL.kernel" ... end
```
To compile the program, please edit paths in the "Makefile" according to the location of Eigen and igraph libraries in your environment.

#### Command-line arguments

  `-k <kernel_name>` : An abbreviated kernel name (see the list below)  
  `-p <parameter>` : Parameter(s) in a kernel (if there are more than two, they should be comma-separated)  
  `-i <input_file_list>` : A file describing the list of graph file names  
  `-g <input_file_path>` : A path to the directory of graph files (the GraphML format is supported)  
  `-o <output_file>` : An output file of the full kernel matrix




## List of graph kernels
The following graph kernels are supported:  
The second column (Abbrev.) is used for the third argument of `graphKernelMatrix`.

Kernels                                            | Abbrev. | Parameter
-------------------------------------------------- | ------- | ---------
Linear kernel between vertex histograms            |       V | None
Linear kernel between edge histograms              |       E | None
Linear kernel between vertex-edge histograms       |      VE | None
Linear kernel combination (V + &#955;VE)           |       H | &#955;
Gaussian RBF kernel between vertex histograms      |      VG | &#963;
Gaussian RBF kernel between edge histograms        |      EG | &#963;
Gaussian RBF kernel between vertex-edge histograms |     VEG | &#963;
Geometric random walk kernel                       |      GR | &#955;
Exponential random walk kernel                     |      ER | &#946;
k-step random walk kernel                          |      kR | &#955;0, &#955;1, ..., &#955;k
Weisfeiler-Lehman subtree kernel                   |      WL | h

## Contact
Author: Mahito Sugiyama  
Affiliation: ISIR, Osaka University, Japan  
E-mail: mahito@ar.sanken.osaka-u.ac.jp
