# Graph Kernels
C++ implementation of graph kernels including:
* simple kernels between label histograms,
* random walk kernels, and
* Weisfeiler-Lehman graph kernel (kwon to be the state-of-the-art).

Please see the following papers for more detail:
* Sugiyama, M., Borgwardt, K. M.: **Halting in Random Walk Kernels,**, *Advances in Neural Information Processing Systems (NIPS 2015)*, 2015

## Usage
### From your program
You can compute the full kernel matrix by include the solo header file "graphKernels.h" and run the function "graphKernelMatrix" in your program.
To compile it, the [Eigen](http://eigen.tuxfamily.org) and [igraph](http://igraph.org/c/) libraries are needed.

The main function is "graphKernelMatrix", defined as:
```
void graphKernelMatrix(vector<igraph_t>& g, vector<double>& par, string& kernel_name, MatrixXd& K);
```
* `g`: a vector of input graphs
* `par`: a vector of parameters
* `kernel_name`: a string to specify a kernel (see the list below)
* `K`: the full kernel matrix will be returned here

### From terminal
To try the code, we also provide the test code "graphKernels_test.h", which includes input and output interface for graph files.

For example:
```
$ make
$ ./gkernel -i graphs/mutag.list -g graphs/mutag/ -k kR -p 1,2,1 -o mutag_kR.kernel
>> Reading files ... end
>> Information:
   Kernel:    k-step random walk kernel
   Parameter: k = 2, lambda = (1, 2, 1)
   Number of graphs: 188
   The average number of vertices: 17.9309
   The maximum number of vertices: 28
   The average number of edges:    19.7926
   The maximum number of edges:    33
>> Computing the kernel matrix ... end
   Runtime for the kernel matrix computation: 2.89738 (sec.)
>> Writing the kernel matrix to "mutag_kR.kernel" ... end
$ ./gkernel -i graphs/mutag.list -g graphs/mutag/ -k WL -p 5 -o mutag_WL.kernel
>> Reading files ... end
>> Information:
   Kernel:    Weisfeiler-Lehman kernel
   Parameter: h = 5
   Number of graphs: 188
   The average number of vertices: 17.9309
   The maximum number of vertices: 28
   The average number of edges:    19.7926
   The maximum number of edges:    33
>> Computing the kernel matrix ... end
   Runtime for the kernel matrix computation: 0.008043 (sec.)
>> Writing the kernel matrix to "mutag_WL.kernel" ... end
```
In compilation, please edit paths in the "Makefile" according to the location of Eigen and igraph libraries in your environment.

#### Argument list

  `-i <input_file_list>` : A file describing graph file names  
  `-i <input_file_path>` : A path to the directory of graphfiles  
  `-k <kernel_name>` : The abbreviated kernel name (see the list below)  
  `-p <parameter>` : Parameter(s) (if there are more than two, they should be comma-separated)  
  `-o <output_file>` : Output file of the full kernel matrix




## List of graph kernels
The following graph kernels are implemented:  
The second column (Abbrev.) is used for the third argument of `graphKernelMatrix`.

Kernels                                            | Abbrev. | Parameter
-------------------------------------------------- | ------- | ---------
Linear kernel between vertex histograms            |       V | None
Linear kernel between edge histograms              |       E | None
Linear kernel between vergex-edge histograms       |      VE | None
Linear kernel combination (V + &#955;VE)           |       H | &#955;
Gaussian RBF kernel between vertex histograms      |      VG | &#963;
Gaussian RBF kernel between edge histograms        |      EG | &#963;
Gaussian RBF kernel between vergex-edge histograms |     VEG | &#963;
Geometric random walk kernel                       |      GR | &#955;
Exponential random walk kernel                     |      ER | &#946;
k-step random walk kernel                          |      kR | &#955;0, &#955;1, ..., &#955;k
Weisfeiler-Lehman subtree kernel                   |      WL | h

## Contact
Author: Mahito Sugiyama  
Affiliation: ISIR, Osaka University, Japan  
E-mail: mahito@ar.sanken.osaka-u.ac.jp
