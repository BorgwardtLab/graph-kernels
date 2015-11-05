// ./gkernel -i /Users/mahito/sync/data/graph/graphml/mutag.list -g /Users/mahito/sync/data/graph/graphml/mutag/ -k WL -p 5 -o output
// ./gkernel -i /home/mahito/sync/data/graph/graphml/mutag.list -g /home/mahito/sync/data/graph/graphml/mutag/ -k WL -p 5 -o output
#include "graphKernels.h"
#include <unistd.h>

#define BUF 256

using namespace std;

// =========================================================== //
// ==================== Utility functions ==================== //
// =========================================================== //
// Measure running time
double measureTime() {
  struct rusage t;
  struct timeval tv;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}
// compute simple statistics of graphs
void computeStatistics(vector<igraph_t>& g) {
  int n = (int)g.size();
  // initialization
  vector<MatrixXi> E(n);
  vector<vector<int> > V_label(n);
  vector<int> V_count(n);
  vector<int> E_count(n);
  vector<int> D_max(n);

  // from igraph to Eigen
  for (int i = 0; i < n; i++) {
    igraphToEigen(g[i], E[i], V_label[i], &V_count[i], &E_count[i], &D_max[i]);
  }

  // compute simple statistics
  double vmean = 0, vmax = 0, emean = 0, emax = 0;
  for (int i = 0; i < n; i++) {
    vmean += (double)V_count[i];
    emean += (double)E_count[i];
    if (vmax < (double)V_count[i]) {
      vmax = (double)V_count[i];
    }
    if (emax < (double)E_count[i]) {
      emax = (double)E_count[i];
    }
  }
  vmean = vmean / (double)n;
  emean = emean / (double)n;

  // cout << ">> Statistics:" << endl;
  cout << "   Number of graphs: " << n << endl;
  cout << "   The average number of vertices: " << vmean << endl;
  cout << "   The maximum number of vertices: " << vmax << endl;
  cout << "   The average number of edges:    " << emean << endl;
  cout << "   The maximum number of edges:    " << emax << endl;
}
// get kernel name
string getKernelName(const char* kernel_name) {
  if (strcmp(kernel_name, "E")   == 0) return "edge label histogram";
  else if (strcmp(kernel_name, "V")   == 0) return "vertex label histogram";
  else if (strcmp(kernel_name, "VE")  == 0) return "vertex-edge label histogram";
  else if (strcmp(kernel_name, "H")   == 0) return "vertex-vertex-edge label histogram";
  else if (strcmp(kernel_name, "EG")  == 0) return "vertex label histogram (Gaussian)";
  else if (strcmp(kernel_name, "VG")  == 0) return "vertex-edge label histogram (Gaussian)";
  else if (strcmp(kernel_name, "VEG") == 0) return "vertex-vertex-edge label histogram (Gaussian)";
  else if (strcmp(kernel_name, "GR")  == 0) return "geometric random walk";
  else if (strcmp(kernel_name, "ER")  == 0) return "exponential random walk";
  else if (strcmp(kernel_name, "kR")  == 0) return "k-step random walk";
  else if (strcmp(kernel_name, "WL")  == 0) return "Weisfeiler-Lehman";
  else return "edge label histogram";
}
// print kernel name and parameters
void printKernel(const char* kernel_name, vector<double>& par) {
  cout << ">> Information:" << endl;
  cout << "   Kernel:    " << getKernelName(kernel_name) << " kernel" << endl;
  cout << "   Parameter: ";
  if (kernelTable(kernel_name) == 1 || kernelTable(kernel_name) == 2 || kernelTable(kernel_name) == 3) {
    cout << "NONE" << endl;
  } else if (kernelTable(kernel_name) == 4 || kernelTable(kernel_name) == 8) {
    cout << "lambda = " << par[0] << endl;
  } else if (kernelTable(kernel_name) == 5 || kernelTable(kernel_name) == 6 || kernelTable(kernel_name) == 7) {
    cout << "sigma = " << par[0] << endl;
  } else if (kernelTable(kernel_name) == 9) {
    cout << "beta = " << par[0] << endl;
  } else if (kernelTable(kernel_name) == 10) {
    cout << "k = " << par.size() << ", lambda = (";
    for (int i = 0; i < (int)par.size() - 1; i++) {
      cout << par[i] << ", ";
    }
    cout << par[par.size() - 1] << ")" << endl;
  } else if (kernelTable(kernel_name) == 11) {
    cout << "h = " << (int)par[0] << endl;
  } else {
    cout << "NONE" << endl;
  }
}

// ========================================================== //
// ==================== Read graph files ==================== //
// ========================================================== //
// Currently GraphML format is used
// Please change the function "igraph_read_graph_graphml" if you want to use other format
// The list of formats supported by iGraph:
// http://igraph.org/c/doc/igraph-Foreign.html
igraph_t readiGraph(const char *input) {
  igraph_t g;
  FILE *fp;
  int res;

  // GraphML
  fp = fopen(input, "r");
  if (fp == 0) {
    cerr << "ERROR: cannot open " << input << endl;
    exit(1);
  }
  res = igraph_read_graph_graphml(&g, fp, 0);
  if (res == IGRAPH_PARSEERROR) {
    cerr << "ERROR: There is a problem reading the file, or the file is syntactically incorrect." << endl;
    exit(1);
  } else if (res == IGRAPH_UNIMPLEMENTED) {
    cerr << "ERROR: The GraphML functionality was disabled at compile-time." << endl;
    exit(1);
  }

  fclose(fp);
  return g;
}


// ======================================================= //
// ==================== Main function ==================== //
// ======================================================= //
int main(int argc, char *argv[]) {
  // use attributes
  igraph_i_set_attribute_table(&igraph_cattribute_table);

  vector<double> par;
  vector<igraph_t> g;
  char *val, *cpar = NULL;
  const char *glist, *gpath, *kernel_name, *output;
  const char *kernel_default = "edge";
  kernel_name = kernel_default;

  // get arguments
  char opt;
  while ( ( opt = getopt ( argc, argv, "i:g:p:k:o:" ) ) != -1 ) {
    switch ( opt ) {
    case 'i': glist = optarg; break;
    case 'g': gpath = optarg; break;
    case 'p': cpar = optarg; break;
    case 'k': kernel_name = optarg; break;
    case 'o': output = optarg; break;
    }
  }

  // read a graph file
  cout << ">> Reading files ... ";
  ifstream ifs(glist);
  char str[BUF], pstr[BUF];
  if (ifs.fail()) {
    cerr << "ERROR: cannot open the list file " << glist << endl;
    exit(1);
  }
  while (ifs.getline(str, BUF - 1)) {
    sprintf(pstr, "%s%s", gpath, str);
    g.push_back(readiGraph(pstr));
  }
  cout << "end" << endl;

  // get parameters
  if (cpar == NULL) {
    par.push_back(1.0);
  } else {
    val = strtok(cpar, ",");
    par.push_back(atof(val));
    while (val != NULL) {
      val = strtok(NULL, ",");
      if (val != NULL ) par.push_back(atof(val));
    }
  }

  // print kernel information
  printKernel(kernel_name, par);
  // check statistics
  computeStatistics(g);

  // compute the kernel value
  MatrixXd K(g.size(), g.size());
  cout << ">> Computing the kernel matrix ... " << flush;
  double t1 = measureTime();
  graphKernelMatrix(g, par, kernel_name, K);
  double t2 = measureTime();
  cout << "end" << endl;
  cout << "   Runtime for the kernel matrix computation: " << t2 - t1 << " (sec.)" << endl;

  // write the kernel to a file
  cout << ">> Writing the kernel matrix to \"" << output << "\" ... ";
  ofstream ofs(output);
  for (int i = 0; i < (int)g.size(); i++) {
    for (int j = 0; j < (int)g.size(); j++) {
      ofs << K(i, j) << ",";
    }
    ofs << endl;
  }
  ofs.close();
  cout << "end" << endl;

  // free
  for (int i = 0; i < (int)g.size(); i++) {
    igraph_destroy(&g[i]);
  }
  return 0;
}
