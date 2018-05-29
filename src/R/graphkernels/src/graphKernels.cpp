#include <RcppEigen.h>
#define Int int32_t

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::SparseMatrix;
using Eigen::LLT;

typedef Eigen::Triplet<double> T;

using namespace Rcpp;
using namespace Eigen;
using namespace std;

template <typename S>
ostream &operator<<(ostream& out, const vector<S>& vec) {
  for (Int i = 0; i < vec.size() - 1; ++i) { cout << vec[i] << " "; }
  cout << vec[vec.size() - 1];
  return out;
}

// =================================================================== //
// ==================== Functions used in kernels ==================== //
// =================================================================== //
// select linear kernel or Gaussian kernel in histogram kernels
double selectLinearGaussian(vector<Int>& h1, vector<Int>& h2, double sigma) {
  double K = 0;
  if (sigma < 0) {
    // linear kernel
    for (Int i = 0; i < (Int)h1.size(); i++) {
      K += (double)h1[i] * (double)h2[i];
    }
  } else {
    // Gaussian kernel
    for (Int i = 0; i < (Int)h1.size(); i++) {
      K += ((double)h1[i] - (double)h2[i]) * ((double)h1[i] - (double)h2[i]);
    }
    K = exp(-1.0 * K / (2.0 * sigma * sigma));
  }
  return K;
}
// map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
Int productMapping(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, MatrixXi& H) {
  Int n_vx = 0;
  for (Int i = 0; i < (Int)v1_label.size(); i++) {
    for (Int j = 0; j < (Int)v2_label.size(); j++) {
      if (v1_label[i] == v2_label[j]) {
	H(i, j) = n_vx;
	n_vx++;
      }
    }
  }
  return n_vx;
}
// compute the adjacency matrix Ax of the direct product graph (sparse)
void productAdjacency(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, MatrixXi& H, SparseMatrix<double>& Ax) {
  vector<T> v;
  for (Int i = 0; i < e1.rows(); i++) {
    for (Int j = 0; j < e2.rows(); j++) {      
      if (   v1_label[e1(i, 0)] == v2_label[e2(j, 0)]
	  && v1_label[e1(i, 1)] == v2_label[e2(j, 1)]
	  && e1(i, 2) == e2(j, 2)) {
	v.push_back(T(H(e1(i, 0), e2(j, 0)), H(e1(i, 1), e2(j, 1)), 1.0));
	v.push_back(T(H(e1(i, 1), e2(j, 1)), H(e1(i, 0), e2(j, 0)), 1.0));
      }
      if (   v1_label[e1(i, 0)] == v2_label[e2(j, 1)]
	  && v1_label[e1(i, 1)] == v2_label[e2(j, 0)]
	  && e1(i, 2) == e2(j, 2)) {
	v.push_back(T(H(e1(i, 0), e2(j, 1)), H(e1(i, 1), e2(j, 0)), 1.0));
	v.push_back(T(H(e1(i, 1), e2(j, 0)), H(e1(i, 0), e2(j, 1)), 1.0));
      }
    }
  }
  Ax.setFromTriplets(v.begin(), v.end());
}
// bucket sort used in Weisfeiler-Leiman graph kernel
void bucketsort(vector<Int>& x, vector<Int>& index, Int label_max) {
  vector<vector<Int>> buckets;
  buckets.resize(label_max + 1);

  for (vector<Int>::iterator itr = index.begin(), end = index.end(); itr != end; ++itr) {
    buckets[ x[*itr] ].push_back(*itr);
  }

  Int counter = 0;
  for (vector<vector<Int>>::iterator itr = buckets.begin(), end = buckets.end(); itr != end; ++itr) {
    for (vector<Int>::iterator itr2 = (*itr).begin(), end2 = (*itr).end(); itr2 != end2; ++itr2) {
      index[counter] = *itr2;
      counter++;
    }
  }
}

// ===================================================== //
// ==================== Each kernel ==================== //
// ===================================================== //
// edge histogram karnel
double edgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, double sigma) {
  Int e_label_max = 0;
  for (Int i = 0; i < e1.rows(); i++) {
    if (e1(i, 2) > e_label_max) e_label_max = e1(i, 2);
  }
  for (Int i = 0; i < e2.rows(); i++) {
    if (e2(i, 2) > e_label_max) e_label_max = e2(i, 2);
  }

  vector<Int> h1(e_label_max + 1, 0);
  vector<Int> h2(e_label_max + 1, 0);

  for (Int i = 0; i < e1.rows(); i++) {
    (h1[e1(i, 2)])++;
  }
  for (Int i = 0; i < e2.rows(); i++) {
    (h2[e2(i, 2)])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}
// vertex histogram karnel
double vertexHistogramKernel(vector<Int>& v1_label, vector<Int>& v2_label, double sigma) {
  Int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  Int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  Int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;

  vector<Int> h1(v_label_max + 1, 0);
  vector<Int> h2(v_label_max + 1, 0);

  for (Int i = 0; i < (Int)v1_label.size(); i++) {
    (h1[v1_label[i]])++;
  }
  for (Int i = 0; i < (Int)v2_label.size(); i++) {
    (h2[v2_label[i]])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}
// vertex-edge histogram karnel
double vertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, double sigma) {
  Int e_label_max = 0;
  for (Int i = 0; i < e1.rows(); i++) {
    if (e1(i, 2) > e_label_max) e_label_max = e1(i, 2);
  }
  for (Int i = 0; i < e2.rows(); i++) {
    if (e2(i, 2) > e_label_max) e_label_max = e2(i, 2);
  }
  e_label_max++;

  Int v1_label_max = *max_element(v1_label.begin(), v1_label.end());
  Int v2_label_max = *max_element(v2_label.begin(), v2_label.end());
  Int v_label_max = v1_label_max > v2_label_max ? v1_label_max : v2_label_max;
  v_label_max++;

  vector<Int> h1(v_label_max * v_label_max * e_label_max, 0);
  vector<Int> h2(v_label_max * v_label_max * e_label_max, 0);

  Int v1, v2;
  for (Int i = 0; i < e1.rows(); i++) {
    v1 = e1(i, 0);
    v2 = e1(i, 1);
    if (v2 > v1) { Int v_tmp = v1; v1 = v2; v2 = v_tmp; }
    (h1[v1_label[v1] + v1_label[v2] * v_label_max + e1(i, 2) * v_label_max * v_label_max])++;
  }
  for (Int i = 0; i < e2.rows(); i++) {
    v1 = e2(i, 0);
    v2 = e2(i, 1);
    if (v2 > v1) { Int v_tmp = v1; v1 = v2; v2 = v_tmp; }
    (h2[v2_label[v1] + v2_label[v2] * v_label_max + e2(i, 2) * v_label_max * v_label_max])++;
  }

  return selectLinearGaussian(h1, h2, sigma);
}
// vertex-vertex-edge histogram karnel
double vertexVertexEdgeHistogramKernel(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, double lambda) {
  return vertexHistogramKernel(v1_label, v2_label, -1.0) + lambda * vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, -1.0);
}
// geometric random walk karnel
double geometricRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, double lambda) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());  
  Int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  productAdjacency(e1, e2, v1_label, v2_label, H, Ax);

  // inverse of I - lambda * Ax by fixed-poInt iterations
  VectorXd I_vec(n_vx);
  for (Int i  = 0; i < n_vx; i++) I_vec[i] = 1;
  VectorXd x = I_vec;
  VectorXd x_pre(n_vx); x_pre.setZero();

  double eps = pow(10, -10);
  Int count = 0;
  while ((x - x_pre).squaredNorm() > eps) {
    if (count > 100) {
      // cout << "does not converge until " << count - 1 << " iterations" << endl;
      break;
    }
    x_pre = x;
    x = I_vec + lambda * Ax * x_pre;
    count++;
  }
  return x.sum();
}
// exponential random walk karnel
double exponentialRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, double beta) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());  
  Int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  productAdjacency(e1, e2, v1_label, v2_label, H, Ax);

  // compute e^{beta * Ax}
  SelfAdjointEigenSolver<MatrixXd> es(Ax);
  VectorXd x = (beta * es.eigenvalues()).array().exp();
  MatrixXd D = x.asDiagonal();  
  MatrixXd V = es.eigenvectors();

  MatrixXd I(n_vx, n_vx);
  I.setIdentity();
  FullPivLU<MatrixXd> solver(V);
  MatrixXd V_inv = solver.solve(I);
  MatrixXd Res = V * D * V_inv;

  // compute the total sum
  double K = 0;
  for (Int i = 0; i < Res.rows(); i++) {
    for (Int j = 0; j < Res.cols(); j++) {
      K += Res(i, j);
    }
  }

  return K;
}
// k-step product graph karnel
double kstepRandomWalkKernel(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, vector<double>& lambda_list) {
  // map each product (v_1, v_2) of vertics to a number H(v_1, v_2)
  MatrixXi H(v1_label.size(), v2_label.size());
  Int n_vx = productMapping(e1, e2, v1_label, v2_label, H);

  // prepare identity matrix
  SparseMatrix<double> I(n_vx, n_vx);
  I.setIdentity();

  // compute the adjacency matrix Ax of the direct product graph
  SparseMatrix<double> Ax(n_vx, n_vx);
  productAdjacency(e1, e2, v1_label, v2_label, H, Ax);

  // compute products until k
  Int k_max = (Int)lambda_list.size() - 1;
  SparseMatrix<double> Ax_pow = I;
  SparseMatrix<double> Sum = lambda_list[0] * I;
  for (Int k = 1; k <= k_max; k++) {
    Ax_pow = Ax * Ax_pow;
    Sum += lambda_list[k] * Ax_pow;
  }

  // compute the total sum
  double K = 0;
  for (Int i = 0; i < Sum.outerSize(); ++i) {
    for (SparseMatrix<double>::InnerIterator it(Sum, i); it; ++it) {
      K += it.value();
    }
  }

  return K;
}
// Weisfeiler-Leiman graph kernel
void WLKernelMatrix(vector<MatrixXi>& E, vector<vector<Int>>& V_label, vector<Int>& num_v, vector<Int>& num_e, vector<Int>& degree_max, Int h_max, NumericMatrix& K_mat) {
  // K_mat.setZero();
  Int n = (Int)E.size();
  Int v_all = accumulate(num_v.begin(), num_v.end(), 0);
  Int degree_max_all = *max_element(degree_max.begin(), degree_max.end());
  vector<Int> label_max_vec(n);
  for (Int i = 0; i < n; i++) {
    label_max_vec[i] = *max_element(V_label[i].begin(), V_label[i].end());
  }
  Int label_max = *max_element(label_max_vec.begin(), label_max_vec.end());

  Int raise = 0;
  vector<Int> counter(*max_element(num_v.begin(), num_v.end()));
  MatrixXi nei_list(v_all, degree_max_all + 1);
  MatrixXi label_list(v_all, h_max + 1);
  vector<Int> x(v_all);
  vector<Int> index(v_all);
  vector<Int> index_org(v_all);
  vector<Int> graph_index(v_all);

  label_list.setZero();

  for (Int i = 0; i < n; i++) {
    for (Int j = 0; j < num_v[i]; j++) {
      label_list(j + raise, 0) = V_label[i][j];
      graph_index[j + raise] = i;
    }
    raise += num_v[i];
  }

  // ===== Increment kernel values using the initial vertex labels =====
  // radix sort
  for (Int i = 0; i < v_all; i++) {
    index[i] = i;
    index_org[i] = i;
    x[i] = label_list(i, 0);
  }
  bucketsort(x, index, label_max);
  // add kernel values
  vector<Int> count(n);
  set<Int> count_index;
  Int k_value;
  for (Int i = 0; i < v_all; i++) {
    count_index.insert(graph_index[index_org[index[i]]]);
    count[graph_index[index_org[index[i]]]]++;
    if (i == v_all - 1 || label_list(index[i], 0) != label_list(index[i + 1], 0)) {
      for (set<Int>::iterator itr = count_index.begin(), end = count_index.end(); itr != end; ++itr) {
	for (set<Int>::iterator itr2 = itr, end2 = count_index.end(); itr2 != end2; ++itr2) {
	  k_value = count[*itr] * count[*itr2];
	  K_mat(*itr, *itr2) += k_value;
	  K_mat(*itr2, *itr) += k_value;
	}
	count[*itr] = 0;
      }
      count_index.clear();
    }
  }

  Int v_raised_1, v_raised_2;
  for (Int h = 0; h < h_max; h++) {
    nei_list.setZero();

    // first put vertex label
    nei_list.col(0) = label_list.col(h);
    // second put neibor labels
    raise = 0;
    for (Int i = 0; i < n; i++) {
      fill(counter.begin(), counter.end(), 1);
      for (Int j = 0; j < num_e[i]; j++) {
	v_raised_1 = E[i](j, 0) + raise;
	v_raised_2 = E[i](j, 1) + raise;
	nei_list(v_raised_1, counter[E[i](j, 0)]) = label_list(v_raised_2, h);
	nei_list(v_raised_2, counter[E[i](j, 1)]) = label_list(v_raised_1, h);
	counter[E[i](j, 0)]++;
	counter[E[i](j, 1)]++;
      }
      raise += num_v[i];
    }

    // sort each row w.r.t. neighbors
    vector<int> y(nei_list.cols() - 1);
    for (int i = 0; i < v_all; i++) {
      for (int j = 1; j < nei_list.cols(); ++j) {
	y[j - 1] = nei_list(i, j);
      }
      sort(y.begin(), y.end(), greater<int>());
      for (int j = 1; j < nei_list.cols(); ++j) {
	nei_list(i, j) = y[j - 1];
      }
    }
    // radix sort
    for (Int i = 0; i < v_all; i++) {
      index[i] = i;
      index_org[i] = i;
    }
    for (Int k = nei_list.cols() - 1; k >= 0; k--) {
      for (Int i = 0; i < v_all; i++) {
	x[i] = nei_list(i, k);
      }
      bucketsort(x, index, label_max);
    }
    // re-labeling and increment kernel values
    label_max++;
    for (Int i = 0; i < v_all; i++) {
      label_list(index_org[index[i]], h + 1) = label_max;
      count_index.insert(graph_index[index_org[index[i]]]);
      count[graph_index[index_org[index[i]]]]++;
      if (i == v_all - 1 ||
	  (nei_list.row(index[i]) - nei_list.row(index[i + 1])).array().abs().sum() != 0) {
	for (set<Int>::iterator itr = count_index.begin(), end = count_index.end(); itr != end; ++itr) {
	  for (set<Int>::iterator itr2 = itr, end2 = count_index.end(); itr2 != end2; ++itr2) {
	    k_value = count[*itr] * count[*itr2];
	    K_mat(*itr, *itr2) += k_value;
	    K_mat(*itr2, *itr) += k_value;
	  }
	  count[*itr] = 0;
	}
	count_index.clear();	
	label_max++;
      }
    }
  }
  for (Int i = 0; i < n; ++i) {
    K_mat(i, i) /= 2;
  }
}

// compute a kernel value of a pair of graphs
double computeKernelValue(MatrixXi& e1, MatrixXi& e2, vector<Int>& v1_label, vector<Int>& v2_label, vector<double>& par, Int kernel_type) {
  double Kval;
  switch (kernel_type) {
  case 1: // edge histogram kernel
    Kval = edgeHistogramKernel(e1, e2, -1.0);
    break;
  case 2: // vertex histogram kernel
    Kval = vertexHistogramKernel(v1_label, v2_label, -1.0);
    break;
  case 3: // vertex-edge histogram kernel
    Kval = vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, -1.0);
    break;
  case 4: // vertex-vertex-edge histogram kernel
    Kval = vertexVertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 5: // edge histogram kernel (Gaussian)
    Kval = edgeHistogramKernel(e1, e2, par[0]);
    break;
  case 6: // vertex histogram kernel (Gaussian)
    Kval = vertexHistogramKernel(v1_label, v2_label, par[0]);
    break;
  case 7: // vertex-edge histogram kernel (Gaussian)
    Kval = vertexEdgeHistogramKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 8: // geometric random walk kernel
    Kval = geometricRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 9: // exponential random walk kernel
    Kval = exponentialRandomWalkKernel(e1, e2, v1_label, v2_label, par[0]);
    break;
  case 10: // k-step random walk kernel
    Kval = kstepRandomWalkKernel(e1, e2, v1_label, v2_label, par);
    break;
  default:
    Kval = 0;
    break;
  }
  return Kval;
}

// get information of graphs from an input R list
void getGraphInfo(List& graph_info_list, vector<MatrixXi>& E, vector<vector<Int>>& V_label, vector<Int>& V_count, vector<Int>& E_count, vector<Int>& D_max) {
  List tmpl;
  SEXP tmpr;
  for (Int i = 0; i < (Int)graph_info_list.size(); ++i) {
    tmpl = graph_info_list[i];
    // an edge matrix
    tmpr = tmpl[0];
    IntegerMatrix E_each(tmpr);
    E.push_back(as<Map<MatrixXi>>(E_each).array() - 1);
    // a vector of a vertex attribute
    tmpr = tmpl[1];
    IntegerVector V_label_each(tmpr);
    V_label.push_back(as<vector<Int>>(V_label_each));
    // the number of vertices
    V_count.push_back(tmpl[2]);
    // the number of edges
    E_count.push_back(tmpl[3]);
    // the maximum degrees
    D_max.push_back(tmpl[4]);
  }
}
// [[Rcpp::export]]
NumericMatrix CalculateKernelCpp(List graph_info_list, NumericVector par_r, Int kernel_type) {
  vector<MatrixXi> E;
  vector<vector<Int>> V_label;
  vector<Int> V_count;
  vector<Int> E_count;
  vector<Int> D_max;
  NumericMatrix K(graph_info_list.size(), graph_info_list.size());
  getGraphInfo(graph_info_list, E, V_label, V_count, E_count, D_max);
  vector<double> par(par_r.begin(), par_r.end());
  if (kernel_type == 11) {
    WLKernelMatrix(E, V_label, V_count, E_count, D_max, (Int)par[0], K);
  } else {
    vector<Int> idx(graph_info_list.size());
    iota(idx.begin(), idx.end(), 0);
    for (auto&& i : idx) {
      for (auto&& j : idx) {
	K(i, j) = computeKernelValue(E[i], E[j], V_label[i], V_label[j], par, kernel_type);
	K(j, i) = K(i, j);
      }
    }
  }
  return K;
}

// ========================================================= //
// ==================== Graphlet kernel ==================== //
// ========================================================= //
// ===== graphlet kernel for k = 4 ===== //
Int find_min(Int a, Int b, Int c) {
  Int m;
  Int mini = a;
  if (b < mini) mini = b;
  if (c < mini) mini = c;
  if (mini == a) {
    if (mini == b) {
      if (mini == c) {
	m = 7;
      } else {
	m = 4;
      }
    } else {
      if (mini == c) {
	m = 5;
      } else {
	m = 1;
      }
    }
  } else {
    if (mini == b) {
      if (mini == c) {
	m = 6;
      } else {
	m = 2;
      }
    } else {
      m = 3;
    }
  }
  return m;
}
void card_ThreeInter(vector<Int>& L1, vector<Int>& L2, vector<Int>& L3, vector<Int>& card) {
  card.resize(7);
  fill(card.begin(), card.end(), 0);
  Int i = 0, j = 0, k = 0;

  while (i < (Int)L1.size() && j < (Int)L2.size() && k < (Int)L3.size()) {
    Int m = find_min(L1[i], L2[j], L3[k]);
    card[m - 1] += 1;
    switch(m) {
    case 1:
      i++; break;
    case 2:
      j++; break;
    case 3:
      k++; break;
    case 4:
      i++; j++; break;
    case 5:
      i++; k++; break;
    case 6:
      j++; k++; break;
    case 7:
      i++; j++; k++; break;
    }
  }

  if (i < (Int)L1.size() || j < (Int)L2.size() || k < (Int)L3.size()) {
    if (i >= (Int)L1.size() && j >= (Int)L2.size()) {
      card[2] += (Int)L3.size() - k;
      k = (Int)L3.size();
    } else {
      if (i >= (Int)L1.size() && k >= (Int)L3.size()) {
	card[1] += (Int)L2.size() - j;
	j = (Int)L2.size();
      } else {
	if (j >= (Int)L2.size() && k >= (Int)L3.size()) {
	  card[0] += (Int)L1.size() - i;
	  i = (Int)L1.size();
	} else {
	  if (i >= (Int)L1.size()) {
	    while (j < (Int)L2.size() && k < (Int)L3.size()) {
	      if (L2[j] < L3[k]) {
		card[1]++;
		j++;
	      } else {
		if (L2[j] > L3[k]) {
		  card[2]++;
		  k++;
		} else {
		  card[5]++;
		  j++;
		  k++;
		}
	      }
	    }
	  } else {
	    if (j >= (Int)L2.size()) {
	      while (i < (Int)L1.size() && k < (Int)L3.size()) {
		if (L1[i] < L3[k]) {
		  card[0]++;
		  i++;
		} else {
		  if (L1[i] > L3[k]) {
		    card[2]++;
		    k++;
		  } else {
		    card[4]++;
		    i++;
		    k++;
		  }
		}
	      }
	    } else {
	      if (k >= (Int)L3.size()) {
		while (i < (Int)L1.size() && j < (Int)L2.size()) {
		  if (L1[i] < L2[j]) {
		    card[0]++;
		    i++;
		  } else {
		    if (L1[i] > L2[j]) {
		      card[1]++;
		      j++;
		    } else {
		      card[3]++;
		      i++;
		      j++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (i < (Int)L1.size() || j < (Int)L2.size() || k < (Int)L3.size()) {
    if (i >= (Int)L1.size() && j >= (Int)L2.size()) {
      card[2] += (Int)L3.size() - k;
    } else if (i >= (Int)L1.size() && k >= (Int)L3.size()) {
      card[1] += (Int)L2.size() - j;
    } else if (j >= (Int)L2.size() && k >= (Int)L3.size()) {
      card[0] += (Int)L1.size() - i;
    }
  }
}
void getIndices(vector<Int>& o_set1, vector<Int>& o_set2, vector<Int>& inter, vector<Int>& diff1, vector<Int>& diff2) {
  vector<Int> inter_(min(o_set1.size(), o_set2.size()), -1);
  vector<Int> diff1_(max(o_set1.size(), o_set2.size()), -1);
  vector<Int> diff2_(max(o_set1.size(), o_set2.size()), -1);

  Int i = 0, j = 0;
  while (i < (Int)o_set1.size() && j < (Int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      diff1_[i] = o_set1[i];
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      diff2_[j] = o_set2[j];
      j++;
    } else {
      inter_[i] = o_set1[i];
      i++;
      j++;
    }
  }

  if (i < (Int)o_set1.size()) {
    for (Int k = i; k < (Int)o_set1.size(); ++k) {
      diff1_[k] = o_set1[k];
    }
  } else if (j < (Int)o_set2.size()) {
    for (Int k = j; k < (Int)o_set2.size(); ++k) {
      diff2_[k] = o_set2[k];
    }
  }

  inter.clear();
  for (auto&& x : inter_) {
    if (x >= 0) inter.push_back(x);
  }
  diff1.clear();
  for (auto&& x : diff1_) {
    if (x >= 0) diff1.push_back(x);
  }
  diff2.clear();
  for (auto&& x : diff2_) {
    if (x >= 0) diff2.push_back(x);
  }
}
template<typename V>
void countGraphletsFour(vector<vector<Int>>& al, V&& count) {
  double n = (double)al.size();
  vector<double> w = {1.0/12.0, 1.0/10.0, 1.0/8.0, 1.0/6.0, 1.0/8.0, 1.0/6.0, 1.0/6.0, 1.0/4.0, 1.0/4.0, 1.0/2.0, 0};
  vector<Int> inter, diff1, diff2, card;
  vector<double> inter_count(11);
  vector<Int> v;
  vector<Int>::iterator it;

  double m = 0.0;
  for (auto&& vec : al) {
    m += (double)vec.size();
  }
  m /= 2.0;

  vector<Int> v1(al.size());
  iota(v1.begin(), v1.end(), 0);
  for (auto&& i : v1) {
    for (auto&& j : al[i]) {
      double K = 0.0;
      fill(inter_count.begin(), inter_count.end(), 0.0);
      getIndices(al[i], al[j], inter, diff1, diff2);
      for (auto&& k : inter) {
	card_ThreeInter(al[i], al[j], al[k], card);
	inter_count[0] += 0.5 * (double)card[6];
	inter_count[1] += 0.5 * (double)(card[3] - 1.0);
	inter_count[1] += 0.5 * (double)(card[4] - 1.0);
	inter_count[1] += 0.5 * (double)(card[5] - 1.0);
	inter_count[2] += 0.5 * (double)card[0];
	inter_count[2] += 0.5 * (double)card[1];
	inter_count[2] += (double)card[2];
	inter_count[6] += n - (double)accumulate(card.begin(), card.end(), 0);
	K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) + 0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff1.size());
      sort(diff1.begin(), diff1.end());
      sort(al[i].begin(), al[i].end());
      it = set_difference(diff1.begin(), diff1.end(), al[i].begin(), al[i].end(), v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
	card_ThreeInter(al[i], al[j], al[k], card);
	inter_count[1] += 0.5 * (double)card[6];
	inter_count[2] += 0.5 * (double)card[3];
	inter_count[2] += 0.5 * (double)card[4];
	inter_count[4] += 0.5 * (double)(card[5] - 1.0);
	inter_count[3] += 0.5 * (double)(card[0] - 2.0);
	inter_count[5] += 0.5 * (double)card[1];
	inter_count[5] += (double)card[2];
	inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
	K += 0.5 * (double)card[6] + 0.5 * (double)card[4] + 0.5 * (double)(card[5] - 1.0) + card[2];
      }
      v.clear();
      v.resize(diff2.size());
      sort(diff2.begin(), diff2.end());
      it = set_difference(diff2.begin(), diff2.end(), v1.begin(), v1.end(), v.begin());
      v.resize(it - v.begin());
      for (auto&& k : v) {
	card_ThreeInter(al[i], al[j], al[k], card);
	inter_count[1] += 0.5 * (double)card[6];
	inter_count[2] += 0.5 * (double)card[3];
	inter_count[4] += 0.5 * (double)(card[4] - 1.0);
	inter_count[2] += 0.5 * (double)card[5];
	inter_count[5] += 0.5 * (double)card[0];
	inter_count[3] += 0.5 * (double)(card[1] - 2.0);
	inter_count[5] += (double)card[2];
	inter_count[7] += n - (double)accumulate(card.begin(), card.end(), 0);
	K += 0.5 * (double)card[6] + 0.5 * (double)(card[4] - 1.0) + 0.5 * (double)card[5] + card[2];
      }
      inter_count[8] += m + 1.0 - (double)v1.size() - (double)al[i].size() - K;
      inter_count[9] += (n - (double)inter.size() - (double)diff1.size() - (double)diff2.size())
	* (n - (double)inter.size() - (double)diff1.size() - (double)diff2.size() - 1.0) / 2
	- (m + 1.0 - (double)v1.size() - (double)al[i].size() - K);

      for (Int k = 0; k < (Int)count.size(); ++k) {
	count(k) += inter_count[k] * w[k];
      }
    }
  }

  count(10) = n * (n - 1.0) * (n - 2.0) * (n - 3.0) / (4.0 * 3.0 * 2.0) - count.head(10).sum();
}
// ===== graphlet kernel for k = 3 ===== //
void getCardinality(vector<Int>& o_set1, vector<Int>& o_set2, vector<double>& card) {
  card.resize(3);
  fill(card.begin(), card.end(), 0.0);
  Int i = 0, j = 0;
  while (i < (Int)o_set1.size() && j < (Int)o_set2.size()) {
    if (o_set1[i] < o_set2[j]) {
      card[0] += 1.0;
      i++;
    } else if (o_set1[i] > o_set2[j]) {
      card[1] += 1.0;
      j++;
    } else {
      i++;
      j++;
      card[2] += 1.0;
    }
  }
  card[0] += (double)((Int)o_set1.size() - i);
  card[1] += (double)((Int)o_set2.size() - j);
}
template<typename V>
void countGraphletsThree(vector<vector<Int>>& al, V&& count) {
  double n = (double)al.size();
  vector<double> w = {1.0/6.0, 1.0/4.0, 1.0/2.0};
  vector<double> card(3);

  vector<Int> L1(al.size());
  iota(L1.begin(), L1.end(), 0);
  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      getCardinality(al[i], al[j], card);
      count(0) += w[0] * card[2];
      count(1) += w[1] * (card[0] + card[1] - 2.0);
      count(2) += w[2] * (n - accumulate(card.begin(), card.end(), 0.0));
    }
  }
  count(3) = n * (n - 1.0) * (n - 2.0) / 6.0 - (count(0) + count(1) + count(2));
}
// ===== connected graphlet kernel for k = 5 ===== //
void getMinValue(SparseMatrix<Int>& am, vector<Int>& idx, vector<Int>& sums) {
  sums.clear();
  sums.resize(idx.size());
  fill(sums.begin(), sums.end(), 0);
  for (Int i = 0; i < (Int)idx.size(); ++i) {
    Int k = idx[i];
    for (SparseMatrix<Int>::InnerIterator it(am, k); it; ++it) {
      if(find(idx.begin(), idx.end(), it.row()) != idx.end()) {
	sums[i] += it.value();
      }
    }
  }
  sums.push_back(1);
}
template<typename V>
void countConnectedGraphletsFive(SparseMatrix<Int>& am, vector<vector<Int>>& al, V&& count) {
  vector<double> w = {1.0/120.0, 1.0/72.0, 1.0/48.0, 1.0/36.0, 1.0/28.0, 1.0/20.0, 1.0/14.0, 1.0/10.0, 1.0/12.0, 1.0/8.0, 1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0/12.0, 1.0/12.0, 1.0/4.0, 1.0/4.0, 1.0/2.0, 0.0, 0.0, 0.0};
  Int n = (Int)am.rows();
  vector<Int> L1(n);
  iota(L1.begin(), L1.end(), 0);
  vector<Int> idx(5);
  vector<Int> sums;

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
	if (k != i) {
	  for (auto&& l : al[k]) {
	    if (l != i && l != j) {
	      for (auto&& m : al[l]) {
		if (m != i && m != j && m != k) {
		  Int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(i, m) + am.coeff(j, l) + am.coeff(j, m) + am.coeff(k, m);
		  if (aux == 6) {
		    count[0] += w[0];
		  } else if (aux == 5) {
		    count[1] += w[1];
		  } else if (aux == 4) {
		    idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
		    getMinValue(am, idx, sums);
		    Int aux1 = *min_element(sums.begin(), sums.end());
		    if (aux1 == 2) {
		      count[3] += w[3];
		    } else {
		      count[2] += w[2];
		    }
		  } else if (aux == 3) {
		    idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
		    getMinValue(am, idx, sums);
		    sort(sums.begin(), sums.end());
		    if (sums[0] == 1) {
		      count[8] += w[8];
		    } else if (sums[1] == 3) {
		      count[4] += w[4];
		    } else if (sums[2]== 2) {
		      count[13] += w[13];
		    } else {
		      count[5] += w[5];
		    }
		  } else if (aux == 2) {
		    idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
		    getMinValue(am, idx, sums);
		    vector<Int> aux1;
		    copy(sums.begin(), sums.end(), back_inserter(aux1));
		    sort(aux1.begin(), aux1.end());
		    if (aux1[0] == 1) {
		      if (aux1[2] == 2) {
			count[15] += w[15];
		      } else {
			count[9] += w[9];
		      }
		    } else {
		      if (aux1[3] == 2) {
			count[10] += w[10];
		      } else {
			vector<Int> ind;
			for (Int ii = 0; ii < (Int)sums.size(); ++ii) {
			  if (sums[ii] == 3) ind.push_back(ii);
			}
			// idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
			if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
			  count[6] += w[6];
			} else {
			  count[14] += w[14];
			}
		      }
		    }
		  } else if (aux == 1) {
		    idx[0] = i; idx[1] = j; idx[2] = k; idx[3] = l; idx[4] = m;
		    getMinValue(am, idx, sums);
		    vector<Int> aux1;
		    copy(sums.begin(), sums.end(), back_inserter(aux1));
		    sort(aux1.begin(), aux1.end());
		    if (aux1[0] == 2) {
		      count[7] += w[7];
		    } else if (aux1[1] == 1) {
		      count[17] += w[17];
		    } else {
		      vector<Int> ind;
		      for (Int ii = 0; ii < (Int)sums.size(); ++ii) {
			if (sums[ii] == 3) ind.push_back(ii);
		      }
		      for (Int ii = 0; ii < (Int)sums.size(); ++ii) {
			if (sums[ii] == 1) ind.push_back(ii);
		      }
		      if (am.coeff(idx[ind[0]], idx[ind[1]]) == 1) {
			count[16] += w[16];
		      } else {
			count[11] += w[11];
		      }
		    }
		  } else {
		    count[12] += w[12];
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    // count graphlets of type 20
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
	if (k != i && am.coeff(i, k) == 0) {
	  for (auto&& l : al[k]) {
	    if (l != i && l != j && am.coeff(i, l) == 0 && am.coeff(j, l) == 0) {
	      for (auto&& m : al[k]) {
		if (m != i && m != j && m != l && am.coeff(i, m) == 0 && am.coeff(j, m) == 0 && am.coeff(l, m) == 0) {
		  count[19] += w[19];
		}
	      }
	    }
	  }
	}
      }
    }
    // count graphlets of type 19 and 21
    for (Int j = 0; j < (Int)al[i].size() - 3; ++j) {
      for (Int k = j + 1; k < (Int)al[i].size() - 2; ++k) {
	for (Int l = k + 1; l < (Int)al[i].size() - 1; ++l) {
	  for (Int m = l + 1; m < (Int)al[i].size(); ++m) {
	    Int aux = am.coeff(al[i][j], al[i][k]) + am.coeff(al[i][j], al[i][l])
	      + am.coeff(al[i][j], al[i][m]) + am.coeff(al[i][k], al[i][l])
	      + am.coeff(al[i][k], al[i][m]) + am.coeff(al[i][l], al[i][m]);
	    if (aux == 1) {
	      count[18]++;
	    } else if (aux == 0) {
	      count[20]++;
	    }
	  }
	}
      }
    }
  }
}
template<typename V>
// ===== connected graphlet kernel for k = 4 ===== //
void countConnectedGraphletsFour(SparseMatrix<Int>& am, vector<vector<Int>>& al, V&& count) {
  vector<double> w = {1.0/24.0, 1.0/12.0, 1.0/4.0, 0.0, 1.0/8.0, 1.0/2.0};
  Int n = (Int)am.rows();
  vector<Int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
	if (k != i) {
	  for (auto&& l : al[k]) {
	    if (l != i && l != j){
	      Int aux = am.coeff(i, k) + am.coeff(i, l) + am.coeff(j, l);
	      if (aux == 3) {
		count[0] += w[0];
	      } else if (aux == 2) {
		count[1] += w[1];
	      } else if (aux == 1) {
		if (am.coeff(i, l) == 1) {
		  count[4] += w[4];
		} else {
		  count[2] += w[2];
		}
	      } else {
		count[5] += w[5];
	      }
	    }
	  }
	}
      }
    }

    // count "stars"
    for (Int j = 0; j < (Int)al[i].size() - 2; ++j) {
      for (Int k = j + 1; k < (Int)al[i].size() - 1; ++k) {
	for (Int l = k + 1; l < (Int)al[i].size(); ++l) {
	  if (am.coeff(al[i][j], al[i][k]) == 0 && am.coeff(al[i][j], al[i][l]) == 0 && am.coeff(al[i][k], al[i][l]) == 0) {
	    count[3]++;
	  }
	}
      }
    }
  }
}
// ===== connected graphlet kernel for k = 3 ===== //
template<typename V>
void countConnectedGraphletsThree(SparseMatrix<Int>& am, vector<vector<Int>>& al, V&& count) {
  vector<double> w = {1.0/2.0, 1.0/6.0};
  Int n = (Int)am.rows();
  vector<Int> L1(n);
  iota(L1.begin(), L1.end(), 0);

  for (auto&& i : L1) {
    for (auto&& j : al[i]) {
      for (auto&& k : al[j]) {
	if (k != i) {
	  if (am.coeff(i, k) == 1) {
	    count[1] += w[1];
	  } else {
	    count[0] += w[0];
	  }
	}
      }
    }
  }
}
// [[Rcpp::export]]
NumericMatrix CalculateGraphletKernelCpp(vector<SparseMatrix<Int>>& graph_adj_all, vector<vector<vector<Int>>>& graph_adjlist_all, Int k, bool connected) {
  // decrement one to start indices from zero
  for (auto&& X : graph_adjlist_all)
    for (auto&& vec : X)
      for (auto&& x : vec) x--;

  MatrixXd freq;
  Int freq_size;
  if (connected) {
    switch (k) {
    case 3: freq_size = 2; break;
    case 4: freq_size = 6; break;
    case 5: freq_size = 21; break;
    }
  } else {
    switch (k) {
    case 3: freq_size = 4; break;
    case 4: freq_size = 11; break;
    }
  }
  freq = MatrixXd::Zero(graph_adjlist_all.size(), freq_size);

  vector<Int> idx_graph(graph_adjlist_all.size());
  iota(idx_graph.begin(), idx_graph.end(), 0);
  for (auto&& i : idx_graph) {
    if (k == 3) {
      if (connected) {
	countConnectedGraphletsThree(graph_adj_all[i], graph_adjlist_all[i], freq.row(i));
      } else {
	countGraphletsThree(graph_adjlist_all[i], freq.row(i));
      }
    } else if (k == 4) {
      if (connected) {
	countConnectedGraphletsFour(graph_adj_all[i], graph_adjlist_all[i], freq.row(i));
      } else {
	countGraphletsFour(graph_adjlist_all[i], freq.row(i));
      }
    } else if (k == 5) {
      if (connected) {
	countConnectedGraphletsFive(graph_adj_all[i], graph_adjlist_all[i], freq.row(i));
      }
    }
    if (freq.row(i).sum() != 0) {
      freq.row(i) /= freq.row(i).sum();
    }
  }
  MatrixXd K = freq * freq.transpose();
  return wrap(K);
}

SEXP graphkernels_CalculateKernelCpp(SEXP graph_info_listSEXP, SEXP par_rSEXP, SEXP kernel_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type graph_info_list(graph_info_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type par_r(par_rSEXP);
    Rcpp::traits::input_parameter< Int >::type kernel_type(kernel_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(CalculateKernelCpp(graph_info_list, par_r, kernel_type));
    return rcpp_result_gen;
END_RCPP
}
SEXP graphkernels_CalculateGraphletKernelCpp(SEXP graph_adj_allSEXP, SEXP graph_adjlist_allSEXP, SEXP kSEXP, SEXP connectedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< vector<SparseMatrix<Int>>& >::type graph_adj_all(graph_adj_allSEXP);
    Rcpp::traits::input_parameter< vector<vector<vector<Int>>>& >::type graph_adjlist_all(graph_adjlist_allSEXP);
    Rcpp::traits::input_parameter< Int >::type k(kSEXP);
    Rcpp::traits::input_parameter< bool >::type connected(connectedSEXP);
    rcpp_result_gen = Rcpp::wrap(CalculateGraphletKernelCpp(graph_adj_all, graph_adjlist_all, k, connected));
    return rcpp_result_gen;
END_RCPP
}
