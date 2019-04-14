#include <mex.h>
#include <matrix.h>
#include <cmath>
#include <limits>

/**
 * [mn, mx] = computeDists(l2, d, vec, dists);
 */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {

  // input [mn, mx, scales, dists]
  double mn = (double) *mxGetPr(prhs[0]);
  double mx = (double) *mxGetPr(prhs[1]);
  int scales = (int) *mxGetPr(prhs[2]);
  float* dists = (float*) mxGetPr(prhs[3]);

  const mwSize *dims = mxGetDimensions(prhs[3]);
  size_t l2 = (size_t) dims[0];

  // output [ratio, n]
  const size_t two = 2;
  plhs[0]= mxCreateDoubleMatrix(two, scales, mxREAL);
  double* ratio = (double*) mxGetPr(plhs[0]);

  int n = 1;
  double epsilon = std::pow(2.0, -n);
  double mx2 = epsilon * mx;

  while (mx2 > 2*mn && n < scales) {
    size_t count = 0;
    for (size_t i = 0; i < l2; ++i) {
	for (size_t j = i+1; j < l2; ++j) {
      if (dists[i*(l2-1)-(i*(i+1))/2+j-1] < mx2) {
        count++;
	}
      }
    }
    if (count > 0) {
        int nm1 = n - 1;
        ratio[0 + nm1*two] = epsilon;
        ratio[1 + nm1*two] = count;
        n++;
    }
    epsilon = std::pow(2.0, -n);
    mx2 = epsilon * mx;
  }

  plhs[1] = mxCreateDoubleScalar(n);
}

