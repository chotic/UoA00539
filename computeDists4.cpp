#include <mex.h>
#include <matrix.h>
#include <cmath>
#include <limits>

/**
 * [mn, mx] = computeDists(l2, d, vec, dists);
 */

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[]) {

  // input [l2, d, vec]
  int l2 = (int) *mxGetPr(prhs[0]);
  int d = (int) *mxGetPr(prhs[1]);
  double* vec = (double*) mxGetPr(prhs[2]);
  float* dists = (float*) mxGetPr(prhs[3]);

  double mn = 0;
  double mx = 0;
  for( int i = 0; i < l2; ++i) {
    for (int j = i+1; j < l2; ++j) {
      double sum = 0;
        for (int k = 0; k < d; ++k) {
        double vki = vec[k + i*d];
        double vkj = vec[k + j*d];
        double dv = vki - vkj;
        sum += dv*dv;
      }
      sum = std::sqrt(sum);

      if (sum > mx) {
        mx = sum;
      }

      if (i == 1 && j == 2) {
        mn = sum;
      }
      else {
        if( (sum < mn) && (sum > 0) ) {
            mn = sum;
        }
      }

      dists[i*(l2-1)-(i*(i+1))/2+j-1] = sum;
    }

    //dists[i + i*l2] = std::numeric_limits<double>::max(); //1.e10;
  }

  // output [mn, mx]
  plhs[0] = mxCreateDoubleScalar(mn);
  plhs[1] = mxCreateDoubleScalar(mx);
}
