#include <math.h>
#include "linear.h"

/*
 * solveLinearSystem()
 * Solve linear system of equations using Gaussian row
 * reduction to place system in upper-echelon form, and then
 * backsolve.
 *
 * AX = B
 * A is RxR, X is RxS, B is RxS
 * M = [A|B] (augmented Rx(R+S) matrix)
 * system is assumed to *not* be singular
 */
int solveLinearSystem(int R, int S, double *M[], double *X[]) {
  int r, c, s;
  int C = R + S;

  for (r = 0; r < R; r++) {
    double max = fabs(M[r][r]);
    int rr, maxr = r;

    for (rr = r+1; rr < R; rr++) {     /* find pivot */
      double val = fabs(M[rr][r]);
      if (val > max) {
	max = val;
	maxr = rr;
      }
    }

    if (maxr != r) {                   /* swap rows if necessary */
      double *tmp = M[r];
      M[r] = M[maxr];
      M[maxr] = tmp;
    }

    if (max <= 0.0) return 0;          /* 0 pivot -- abort! */

    for (rr = r+1; rr < R; rr++) {     /* get 0's below pivot */
      double scale = M[rr][r]/M[r][r];
      for (c = r+1; c < C; c++)
	M[rr][c] -= scale*M[r][c];
    }
  }

  for (s = 0; s < S; s++)              /* back solve for each column of X */
    for (r = R-1; r >= 0; r--) {
      double sum = 0.0;
      for (c = r+1; c < R; c++)
	sum += M[r][c]*X[s][c];
      X[s][r] = (M[r][R+s] - sum)/M[r][r];
    }

  return 1;
}  
