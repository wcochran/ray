/*
 * This module is for generating single octave noise
 * and multiple octave turbulance as described in:
 *
 *     Hypertexture
 *     Ken Perlin and Eric M. Hoffert
 *     Computer Graphics, Siggraph Proceedings 1989
 *     Vol. 23, No 3
 *     pp 253--261.
 */
#include <stdlib.h>
#include <math.h>
#include "noise.h"

#define LOGN 8
#define N (1<<LOGN)    /* length of P and G arrays */

static double G3[N][3];  /* N random unit vectors in 3D */
static int P[N];         /* random permutation of first N integers */

#define PHI(i) P[(i)&(N-1)]  /* swizzler : picks random index */

#define TAU3(i,j,k) G3[PHI((i)+PHI((j)+(PHI(k))))]

static int initialized = 0;  /* have we done our precomputations yet? */

static void init(void) {
  int i,j;

  /*
   * Initialize G3.
   * Each vector in G is a uniformly distributed vector
   * on the unit sphere. We pick samples from the unit cube,
   * and discard the samples that are not also w/in the unit sphere.
   */
  for (i = 0; i < N; i++) {
    double v[3];
    double mag;
    do {
      v[0] = 2*drand48() - 1;
      v[1] = 2*drand48() - 1;
      v[2] = 2*drand48() - 1;
      mag = fabs(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    } while (mag < 0.0001 || mag > 1.0);
    mag = sqrt(mag);
    G3[i][0] = v[0]/mag;
    G3[i][1] = v[1]/mag;
    G3[i][2] = v[2]/mag;
  }

  /* 
   * Initialize P.
   * Create random permutation of first N ints.
   */
  for (i = 0; i < N; i++)
    P[i] = i;

  for (i = 0; i < N; i++) {  
    int tmp = P[i];
    j = rand() & (N-1);
    P[i] = P[j];
    P[j] = tmp;
  }

  initialized = 1;  /* initialization is complete */
}

/*             /
 * omega(t) = | 2|t|^3 - 3|t|^2 + 1,  if |t| <= 1
 *            | 0,                    otherwise
 *             \
 */
static double omega(double t) {
  t = fabs(t);
  if (t > 1)
    return 0.0;
  return t*t*(2*t - 3) + 1;
}
  
double noise3(double x, double y, double z) {
  int i,j,k;
  int xfloor = (int) floor(x);
  int yfloor = (int) floor(y);
  int zfloor = (int) floor(z);
  double sum = 0.0;

  if (!initialized)
    init();

  for (i = xfloor+1; i >= xfloor; i--) {
    double u = x - i;
    for (j = yfloor+1; j >= yfloor; j--) {
      double v = y - j;
      for (k = zfloor+1; k >= zfloor; k--) {
	double w = z - k;
	double scale = omega(u)*omega(v)*omega(w);
	double *grad = TAU3(i,j,k);
	double dot = grad[0]*u + grad[1]*v + grad[2]*w;
	sum += scale*dot;
      }
    }
  }

  return sum;
}

double turbulence3(double x, double y, double z, double freq) {
  double t, v[3];
  
  for (t = 0.0; freq >= 1.0; freq *= 0.5) {
    v[0] = freq*x;
    v[1] = freq*y;
    v[2] = freq*z;
    t += fabs(noise3(v[0], v[1], v[2])) / freq;
  }
  
  return t;
}
