#ifndef RIS_H
#define RIS_H

#include "func2.h"

/*
 * Transformation for Recurrent Interpolation Surface (RIS).
 * 
 * w(x,y,z) --> (x',y',z')
 * x' = a*x + b
 * y' = c*y + d
 * z' = q(x',y') + alpha(x',y')*z
 *
 * alpha() is user defined
 * q() defines the RIS "flavor" (e.g., bilinear or hermite).
 */
typedef struct {
  double a, b;
  double c, d;
  Func2 *alpha;
  Func2 *q;
} RisMap;

/*
 * dst = w(src);
 * Assumes src is over the domain cell support.
 */
double risMapTransform(RisMap *map, double src[3], double dst[3]);

/*
 * Transform src = (x, y, z, zx, zy, zxy)
 * into dst = (x', y', z', zx', zy', zxy')
 * where
 * x' = a*x + b
 * y' = c*y + d
 * z' = q(x',y') + alpha(x',y')*z
 * zx' = qx(x',y') + Sx*alpha(x',y')*zx + alpha_x(x',y')*z
 * zy' = qy(x',y') + Sy*alpha(x',y')*zy + alpha_y(x',y')*z
 * zxy' = qxy(x',y') + Sx*Sy*alpha(x',y')*zxy +
 *           Sx*alpha_y(x',y')*zx + Sy*alpha_x(x',y')*zy
 *           alpha_xy(x',y')*z
 */
void risMapTransformWithDerivs(RisMap *map, int Sx, int Sy,
                               double src[6], double dst[6]);

/*
 * Recurrent Interpolation Surface (RIS) z = f(x,y).
 * Interpolates (M+1)*(N+1) grid points (x[i], y[j], z[j][i]).
 * Surface is considered topologically equivalent to
 * a torus so left edge = right edge and bottom edge = top edge;
 * This constraint is assumed to be enforced elsewhere (i.e., we
 * blindly proceed as if this condition is true).
 */
typedef struct {
  int M, N;                  /* (M+1)*(N+1) grid points */
  double *x;                 /* x0 < x1 < ... < xM */
  double *y;                 /* y0 < y1 < ... < yN */
  double **z;                /* original grid points */
  double **zx, **zy, **zxy;  /* derivatives at grid points */
  int Sx, Sy;                /* size of domain cells */
  RisMap **maps;             /* M*N transformation */
} Ris;

/*
 * Allocates and initializes bilinear RIS.
 * M, N : grid size
 * Sx, Sy: size of domain cell
 * x0, y0 : origin
 * W, H : size of support
 * alpha : "roughness" function (note: we only store alpha reference
 *   so don't go destroyin' it while its being used here).
 * All z[j][i] are initialized to 0.
 * zx = zy = zxy = NULL since we ignore them.
 */
Ris *createBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
                       double W, double H, Func2 *alpha);
void destroyBilinearRis(Ris *ris);

/*
 * Anytime and of the z[j][i] values change or alpha() is modified
 * externally, this needs to be called.
 */
void computeBilinearRisMaps(Ris *ris);

/*
 * Allocates an initializes hermite RIS.
 * See createBilinearRis() comments above.
 * zx, zy, zxy are now allocated.
 */
Ris *createHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                      double W, double H, Func2 *alpha);
void destroyHermiteRis(Ris *ris);

/*
 * Anytime and of the z[j][i], zx[j][i], zy[j][i], zxy[j][i] values change 
 * or alpha() is modified externally, this needs to be called.
 */
void computeHermiteRisMaps(Ris *ris);

/*
 * We can evaluate a RIS z = f(x,y) by finding a sequence of transformations
 * that map some initial point to the point (x,y,f(x,y)). 
 * We find this sequence by observing the effect of each map
 * on the x & y axes independently.
 */
void risMapSeqX(Ris *ris, double x, int n, int seqx[]);
void risMapSeqY(Ris *ris, double y, int n, int seqy[]);
double risEvalSeq(Ris *ris, int n, int seqx[], int seqy[]);

/*
 * Combines the the process of finding the map sequences and
 * then applying the maps.
 */
double risEval(Ris *ris, double x, double y);

/*
 * Given map sequence, computes
 * z = f(x,y)
 * zx = fx(x,y)
 * zy = fy(x,y)
 * zxy = fxy(x,y)
 */
void risEvalSeqWithDerivs(Ris *ris, int n, int seqx[], int seqy[],
                          double *z, double *zx, double *zy, double *zxy);

/*
 * Treat an RIS like a function.
 */
Func2 *risToFunc2(Ris *ris);

#endif /* RIS_H */
