#ifndef RECTFIS_H
#define RECTFIS_H

#include "func2.h"

/*
 * Transformation parameters that map the source point (x,y,z) to
 * the point (x',y',z'). This actually defines a whole class
 * of maps depending on the functions q and alpha.
 *  x' = a*x + b
 *  y' = b*y + d
 *  z' = q(x',y') + alpha(x',y')*z
 */
typedef struct {
  double a, b;
  double c, d;
  Func2 *q;
  Func2 *alpha;
} RectFisMap;

/*
 *  (x',y',z') <-- map(x,y,z)
 */
double rectFisTransform(RectFisMap *map, double src[3], double dst[3]);

/*
 * The following structure is an abstract base class for defining
 * Fractal Interpolation Surfaces (FIS) z = f(x,y) supported by a 
 * regular grid. The FIS is an Iterated Function System (IFS) defined
 * by MxN maps. The surface interpolates the points Z[j][i] = f(X[i],Y[j]).
 *
 *  (N=3)
 *   Y3  *-----*-----*-----*
 *       |     |     |     |
 *   Y2  *-----*-----*-----*
 *       |     |     |     |
 *   Y1  *-----*-----*-----*
 *       |     |     |     |
 *   Y0  *-----*-----*-----*
 *       X0    X1    X2    X3  (M=3)
 *
 */
typedef struct {
  int M, N;
  double *X, *Y, **Z;
  RectFisMap **maps;
} RectFis;

/*
 * Virtual constructor that allocates and initiazes fields
 * common to all derived classes. The q and alpha functions
 * are really what is unique to each derived concrete class --
 * these are set to NULL is each map (i.e., pure virtual methods).
 */
RectFis *rectFisCreate(int M, int N, double *X, double *Y);

/*
 * Given (x,y), these functions determine the sequences of transformations
 * that map some arbirary point (x*,y*,z*) to (x,y,z) (as close as
 * possible anyway).
 */
void rectFisMapSequenceX(RectFis *fis, double x, int len, int seq[]);
void rectFisMapSequenceY(RectFis *fis, double y, int len, int seq[]);

/*
 * This applies the composite sequence of transformations to some 
 * given initial point (x0,y0,z0). The resulting z-value is returned.
 */
double rectFisEvalAux(RectFis *fis, double x0, double y0, double z0,
		      int len, int seqX[], int seqY[]);

/*
 * Approximates z = f(x,y) for FIS using
 * rectFisEvalAux() to do the dirty work.
 */
double rectFisEval(RectFis *fis, double x, double y);

/*
 * Convert a (concrete) FIS to a function.
 */
Func2 *rectFisToFunc2(RectFis *fis);

#endif /* RECTFIS_H */
