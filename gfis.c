/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "raytrace.h"
#include "gfis.h"
#include "linear.h"
#include "bbox.h"
#include "rectgrid.h"

#define EVEN(n) (((n)&1)==0)
#define SEQUENCE_LEN 16         /* for normal computation */
#ifdef COMMENT_OUT_BLAH
#define MAX_SEQUENCE_LEN 14     /* for ray intersection */
#define HORZ_EPSILON 1e-6
#define VERT_EPSILON 1e-6
#endif
#define MAX_SEQUENCE_LEN 16     /* for ray intersection */
#define HORZ_EPSILON 1e-7
#define VERT_EPSILON 1e-7
#define EVAL_BILINEAR(q, x, y) \
  (q[0] + (x)*(q[1] + q[3]*(y)) + q[2]*(y))

/*
 * x' = a*x + b
 * y' = c*y + d
 * z' = q(x',y') + alpha(x',y')*z
 *
 * Both q(x,y) and alpha(x,y) are bilinear functions.
 *   q(x,y) = q0 + q1*x + q2*y + q3*x*y
 */
typedef struct {        /* transformation */
  double a, b;          /* x-contraction coefficients */
  double c, d;          /* y-contraction coefficients */
  double q[4];          /* q(x,y) coefficients */
  double alpha[4];      /* alpha(x,y) coefficients */
} MAP;

/*
 * Our Fractal Interpolation Surfaces (FIS) are height fields
 * defined by Fractal Interpolation Functions (FIF) each with a
 * rectangular support S = [x0,xM] X [y0,yN] that interpolates
 * (M+1)x(N+1) knots defined over a regular grid:
 *
 * y3 *------*------*------*
 *    |      |      |      |
 *    |      |      |      |
 * y2 *------*------*------*
 *    |      |      |      |
 *    |      |      |      |
 * y1 *------*------*------*
 *    |      |      |      |
 *    |      |      |      |
 * y0 *------*------*------*
 *    x0    x1     x2     x3      (M = 3, N = 3 in this figure)
 *
 *  At each point in the grid we have a z-value:
 *    z(i,j) = f(x[i], y[j])
 *  The fractal surface is defined by MxN transformations.
 */
typedef struct {         /* fractal interpolation surface (FIS) */
  int M, N;              /* num of knots - 1 in x & y direction */  
  POINT3 **knots;        /* (M+1) x (N+1) grid of knots */
  double **alpha;        /* alpha value corresponding to each grid point */
  MAP **maps;            /* M x N transformations */
  BBOX bbox;             /* bounding box */
  double color[3];       /* rgb color of surface */
} FIS;

/*
 * mapPoint()
 * Map source point to destination point with given transform.
 * Note that the transformation can be done in place (i.e. src == dst
 * is ok).
 */
static
void mapPoint(MAP *map, POINT3 *src, POINT3 *dst) {
  double x = map->a*src->x + map->b;
  double y = map->c*src->y + map->d;
  dst->z = EVAL_BILINEAR(map->q, x, y) + 
    EVAL_BILINEAR(map->alpha, x, y)*src->z;
  dst->x = x;
  dst->y = y; 
}

/*
 * color()
 * Method for retrieving solid color of FIS object at hit point.
 */
static
void color(OBJECT *this,double hitPoint[3], HIT_INFO *info, double color[3]) {
  FIS *fis = (FIS *) this->data;
  color[0] = fis->color[0];
  color[1] = fis->color[1];
  color[2] = fis->color[2];
}

/*
 * setColor()
 * Method for setting solid color of FIS object.
 */
static
void setColor(OBJECT *this, double color[3]) {
  FIS *fis = (FIS *) this->data;
  fis->color[0] = color[0];
  fis->color[1] = color[1];
  fis->color[2] = color[2];
}

/*
 * bilinearExtrema()
 * Computes the minimum and maximum of a bilinear function f(x,y)
 * over a rectangular region. Since their are no local extrema
 * for a bilinear function, the extrema will occur at the corners
 * of the support.
 */
static
void bilinearExtrema(double f[4], 
		     double minx, double maxx, double miny, double maxy,
		     double *minf, double *maxf) {
  double minz, maxz;
  double z1, z2, z3;
  
  minz = maxz = EVAL_BILINEAR(f, minx, miny);
  z1 = EVAL_BILINEAR(f, maxx, miny);
  z2 = EVAL_BILINEAR(f, minx, maxy);
  z3 = EVAL_BILINEAR(f, maxx, maxy);

  if (z1 > maxz) maxz = z1;
  else if (z1 < minz) minz = z1;
  if (z2 > maxz) maxz = z2;
  else if (z2 < minz) minz = z2;
  if (z3 > maxz) maxz = z3;
  else if (z3 < minz) minz = z3;

  *minf = minz;
  *maxf = maxz;
}

/*
 * subBoundingBox();
 * Given a bounding box (axes alligned) and transformation, 
 * find the bounding box for the transformed input bounding box.
 * Returns the maximum magnitude of alpha(x,y) over the sub-bounding
 * box's support.
 */
static
double subBoundingBox(BBOX *bbox, MAP *map, BBOX *subBbox) {
  double alphaMin, alphaMax;
  double qMin, qMax;
  double zmin, zmax;
  double z1, z2, z3;

  if (map->a >= 0) {
    subBbox->min.x = map->a*bbox->min.x + map->b;
    subBbox->max.x = map->a*bbox->max.x + map->b;
  } else {
    subBbox->min.x = map->a*bbox->max.x + map->b;
    subBbox->max.x = map->a*bbox->min.x + map->b;
  }

  if (map->c >= 0) {
    subBbox->min.y = map->c*bbox->min.y + map->d;
    subBbox->max.y = map->c*bbox->max.y + map->d;
  } else {
    subBbox->min.y = map->c*bbox->max.y + map->d;
    subBbox->max.y = map->c*bbox->min.y + map->d;
  }
  
  bilinearExtrema(map->q, 
		  subBbox->min.x, subBbox->max.x,
		  subBbox->min.y, subBbox->max.y,
		  &qMin, &qMax);

  bilinearExtrema(map->alpha, 
		  subBbox->min.x, subBbox->max.x,
		  subBbox->min.y, subBbox->max.y,
		  &alphaMin, &alphaMax);

  zmin = zmax = alphaMin*bbox->min.z;
  z1 = alphaMin*bbox->max.z;
  z2 = alphaMax*bbox->min.z;
  z3 = alphaMax*bbox->max.z;

  if (z1 < zmin) zmin = z1;
  else if (z1 > zmax) zmax = z1;
  if (z2 < zmin) zmin = z2;
  else if (z2 > zmax) zmax = z2;
  if (z3 < zmin) zmin = z3;
  else if (z3 > zmax) zmax = z3;

  subBbox->min.z = qMin + zmin;
  subBbox->max.z = qMax + zmax;

  alphaMin = fabs(alphaMin);
  alphaMax = fabs(alphaMax);
  if (alphaMin > alphaMax)
    alphaMax = alphaMin;

  return alphaMax;
}

/*
 * rayHitAux()
 * Auxilary recursive routine for determining if a ray intersects a FIS.
 * We acomplish this by recursively generating smaller and smaller
 * bounding boxes that envelope a piece of the surface as well as
 * the ray. See comments below for details.
 */
static 
double rayHitAux(FIS *fis,
		 double rayOrg[3], double rayDir[3],
		 int len, int sx[], int sy[]) {
  int k;
  double a,b, c,d, alphaMax;
  BBOX bbox = fis->bbox;
  double trange[2];
  double rr, sum;

  /*
   * Contruct bounding box from sequence of transformations
   * determined so far. We track the xy portion of the
   * resulting composite transformation.
   */
  alphaMax = a = c = 1.0;
  b = d = 0.0;

  for (k = len-1; k >= 0; k--) {
    MAP *map = &fis->maps[sy[k]][sx[k]];
    BBOX newbbox;
    
    b += a*map->b;   /* a*(a[i] + b[i]) + b */
    a *= map->a;
		 
    d += c*map->d;   /* c*(c[j] + d[j]) + d */
    c *= map->c;

    alphaMax = subBoundingBox(&bbox, map, &newbbox);
    bbox = newbbox;
  }

  /*
   * If ray misses sub-surface's bounding box then return "miss".
   */
  if (!rayHitsBoundingBox(&bbox, rayOrg, rayDir, trange))
    return -1.0;

  /*
   * If the subsurface bounding box is sufficiently small, then
   * we'll just return the average value of the entry and
   * exit t-values for where the ray intersected the bounding box.
   */
  if (len >= MAX_SEQUENCE_LEN || 
      (fabs(a) <= HORZ_EPSILON && fabs(c) <= HORZ_EPSILON &&
       alphaMax <= VERT_EPSILON))
    return 0.5*(trange[0] + trange[1]);  /* return average t-value */
  
  /*
   * We traverse through the fif's maps in an order that depends
   * on the direction of the ray as it is projected onto
   * the xy-plane. We track the next map in the sequence of transformations
   * that may contain the ray. Because of the ordering, 
   * the first ray intersection we encounter is the closest.
   */
  for (sum = 0; sum <= fis->M + fis->N - 2; sum++)
    for (rr = 0; rr < fis->N; rr++) {
      int i,j;
      double t;

      i = sum - rr;
      if (i < 0 || i >= fis->M) continue;

      j = rr;

      if (rayDir[0] < 0.0) i = fis->M - i - 1;
      if (rayDir[1] < 0.0) j = fis->N - j - 1;
      
      if (a < 0.0) i = fis->M - i - 1;
      if (c < 0.0) j = fis->N - j - 1;

      sx[len] = i;
      sy[len] = j;

      if ((t = rayHitAux(fis, rayOrg, rayDir, len+1, sx, sy)) > EPSILON)
	return t;  /* hit! -- unwind recursion */
    }

  return -1.0;  /* ray missed surface */
}

/*
 * rayHit()
 * Method for determining ray intersection with fractal surface.
 */
static
double rayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
	      HIT_INFO *hitInfo) {
  FIS *fis = (FIS *) this->data;
  int sx[MAX_SEQUENCE_LEN], sy[MAX_SEQUENCE_LEN];

  return rayHitAux(fis, rayOrg, rayDir, 0, sx, sy);
} 

/*
 * mapSequenceX(), mapSequenceY()
 * Determine sequence of transformations that map any
 * point in the fif's domain to (x,y) (approximately).
 */

static
void mapSequenceX(FIS *fis, double x, int len, int sx[]) {
  double a,b;
  int M = fis->M;
  int k, maxk = len-1;

  a = 1.0;    /* identity */
  b = 0.0;

  for (k = 0; k <= maxk; k++) {
    int i;

    if (a >= 0.0) {
      i = 0;
      while (i < M-1 && x > a*fis->knots[0][i+1].x + b)
	i++;
    } else {
      i = M-1;
      while (i > 0 && x > a*fis->knots[0][i].x + b)
	i--;
    }

    b += a*fis->maps[0][i].b;
    a *= fis->maps[0][i].a;

    sx[maxk - k] = i;
  }
}

static
void mapSequenceY(FIS *fis, double y, int len, int sy[]) {
  double c,d;
  int N = fis->N;
  int k, maxk = len-1;

  c = 1.0;  /* identity */
  d = 0.0;

  for (k = 0; k <= maxk; k++) {
    int j;

    if (c >= 0.0) {
      j = 0;
      while (j < N-1 && y > c*fis->knots[j+1][0].y + d)
	j++;
    } else {
      j = N-1;
      while (j > 0 && y > c*fis->knots[j][0].y + d)
	j--;
    }

    d += c*fis->maps[j][0].d;
    c *= fis->maps[j][0].c;

    sy[maxk - k] = j;
  }
}

/*
 * evalFis()
 * Given that the appropriate sequence of transformations have already
 * be computed (via mapSequenceX() and mapSequenceY()), this repeatedly
 * transforms a point (initially set to one of the middle knot values)
 * with the given sequence of maps. The resulting z-value gives us the
 * (approximate) value z = f(x,y).
 */
static
double evalFis(FIS *fis, int len, int sx[], int sy[]) {
  POINT3 p;
  int r = fis->N/2, c = fis->M/2;
  int k;

  p.x = fis->knots[r][c].x;
  p.y = fis->knots[r][c].y;
  p.z = fis->knots[r][c].z;

  for (k = 0; k < len; k++)
    mapPoint(&fis->maps[sy[k]][sx[k]], &p, &p);

  return p.z;
}

#ifdef COMMENT_OUT_OLD_MAP_SEQUENCE_XXX

/*
 * mapSequence()
 * Determine sequence of transformations that map any
 * point in the fif's domain to (x,y) (approximately).
 */
static
void mapSequence(FIS *fis, double x, double y, int len, int sx[], int sy[]) {
  double a,b, c,d;
  int M = fis->M, N = fis->N;
  int k, maxk = len-1;

  a = c = 1.0;  /* identity */
  b = d = 0.0;

  for (k = 0; k <= maxk; k++) {
    int i,j;

    if (a >= 0.0) {
      i = 0;
      while (i < M-1 && x > a*fis->knots[0][i+1].x + b)
	i++;
    } else {
      i = M-1;
      while (i > 0 && x > a*fis->knots[0][i].x + b)
	i--;
    }

    if (c >= 0.0) {
      j = 0;
      while (j < N-1 && y > c*fis->knots[j+1][0].y + d)
	j++;
    } else {
      j = N-1;
      while (j > 0 && y > c*fis->knots[j][0].y + d)
	j--;
    }

    b += a*fis->maps[j][i].b;
    a *= fis->maps[j][i].a;

    d += c*fis->maps[j][i].d;
    c *= fis->maps[j][i].c;

    sx[maxk - k] = i;
    sy[maxk - k] = j;
  }
}

#endif

/*
 *  normal()
 *  Computes normal to surface at hit point (x,y,z).
 *
 *  The key recursive relationship that defined f(x,y) is
 *   f(a*x + b, c*y + d) = q(a*x + b, c*y + d) + alpha(a*x + b, c*y + d)*f(x,y)
 *  If we take the partial derivative wrt x on each side we get
 *
 *   fx(a*x + b, c*y + d)*a = 
 *           qx(a*x + b, c*y + d)*a +
 *           alpha(a*x + b, c*y + d)*fx(x,y) +
 *           a*alpha_x(a*x + b, c*y + d)*f(x,y)
 *
 *   fx(a*x + b, c*y + d) = 
 *           qx(a*x + b, c*y + d) +
 *           (1/a)*alpha(a*x + b, c*y + d)*fx(x,y) +
 *           alpha_x(a*x + b, c*y + d)*f(x,y)
 *
 *  which tells us that if we know the derivative fx(x,y) and the value
 *  f(x,y) then we know the derivative fx(a*x + b, c*y + d).
 *  We arrive at a similar equation for the partial deriv wrt y:
 *
 *   fy(a*x + b, c*y + d) = 
 *           qy(a*x + b, c*y + d) +
 *           (1/c)*alpha(a*x + b, c*y + d)*fy(x,y) +
 *           alpha_y(a*x + b, c*y + d)*f(x,y).
 *
 *  The update function for f(x,y) is given by the key recursive equation.
 *
 *  We determine a sequence of transformations that map any initial
 *  point (x0,y0) arbitrarily close to (x,y) -- the point at which
 *  we want to dertermine the normal. We just assume that the
 *  initial partial derivatives are zero:
 *     fx(x0,y0) = fy(x0,y0) = 0
 *  We do know the value of f(x0,y0) which is just one of the knots.
 *  The normal vector N is then
 *     N = (-fx, -fy, 1).
 *  The derivatives of a bilinear function are computed as follows:
 *    q(x,y) = q0 + q1*x + q2*y + q3*x*y
 *    qx(x,y) = q1 + q3*y
 *    qy(x,y) = q2 + q3*x
 */
static
void normal(OBJECT *this, double hitPoint[3], HIT_INFO *info, 
	    double normal[3]) {
  FIS *fis = (FIS *) this->data;
  int k, sx[SEQUENCE_LEN], sy[SEQUENCE_LEN];
  double x,y,z, dfdx,dfdy;

  x = fis->knots[0][0].x;
  y = fis->knots[0][0].y;
  z = fis->knots[0][0].z;
  dfdx = dfdy = 0.0;  /* initial guess at derivatives */
 
#ifdef COMMENT_OUT_OLD_MAP_SEQUENCE_XXX
  mapSequence(fis, hitPoint[0], hitPoint[1], SEQUENCE_LEN, sx, sy);
#endif
  mapSequenceX(fis, hitPoint[0], SEQUENCE_LEN, sx);
  mapSequenceY(fis, hitPoint[1], SEQUENCE_LEN, sy);  
  
  for (k = 0; k < SEQUENCE_LEN; k++) {
    MAP *map = &fis->maps[sy[k]][sx[k]];
    double xx = map->a*x + map->b;
    double yy = map->c*y + map->d;
    double alpha = EVAL_BILINEAR(map->alpha, xx, yy);
    dfdx = 
      map->q[1] + map->q[3]*yy +
      alpha*dfdx/map->a +
      (map->alpha[1] + map->alpha[3]*yy)*z;
    dfdy = 
      map->q[2] + map->q[3]*xx +
      alpha*dfdy/map->c +
      (map->alpha[2] + map->alpha[3]*xx)*z;
    z = EVAL_BILINEAR(map->q, xx, yy) + alpha*z;
    x = xx;
    y = yy;
  }

  x = 1.0/sqrt(dfdx*dfdx + dfdy*dfdy + 1);
  normal[0] = -x*dfdx;
  normal[1] = -x*dfdy;
  normal[2] = x;
}

OBJECT *createGeneralBilinearFisObject(int M, int N, 
				       POINT3 **knots, double **alpha) {
  FIS *fis;
  OBJECT *obj;
  int i,j;
  double dz;
  POINT3 p;

  if ((fis = (FIS *) malloc(sizeof(FIS))) == NULL ||
      (fis->knots = (POINT3 **) malloc((N+1)*sizeof(POINT3 *))) == NULL ||
      (fis->alpha = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (fis->maps = (MAP **) malloc(N*sizeof(MAP *))) == NULL ||
      (obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating bilinear fis");
    exit(-1);
  }
  
  for (j = 0; j <= N; j++)
    if ((fis->knots[j] = (POINT3 *) malloc((M+1)*sizeof(POINT3))) == NULL ||
	(fis->alpha[j] = (double *) malloc((M+1)*sizeof(double))) == NULL) {
      perror("allocating row of knots and alphas for fis");
      exit(-1);
    }

  for (j = 0; j < N; j++)
    if ((fis->maps[j] = (MAP *) malloc(M*sizeof(MAP))) == NULL) {
      perror("allocating row of maps for fis");
      exit(-1);
    }

  fis->M = M;
  fis->N = N;

  for (j = 0; j <= N; j++)
    for (i = 0; i <= M; i++) {
      fis->knots[j][i] = knots[j][i];
      fis->alpha[j][i] = alpha[j][i];
    }

  /*
   * For each map i,j
   *
   * a = (-1)^i/M,
   * b = x[i] - a*x[0],   i even
   *     x[i+1] - a*x[0], i odd
   *
   * c = (-1)^j/N;
   * d = y[j] - c*y[0],   j even
   *   = y[j+1] - c*y[0], j odd
   *
   * The bilinear function alpha(x,y) interpolates the four alpha values
   * defined at the map's destination rectangle's corners. We compute
   * the coefficients for alpha(x,y) by solving the system:
   *   _                                _  _        _      _               _
   *  | 1  x[i]    y[j]    x[i]*y[j]     || alpha[0] |    | alpha[j][i]     | 
   *  | 1  x[i+1]  y[j]    x[i+1]*y[j]   || alpha[1] | =  | alpha[j][i+1]   |
   *  | 1  x[i]    y[j+1]  x[i]*y[j+1]   || alpha[2] |    | alpha[j+1][i]   |
   *  | 1  x[i+1]  y[j+1]  x[I+1]*y[j+1] || alpha[3] |    | alpha[j+1][i+1] |
   *   -                                -  -        -      -               -
   *                                           |                  |
   *             alpha(x,y) coefficients ------                   |
   *                                    knot alpha values --------
   *
   * We solve the following linear system for the coefficients of q(x,y):
   *  let x0' = a*x0 + b,  xM' = a*xM + b,
   *      y0' = c*y0 + d,  yN' = c*yN + d
   *   _                    _  _  _     _                                    _
   *  | 1  x0'  y0'  x0'*y0' || q0 |   | f(x0',y0') - alpha(x0',y0')*f(x0,y0) |
   *  | 1  xM'  y0'  xM'*y0' || q1 | = | f(xM',y0') - alpha(xM',y0')*f(xM,y0) |
   *  | 1  x0'  yN'  x0'*yN' || q2 |   | f(x0',yN') - alpha(x0',yN')*f(x0,yN) |
   *  | 1  xM'  yN'  xM'*yN' || q2 |   | f(xM',yN') - alpha(xM',yN')*f(xM,yN) |
   *   -                    -  -  -     -                                    -
   * 
   * All the above values for x, y, and f(x,y) are just knots values.
   *
   *  i even:  x0' = x[i],    xM' = x[i+1]
   *  i odd:   x0' = x[i+1],  xM' = x[i]
   *
   *  j even:  y0' = y[j],    yN' = x[j+1]
   *  j odd:   y0' = y[j+1],  yN' = x[j]
   */
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++) {
      MAP *map = &fis->maps[j][i];
      int r, ii[2], jj[2];
      double A[4][5], *Ap[4], *X[1];
      double x0_, xM_, y0_, yN_;

      if (EVEN(i)) {
	map->a = 1.0/M;
	map->b = fis->knots[0][i].x - map->a*fis->knots[0][0].x;
	ii[0] = i; ii[1] = i+1;
      } else {
	map->a = -1.0/M;
	map->b = fis->knots[0][i+1].x - map->a*fis->knots[0][0].x;
	ii[0] = i+1; ii[1] = i;
      }

      if (EVEN(j)) {
	map->c = 1.0/N;
	map->d = fis->knots[j][0].y - map->c*fis->knots[0][0].y;
	jj[0] = j; jj[1] = j+1;
      } else {
	map->c = -1.0/N;
	map->d = fis->knots[j+1][0].y - map->c*fis->knots[0][0].y;
	jj[0] = j+1; jj[1] = j;
      }

      /*
       * Solve for coefficients of alpha(x,y).
       */
      A[0][1] = A[2][1] = fis->knots[0][i].x;
      A[1][1] = A[3][1] = fis->knots[0][i+1].x;
      A[0][2] = A[1][2] = fis->knots[j][0].y;
      A[2][2] = A[3][2] = fis->knots[j+1][0].y;
      A[0][4] = alpha[j][i];
      A[1][4] = alpha[j][i+1];
      A[2][4] = alpha[j+1][i];
      A[3][4] = alpha[j+1][i+1];

      for (r = 0; r < 4; r++) {
	A[r][0] = 1.0;
	A[r][3] = A[r][1]*A[r][2];
	Ap[r] = A[r];
      }

      X[0] = map->alpha;

      solveLinearSystem(4, 1, Ap, X);

      /*
       * Solve for coefficients of q(x,y).
       */
      x0_ = fis->knots[0][ii[0]].x;
      xM_ = fis->knots[0][ii[1]].x;
      y0_ = fis->knots[jj[0]][0].y;
      yN_ = fis->knots[jj[1]][0].y;

      A[0][1] = A[2][1] = x0_;
      A[1][1] = A[3][1] = xM_;
      A[0][2] = A[1][2] = y0_;
      A[2][2] = A[3][2] = yN_;

      for (r = 0; r < 4; r++) {
	A[r][0] = 1.0;
	A[r][3] = A[r][1]*A[r][2];
	Ap[r] = A[r];	
      }

      A[0][4] = knots[jj[0]][ii[0]].z - 
	EVAL_BILINEAR(map->alpha, x0_, y0_)*knots[0][0].z;
      A[1][4] = knots[jj[0]][ii[1]].z - 
	EVAL_BILINEAR(map->alpha, xM_, y0_)*knots[0][M].z;
      A[2][4] = knots[jj[1]][ii[0]].z - 
	EVAL_BILINEAR(map->alpha, x0_, yN_)*knots[N][0].z;
      A[3][4] = knots[jj[1]][ii[1]].z - 
	EVAL_BILINEAR(map->alpha, xM_, yN_)*knots[N][M].z;

      X[0] = map->q;
      
      solveLinearSystem(4, 1, Ap, X);
    }

  /*
   * Play chaos game to get bounding box.
   */
  p.x = fis->knots[0][0].x;
  p.y = fis->knots[0][0].y;
  p.z = fis->knots[0][0].z;

  fis->bbox.min = fis->knots[0][0];
  fis->bbox.max = fis->knots[N][M];
  if (fis->bbox.min.z > fis->bbox.max.z) {
    double tmp = fis->bbox.min.z;
    fis->bbox.min.z = fis->bbox.max.z;
    fis->bbox.max.z = tmp;
  }

  for (i = 0; i < 100000; i++) {
    int ii = rand() % M;
    int jj = rand() % N;
    MAP *map = &fis->maps[jj][ii];
    mapPoint(map, &p, &p);
    if (p.z < fis->bbox.min.z)
      fis->bbox.min.z = p.z;
    if (p.z > fis->bbox.max.z)
      fis->bbox.max.z = p.z;
  }

  /*
   * inflate bounding-box by 2.5% in z-direction to be safe.
   */
  dz = 0.025*(fis->bbox.max.z - fis->bbox.min.z);
  if (dz <= 0.0) dz = 1.0;  /* rare case (flat plane) */
  fis->bbox.min.z -= dz;
  fis->bbox.max.z += dz;

  fis->color[0] = 0.2;   /* nice blue surface */
  fis->color[1] = 0.2;
  fis->color[2] = 0.6;

  obj->data = fis;       /* hook fis data to ray trace object */
  
  obj->ka = 0.15;        /* some nice material coefficients */
  obj->kd = 0.40;
  obj->ks = 0.45;
  obj->kt = 0.55;
  obj->ni = 1.52;
  obj->phong = 6.0;

  obj->rayHit = rayHit;  /* set object's polymorphic method operators */
  obj->normal = normal;
  obj->color = color;
  obj->setColor = setColor;

  return obj;
}

/*
 * createGeneralBilinearFisObjectForQuadratic()
 * Quadratic equation:
 *
 *  f(x,y) = A x^2 + B y^2 + C xy + D x + E y + F    (1)
 *
 *  f(ax+b,cy+d) = q(ax+b,cy+d) + alpha f(x,y)       (2)
 *
 *  (ax + b)^2 = a^2 x^2 + 2 ab x + b^2
 *  (cy + d)^2 = c^2 y^2 + 2 cd y + d^2
 *  (ax + b)*(cy + d) = ac xy + ad x + bc y + bd
 *
 *  Plug Equation 1 into Equation 2:
 *
 *  A*(a^2 x^2 + 2 ab x + b^2) + B*(c^2 y^2 + 2 cd y + d^2) + 
 *  C*(ac xy + ad x + bc y + bd) + D*(ax + b) + E(cy + d) + F
 *  =
 *  q0 + q1*(ax + b) + q2*(cy + d) + q3*(ac xy + ad x + bc y + bd) +
 *  alpha*(A x^2 + B y^2 + C xy + D x + E y + F)
 *
 *  Group like terms on each side:
 *
 *  A a^2 = alpha A    (x^2 term)
 *  B c^2 = alpha B    (y^2 term)
 *
 *  therefore alpha = a^2 = c^2
 *
 *  C*ac = q3*ac + alpha*C   (xy term)
 *
 *  therefore q3 = (C*ac - alpha*C)/(ac) = 0
 *
 *  A*2ab + C*ad + D*a = q1*a + alpha*D  (x term)
 *
 *  therefore q1 = (A*2ab + C*ad + D*a - alpha*D)/a
 *                = 2Ab + Cd + D*(1 - a)
 *
 *  B*2cd + C*bc + E*c = q2*c + alpha*E  (y term)
 *
 *  therefore q2 = (B*2cd + C*bc + E*c - alpha*E)/c
 *               = 2Bd + Cb + E*(1 - c)
 *
 *  A*b^2 + B*d^2 + C*bd + D*b + E*d + F =    (constant term)
 *       q0 + q1*b + q2*d + alpha*F
 *
 *  therefore q0 = f(b,d) - (q1*b + q2*d + alpha*F)
 *
 *  Since q3 = 0, we actually have an affine transformation:
 *   _  _     _                 _   _ _     _                _
 *  | x' |   | a      0     0    | | x |   | b                |
 *  | y' | = | 0      c     0    | | y | + | d                |
 *  | z' |   | a*q1  c*q2  alpha | | z |   | q0 + q1*b + q2*d |
 *   -  -     -                 -   - -     -                -
 */
OBJECT *createGeneralBilinearFisObjectForQuadratic(double coeffs[6],
						   int N, /* M = N */
						   double orgx, double orgy,
						   double w, double h) {
  FIS *fis;
  OBJECT *obj;
  int i,j;
  double x,y, dx,dy, a, alpha, dz;
  POINT3 p;

  if ((fis = (FIS *) malloc(sizeof(FIS))) == NULL ||
      (fis->knots = (POINT3 **) malloc((N+1)*sizeof(POINT3 *))) == NULL ||
      (fis->alpha = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (fis->maps = (MAP **) malloc(N*sizeof(MAP *))) == NULL ||
      (obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating bilinear fis");
    exit(-1);
  }

  for (j = 0; j <= N; j++)
    if ((fis->knots[j] = (POINT3 *) malloc((N+1)*sizeof(POINT3))) == NULL ||
	(fis->alpha[j] = (double *) malloc((N+1)*sizeof(double))) == NULL) {
      perror("allocating row of knots and alphas for fis");
      exit(-1);
    }

  for (j = 0; j < N; j++)
    if ((fis->maps[j] = (MAP *) malloc(N*sizeof(MAP))) == NULL) {
      perror("allocating row of maps for fis");
      exit(-1);
    }

  fis->M = fis->N = N;

  /*
   * f(x,y) = A x^2 + B y^2 + C xy + D x + E y + F
   */
#define _A_ coeffs[0]
#define _B_ coeffs[1]
#define _C_ coeffs[2]
#define _D_ coeffs[3]
#define _E_ coeffs[4]
#define _F_ coeffs[5]

#define EVAL_QUAD(x,y) \
  ((x)*(_A_*(x) + _C_*(y) + _D_) + (y)*(_B_*(y) + _E_) +_F_)

  dx = w/N;
  dy = h/N;

  alpha = 1.0/(N*N);

  for (j = 0, y = orgx; j <= N; j++, y += dy)
    for (i = 0, x = orgx; i <= N; i++, x += dx) {
      fis->knots[j][i].x = x;
      fis->knots[j][i].y = y;
      fis->knots[j][i].z = EVAL_QUAD(x,y);
      fis->alpha[j][i] = alpha;
    }

  a = 1.0/N;

  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++) {
      MAP *map = &fis->maps[j][i];
  
      if (EVEN(i)) {
	map->a = a;
	map->b = fis->knots[0][i].x - a*fis->knots[0][0].x;
      } else {
	map->a = -a;
	map->b = fis->knots[0][i+1].x + a*fis->knots[0][0].x;
      }

      if (EVEN(j)) {
	map->c = a;
	map->d = fis->knots[j][0].y - a*fis->knots[0][0].y;
      } else {
	map->c = -a;
	map->d = fis->knots[j+1][0].y + a*fis->knots[0][0].y;
      }

      map->alpha[0] = alpha;
      map->alpha[1] = map->alpha[2] = map->alpha[3] = 0.0;

      map->q[3] = 0.0;
      map->q[2] = 2*_B_*map->d + _C_*map->b + _E_*(1.0 - map->a);
      map->q[1] = 2*_A_*map->b + _C_*map->d + _D_*(1.0 - map->a);
      map->q[0] = 
	EVAL_QUAD(map->b, map->d) -
	(map->q[1]*map->b + map->q[2]*map->d + alpha*_F_);
    }
	
  /*
   * Play chaos game to get bounding box.
   */
  p.x = fis->knots[0][0].x;
  p.y = fis->knots[0][0].y;
  p.z = fis->knots[0][0].z;

  fis->bbox.min = fis->knots[0][0];
  fis->bbox.max = fis->knots[N][N];
  if (fis->bbox.min.z > fis->bbox.max.z) {
    double tmp = fis->bbox.min.z;
    fis->bbox.min.z = fis->bbox.max.z;
    fis->bbox.max.z = tmp;
  }

  for (i = 0; i < 100000; i++) {
    int ii = rand() % N;
    int jj = rand() % N;
    MAP *map = &fis->maps[jj][ii];
    mapPoint(map, &p, &p);
    if (p.z < fis->bbox.min.z)
      fis->bbox.min.z = p.z;
    if (p.z > fis->bbox.max.z)
      fis->bbox.max.z = p.z;
  }

  /*
   * inflate bounding-box by 2.5% in z-direction to be safe.
   */
  dz = 0.025*(fis->bbox.max.z - fis->bbox.min.z);
  if (dz <= 0.0) dz = 1.0;  /* rare case (flat plane) */
  fis->bbox.min.z -= dz;
  fis->bbox.max.z += dz;

  fis->color[0] = 0.2;   /* nice blue surface */
  fis->color[1] = 0.2;
  fis->color[2] = 0.6;

  obj->data = fis;       /* hook fis data to ray trace object */
  
  obj->ka = 0.15;        /* some nice material coefficients */
  obj->kd = 0.40;
  obj->ks = 0.45;
  obj->kt = 0.55;
  obj->ni = 1.52;
  obj->phong = 6.0;

  obj->rayHit = rayHit;  /* set object's polymorphic method operators */
  obj->normal = normal;
  obj->color = color;
  obj->setColor = setColor;

  return obj;
}
  

/*
 * createGridObjectFromGeneralBilinearFisObject()
 * Since ray tracing our FIS's can take forever, we can
 * first "tesselate the surface" into into a buffer
 * and then ray trace that.
 */
OBJECT *createGridObjectFromGeneralBilinearFisObject(OBJECT *obj,
						     int W, int H) {
  FIS *fis = (FIS *) obj->data;
  double xmin,ymin, xmax,ymax;
  double *z;
  double x,y, dx,dy;
  int i,j;

  /*
   * Allocate buffer of samples; Note that it may be difficult
   * to deallocate this memory since you will have lost access to it
   * due to the fact that it will be hidden in a private field of another
   * object.
   */
  if ((z = (double *) malloc(W*H*sizeof(double))) == NULL) {
    perror("z[]:createGridObjectFrom...()");
    exit(-1);
  }

  xmin = fis->knots[0][0].x;
  ymin = fis->knots[0][0].y;

  xmax = fis->knots[fis->N][fis->M].x;
  ymax = fis->knots[fis->N][fis->M].y;

  dx = (xmax - xmin)/(W-1);
  dy = (ymax - ymin)/(H-1);

#define GRID_SEQUENCE_LEN 20

  for (y = ymin, j = 0; j < H; j++, y += dy) {
    int sy[GRID_SEQUENCE_LEN];

    mapSequenceY(fis,y, GRID_SEQUENCE_LEN, sy);

    for (x = xmin, i = 0; i < W; i++, x += dx) {
      int sx[GRID_SEQUENCE_LEN];

      mapSequenceX(fis,x, GRID_SEQUENCE_LEN, sx);

      z[j*W + i] = evalFis(fis, GRID_SEQUENCE_LEN, sx, sy);
    }
  }

  return createSurfaceSupportedByRectangularGrid(W, H,
						 xmin, ymin,
						 xmax - xmin, ymax - ymin,
						 z, W);
}







