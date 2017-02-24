/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "raytrace.h"
#include "fis.h"
#include "linear.h"
#include "bbox.h"

#define EVEN(n) (((n)&1)==0)
#define HORZ_EPSILON 1e-7
#define VERT_EPSILON 1e-7
#define SEQUENCE_LEN 16
#define EVAL_BILINEAR(q, x, y) \
  (q[0] + (x)*(q[1] + q[3]*(y)) + q[2]*(y))

/*
 * x' = a*x + b
 * y' = c*y + d
 * z' = q(x,y) + alpha*z
 *
 * q(x,y) = q0 + q1*x + q2*y + q3*x*y
 *
 * Note the same alpha values is shared for all maps, but we give
 * each map its own copy.
 */
typedef struct {        /* transformation */
  double a, b;          /* x-contraction coefficients */
  double c, d;          /* y-contraction coefficients */
  double q[4];          /* bilinear coefficients */
  double alpha;         /* z-contractivity constant */
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
  int M, N;              /* num of knots in x & y direction */
  POINT3 **knots;        /* (M+1) x (N+1) grid of knots */
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
  dst->z = EVAL_BILINEAR(map->q, src->x, src->y) + map->alpha*src->z;
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
 * composeMaps()
 * Build composite map C(x,y,z) = (M2 o M1)(x,y,z)
 */
static
void composeMaps(MAP *map2, MAP *map1, MAP *composite) {
  composite->a = map2->a*map1->a;
  composite->b = map2->a*map1->b + map2->b;
  composite->c = map2->c*map1->c;
  composite->d = map2->c*map1->d + map2->d;

  composite->q[0] = 
    map2->q[0] + 
    map2->q[1]*map1->b + 
    map2->q[2]*map1->d +
    map2->q[3]*map1->b*map1->d +
    map2->alpha*map1->q[0];
  composite->q[1] = 
    map2->q[1]*map1->a +
    map2->q[3]*map1->a*map1->d +
    map2->alpha*map1->q[1];
  composite->q[2] = 
    map2->q[2]*map1->c +
    map2->q[3]*map1->b*map1->c +
    map2->alpha*map1->q[2];
  composite->q[3] = 
    map2->q[3]*map1->a*map1->c +
    map2->alpha*map1->q[3];

  composite->alpha = map2->alpha*map1->alpha;
}

/*
 * subBoundingBox();
 * Given a bounding box (axes alligned) and transformation, 
 * find the bounding box for the transformed input bounding box.
 * Since the map uses a bilinear function for q(x,y), we know
 * the extrema for q(x,y) are at the corners of the rectangle
 * projected onto the xy-plane.
 */
static
void subBoundingBox(BBOX *bbox, MAP *map, BBOX *subBbox) {
  double q, minq, maxq;

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

  minq = maxq = EVAL_BILINEAR(map->q, bbox->min.x, bbox->min.y);
  q = EVAL_BILINEAR(map->q, bbox->max.x, bbox->min.y);
  if (q < minq) minq = q;
  else if (q > maxq) maxq = q;
  q = EVAL_BILINEAR(map->q, bbox->min.x, bbox->max.y);
  if (q < minq) minq = q;
  else if (q > maxq) maxq = q;
  q = EVAL_BILINEAR(map->q, bbox->max.x, bbox->max.y);
  if (q < minq) minq = q;
  else if (q > maxq) maxq = q;
    
  if (map->alpha >= 0.0) {
    subBbox->min.z = minq + map->alpha*bbox->min.z;
    subBbox->max.z = maxq + map->alpha*bbox->max.z;
  } else {
    subBbox->min.z = minq + map->alpha*bbox->max.z;
    subBbox->max.z = maxq + map->alpha*bbox->min.z;
  }
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
		 MAP *compositeMap) {
  double trange[2];
  BBOX bbox;
  int r,c, sum;

  /*
   * Compute the bounding box for piece of the surface
   * corresponding to the given composite transformation.
   * We accomplish this by transforming the bounding box 
   * that bounds the entire FIS with the given composite transformation, 
   * and then computing the resulting smaller bounding box.
   */
  subBoundingBox(&fis->bbox, compositeMap, &bbox);
  
  /*
   * If ray misses sub-surface's bounding box then return "miss".
   */
  if (!rayHitsBoundingBox(&bbox, rayOrg, rayDir, trange))
    return -1.0;

  /*
   * If the composite transformation is sufficiently contractive
   * (i.e. the sub-surface bounding box is sufficiently small) then
   * we'll just return the average value of the entry and
   * exit t-values for where the ray intersected the bounding box.
   */
  if (fabs(compositeMap->a) <= HORZ_EPSILON &&
      fabs(compositeMap->c) <= HORZ_EPSILON &&
      fabs(compositeMap->alpha) <= VERT_EPSILON)
    return 0.5*(trange[0] + trange[1]);  /* return average t-value */

  /*
   * We traverse through the fif's maps in an order that depends
   * on the direction of the ray as it is projected onto
   * the xy-plane. We build a new composite transformation for the
   * corresponding sub-bounding box, and call ourseleves recursively.
   * Because of the ordering, the first ray intersection we encounter
   * is the closest.
   */
  for (sum = 0; sum <= fis->M + fis->N - 2; sum++)
    for (r = 0; r < fis->N; r++) {
      int i, j;
      MAP map;
      double t;

      c = sum - r;
      if (c < 0 || c >= fis->M) continue;

      i = (rayDir[0] >= 0.0) ? c : fis->M - c - 1;
      j = (rayDir[1] >= 0.0) ? r : fis->N - r - 1;

      if (compositeMap->a < 0) i = fis->M - i - 1;
      if (compositeMap->c < 0) j = fis->N - j - 1;

      /* map(x,y,z) = (compositeMap o fis->maps[j][i])(x,y,z) */
      composeMaps(compositeMap, &fis->maps[j][i], &map);

      if ((t = rayHitAux(fis, rayOrg, rayDir, &map)) > EPSILON)
	return t;
    }

  return -1.0;
}

/*
 * rayHit()
 * Method for determining ray intersection with fractal surface.
 */
static
double rayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
	      HIT_INFO *hitInfo) {
  FIS *fis = (FIS *) this->data;
  MAP identity;

  identity.a = identity.c = 1.0;
  identity.b = identity.d = 0.0;
  identity.q[0] = identity.q[1] = identity.q[2] = identity.q[3] = 0.0;
  identity.alpha = 1.0;
  
  return rayHitAux(fis, rayOrg, rayDir, &identity);
}

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

/*
 *  normal()
 *  Computes normal to surface at hit point (x,y,z).
 *  
 *  For any given map, the key recursive relationship that defines f is
 *     f(a*x + b, c*y + d) = q(x,y) + alpha*f(x,y)
 *  If we take the partial derivatives wrt x & y on each side we get
 *     fx(a*x + b, c*y + d)*a = qx(x,y) + alpha*fx(x,y)
 *     fy(a*x + b, c*y + d)*c = qy(x,y) + alpha*fy(x,y)
 *  This tells us that if we know the partial derivatives at (x,y)
 *  then we can compute the partial derivatives at (a*x + b, c*y + d).
 *  We determine a sequence of transformations that map any initial
 *  point (x0,y0) arbitrarily close to (x,y) -- the point at which
 *  we want to dertermine the normal. We just assume that the
 *  initial partial derivatives are zero:
 *     fx(x0,y0) = fy(x0,y0) = 0
 *  and then use the map sequence to approximate the target partial
 *  derivatives using the following formulas at each step:
 *     fx(a*x + b, c*y + d) = (qx(x,y) + alpha*fx(x,y))/a
 *     fy(a*x + b, c*y + d) = (qy(x,y) + alpha*fy(x,y))/c
 *  The normal vector N is then
 *     N = (-fx, -fy, 1).
 *
 *  q(x,y) = q0 + q1*x + q2*y + q3*x*y
 *  qx(x,y) = q1 + q3*y
 *  qy(x,y) = q2 + q3*x
 */
static
void normal(OBJECT *this, double hitPoint[3], HIT_INFO *info, 
	    double normal[3]) {
  FIS *fis = (FIS *) this->data;
  int k, sx[SEQUENCE_LEN], sy[SEQUENCE_LEN];
  double x,y, dfdx,dfdy;
  
  mapSequence(fis, hitPoint[0], hitPoint[1], SEQUENCE_LEN, sx, sy);

  x = fis->knots[0][0].x;
  y = fis->knots[0][0].y;
  dfdx = dfdy = 0.0;  /* initial guess at derivatives */

  for (k = 0; k < SEQUENCE_LEN; k++) {
    MAP *map = &fis->maps[sy[k]][sx[k]];
    dfdx = (map->q[1] + map->q[3]*y + map->alpha*dfdx)/map->a;
    dfdy = (map->q[2] + map->q[3]*x + map->alpha*dfdy)/map->c;
    x = map->a*x + map->b;
    y = map->c*y + map->d;
  }

  x = 1.0/sqrt(dfdx*dfdx + dfdy*dfdy + 1);
  normal[0] = -x*dfdx;
  normal[1] = -x*dfdy;
  normal[2] = x;
}

OBJECT *createBilinearFisObject(int M, int N, POINT3 **knots, double alpha) {
  FIS *fis;
  OBJECT *obj;
  int i,j;
  double dz;
  POINT3 p;

  if ((fis = (FIS *) malloc(sizeof(FIS))) == NULL ||
      (fis->knots = (POINT3 **) malloc((N+1)*sizeof(POINT3 *))) == NULL ||
      (fis->maps = (MAP **) malloc(N*sizeof(MAP *))) == NULL ||
      (obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating bilinear fis");
    exit(-1);
  }
  
  for (j = 0; j <= N; j++)
    if ((fis->knots[j] = (POINT3 *) malloc((M+1)*sizeof(POINT3))) == NULL) {
      perror("allocating row of knots for fis");
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
      fis->knots[j][i].x = knots[j][i].x;
      fis->knots[j][i].y = knots[j][i].y;
      fis->knots[j][i].z = knots[j][i].z;
    }

  /*
   * for each map i,j
   * a = (-1)^i/M,
   * b = x[i] - a*x[0],   i even
   *     x[i+1] - a*x[0], i odd
   *
   * c = (-1)^j/N;
   * d = y[j] - c*y[0],   j even
   *   = y[j+1] - c*y[0], j odd
   *
   * alpha(i,j) = alpha  (i.e. constant)
   *
   * We solve the following linear system for the coefficients of q(x,y):
   *   _                _  _  _     _                                      _
   *  | 1  x0  y0  x0*y0 || q0 |   | f(a*x0 + b, c*y0 + d) - alpha*f(x0,y0) |
   *  | 1  xM  y0  xM*y0 || q1 | = | f(a*xM + b, c*y0 + d) - alpha*f(xM,y0) |
   *  | 1  x0  yN  x0*yN || q2 |   | f(a*x0 + b, c*yN + d) - alpha*f(x0,yN) |
   *  | 1  xM  yN  xM*yN || q3 |   | f(a*xM + b, c*yN + d) - alpha*f(xM,yN) |
   *   -                -  -  -     -                                      -
   *
   * The values for f() on the right hand side of the above equation are
   * just the z-values at the knots as follows:
   *
   * f(xi,yj) = z(i,j);
   *
   *  i even & j even  (i.e. a > 0 & c > 0)
   *     f(a*x0 + b, c*y0 + d) = z(i,j)
   *     f(a*xM + b, c*y0 + d) = z(i+1,j)
   *     f(a*x0 + b, c*yN + d) = z(i,j+1)
   *     f(a*xM + b, c*yN + d) = z(i+1,j+1)
   *  i odd & j even (i.e. a < 0 & c > 0)
   *     f(a*x0 + b, c*y0 + d) = z(i+1,j)
   *     f(a*xM + b, c*y0 + d) = z(i,j)
   *     f(a*x0 + b, c*yN + d) = z(i+1,j+1)
   *     f(a*xM + b, c*yN + d) = z(i,j+1)
   *  i even & j odd (i.e. a > 0 & c < 0)
   *     f(a*x0 + b, c*y0 + d) = z(i,j+1)
   *     f(a*xM + b, c*y0 + d) = z(i+1,j+1)
   *     f(a*x0 + b, c*yN + d) = z(i,j)
   *     f(a*xM + b, c*yN + d) = z(i+1,j)
   *  i odd & j odd  (i.e. a < 0 && j < 0)
   *     f(a*x0 + b, c*y0 + d) = z(i+1,j+1)
   *     f(a*xM + b, c*y0 + d) = z(i,j+1)
   *     f(a*x0 + b, c*yN + d) = z(i+1,j)
   *     f(a*xM + b, c*yN + d) = z(i,j)
   */
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++) {
      MAP *map = &fis->maps[j][i];
      int r, ii[2], jj[2];
      double A[4][5], *Ap[4], *X[1];

      map->alpha = alpha;

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

      A[0][1] = A[2][1] = fis->knots[0][0].x;
      A[1][1] = A[3][1] = fis->knots[0][M].x;
      A[0][2] = A[1][2] = fis->knots[0][0].y;
      A[2][2] = A[3][2] = fis->knots[N][0].y;

      for (r = 0; r < 4; r++) {
	A[r][0] = 1.0;
	A[r][3] = A[r][1]*A[r][2];
	Ap[r] = A[r];	
      }

      A[0][4] = knots[jj[0]][ii[0]].z - alpha*knots[0][0].z;
      A[1][4] = knots[jj[0]][ii[1]].z - alpha*knots[0][M].z;
      A[2][4] = knots[jj[1]][ii[0]].z - alpha*knots[N][0].z;
      A[3][4] = knots[jj[1]][ii[1]].z - alpha*knots[N][M].z;

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

  for (i = 0; i < 10000; i++) {
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
   * inflate bounding-box by 5% in z-direction to be safe.
   */
  dz = 0.05*(fis->bbox.max.z - fis->bbox.min.z);
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
