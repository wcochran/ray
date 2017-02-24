#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bbox.h"
#include "superellipsoid.h"
#include "unimodalroot.h"

#define DOT(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])
#define CROSS(U,V,UxV) ((UxV)[0] = (U)[1]*(V)[2] - (U)[2]*(V)[1],\
                        (UxV)[1] = (U)[2]*(V)[0] - (U)[0]*(V)[2],\
                        (UxV)[2] = (U)[0]*(V)[1] - (U)[1]*(V)[0])

static
void normalize(double U[3]) {
  double s = DOT(U,U);
  if (s >= 0.0) {
    s = 1/sqrt(s);
    U[0] *= s;
    U[1] *= s;
    U[2] *= s;
  }
}

/*
 * Data type for state info used by evalUnitSuperellipsoid() below.
 */
typedef struct {
  double n,m;                /* bulge factors */
  double *rayOrg, *rayDir;   /* ray in coord system of unit superellipsoid */
} G_DATA;

/*
 * Function g(t) passed to unimodalRoot() solver.
 *    f(x,y,z) = (|x|^n + |y|^n)^(m/n) + |z|^m - 1
 *    g(t) = f(R(t)) where R is the ray
 */
static
double evalUnitSuperellipsoid(void *data, double t) {
  G_DATA *gdata = (G_DATA *) data;
  double P[3];
  double xn, yn, zm, s;
  P[0] = fabs(gdata->rayOrg[0] + t*gdata->rayDir[0]);
  P[1] = fabs(gdata->rayOrg[1] + t*gdata->rayDir[1]);
  P[2] = fabs(gdata->rayOrg[2] + t*gdata->rayDir[2]);
  xn = (P[0] == 0.0) ? 0.0 : pow(P[0],gdata->n);
  yn = (P[1] == 0.0) ? 0.0 : pow(P[1],gdata->n);
  zm = (P[2] == 0.0) ? 0.0 : pow(P[2],gdata->m);
  s = xn + yn;
  s = (s == 0.0) ? 0.0 : pow(s, gdata->m/gdata->n);
  return s + zm - 1;
}

/*
 * We transform a point P on the superellipsoid by
 * the matrix M yielding the point P'
 *
 *        P' = M P.
 *
 * So P is in "modeling space" (i.e., our "unit superellipsoid")
 * and P' is in "world space" (i.e., the instance the user wants
 * to render). If the ray R(t) intersects P' we have
 *
 *       R(t) = P' = M P
 *
 * Instead of transforming the inifinte number of points on our
 * solid's surface, we'll transform the ray instead.
 *
 *      M^-1 R(t) = P
 *
 * Say R(t) has origin A = (Ax,Ay,Az,1) and direction D = (Dx,Dy,Dz,0).
 * Say M = T(cx,cy,cx)*R*S(w,h,d) which defines a scale, followed
 * by a rotation, followed by a translation. The inverse of M is then
 *
 *     M^-1 = S(1/w,1/h,1/z) * R^T * T(-cx,-cy,-cz)
 *
 * Therefore, we transform A and D as follows
 *
 *     M^-1 A = S(1/w,1/h,1/z) * R^T * [Ax-cx, Ay-cy, Az-cz, 0]^T
 *    
 *     M^-1 D = S(1/w,1/h,1/z) * R^T * [Dx, Dy, Dz].
 *
 * Then we find the intersection of this transformed ray with
 * our unit superellipsoid.
 */
static
double myRayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
		HIT_INFO *hitInfo) {
  SUPERELLIPSOID_DATA *data = (SUPERELLIPSOID_DATA *) this->data;
  double xRayOrg[3], xRayDir[3];
  int i,j;
  BBOX bbox;
  int numt;
  double t[5];
  G_DATA gdata;
  double ga, gb;

  /*
   * xRayOrg = S^-1 R^-1 T^-1 rayOrg
   *         = S^-1 R^-1 (rayOrg - center)
   */
  for (i = 0; i < 3; i++) {
    double v = 0.0;
    for (j = 0; j < 3; j++)
      v += data->orientation[j][i]*(rayOrg[j] - data->center[j]);
    xRayOrg[i] = v/data->size[i];
  }

  /*
   * xRayDir = S^-1 R^-1 rayDir
   */
  for (i = 0; i < 3; i++) {
    double v = 0.0;
    for (j = 0; j < 3; j++)
      v += data->orientation[j][i]*rayDir[j];
    xRayDir[i] = v/data->size[i];
  }

  /*
   * Determine where transformated ray intersects
   * the surface's slightly inflated bounding box.
   * If the ray misses the box, then exit.
   * The entry and exit t-values are inserted into t list.
   */
#define INFLATE 0.01
  bbox.min.x = bbox.min.y = bbox.min.z = -1 - INFLATE;
  bbox.max.x = bbox.max.y = bbox.max.z = +1 + INFLATE;
  if (!rayHitsBoundingBox(&bbox, xRayOrg, xRayDir, t))
    return -1;
  numt = 2;

  if (t[1] < EPSILON)  /* surface behind ray */
    return -1;

  /*
   * Add to t-list the t-values where the transformed
   * ray intersects the x=0, y=0, and z=0 planes.
   */
  for (i = 0; i < 3; i++)
    if (xRayDir[i] != 0.0) {
      double s = -xRayOrg[i]/xRayDir[i];
      if (s > t[0] && s < t[1])
	t[numt++] = s;
    }

  /*
   * Bubble sort t list.
   */
  for (i = 0; i < numt-1; i++)
    for (j = i+1; j < numt; j++) 
      if (t[i] > t[j]) {
	double tmp = t[i];
	t[i] = t[j];
	t[j] = tmp;
      }

  /*
   * Stuff gdata with necessary state information
   * needed to evaluate evalUnitSuperellipsoid().
   */
  gdata.n = data->n;
  gdata.m = data->m;
  gdata.rayOrg = xRayOrg;
  gdata.rayDir = xRayDir;

  /*
   * Search each interval [t[i],t[i+1]] for smallest t-intersection
   * value >= EPSILON.
   */
  for (i = 0; i < numt-1 && t[i+1] < EPSILON; i++)
    ;
  gb = evalUnitSuperellipsoid(&gdata, t[i]);
  for (; i < numt-1; i++) {
    double r;
    ga = gb;
    gb = evalUnitSuperellipsoid(&gdata, t[i+1]);
    r = unimodalRoot(&gdata, evalUnitSuperellipsoid,
		     t[i], t[i+1], ga, gb, EPSILON);
    if (r >= EPSILON) {
      hitInfo->superellipsoid.h[0] = xRayOrg[0] + r*xRayDir[0];
      hitInfo->superellipsoid.h[1] = xRayOrg[1] + r*xRayDir[1];
      hitInfo->superellipsoid.h[2] = xRayOrg[2] + r*xRayDir[2];
      return r;
    }
  }

  return -1;  /* no ray intersection */
}

/*
 * f(x,y,z) = (|x|^n + |y|^n)^(m/n) + |z|^m - 1
 *
 * grad f = (fx, fy, fz)
 *
 * We need to take care with derivatives involving absolutes values.
 * Let p(x,n) = |x|^n. By the chain rule we have
 * p'(x,n) = n*|x|^(n-1) * (d/dx)|x|
 *          / n*|x|^(n-1) if x > 0
 *       = |  undef     if x = 0
 *          \ -n*|x|^(n-1) if x < 0
 *
 * fx = (m/n)*(|x|^n + |y|^n)^(m/n - 1)*p'(x,n)
 * fy = (m/n)*(|x|^n + |y|^n)^(m/n - 1)*p'(y,n)
 * fz = p'(z,m)
 *                               / |x|^(n-1) if x > 0
 * Let q(x,n) = (1/n)*p'(x,n) = |  undef   if x = 0
 *                               \ -|x|^(n-1) if x < 0
 * Then we have
 * fx = m*(|x|^n + |y|^n)^(m/n - 1)*q(x,n)
 * fy = m*(|x|^n + |y|^n)^(m/n - 1)*q(y,n)
 * fz = m*q(z,m)
 *
 * normal N = grad f / || grad f ||
 *
 * We can factor the m out and use N = (Nx,Ny,Nz) as follows
 * Nx = sign(x)*(|x|^n + |y|^n)^(m/n - 1)*|x|^(n-1)
 * Ny = sign(y)*(|x|^n + |y|^n)^(m/n - 1)*|y|^(n-1)
 * Nz = sign(z)*|z|^(m-1)
 * and then normalize.
 *
 * We stored the hit point h in modelling coordinates in ray intersection
 * routine. We use that to compute the normal N in modelling coordinates.
 * We need to transform the normal N' into world coordinates.
 * As we all know, normals are transformed via the inverse transpose of M:
 * 
 *     N' = M^-T N
 *        = (R*S(w,h,d))^-T N
 *        = (S(w,h,d)*R^T)^-1 N
 *        = R * S(1/w,1/h,1/d) N
 */
static
void myNormal (struct OBJECT *this, double hit[3], 
		   HIT_INFO *hitInfo, double normal[3]) {
  SUPERELLIPSOID_DATA *data = (SUPERELLIPSOID_DATA *) this->data;
  int i,j;
  double grad[3];
  double x,y,z, xn1,yn1,zm1, s;

  /*
   * Compute gradient of surface at hit point in modelling coordinates.
   */
  x = hitInfo->superellipsoid.h[0];
  y = hitInfo->superellipsoid.h[1];
  z = hitInfo->superellipsoid.h[2];
  xn1 = (x == 0.0) ? 0.0 : pow(fabs(x),data->n-1);
  yn1 = (y == 0.0) ? 0.0 : pow(fabs(y),data->n-1);
  zm1 = (z == 0.0) ? 0.0 : pow(fabs(z),data->m-1);
  s = fabs(x)*xn1 + fabs(y)*yn1;
  s = (s == 0.0) ? 0.0 : pow(s, data->m/data->n - 1);
  grad[0] = s*xn1;
  grad[1] = s*yn1;
  grad[2] = zm1;
  if (x < 0) grad[0] = -grad[0];
  if (y < 0) grad[1] = -grad[1];
  if (z < 0) grad[2] = -grad[2];

  /*
   * Transform gradient back into world coordinates
   * using the M^-T (ignoring translation) which is 
   * (R S)^-T = (S^-1 R^-1)^T = R S^-1
   */
  for (i = 0; i < 3; i++)
    grad[i] /= data->size[i];
  for (i = 0; i < 3; i++) {
    normal[i] = 0.0;
    for (j = 0; j < 3; j++)
      normal[i] += data->orientation[i][j]*grad[j];
  }

  /*
   * Normalize surface normal.
   */
  normalize(normal);
}

static
void myColor (struct OBJECT *this, double hit[3], HIT_INFO *hitInfo,
		  double color[3]) {
  SUPERELLIPSOID_DATA *data = (SUPERELLIPSOID_DATA *) this->data;
  color[0] = data->color[0];
  color[1] = data->color[1];
  color[2] = data->color[2];
}

static
void mySetColor(OBJECT *this, double color[3]) {
  SUPERELLIPSOID_DATA *data = (SUPERELLIPSOID_DATA *) this->data;
  data->color[0] = color[0];
  data->color[1] = color[1];
  data->color[2] = color[2];
}

OBJECT *createSuperellipsoidObject(double n, double m, 
				   double size[3], 
				   double zdir[3], double xdir[3],
				   double center[3]) {
  SUPERELLIPSOID_DATA *data;
  OBJECT *object;

  if ((data = (SUPERELLIPSOID_DATA *)
       malloc(sizeof(SUPERELLIPSOID_DATA))) == NULL ||
      (object = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating Superellipsoid object");
    exit(-1);
  }
      
  /*
   * Set private data.
   */
  data->n = n;
  data->m = m;

  data->size[0] = size[0];
  data->size[1] = size[1];
  data->size[2] = size[2];

  data->center[0] = center[0];
  data->center[1] = center[1];
  data->center[2] = center[2];

  /*
   * Create orientation rotation matrix from
   * orthonormal vectors forming basis for the shape.
   */
  CROSS(zdir, xdir, data->orientation[1]);
  CROSS(data->orientation[1], zdir, data->orientation[0]);
  data->orientation[2][0] = zdir[0];
  data->orientation[2][1] = zdir[1];
  data->orientation[2][2] = zdir[2];
  normalize(data->orientation[0]);
  normalize(data->orientation[1]);
  normalize(data->orientation[2]);

  /*
   * Assign appropriate private data and methods
   */
  object->data = data;
  object->rayHit = myRayHit;
  object->normal = myNormal;
  object->color = myColor;
  object->setColor = mySetColor;

  /*
   * Assign some nice defaults for the rest.
   */
  data->color[0] = 0.2;
  data->color[1] = 0.8;
  data->color[2] = 0.2;
  object->ka = 0.15;
  object->kd = 0.40;
  object->ks = 0.45;
  object->kt = 0.55;
  object->ni = 1.52;
  object->phong = 6.0;

  return object;
}

void destroySuperellipsoidObject(OBJECT *object) {
  free(object->data);
  free(object);
}

