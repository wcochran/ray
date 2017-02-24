#include <stdlib.h>
#include <math.h>
#include "noise.h"
#include "marble.h"

/*
 * Hermite curve:
 *                             _     _
 * f(u) = [u^3 u^2 u 1] * M * | f(0)  | = U * M * G
 *                            | f(1)  |
 *                            | f'(0) |
 *                            | f'(1) |
 *                             -     -
 *  u = (t - tl)/(tr - tl)
 *  g(t) = f(u)
 *  g'(t) = f'(u)/(tr - tl)
 */
static double M[4][4] = { /* Hermite Basis Matrix */
  { 2, -2,  1,  1},
  {-3,  3, -2, -1},
  { 0,  0,  1,  0},
  { 1,  0,  0,  0}
};

/*
 * input: G = [f(0) f(1) f'(0) f'(1)]
 * output: [A0 A1 A2 A3]^T = M*G
 */
static void hermite(double G[4], double A[4]) {
  int i,j;
  for (j = 0; j < 4; j++) {
    A[j] = 0.0;
    for (i = 0; i < 4; i++)
      A[j] += M[j][i]*G[i];
  }
}

typedef struct {      /* collection of cubic functions */
  int n;              /* number of cubics */
  double *t;          /* n+1 domain partition values: t0 < t1 < ... < tn */
  double (*A)[4];     /* n sets of coefficients */
} SPLINE;

/*
 * Input:
 *   n: number of intervals;
 *   t: n+1 partitions points
 *   g: n+1 function values: g[i] = g(t[i])
 * Returns:
 *   Smooth interpolating spline
 */
static SPLINE *createSpline(int n, double t[], double g[]) {
  SPLINE *spline;
  int i;
  double dgdt;

  if ((spline = (SPLINE *) malloc(sizeof(SPLINE))) == NULL ||
      (spline->t = (double *) malloc((n+1)*sizeof(double))) == NULL ||
      (spline->A = (double (*)[4]) malloc(n*4*sizeof(double))) == NULL) {
    perror("createSpline():malloc()");
    exit(-1);
  }

  spline->n = n;
  
  for (i = 0; i <= n; i++)
    spline->t[i] = t[i];

  dgdt = (g[1] - g[0])/(t[1] - t[0]);  /* derivative @ 1st knot */

  for (i = 0; i < n; i++) {
    double G[4];
    double dtdu = (t[i+1] - t[i]);  /* dt/du */
    
    G[0] = g[i];
    G[1] = g[i+1];
    G[2] = dgdt*dtdu;  /* dg/du */

    if (i == n-1)
      dgdt = (g[i+1] - g[i])/(t[i+1] - t[i]);  /* last knot */
    else
      dgdt = (g[i+2] - g[i])/(t[i+2] - t[i]); /* interior knot */

    G[3] = dgdt*dtdu;

    hermite(G, spline->A[i]);
  }   

  return spline;
} 

/*
 * Evaluate Spline: g(t)
 */
static double evalSpline(SPLINE *spline, double t) {
  int i;
  double u;
  double *A;

  for (i = spline->n-1; i > 0 && t < spline->t[i]; i--)
    ;

  A = spline->A[i];

  u = (t - spline->t[i])/(spline->t[i+1] - spline->t[i]);

  /* A0*u^3 + A1*u^2 + A2*u + A3 = u*(u*(u*A0 + A1) + A2) + A3 */
  return u*(u*(u*A[0] + A[1]) + A[2]) + A[3];
} 

typedef struct {        /* new data for "marblized object" */
  OBJECT *old;          /* ptr to original object */
  double veinDir[3];    /* direction of marble "grain" */
  double freq;          /* top frequency for turbulence */
  SPLINE *R, *G, *B;    /* color splines */
} MARBLE_DATA;

static double rayHit(struct OBJECT *this, 
		     double rayOrg[3],
		     double rayDir[3],
		     HIT_INFO *info) {
  MARBLE_DATA *data = (MARBLE_DATA *) this->data;
  return (*data->old->rayHit)(data->old, rayOrg, rayDir, info);
}

void normal(struct OBJECT *this,
	    double hitPoint[3],
	    HIT_INFO *info,
	    double normal[3]) {
  MARBLE_DATA *data = (MARBLE_DATA *) this->data;
  (*data->old->normal)(data->old, hitPoint, info, normal);
}

static void color(struct OBJECT *this,
		  double hitPoint[3],
		  HIT_INFO *info,
		  double color[3]) {
  MARBLE_DATA *data = (MARBLE_DATA *) this->data;
  double x = 
    hitPoint[0]*data->veinDir[0] +
    hitPoint[1]*data->veinDir[1] +
    hitPoint[2]*data->veinDir[2];
  double t = 
    sin(x 
	+
	turbulence3(hitPoint[0], 
		    hitPoint[1], 
		    hitPoint[2], data->freq)
	);

  color[0] = evalSpline(data->R, t);
  color[1] = evalSpline(data->G, t);
  color[2] = evalSpline(data->B, t);

  if (color[0] < 0.0) color[0] = 0.0;    /* clamp color */
  else if (color[0] > 1.0) color[0] = 1.0;

  if (color[1] < 0.0) color[1] = 0.0;
  else if (color[1] > 1.0) color[1] = 1.0;

  if (color[2] < 0.0) color[2] = 0.0;
  else if (color[2] > 1.0) color[2] = 1.0;
}

void setColor(struct OBJECT *this,
	      double color[3]) {
  MARBLE_DATA *data = (MARBLE_DATA *) this->data;
  (*data->old->setColor)(data->old, color);
}

OBJECT *marbleObject(OBJECT *old, int n, double rgbSpline[][3], 
		     double veinDir[3], int octaves) {
  OBJECT *obj;
  MARBLE_DATA *data;
  double scale;
  int i;
  double v, *t, dt, *g;
  
  if ((obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL ||
      (data = (MARBLE_DATA *) malloc(sizeof(MARBLE_DATA))) == NULL ||
      (t = (double *) malloc(n*sizeof(double))) == NULL ||
      (g = (double *) malloc(n*sizeof(double))) == NULL) {
    perror("marbleObject():malloc()");
    exit(-1);
  }

  *obj = *old;
  data->old = old;

  scale = 1.0/sqrt(veinDir[0]*veinDir[0] + 
		   veinDir[1]*veinDir[1] +
		   veinDir[2]*veinDir[2]);
  data->veinDir[0] = scale*veinDir[0];
  data->veinDir[1] = scale*veinDir[1];
  data->veinDir[2] = scale*veinDir[2];

  data->freq = (double) (1 << octaves);

  dt = 2.0/(n-1);
  for (i = 0, v = -1.0; i < n; i++, v += dt)
    t[i] = v;

  for (i = 0; i < n; i++)
    g[i] = rgbSpline[i][0];
  data->R = createSpline(n-1, t, g);

  for (i = 0; i < n; i++)
    g[i] = rgbSpline[i][1];
  data->G = createSpline(n-1, t, g);

  for (i = 0; i < n; i++)
    g[i] = rgbSpline[i][2];
  data->B = createSpline(n-1, t, g);

  free(g);
  free(t);

  obj->data = data;
  obj->rayHit = rayHit;
  obj->normal = normal;
  obj->color = color;
  obj->setColor = setColor;

  return obj;
}
