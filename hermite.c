#include <stdio.h>
#include <stdlib.h>
#include "func2.h"
#include "hermite.h"

typedef struct {
  double A[4][4];
} Data;

/*
 * f(x) = a0 x^3 + a1 x^2 + a2 x + a3
 *      = x(x(x*a0 + a1) + a2) + a3
 */
#define HORNER(a, x)\
  ((x)*((x)*((x)*a[0] + a[1]) + a[2]) + a[3])

/*
 * f'(x) = 3 a0 x^2 + 2 a1 x + a2
 *       = x(x*3*a0 + 2*a1) + a2
 */
#define HORNERDERIV(a, x)\
  ((x)*((x)*3*a[0] + 2*a[1]) + a[2])

/*                         _               _   _   _
 * f(u,v) = [u^3 u^2 u 1] | A00 A01 A02 A03 | | v^3 |
 *                        | A10 A11 A12 A13 | | v^2 |
 *                        | A20 A21 A22 A23 | |  v  |
 *                        | A30 A31 A32 A33 | |  1  |
 *                         -               -   -   -
 */
static double eval(struct Func2 *this, double x, double y) {
  double (*A)[4] = ((Data *) this->data)->A;
  double b[4];  
  int i;
  for (i = 0; i < 4; i++)
    b[i] = HORNER(A[i], y);
  return HORNER(b, x);
}

static double dx(struct Func2 *this, double x, double y) {
  double (*A)[4] = ((Data *) this->data)->A;
  double b[4];  
  int i;
  for (i = 0; i < 3; i++)
    b[i] = HORNER(A[i], y);
  return HORNERDERIV(b, x);
}

static double dy(struct Func2 *this, double x, double y) {
  double (*A)[4] = ((Data *) this->data)->A;
  double b[4];  
  int i;
  for (i = 0; i < 4; i++)
    b[i] = HORNERDERIV(A[i], y);
  return HORNER(b, x);
}

static double dxdy(struct Func2 *this, double x, double y) {
  double (*A)[4] = ((Data *) this->data)->A;
  double b[4];  
  int i;
  for (i = 0; i < 3; i++)
    b[i] = HORNERDERIV(A[i], y);
  return HORNERDERIV(b, x);
}

Func2 *createHermiteFunc2(void) {
  Func2 *f;
  Data *data;
  int i,j;

  if ((f = (Func2 *) malloc(sizeof(Func2))) == NULL ||
      (data = (Data *) malloc(sizeof(Data))) == NULL) {
    perror("createHermiteFunc2():malloc()");
    exit(-1);
  }

  f->eval = eval;
  f->dx = dx;
  f->dy = dy;
  f->dxdy = dxdy;
  f->data = data;

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
      data->A[j][i] = 0.0;

  return f;
}

void setHermiteFunc2(Func2 *func,
                     double f[2][2],
                     double fu[2][2], double fv[2][2], double fuv[2][2]) {
  Data *data = (Data *) func->data;
  static double M[4][4] = { /* Hermite Basis Matrix (Foley, van Dam. p 484) */
    { 2, -2,  1,  1},
    {-3,  3, -2, -1},
    { 0,  0,  1,  0},
    { 1,  0,  0,  0}
  };
  double G[4][4];
  double B[4][4];
  int i,j,k;

  /*
   * Build Hermite surface geometry matrix (Foley, van Dam. p 519) 
   * Note array indexing: f(0,1) = f[1][0]
   */
  G[0][0] = f[0][0]; G[0][1] = f[1][0];
  G[1][0] = f[0][1]; G[1][1] = f[1][1];

  G[0][2] = fv[0][0]; G[0][3] = fv[1][0];
  G[1][2] = fv[0][1]; G[1][3] = fv[1][1];

  G[2][0] = fu[0][0]; G[2][1] = fu[1][0];
  G[3][0] = fu[0][1]; G[3][1] = fu[1][1];

  G[2][2] = fuv[0][0]; G[2][3] = fuv[1][0];
  G[3][2] = fuv[0][1]; G[3][3] = fuv[1][1];

  /*
   * A = M G M^T
   *                         _               _   _   _
   * f(u,v) = [u^3 u^2 u 1] | A00 A01 A02 A03 | | v^3 |
   *                        | A10 A11 A12 A13 | | v^2 |
   *                        | A20 A21 A22 A23 | |  v  |
   *                        | A30 A31 A32 A33 | |  1  |
   *                         -               -   -   -
   *
   *          ___3    ___3
   * f(u,v) = \       \       A[i][j] u^(3-i) v^(3-j)
   *          /__j=0  /__i=0
   *
   */

  for (i = 0; i < 4; i++)          /* B = M G */
    for (j = 0; j < 4; j++) {
      B[i][j] = 0.0;
      for (k = 0; k < 4; k++)
        B[i][j] += M[i][k]*G[k][j];
    }

  for (i = 0; i < 4; i++)          /* A = M G M^T = B M^T */
    for (j = 0; j < 4; j++) {
      data->A[i][j] = 0.0;
      for (k = 0; k < 4; k++)
        data->A[i][j] += B[i][k]*M[j][k];
    }
}

/*
 * Inherited from Data.
 */
typedef struct {
  double A[4][4];  /* must be first field to be compatible with Data struct */
  double a, b;
  double c, d;
} GData;

static double geval(struct Func2 *this, double x, double y) {
  GData *data = (GData *) this->data;
  return eval(this, data->a*x + data->b, data->c*y + data->d);
}

static double gdx(struct Func2 *this, double x, double y) {
  GData *data = (GData *) this->data;
  return data->a*dx(this, data->a*x + data->b, data->c*y + data->d);
}

static double gdy(struct Func2 *this, double x, double y) {
  GData *data = (GData *) this->data;
  return data->c*dy(this, data->a*x + data->b, data->c*y + data->d);
}

static double gdxdy(struct Func2 *this, double x, double y) {
  GData *data = (GData *) this->data;
  return data->a*data->c*dxdy(this, data->a*x + data->b, data->c*y + data->d);
}

Func2 *createGeneralHermiteFunc2(double x0, double x1, double y0, double y1) {
  Func2 *g;
  GData *data;
  int i,j;

  if ((g = (Func2 *) malloc(sizeof(Func2))) == NULL ||
      (data = (GData *) malloc(sizeof(GData))) == NULL) {
    perror("createGeneralHermiteFunc2():malloc()");
    exit(-1);
  }

  g->eval = geval;
  g->dx = gdx;
  g->dy = gdy;
  g->dxdy = gdxdy;
  g->data = data;

  for (j = 0; j < 4; j++)
    for (i = 0; i < 4; i++)
      data->A[j][i] = 0.0;

  data->a = 1.0/(x1 - x0);
  data->b = -data->a*x0;
  data->c = 1.0/(y1 - y0);
  data->d = -data->c*y0;

  return g;
}

void destroyGeneralHermiteFunc2(Func2 *f) {
  free(f->data);
  free(f);
}

void setGeneralHermiteFunc2(Func2 *func,
                            double g[2][2],
                            double gx[2][2], double gy[2][2], 
                            double gxy[2][2]) {
  GData *data = (GData *) func->data;
  double inva = 1.0/data->a, invc = 1.0/data->c;
  double invac = inva*invc;
  double fu[2][2];
  double fv[2][2];
  double fuv[2][2];  

  fu[0][0] = gx[0][0]*inva;  fu[0][1] = gx[0][1]*inva;
  fu[1][0] = gx[1][0]*inva;  fu[1][1] = gx[1][1]*inva;

  fv[0][0] = gy[0][0]*invc;  fv[0][1] = gy[0][1]*invc;
  fv[1][0] = gy[1][0]*invc;  fv[1][1] = gy[1][1]*invc;

  fuv[0][0] = gxy[0][0]*invac;  fuv[0][1] = gxy[0][1]*invac;
  fuv[1][0] = gxy[1][0]*invac;  fuv[1][1] = gxy[1][1]*invac;

  setHermiteFunc2(func, g, fu, fv, fuv);
}
