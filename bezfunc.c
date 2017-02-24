#include <stdlib.h>
#include <stdio.h>
#include "bezfunc.h"

/*
 * Cubic Bezier function z = f(x,y);
 * Control mesh: (x[i], y[j], z[j][i]), 0 <= i,j <= 3
 * Assume evenly spaced mesh in x & y.
 *   dx = x[3] - x[0], dy = y[3] - y[0];
 */
typedef struct {
  double x0, y0;
  double dx, dy;
  double z[4][4];
} Data;

#define LERP(a, b, t) (((b) - (a))*(t) + (a))

static double deCasteljau1(double z[4], int n, double u) {
  int i,k;
  double q[4];
  for (i = 0; i <= n; i++)
    q[i] = z[i];
  for (k = 1; k <= n; k++)
    for (i = 0; i <= n-k; i++)
      q[i] = LERP(q[i], q[i+1], u);
  return q[0];
}

static double deCasteljau2(double z[4][4], int n, int m, double u, double v) {
  int j;
  double q[4];
  for (j = 0; j <= m; j++)
    q[j] = deCasteljau1(z[j], n, u);
  return deCasteljau1(q, m, v);
}

static double eval(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  double u = (x - data->x0)/data->dx;
  double v = (y - data->y0)/data->dy;
  return deCasteljau2(data->z, 3, 3, u, v);
}

/*
 *  df/dx = df/du * du/dx
 */
static double dx(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  int i,j;
  double z[4][4];
  double u = (x - data->x0)/data->dx;
  double v = (y - data->y0)/data->dy;
  for (j = 0; j <= 3; j++)
    for (i = 0; i <= 2; i++)
      z[j][i] = data->z[j][i+1] - data->z[j][i];
  return deCasteljau2(z, 2, 3, u, v)/data->dx;
}

/*
 *  df/dy = df/dv * dv/dy
 */
static double dy(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  int i,j;
  double z[4][4];
  double u = (x - data->x0)/data->dx;
  double v = (y - data->y0)/data->dy;
  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 2; j++)
      z[j][i] = data->z[j+1][i] - data->z[j][i];
  return deCasteljau2(z, 3, 2, u, v)/data->dy;
}

static double dxdy(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  int i,j;
  double z[4][4];
  double u = (x - data->x0)/data->dx;
  double v = (y - data->y0)/data->dy;
  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 2; j++)
      z[j][i] = data->z[j+1][i] - data->z[j][i];
  for (j = 0; j <= 2; j++)
    for (i = 0; i <= 2; i++)
      z[j][i] = z[j][i+1] - z[j][i];
  return deCasteljau2(z, 2, 2, u, v)/(data->dx*data->dy);
}

Func2 *createBezFunc2(double x0, double x3, double y0, double y3) {
  Func2 *f;
  Data *data;
  int i,j;

  if ((f = (Func2 *) malloc(sizeof(Func2))) == NULL ||
      (data = (Data *) malloc(sizeof(Data))) == NULL) {
    perror("createBezFunc()");
    exit(-1);
  }

  data->x0 = x0;
  data->dx = x3 - x0;

  data->y0 = y0;
  data->dy = y3 - y0;

  for (j = 0; j <= 3; j++)
    for (i = 0; i <= 3; i++)
      data->z[j][i] = 0.0;

  f->eval = eval;
  f->dx = dx;
  f->dy = dy;
  f->dxdy = dxdy;
  f->data = data;

  return f;
}

void destroyBezFunc2(Func2 *f) {
  free(f->data);
  free(f);
}

void setBezFunc2(Func2 *f, double z[4][4]) {
  Data *data = (Data *) f->data;
  int i,j;
  for (j = 0; j <= 3; j++)
    for (i = 0; i <= 3; i++)
      data->z[j][i] = z[j][i];
}

/*
 *  z[0][1] - z[0][0] = (d/dx) z[0][0] = gx[0][0]
 *  z[0][3] - z[0][2] = (d/dx) z[0][3] = gx[0][1]
 *  z[3][1] - z[3][0] = (d/dx) z[3][0] = gx[1][0]
 *  z[3][3] - z[3][2] = (d/dx) z[3][3] = gx[1][1]
 *
 *  z[1][0] - z[0][0] = (d/dy) z[0][0] = gy[0][0]
 *  z[3][0] - z[2][0] = (d/dy) z[3][0] = gy[1][0]
 *  z[1][3] - z[0][3] = (d/dy) z[0][3] = gy[0][1]
 *  z[3][3] - z[2][3] = (d/dy) z[3][3] = gy[1][1]
 *
 *  (z[1][1] - z[1][0]) - (z[0][1] - z[0][0]) = (d2/dxdy) z[0][0] = gxy[0][0]
 *  (z[3][1] - z[3][0]) - (z[2][1] - z[2][0]) = (d2/dxdy) z[3][0] = gxy[1][0]
 *  (z[1][3] - z[1][2]) - (z[0][3] - z[0][2]) = (d2/dxdy) z[0][3] = gxy[0][1]
 *  (z[3][3] - z[3][2]) - (z[2][3] - z[2][2]) = (d2/dxdy) z[3][3] = gxy[1][1]
 */
void setBezFunc2FromDerivs(Func2 *func,
			   double g[2][2],
			   double gx[2][2], double gy[2][2], 
			   double gxy[2][2]) {
  Data *data = (Data *) func->data;
  double z[4][4];
  double fu[2][2], fv[2][2], fuv[2][2];
  double dxdy = data->dx*data->dy;

  fu[0][0] = gx[0][0]*data->dx;
  fu[0][1] = gx[0][1]*data->dx;
  fu[1][0] = gx[1][0]*data->dx;
  fu[1][1] = gx[1][1]*data->dx;

  fv[0][0] = gy[0][0]*data->dy;
  fv[0][1] = gy[0][1]*data->dy;
  fv[1][0] = gy[1][0]*data->dy;
  fv[1][1] = gy[1][1]*data->dy;

  fuv[0][0] = gxy[0][0]*dxdy;
  fuv[0][1] = gxy[0][1]*dxdy;
  fuv[1][0] = gxy[1][0]*dxdy;
  fuv[1][1] = gxy[1][1]*dxdy;

  z[0][0] = g[0][0];
  z[0][3] = g[0][1];
  z[3][0] = g[1][0];
  z[3][3] = g[1][1];

  z[0][1] = fu[0][0] + z[0][0];
  z[0][2] = z[0][3] - fu[0][1];

  z[3][1] = fu[1][0] + z[3][0];
  z[3][2] = z[3][3] - fu[1][1];

  z[1][0] = fv[0][0] + z[0][0];
  z[2][0] = z[3][0] - fv[1][0];

  z[1][3] = fv[0][1] + z[0][3];
  z[2][3] = z[3][3] - fv[1][1];

  z[1][1] = fuv[0][0] + z[0][1] - z[0][0] + z[1][0];
  z[2][1] = z[3][1] - z[3][0] + z[2][0] - fuv[1][0];
  z[1][2] = z[1][3] - z[0][3] + z[0][2] - fuv[0][1];
  z[2][2] = fuv[1][1] + z[3][2] - z[3][3] + z[2][3];
  
  setBezFunc2(func, z);
}
