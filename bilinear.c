#include <stdio.h>
#include <stdlib.h>
#include "bilinear.h"

static void *myalloc(size_t sz) {
  void *p = malloc(sz);
  if (p == NULL) {
    perror("malloc()");
    exit(-1);
  }
  return p;
}

/*
 *             zt                
 *       *-----*----* y1         u = (x - x0)/(x1 - x0)
 *       |     |    |            v = (y - y0)/(y1 - y0)
 *       |   z |    |
 *       |-----*----|            zb = (z01 - z00)*u + z00
 *       |     |    |            zt = (z11 - z10)*u + z10
 *       *-----*----* y0         z = (zt - zb)*v + zb
 *       x0   zb   x1
 *
 *  z = [(z11 - z10)*u + z10 - (z01 - z00)*u - z00]*v + (z01 - z00)*u + z00
 *    = ((z11 - z10) - (z01 - z00))*u*v + (z10 - z00)*v + (z01 - z00)*u + z00
 *    = A*u*v + B*u + C*v + D
 *
 *    A = (z11 - z10) - (z01 - z00)
 *    B = z01 - z00
 *    C = z10 - z00
 *    D = z00
 *
 *  z = A*(a*x + b)*(c*y + d) + B*(a*x + b) + C*(c*y + d) + D
 *    = A*a*c*x*y + (B*a + A*a*d)*x + (C*c + A*b*c)*y + A*b*d + B*b + C*d + D
 *    = A*a*c*x*y + a*(B + A*d)*x + c*(C + A*b)*y + b*(A*d + B) + C*d + D
 *
 *    a = 1/(x1 - x0),  b = -x0/(x1 - x0)
 *    c = 1/(y1 - y0),  b = -y0/(y1 - y0)
 */
typedef struct {
  double A, B, C, D;  /* f(x,y) = A*x*y + B*x + C*y + D */
  double a, b;        /* u = a*x + b */
  double c, d;        /* v = c*y + d */
} Data;

static double eval(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  return x*(data->A*y + data->B) + data->C*y + data->D;
}

static double dx(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  return data->A*y + data->B;
}

static double dy(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  return data->A*x + data->C;
}

static double dxdy(Func2 *this, double x, double y) {
  Data *data = (Data *) this->data;
  return data->A;
}

Func2 *createBilinearFunc2(double x0, double x1, double y0, double y1) {
  Func2 *f = (Func2 *) myalloc(sizeof(Func2));
  Data *data = (Data *) myalloc(sizeof(Data));
  data->a = 1.0/(x1 - x0);
  data->b = -x0*data->a;
  data->c = 1.0/(y1 - y0);
  data->d = -y0*data->c;
  data->A = data->B = data->C = data->D = 0.0;
  f->eval = eval;
  f->dx = dx;
  f->dy = dy;
  f->dxdy = dxdy;
  f->data = data;
  return f;
}

void destroyBilinearFunc2(Func2 *f) {
  free(f->data);
  free(f);
}

void setBilinearFunc2(Func2 *f, double z[2][2]) {
  Data *data = (Data *) f->data;
  double A_ = z[1][1] - z[1][0] - z[0][1] + z[0][0];
  double B_ = z[0][1] - z[0][0];
  double C_ = z[1][0] - z[0][0];
  data->A = A_*data->a*data->c;
  data->B = data->a*(B_ + A_*data->d);
  data->C = data->c*(C_ + A_*data->b);
  data->D = data->b*(A_*data->d + B_) + C_*data->d + z[0][0];
}
