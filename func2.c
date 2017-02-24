#include <stdio.h>
#include <stdlib.h>
#include "func2.h"

typedef struct {
  double val;
} Data;

static double eval(Func2 *this, double x, double y) {
  return ((Data *) this->data)->val;
}

static double dx(Func2 *this, double x, double y) {return 0.0;}
static double dy(Func2 *this, double x, double y) {return 0.0;}
static double dxdy(Func2 *this, double x, double y) {return 0.0;}

Func2 *createConstantFunc2(double c) {
  Func2 *f;
  Data *data;

  if ((f = (Func2 *) malloc(sizeof(Func2))) == NULL ||
      (data = (Data *) malloc(sizeof(Data))) == NULL) {
    perror("malloc:createConstantFunc2()");
    exit(-1);
  }

  f->eval = eval;
  f->dx = dx;
  f->dy = dy;
  f->dxdy = dxdy;
  data->val = c;
  f->data = data;

  return f;
}

void setConstantFunc2(Func2 *f, double c) {
  ((Data *) f->data)->val = c;
}

void getFunc2Samples(Func2 *f,
		     int W, int H, 
		     double x0, double y0, double w, double h,
		     double *z, int rowStride) {
  double x,y, dx,dy;
  int i,j;

  dx = w/(W-1);
  dy = h/(H-1);

  for (j = 0, y = y0; j < H; j++, y += dy)
    for (i = 0, x = x0; i < W; i++, x += dx)
      z[j*rowStride + i] = (*f->eval)(f, x, y);
}
