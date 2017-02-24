#ifndef FUNC2_H
#define FUNC2_H

/*
 * Abstract base class for height function z = f(x,y).
 * Methods
 *  eval : evaluates f(x,y)
 *  dx, dy : evaluate partial derivatives df/dx and df/dy
 *  dxdy : evaluate mixed second partial derivative df^2/dxdy
 */
typedef struct Func2 {
  double (*eval) (struct Func2 *this, double x, double y);
  double (*dx) (struct Func2 *this, double x, double y);
  double (*dy) (struct Func2 *this, double x, double y);
  double (*dxdy) (struct Func2 *this, double x, double y);
  void *data;  /* data for concrete derived class */
} Func2;

Func2 *createConstantFunc2(double c);
void setConstantFunc2(Func2 *f, double c);

void getFunc2Samples(Func2 *f,
		     int W, int H, 
		     double x0, double y0, double w, double h,
		     double *z, int rowStride);

#endif /* FUNC2_H */
