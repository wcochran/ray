#ifndef BEZFUNC_H
#define BEZFUNC_H

#include "func2.h"

Func2 *createBezFunc2(double x0, double x3, double y0, double y3);
void destroyBezFunc2(Func2 *f);

void setBezFunc2(Func2 *f, double z[4][4]);
void setBezFunc2FromDerivs(Func2 *func,
			   double g[2][2],
			   double gx[2][2], double gy[2][2], 
			   double gxy[2][2]);

#endif /* BEZFUNC_H */
