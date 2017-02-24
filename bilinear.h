#ifndef BILINEAR_H
#define BILINEAR_H

#include "func2.h"

Func2 *createBilinearFunc2(double x0, double x1, double y0, double y1);
void destroyBilinearFunc2(Func2 *f);
void setBilinearFunc2(Func2 *f, double z[2][2]);

#endif /* BILINEAR_H */
