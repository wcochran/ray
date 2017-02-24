#ifndef HERMITE_H
#define HERMITE_H

#include "func2.h"

/*
 *  bicubic Hermite surface: f(u,v),  [u,v] \in [0,1] x [0,1]
 *
 *  Surface specified by function value at the corners, as well
 *  as the first derivatives fu, fv, and mixed second derivatives
 *  fuv. See 11.3.1, pp 517-521 in Foley, van Dam.
 *
 *       f[1][0] ------- f[1][1]
 *          |              |
 *          |              |       
 *          |              |       
 *          |              |       
 *       f[0][0] ------- f[0][1]
 */
Func2 *createHermiteFunc2(void);
void setHermiteFunc2(Func2 *func,
                     double f[2][2],
                     double fx[2][2], double fy[2][2], double fxy[2][2]);

/*
 *  The Hermite function above is defined over the support [0,1]x[0,1].
 *  Here we want to define a Hermite function defined over
 *  the support [x0,x1]x[y0,y1]
 *
 *  g(x,y) = f((x-x0)/(x1-x0), (y - y0)/(y1 - y0))
 *         = f(a*x + b, c*y + d)
 *
 *  gx(x,y) = a*fu(a*x + b, c*y + d)
 *  gy(x,y) = c*fv(a*x + b, c*y + d)
 *  gxy(x,y) = a*c*fuv(a*x + b, c*y + d)
 */
Func2 *createGeneralHermiteFunc2(double x0, double x1, double y0, double y1);
void destroyGeneralHermiteFunc2(Func2 *f);
void setGeneralHermiteFunc2(Func2 *func,
                            double g[2][2],
                            double gx[2][2], double gy[2][2], 
                            double gxy[2][2]);

#endif /* HERMITE_H */
