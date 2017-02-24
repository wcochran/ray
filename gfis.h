/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef GFIS_H
#define GFIS_H

OBJECT *createGeneralBilinearFisObject(int M, int N, 
				       POINT3 **knots, double **alpha);

OBJECT *createGridObjectFromGeneralBilinearFisObject(OBJECT *obj,
						     int W, int H);

/*
 *  f(x,y) = A x^2 + B y^2 + C xy + D x + E y + F
 *  coeffs[6] = {A, B, C, D, E, F}.
 * support origin = (orgx,orgx), size = (w, h).
 */
OBJECT *createGeneralBilinearFisObjectForQuadratic(double coeffs[6],
						   int N, /* M = N */
						   double orgx, double orgy,
						   double w, double h);

#endif /* GFIS_H */
