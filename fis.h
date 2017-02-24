/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef FIS_H
#define FIS_H

#include "raytrace.h"

OBJECT *createBilinearFisObject(int M, int N, POINT3 **knots, double alpha);

#endif /* FIS_H */
