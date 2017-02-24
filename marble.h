#ifndef MARBLE_H
#define MARBLE_H

#include "raytrace.h"

OBJECT *marbleObject(OBJECT *old, int n, double rgbSpline[][3], 
		     double veinDir[3], int octaves);

#endif /* MARBLE_H */
