#ifndef SUPERELLIPSOID_H
#define SUPERELLIPSOID_H

#include "raytrace.h"

/*
 * Unit superellipsoid implicit function
 *    f(x,y,z) = (|x|^n + |y|^n)^(m/n) + |z|^m
 * m,n > 0
 */
typedef struct {
  double color[3];
  double n, m;
  double size[3];
  double orientation[3][3];
  double center[3];
} SUPERELLIPSOID_DATA;

OBJECT *createSuperellipsoidObject(double n, double m, 
				   double size[3], 
				   double zdir[3], double xdir[3],
				   double center[3]);
void destroySuperellipsoidObject(OBJECT *object);

#endif /* SUPERELLIPSOID_H */
