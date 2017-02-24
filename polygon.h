/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef POLYGON_H
#define POLYGON_H

#include "raytrace.h"

/*
 * Creates polygon object out of 3D vertex list. For this constructor,
 * it is assumed that the polygon is convex and that the first
 * two edges are not parallel.
 */
OBJECT *createConvexPolygonObject(int numVerts, double (*verts)[3]);

#endif /* POLYGON_H */
