/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef BBOX_H
#define BBOX_H

#include "raytrace.h"

typedef struct {         /* bounding box */
  POINT3 min, max;
} BBOX;

/*
 * rayHitsBoundingBox()
 * Return true if ray intersects bounding box (sides of box
 * are parallel with the principle axes). Entry and exit t-values
 * are returned through t[2] parameter.
 */
int rayHitsBoundingBox(BBOX *bbox, double rayOrg[3], double rayDir[3],
		       double t[2]);

#endif /* BBOX_H */
