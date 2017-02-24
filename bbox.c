/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include "bbox.h"
#include "raytrace.h"

#define MIN2(a,b) (((a) < (b)) ? (a) : (b))
#define MAX2(a,b) (((a) > (b)) ? (a) : (b))

/*
 * rayHitsBoundingBox()
 * Return true if ray intersects bounding box (sides of box
 * are parallel with the principle axes). Entry and exit t-values
 * are returned through t[2] parameter.
 */
int rayHitsBoundingBox(BBOX *bbox, double rayOrg[3], double rayDir[3],
		       double t[2]) {
  double tnear[3], tfar[3];

  /*
   * Find intersection of ray with sides of the box where x is constant.
   */
  if (rayDir[0] != 0.0) {
    tnear[0] = (bbox->min.x - rayOrg[0])/rayDir[0];
    tfar[0]  = (bbox->max.x - rayOrg[0])/rayDir[0];
    if (tfar[0] < tnear[0]) {
      double tmp = tfar[0];
      tfar[0] = tnear[0];
      tnear[0] = tmp;
    }
  } else {
    tnear[0] = -INFINITY;
    tfar[0] = INFINITY;
  }

  /*
   * Find intersection of ray with sides of the box where y is constant.
   */
  if (rayDir[1] != 0.0) {
    tnear[1] = (bbox->min.y - rayOrg[1])/rayDir[1];
    tfar[1]  = (bbox->max.y - rayOrg[1])/rayDir[1];
    if (tfar[1] < tnear[1]) {
      double tmp = tfar[1];
      tfar[1] = tnear[1];
      tnear[1] = tmp;
    }
  } else {
    tnear[1] = -INFINITY;
    tfar[1] = INFINITY;
  }

  /*
   * Find entry and exit t-values computed so far.
   */
  t[0] = MAX2(tnear[0], tnear[1]);
  t[1] = MIN2(tfar[0], tfar[1]);

  /*
   * See if ray misses the sides of the box.
   * If an intersection is possible, then the pair of intervals
   * [tnear[0], tfar[0]] and [tnear[1], tfar[1]] must overlap
   * and the minimum of tfar[0] and tfar[1] must be non-negative.
   */
  if (t[0] <= t[1] && t[1] >= 0.0) {
    if (t[0] < 0.0) t[0] = 0.0;  /* clip to eye point */
  } else
    return 0;  /* no intersection possible */

  /*
   * Find intersection of ray with top and bottom of the box.
   */
  if (rayDir[2] != 0.0) {
    tnear[2] = (bbox->min.z - rayOrg[2])/rayDir[2];
    tfar[2]  = (bbox->max.z - rayOrg[2])/rayDir[2];
    if (tfar[2] < tnear[2]) {
      double tmp = tfar[2];
      tfar[2] = tnear[2];
      tnear[2] = tmp;
    }
  } else {
    tnear[2] = -INFINITY;
    tfar[2] = INFINITY;
  }

  /*
   * Update entry and exit t-values.
   */
  if (tnear[2] > t[0])
    t[0] = tnear[2];
  if (tfar[2] < t[1])
    t[1] = tfar[2];

  /*
   * If we have an intersection and the exit point is in-front of
   * the eye, then we return true (false otherwise).
   */
  return t[0] <= t[1] && t[1] >= 0.0;
}
