
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef PLANE_H
#define PLANE_H

#include "pnmio.h"
#include "raytrace.h"

typedef struct PLANE_IMAGE_DATA {
  short tile;          /* tile image? */
  short interpolate;   /* use bilinear interpolation for color? */
  double org[3];       /* origin of image */
  double x[3];         /* image vertical axis */
  double y[3];         /* image horizontal axis */
  double width;        /* width of image */
  double height;       /* height of image */
  pnm_image *image;    /* mapped image */
} PLANE_IMAGE_DATA;

typedef struct PLANE_DATA {
  double color[3];
  double normal[3];         /* equation of half-space dividing plane */
  double d;                 /* N[0]*x + N[1]*y + N[2]*z + d = 0 */
  PLANE_IMAGE_DATA *image;  /* image mapped to plane data */
} PLANE_DATA;

/*
 * Equation of plane:
 *    coeff[0]*x + coeff[1]*y + coeff[2]*z + coeff[3] = 0.
 */
OBJECT *createPlaneObject(double coeff[4]);
OBJECT *createPlaneObjectFromNormalandPoint(double normal[3], double point[3]);
void destroyPlaneObject(OBJECT *object);
void mapImageToPlane(OBJECT *object,               /* plane object */
		     int tile,                     /* tile image? */
		     int interpolate,              /* use bilinear interp? */
		     double org[3],                /* image origin */
		     double yaxis[3],              /* image vertical axis */
		     double width, double height,  /* size of mapped image */
		     pnm_image *image);            /* the PPM image */

#endif

