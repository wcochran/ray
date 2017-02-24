
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef SPHERE_H
#define SPHERE_H

#include "pnmio.h"
#include "raytrace.h"
#include "trim.h"

typedef struct {
  double pole[3];        /* points from center to north pole */
  double equator[3];     /* points from center to equator ref. point */
  double pxe[3];         /* cached cross product pole x equator */
} SPHERE_UVMAP;

typedef struct {
  short tile;            /* tile image? */
  short interpolate;     /* use bilinear interpolation for color? */
  SPHERE_UVMAP uvmap;    /* info for mapping sphere point to uv plane */
  double org[2];         /* image origin in (u,v) coords */ 
  double width;          /* width (in 0 <= u <= 1) coords */
  double height;         /* height (in 0 <= v <= 1) coords */
  pnm_image *image;      /* mapped image */
} SPHERE_IMAGE_DATA;

typedef struct {         
  SPHERE_UVMAP uvmap;      /* info for mapping sphere point to uv plane */
  TRIM_OBJECT *trimmer;    /* trim object */
} SPHERE_TRIM_DATA;

typedef struct {     /* occupies private data element of object */
  double color[3];          /* solid color of sphere */
  double center[3];         /* center of sphere */
  double rad;               /* radius of sphere */
  SPHERE_IMAGE_DATA *image; /* image mapping info */
  SPHERE_TRIM_DATA *trim;   /* uv-trim shape info */
  PROCEDURAL_TEXTURE *tex;  /* procedural texture data */
} SPHERE_DATA;

OBJECT *createSphereObject(double center[3], double radius);
void destroySphereObject(OBJECT *object);

void mapImageToSphere(OBJECT *object,              /* sphere object */
		      int tile,                    /* tile image? */
		      int interpolate,             /* use bilinear interp? */
		      double pole[3],              /* center to north pole */
		      double equator[3],           /* center to equator pt */
		      double org[2],               /* image org in uv coords */
		      double width, double height, /* image uv size */
		      pnm_image *image);           /* PPM image */

void proceduralTextureSphere(OBJECT *object,
			     PROCEDURAL_TEXTURE *tex);

void trimSphere(OBJECT *object,              /* sphere object */
		double pole[3],              /* center to north pole */
		double equator[3],           /* center to equator pt */
		TRIM_OBJECT *trimmer);       /* object trimmer */

#endif
