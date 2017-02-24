
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef BEZIER3_H
#define BEZIER3_H

#include "raytrace.h"
#include "trim.h"

/*
 * Default tolerance/bound on error in estimating s,t parameter
 * values for ray/patch intersection. If we really wanted to do
 * things right we would compute this as a function of the current
 * "pixel size" in world space.
 */
#define ST_TOLERANCE 0.0000001

typedef struct {
  short tile;          /* tile uv-map? */
  double org[2];       /* uv-org */
  double u[2];         /* u-axis */
  double v[2];         /* v-axis */
  double width;        /* width of uv-map */
  double height;       /* height of uv-map */
} BEZIER3_UVMAP;

typedef struct {
  BEZIER3_UVMAP uvmap;    /* for mapping (s,t) params to (u,v) */
  TRIM_OBJECT *trimmer;   /* trimming object */
} BEZIER3_TRIM_DATA;

typedef struct {
  double st_tolerance;     /* subdivision tolerance in param. space */
  POINT3 P[4][4];          /* control points for cubic Bezier mesh */
  double color[3];         /* solid color of surface */  
  BEZIER3_TRIM_DATA *trim; /* (s,t) -> (u,v) trim info */ 
} BEZIER3_DATA;

/*
 * Creates a bicubic Bezier patch object given 4x4 mesh of
 * control points.
 */
OBJECT *createBezier3Object(POINT3 P[4][4]);
void destroyBezierObject(OBJECT *obj);

void trimBezier3(OBJECT *object,             /* bezier3 object */
		 TRIM_OBJECT *trimObject,    /* trimming object */
		 int tile,                   /* tile trimming? */
		 double origin[2],           /* u,v origin on st-plane */
		 double uaxis[2],            /* u axis on st-plane */
		 double vaxis[2],            /* v axis on st-plane */
		 double width,               /* width of trim */
		 double height);             /* height of trim */

#endif
