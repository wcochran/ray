/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"
#include "polygon.h"

typedef struct {
  double color[3];            /* color of polygonal facet */
  double plane[4];            /* plane of polygon, plane[0..2] is normal */
  double org[3];              /* origin of local coordinate system */
  double x[3], y[3];          /* orthonormal vectors spanning plan */
  short numVerts;             /* number of vertices */
  double (*projVerts)[2];     /* 2D projected version of vertices */
} POLYGON_DATA;

static
double getRayHit(OBJECT *this, 
		 double rayOrg[3], double rayDir[3], 
		 HIT_INFO *hitInfo) {
  POLYGON_DATA *data = (POLYGON_DATA *) this->data;
  double t, numer, denom;
  double hit[3], proj[2], zcross, u[2], v[2];
  double (*projVerts)[2];
  int i, N;

  /*
   * Find intersection of ray with plane of polygon.
   */
  denom = (data->plane[0]*rayDir[0] +
	   data->plane[1]*rayDir[1] +
	   data->plane[2]*rayDir[2]);
  if (denom == 0.0) return -1.0;   /* ray parallel with plane */
  numer = -(data->plane[0]*rayOrg[0] +
	    data->plane[1]*rayOrg[1] +
	    data->plane[2]*rayOrg[2] + data->plane[3]);
  t = numer/denom;

  if (t < EPSILON) return -1.0;

  /*
   * Compute hit point on plane.
   */
  hit[0] = rayOrg[0] + t*rayDir[0];
  hit[1] = rayOrg[1] + t*rayDir[1];
  hit[2] = rayOrg[2] + t*rayDir[2];

  /*
   * Project point to plane local coordinates.
   */
  hit[0] -= data->org[0];
  hit[1] -= data->org[1];
  hit[2] -= data->org[2];
  proj[0] = hit[0]*data->x[0] + hit[1]*data->x[1] + hit[2]*data->x[2];
  proj[1] = hit[0]*data->y[0] + hit[1]*data->y[1] + hit[2]*data->y[2];

  /*
   * Make sure hit point is inside polygon via cross product testing.
   */
  N = data->numVerts;
  projVerts = data->projVerts;
  u[0] = proj[0] - projVerts[N-1][0];
  u[1] = proj[1] - projVerts[N-1][1];
  v[0] = projVerts[0][0] - projVerts[N-1][0];
  v[1] = projVerts[0][1] - projVerts[N-1][1];
  zcross = u[0]*v[1] - u[1]*v[0];
  for (i = N-2; i >= 0; i--) {
    u[0] = proj[0] - projVerts[i][0];
    u[1] = proj[1] - projVerts[i][1];
    v[0] = projVerts[i+1][0] - projVerts[i][0];
    v[1] = projVerts[i+1][1] - projVerts[i][1];
    if (zcross*(u[0]*v[1] - u[1]*v[0]) < 0.0)
      return -1.0;  /* outside polygon */
  }
      
  return t;  /* inside polygon */
}

static
void getNormal(OBJECT *this, double hit[3], 
	       HIT_INFO *hitInfo, double normal[3]) {
  POLYGON_DATA *data = (POLYGON_DATA *) this->data;
  normal[0] = data->plane[0];
  normal[1] = data->plane[1];
  normal[2] = data->plane[2];
}

static
void getColor(OBJECT *this, double hit[3], 
	      HIT_INFO *info, double color[3]) {
  POLYGON_DATA *data = (POLYGON_DATA *) this->data;
  color[0] = data->color[0];
  color[1] = data->color[1];
  color[2] = data->color[2];
}

static
void setColor(OBJECT *this, double color[3]) {
  POLYGON_DATA *data = (POLYGON_DATA *) this->data;
  data->color[0] = color[0];
  data->color[1] = color[1];
  data->color[2] = color[2];
}

OBJECT *createConvexPolygonObject(int numVerts, double (*verts)[3]) {
  POLYGON_DATA *data;
  OBJECT *object;
  double x[3], y[3], org[3], scale;
  int i;
  
  /*
   * Allocate all the memory we need for the  object.
   */
  if ((data = (POLYGON_DATA *) malloc(sizeof(POLYGON_DATA))) == NULL ||
      (data->projVerts = (double (*)[2]) 
       malloc(2*sizeof(double)*numVerts)) == NULL ||     
      (object = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("createConvexPolygonObject:malloc()");
    exit(-1);
  }

  /*
   * Find origin on poly's plane and orthonormal vectors x & y
   * that span the plane. (Note that if the first two edges of
   * the plane are parallel (or close to it) then we are screwed).
   */
  org[0] = verts[1][0];  
  org[1] = verts[1][1];  
  org[2] = verts[1][2];
  x[0] = verts[0][0] - org[0];
  x[1] = verts[0][1] - org[1];
  x[2] = verts[0][2] - org[2];
  scale = 1.0/sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
  x[0] *= scale;
  x[1] *= scale;
  x[2] *= scale;
  y[0] = verts[2][0] - org[0];
  y[1] = verts[2][1] - org[1];
  y[2] = verts[2][2] - org[2];
  scale = x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
  y[0] -= scale*x[0];
  y[1] -= scale*x[1];
  y[2] -= scale*x[2];
  scale = 1.0/sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
  y[0] *= scale;
  y[1] *= scale;
  y[2] *= scale;

  /*
   * Compute plane coeffs (ax + by + cz + d = 0). 
   * First three are the plane's normal
   * determines via cross product. We'll compute the last one
   * given a point on the plane.
   */
  data->plane[0] = x[1]*y[2] - x[2]*y[1];
  data->plane[1] = x[2]*y[0] - x[0]*y[2];
  data->plane[2] = x[0]*y[1] - x[1]*y[0];
  data->plane[3] = -(org[0]*data->plane[0] + 
		     org[1]*data->plane[1] + 
		     org[2]*data->plane[2]);

  /*
   * Store point on plane and spanning set.
   */
  data->org[0] = org[0];
  data->org[1] = org[1];
  data->org[2] = org[2];
  data->x[0] = x[0];
  data->x[1] = x[1];
  data->x[2] = x[2];
  data->y[0] = y[0];
  data->y[1] = y[1];
  data->y[2] = y[2];

  /*
   * Project vertices on plane.
   */
  data->numVerts = numVerts;
  for (i = 0; i < numVerts; i++) {
    double p[3];
    p[0] = verts[i][0] - org[0];
    p[1] = verts[i][1] - org[1];
    p[2] = verts[i][2] - org[2];
    data->projVerts[i][0] = p[0]*x[0] + p[1]*x[1] + p[2]*x[2];
    data->projVerts[i][1] = p[0]*y[0] + p[1]*y[1] + p[2]*y[2];
  }

  /*
   * Set reference to private data.
   */
  object->data = data;

  /*
   * Assign methods.
   */
  object->rayHit = getRayHit;
  object->normal = getNormal;
  object->color = getColor;
  object->setColor = setColor;

  /*
   * Assign some nice defaults for the rest.
   */
  data->color[0] = 0.9;
  data->color[1] = 0.1;
  data->color[2] = 0.1;
  object->ka = 0.25;
  object->kd = 0.40;
  object->ks = 0.35;
  object->kt = 0.05;
  object->ni = 1.33;
  object->phong = 6.0;

  return object;
}

#ifdef TEST_POLYGON

static double cubeVerts[8][3] = {
  {0.0, 0.0, 0.0},
  {1.0, 0.0, 0.0},
  {1.0, 1.0, 0.0},
  {0.0, 1.0, 0.0},
  {0.0, 0.0, 1.0},
  {1.0, 0.0, 1.0},
  {1.0, 1.0, 1.0},
  {0.0, 1.0, 1.0}
};

static int cubePolys[6][4] = {
  {0, 1, 2, 3},
  {2, 6, 7, 3},
  {0, 4, 5, 1},
  {1, 5, 6, 2},
  {0, 3, 7, 4},
  {4, 5, 6, 7}
};

static double cubeColors[6][3] = {
  {1.0, 0.0, 0.0},  /* red */
  {0.0, 1.0, 0.0},  /* green */
  {0.0, 0.0, 1.0},  /* blue */
  {1.0, 1.0, 0.0},  /* yellow */
  {0.0, 1.0, 1.0},  /* cyan */
  {1.0, 0.0, 1.0}   /* magenta */
};

static double base[4][3] = {
  {-1.0, 0.0, -1.0},
  {-1.0, 0.0,  2.0},
  { 2.0, 0.0,  2.0},
  { 2.0, 0.0, -1.0}
};
    
int main(void) {
  double eyePos[3], eyeDir[3], eyeUp[3];
  double ambient[3], background[3];
  OBJECT *objects[7];
  LIGHT light[2];
  pnm_image *image;
  FILE *f;
  int i;

  if ((f = fopen("polygons.ppm","wb")) == NULL) {
    perror("polygons.ppm");
    exit(-1);
  }

  for (i = 0; i < 6; i++) {
    int j,k;
    double poly[4][3];
    for (j = 0; j < 4; j++)
      for (k = 0; k < 3; k++)
	poly[j][k] = cubeVerts[cubePolys[i][j]][k];
    objects[i] = createConvexPolygonObject(4, poly);
    objects[i]->setColor(objects[i], cubeColors[i]);
    /*
    objects[i]->ka = 0.10;
    objects[i]->kd = 0.25;
    objects[i]->ks = 0.65;
    objects[i]->kt = 0.00;
    */
    objects[i]->ka = 0.10;
    objects[i]->kd = 0.90;
    objects[i]->ks = 0.00;
    objects[i]->kt = 0.00;
  }

  objects[6] = createConvexPolygonObject(4, base);
  objects[6]->ka = 0.00;
  objects[6]->kd = 0.00;
  objects[6]->ks = 1.00;
  objects[6]->kt = 0.00;

  eyePos[0] =  3.0; eyePos[1] = 3.0; eyePos[2] = 3.0; 
  eyeDir[0] = 0.5 - eyePos[0]; 
  eyeDir[1] = 0.5 - eyePos[1]; 
  eyeDir[2] = 0.5 - eyePos[2];
  eyeUp[0] = 0.0; eyeUp[1] = 1.0; eyeUp[2] = 0.0;

  light[0].pos[0] = -5.0; light[0].pos[1] = 5.0, light[0].pos[2] = 5.0;
  light[0].color[0] = light[0].color[1] = light[0].color[2] = 1.0;
  light[0].c[0] = 1.0; light[0].c[1] = 0.01; light[0].c[2] = 0.0;

  light[1].pos[0] = 5.0; light[1].pos[1] = 5.0, light[1].pos[2] = -5.0;
  light[1].color[0] = light[1].color[1] = light[1].color[2] = 1.0;
  light[1].c[0] = 1.0; light[1].c[1] = 0.01; light[1].c[2] = 0.0;

  ambient[0] = 0.0; ambient[1] = 0.15; ambient[2] = 0.15;
  background[0] = 0.1; background[1] = 0.1; background[2] = 0.9;

  image = raytrace(eyePos, eyeDir, eyeUp,
		   1.0, 1.0, 1.0, 
		   512, 512, 2,
		   /* 256, 256, 1, */
		   /* 256, 256, 2, */
		   /* 128, 128, 1, */
		   6,
		   ambient, background,
		   /* 1, light, */
		   2, light,
		   7, objects);

  write_pnm_image(image, f);

  return 0;
}

#endif
