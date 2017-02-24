
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "raytrace.h"
#include "plane.h"

static 
double planeRayHit (OBJECT *this, double rayOrg[3], double rayDir[3], 
		    HIT_INFO *hitInfo) {
  PLANE_DATA *data = (PLANE_DATA *) this->data;
  double t, numer, denom = 
    data->normal[0]*rayDir[0] +
    data->normal[1]*rayDir[1] +
    data->normal[2]*rayDir[2];
  if (denom == 0.0) return -1.0;  /* ray parallel to plane */
  numer = -(data->normal[0]*rayOrg[0] +
	    data->normal[1]*rayOrg[1] +
	    data->normal[2]*rayOrg[2] +
	    data->d);
  t = numer/denom;
  if (t < EPSILON) return -1.0;
  return t;
}

static
void planeNormal (struct OBJECT *this, double hit[3], 
		  HIT_INFO *hitInfo, double normal[3]) {
  PLANE_DATA *data = (PLANE_DATA *) this->data;
  normal[0] = data->normal[0];
  normal[1] = data->normal[1];
  normal[2] = data->normal[2];
}

static
void planeColor (struct OBJECT *this, double hit[3], 
		 HIT_INFO *info, double color[3]) {
  PLANE_DATA *data = (PLANE_DATA *) this->data;
  color[0] = data->color[0];
  color[1] = data->color[1];
  color[2] = data->color[2];
}

static
void setPlaneColor(OBJECT *this, double color[3]) {
  PLANE_DATA *data = (PLANE_DATA *) this->data;
  data->color[0] = color[0];
  data->color[1] = color[1];
  data->color[2] = color[2];
}

OBJECT *createPlaneObject(double coeff[4]) {
  PLANE_DATA *data;
  OBJECT *object;
  double mag, scale;
  
  /*
   * Allocate memory for private plane data and object.
   */
  if ((data = (PLANE_DATA *) malloc(sizeof(PLANE_DATA))) == NULL ||
      (object = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating plane object");
    exit(-1);
  }

  /*
   * Set private specific data.
   */
  mag = sqrt(coeff[0]*coeff[0] + coeff[1]*coeff[1] + coeff[2]*coeff[2]);
  scale = 1.0/mag;
  data->normal[0] = scale*coeff[0];
  data->normal[1] = scale*coeff[1];
  data->normal[2] = scale*coeff[2];
  data->d = scale*coeff[3];
  data->image = NULL;
  
  /*
   * Assign appropriate private data and methods
   */
  object->data = data;
  object->rayHit = planeRayHit;
  object->normal = planeNormal;
  object->color = planeColor;
  object->setColor = setPlaneColor;

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

OBJECT *createPlaneObjectFromNormalandPoint(double normal[3], 
					    double point[3]) {
  double coeff[4];
  coeff[0] = normal[0];
  coeff[1] = normal[1];
  coeff[2] = normal[2];
  coeff[3] = -(point[0]*normal[0] + point[1]*normal[1] + point[2]*normal[2]);
  return createPlaneObject(coeff);
}

void destroyPlaneObject(OBJECT *object) {
  if (((PLANE_DATA *) (object->data))->image != NULL)
    free(((PLANE_DATA *) (object->data))->image);
  free(object->data);
  free(object);
}

static
void planeImageColor (struct OBJECT *this, double hit[3], 
		      HIT_INFO *info, double color[3]) {
  PLANE_DATA *pdata = (PLANE_DATA *) this->data;
  PLANE_IMAGE_DATA *data = pdata->image;
  pnm_image *image = data->image;
  double x,y, p[3], scale;  
  int i,j;

  /*
   * Adjust hit point to the image's origin.
   */
  p[0] = hit[0] - data->org[0];
  p[1] = hit[1] - data->org[1];
  p[2] = hit[2] - data->org[2];

  /*
   * Find the hit point in the image's xy-coordinate system
   * via projection onto the proper coords.
   */
  x = p[0]*data->x[0] + p[1]*data->x[1] + p[2]*data->x[2];
  y = p[0]*data->y[0] + p[1]*data->y[1] + p[2]*data->y[2];

  /*
   * If we are tiling the image and the hit point is outside
   * the image region we'll just return the plane's normal color.
   */
  if (!data->tile && (x < 0.0 || y < 0.0 || 
		      x > data->width || y > data->height)) {
    color[0] = pdata->color[0];
    color[1] = pdata->color[1];
    color[2] = pdata->color[2];
    return;
  }

  /*
   * Scale x,y to unit box.
   */
  x /= data->width;
  y /= data->height;

  /*
   * Subtract off any whole terms.
   */
  x -= (double) ((int) x);
  y -= (double) ((int) y);

  /*
   * Wrap negative terms (negative floats truncate towards zero).
   */
  if (x < 0.0) x += 1.0;
  if (y < 0.0) y += 1.0;

  /*
   * Flip y-axis.
   */
  y = 1.0 - y;

  /*
   * Find corner pixel address.
   */
  x *= PNM_NC(image) - 1;
  y *= PNM_NR(image) - 1;
  j = (int) x;
  i = (int) y;

  /*
   * Compute component normalization factor.
   */
  scale = 1.0/PNM_MAXVAL(image);

  if (data->interpolate) { /* use bilinear interpolation */
    double w[4];

    /*
     * Determine coordinate distance from cornel pixel.
     */
    x -= j;
    y -= i;

    w[0] = (1.0 - x)*(1.0 - y);  /* weights */
    w[1] = x*(1.0 - y);
    w[2] = (1.0 - x)*y;
    w[3] = x*y;

    color[0] = scale*(w[0]*PPM_PIXEL_R(image,i,j) +
		      w[1]*PPM_PIXEL_R(image,i,j+1) +
		      w[2]*PPM_PIXEL_R(image,i+1,j) +
		      w[3]*PPM_PIXEL_R(image,i+1,j+1));
    color[1] = scale*(w[0]*PPM_PIXEL_G(image,i,j) +
		      w[1]*PPM_PIXEL_G(image,i,j+1) +
		      w[2]*PPM_PIXEL_G(image,i+1,j) +
		      w[3]*PPM_PIXEL_G(image,i+1,j+1));
    color[2] = scale*(w[0]*PPM_PIXEL_B(image,i,j) +
		      w[1]*PPM_PIXEL_B(image,i,j+1) +
		      w[2]*PPM_PIXEL_B(image,i+1,j) +
		      w[3]*PPM_PIXEL_B(image,i+1,j+1));
  } else {  /* use closest pixel color */
    color[0] = scale*(PPM_PIXEL_R(image,i,j));
    color[1] = scale*(PPM_PIXEL_G(image,i,j));
    color[2] = scale*(PPM_PIXEL_B(image,i,j));
  }
}

void mapImageToPlane(OBJECT *object,               /* plane object */
		     int tile,                     /* tile image? */
		     int interpolate,              /* use bilinear interp? */
		     double org[3],                /* image origin */
		     double yaxis[3],              /* image vertical axis */
		     double width, double height,  /* size of mapped image */
		     pnm_image *image)             /* the PPM image */
{
  PLANE_DATA *data = (PLANE_DATA *) object->data;
  PLANE_IMAGE_DATA *idata;
  int i;
  double mag, v, p[3];

  /*
   * Allocate memory for image info if not allocated already.
   */
  if (data->image == NULL)
    if ((data->image = (PLANE_IMAGE_DATA *) malloc(sizeof(PLANE_IMAGE_DATA)))
	== NULL) {
      perror("mapImageToPlane:malloc");
      exit(-1);
    }

  idata = (PLANE_IMAGE_DATA *) data->image;

  /*
   * Copy local parameters.
   */
  idata->tile = tile; 
  idata->interpolate = interpolate;
  idata->width = width;
  idata->height = height;
  idata->image = image;

  /*
   * Find a point on the plane.
   */
  i = 0;
  mag = fabs(data->normal[0]);
  if ((v = fabs(data->normal[1])) > mag) {
    i = 1;
    mag = v;
  }
  if ((v = fabs(data->normal[2])) > mag)
    i = 2;
  switch(i) {
  case 0:
    p[1] = p[2] = 1.0;
    p[0] = -(data->d + data->normal[1] + data->normal[2])/data->normal[0];
    break;
  case 1:
    p[0] = p[2] = 1.0;
    p[1] = -(data->d + data->normal[0] + data->normal[2])/data->normal[1];
    break;
  case 2:
    p[0] = p[1] = 1.0;
    p[2] = -(data->d + data->normal[0] + data->normal[1])/data->normal[2];
  }

  /*
   * Project origin onto plane.
   */
  p[0] = org[0] - p[0];
  p[1] = org[1] - p[1];
  p[2] = org[2] - p[2];
  mag = p[0]*data->normal[0] + p[1]*data->normal[1] + p[2]*data->normal[2];
  idata->org[0] = org[0] - mag*data->normal[0];
  idata->org[1] = org[1] - mag*data->normal[1];
  idata->org[2] = org[2] - mag*data->normal[2];

  /*
   * Project y-axis onto plane and normalize.
   */
  mag = 
    yaxis[0]*data->normal[0] + 
    yaxis[1]*data->normal[1] + 
    yaxis[2]*data->normal[2];
  idata->y[0] = yaxis[0] - mag*data->normal[0];
  idata->y[1] = yaxis[1] - mag*data->normal[1];
  idata->y[2] = yaxis[2] - mag*data->normal[2];
  mag = sqrt(idata->y[0]*idata->y[0] + 
	     idata->y[1]*idata->y[1] + 
	     idata->y[2]*idata->y[2]);
  if (mag > 0.0) {
    mag = 1.0/mag;
    idata->y[0] *= mag;
    idata->y[1] *= mag;
    idata->y[2] *= mag;
  }

  /*
   * Compute image x-axis by cross product yaxis x normal.
   */
  idata->x[0] = idata->y[1]*data->normal[2] - idata->y[2]*data->normal[1];
  idata->x[1] = idata->y[2]*data->normal[0] - idata->y[0]*data->normal[2];
  idata->x[2] = idata->y[0]*data->normal[1] - idata->y[1]*data->normal[0];

  /*
   * Change color method.
   */
  object->color = planeImageColor;
}


/* #define TESTPLANE */
#ifdef TESTPLANE

#include "sphere.h"

int main(void) {
  double center[3], plane[4], eyePos[3], eyeDir[3], eyeUp[3];
  double ambient[3], background[3];
  double org[3], yaxis[3];
  double pole[3], equator[3];
  double uvorg[2];    
  OBJECT *objects[6];
  LIGHT light[2];
  pnm_image *image, *planeMap, *sphereMap;
  FILE *f;

  if ((f = fopen("planes.ppm","wb")) == NULL) {
    perror("planes.ppm");
    exit(-1);
  }

  light[0].pos[0] = -5.0; light[0].pos[1] = 5.0, light[0].pos[2] = 5.0;
  light[0].color[0] = light[0].color[1] = light[0].color[2] = 1.0;
  light[0].c[0] = 1.0; light[0].c[1] = 0.01; light[0].c[2] = 0.0;

  light[1].pos[0] = 5.0; light[1].pos[1] = 5.0, light[1].pos[2] = 5.0;
  light[1].color[0] = light[1].color[1] = light[1].color[2] = 1.0;
  light[1].c[0] = 1.0; light[1].c[1] = 0.01; light[1].c[2] = 0.0;
  
  /*
  eyePos[0] =  4.0; eyePos[1] = 0.0; eyePos[2] = 0.0; 
  eyeDir[0] = -1.0; eyeDir[1] = 0.0; eyeDir[2] = 0.0;
  eyeUp[0]  =  0.0; eyeUp[1]  = 1.0;  eyeUp[2] = 0.0;
  */
  eyePos[0] =  4.0; eyePos[1] = 1.2; eyePos[2] = 0.0; 
  eyeDir[0] = -eyePos[0]; eyeDir[1] = -eyePos[1]; eyeDir[2] = -eyePos[2];;
  eyeUp[0]  =  0.0; eyeUp[1]  = 1.0;  eyeUp[2] = 0.0;

  ambient[0] = 0.15; ambient[1] = 0.15; ambient[2] = 0.0;
  background[0] = 0.1; background[1] = 0.1; background[2] = 0.9;

  /* main sphere */
  center[0] = center[1] = center[2] = 0.0;
  objects[0] = createSphereObject(center, 1.0);
  /*
  objects[0]->ka = 0.00;
  objects[0]->kd = 0.00;
  objects[0]->ks = 0.25;
  objects[0]->kt = 1.00;
  objects[0]->ni = 1.52;
  */
  /*
  objects[0]->ka = 0.00;
  objects[0]->kd = 0.00;
  objects[0]->ks = 0.75;
  objects[0]->kt = 0.95;
  objects[0]->ni = 1.52;
  */
  /*
  objects[0]->ka = 0.00;
  objects[0]->kd = 0.00;
  objects[0]->ks = 1.00;
  objects[0]->kt = 0.00;
  objects[0]->ni = 1.52;
  */
  objects[0]->ka = 0.00;
  objects[0]->kd = 0.25;
  objects[0]->ks = 0.70;
  objects[0]->kt = 0.00;
  objects[0]->ni = 1.52;
  /*
  objects[0]->ka = 0.00;
  objects[0]->kd = 1.00;
  objects[0]->ks = 0.00;
  objects[0]->kt = 0.00;
  objects[0]->ni = 1.52;
  */

  /* bottom plane */
  plane[0] = 0.0; plane[1] = 1.0; plane[2] = 0.0; plane[3] = 1.0001;
  objects[1] = createPlaneObject(plane);
  objects[1]->kt = 0.0;

  /* top plane */
  plane[0] = 0.0; plane[1] = -1.0; plane[2] = 0.0; plane[3] = -1.5;
  objects[2] = createPlaneObject(plane);
  objects[2]->kt = 0.0;

  /* back plane */
  plane[0] = 1.0; plane[1] = 0.0; plane[2] = 0.0; plane[3] = 1.5;
  objects[3] = createPlaneObject(plane);
  objects[3]->kt = 0.0;

  /*
  if ((planeMap = read_pnm_image_from_file("spiral.ppm")) == NULL) {
    perror("spiral.ppm");
    exit(-1);
  }

  org[0] = 3.742; org[1] = -1.0; org[2] = 3.0;
  yaxis[0] = -1.0; yaxis[1] = 0.0; yaxis[2] = 0.0;
  mapImageToPlane(objects[1], 1, 1, org, yaxis, 7.485, 6.0, planeMap);
  */

  if ((planeMap = read_pnm_image_from_file("lj.ppm")) == NULL) {
    perror("lj.ppm");
    exit(-1);
  }

  org[0] = .9439; org[1] = -1.0; org[2] = 0.5;
  yaxis[0] = -1.0; yaxis[1] = 0.0; yaxis[2] = 0.0;
  mapImageToPlane(objects[1], 1, 1, org, yaxis, 1.8878, 1.0, planeMap);

  /*
  if ((sphereMap = read_pnm_image_from_file("ssc.ppm")) == NULL) {
    perror("ssc.ppm");
    exit(-1);
  }
  */

  /*
  if ((sphereMap = read_pnm_image_from_file("check.ppm")) == NULL) {
    perror("check.ppm");
    exit(-1);
  }
  */

  if ((sphereMap = read_pnm_image_from_file("logo.ppm")) == NULL) {
    perror("logo.ppm");
    exit(-1);
  }

  /*
  pole[0] = 0.0; pole[1] = 1.0; pole[2] = 0.0;
  equator[0] = 1.0; equator[1] = 0.0; equator[2] = 0.0;
  uvorg[0] = 0.0; uvorg[1] = 0.0;
  */
  pole[0] = 0.57735027; pole[1] = 0.57735027; pole[2] = -0.57735027;
  equator[0] = 0.57735027; equator[1] = 0.0; equator[2] = 0.57735027;
  uvorg[0] = 0.0; uvorg[1] = 0.0;
  /*
  mapImageToSphere(objects[0], 0, 0, pole, equator, uvorg, 
                   0.8, 0.375, sphereMap);
  mapImageToSphere(objects[0], 1, 0, pole, equator, uvorg, 
		   0.25, 0.25, sphereMap);
  */
  mapImageToSphere(objects[0], 1, 1, pole, equator, uvorg, 
		   0.125, 0.125, sphereMap);

  image = raytrace(eyePos, eyeDir, eyeUp,
		   1.0, 1.0, 1.0, 
		   512, 512, 2,
		   /* 256, 256, 1, */
		   /* 128, 128, 1, */
		   12,
		   ambient, background,
		   /* 1, light, */
		   2, light,
		   2, objects);

  write_pnm_image(image, f);

  return 0;
}

#endif  /* TEST_PLANE */
