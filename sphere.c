
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "raytrace.h"
#include "sphere.h"
#include "trim.h"

static
int quadratic(double B, double C, double roots[2]) {
  double det = B*B - 4*C;
  if (det < 0.0) return 0;  /* imaginary */
  det = sqrt(det);
  roots[0] = 0.5*(-B + det);
  roots[1] = 0.5*(-B - det);
  return 1;
}

static
double sphereRayHit (OBJECT *this, double rayOrg[3], double rayDir[3],
		     HIT_INFO *hitInfo) {
  SPHERE_DATA *data = (SPHERE_DATA *) this->data;
  double B = 2*(rayDir[0]*(rayOrg[0] - data->center[0]) + 
		rayDir[1]*(rayOrg[1] - data->center[1]) + 
		rayDir[2]*(rayOrg[2] - data->center[2]));
  double dx = rayOrg[0] - data->center[0];
  double dy = rayOrg[1] - data->center[1];
  double dz = rayOrg[2] - data->center[2];
  double C = dx*dx + dy*dy + dz*dz - data->rad*data->rad;
  double roots[2];
  if (quadratic(B,C,roots)) {
    double t;
    if (roots[0] < EPSILON && roots[1] < EPSILON) return -1.0; /* behind eye */
    else if (roots[0] < EPSILON) t = roots[1];
    else if (roots[1] < EPSILON) t = roots[0];
    else t = (roots[0] < roots[1]) ? roots[0] : roots[1];
    return t;
  }
  return -1.0;  /* ray does not intersect sphere */
}

static
void sphereNormal (struct OBJECT *this, double hit[3], 
		   HIT_INFO *hitInfo, double normal[3]) {
  SPHERE_DATA *data = (SPHERE_DATA *) this->data;
  double invRad = 1.0/data->rad;
  normal[0] = invRad*(hit[0] - data->center[0]);
  normal[1] = invRad*(hit[1] - data->center[1]);
  normal[2] = invRad*(hit[2] - data->center[2]);
}

static
void sphereColor (struct OBJECT *this, double hit[3], HIT_INFO *hitInfo,
		  double color[3]) {
  SPHERE_DATA *data = (SPHERE_DATA *) this->data;
  color[0] = data->color[0];
  color[1] = data->color[1];
  color[2] = data->color[2];
}

static
void setSphereColor(OBJECT *this, double color[3]) {
  SPHERE_DATA *data = (SPHERE_DATA *) this->data;
  data->color[0] = color[0];
  data->color[1] = color[1];
  data->color[2] = color[2];
}

#ifndef M_1_PI
#define M_1_PI      0.31830988618379067154      /* 1/pi */
#endif
#ifndef M_PI
#define M_PI        3.14159265358979323846      /* pi */
#endif

static
void sphereUV (SPHERE_DATA *sdata, SPHERE_UVMAP *uvmap,
	       double hit[3], double uv[2]) {
  double mag;
  double normal[3];
  double phi, u,v;
  
  /*
   * Compute sphere normal at hit point.
   */
  mag = 1.0/sdata->rad;
  normal[0] = mag*(hit[0] - sdata->center[0]);
  normal[1] = mag*(hit[1] - sdata->center[1]);
  normal[2] = mag*(hit[2] - sdata->center[2]);
  
  /*
   * Calculate the latitudinal parameter v.
   */
  mag = -(normal[0]*uvmap->pole[0] +
	  normal[1]*uvmap->pole[1] +
	  normal[2]*uvmap->pole[2]);
  phi = acos(mag);
  v = phi*M_1_PI;

  /*
   * Calculate the longitudinal parameter u.
   */
  if (v <= 0.0 || v >= 1.0)
    u = 0.0;
  else {  
    double theta;
    mag = 
      normal[0]*uvmap->equator[0] +
      normal[1]*uvmap->equator[1] +
      normal[2]*uvmap->equator[2];
    mag /= sin(phi);
    theta = (fabs(mag) >= 1.0) ? 0.0 : acos(mag)*(1.0/(2*M_PI));
    mag = 
      uvmap->pxe[0]*normal[0] + 
      uvmap->pxe[1]*normal[1] + 
      uvmap->pxe[2]*normal[2];
    u = (mag > 0.0) ? theta : 1 - theta;
  }

  uv[0] = u;
  uv[1] = v;
}

static
double sphereTrimmedRayHit (OBJECT *this, double rayOrg[3], double rayDir[3],
			    HIT_INFO *hitInfo) {
  SPHERE_DATA *data = (SPHERE_DATA *) this->data;
  double B = 2*(rayDir[0]*(rayOrg[0] - data->center[0]) + 
		rayDir[1]*(rayOrg[1] - data->center[1]) + 
		rayDir[2]*(rayOrg[2] - data->center[2]));
  double dx = rayOrg[0] - data->center[0];
  double dy = rayOrg[1] - data->center[1];
  double dz = rayOrg[2] - data->center[2];
  double C = dx*dx + dy*dy + dz*dz - data->rad*data->rad;
  double roots[2];
  if (quadratic(B,C,roots)) {
    if (roots[0] < EPSILON && roots[1] < EPSILON) 
      return -1.0; /* behind eye */
    else {
      SPHERE_TRIM_DATA *trimData = data->trim;
      TRIM_OBJECT *trimmer = trimData->trimmer;
      double hit[3], uv[2];
      if (roots[0] < EPSILON) { 
	double t = roots[1];
	hit[0] = rayOrg[0] + rayDir[0]*t;
	hit[1] = rayOrg[1] + rayDir[1]*t;
	hit[2] = rayOrg[2] + rayDir[2]*t;
	sphereUV (data, &trimData->uvmap, hit, uv);
	if (trimmer->trim(trimmer, uv[0], uv[1]))
	  return t;
      } else if (roots[1] < EPSILON) {
	double t = roots[0];
	hit[0] = rayOrg[0] + rayDir[0]*t;
	hit[1] = rayOrg[1] + rayDir[1]*t;
	hit[2] = rayOrg[2] + rayDir[2]*t;
	sphereUV (data, &trimData->uvmap, hit, uv);
	if (trimmer->trim(trimmer, uv[0], uv[1]))
	  return t;
      } else {
	double mint = roots[0], maxt = roots[1];
	if (roots[0] > roots[1]) {
	  mint = roots[1]; maxt = roots[0];
	}
	hit[0] = rayOrg[0] + rayDir[0]*mint;
	hit[1] = rayOrg[1] + rayDir[1]*mint;
	hit[2] = rayOrg[2] + rayDir[2]*mint;
	sphereUV (data, &trimData->uvmap, hit, uv);
	if (trimmer->trim(trimmer, uv[0], uv[1]))
	  return mint;
	hit[0] = rayOrg[0] + rayDir[0]*maxt;
	hit[1] = rayOrg[1] + rayDir[1]*maxt;
	hit[2] = rayOrg[2] + rayDir[2]*maxt;
	sphereUV (data, &trimData->uvmap, hit, uv);
	if (trimmer->trim(trimmer, uv[0], uv[1]))
	  return maxt;
      }
    }
  }
  return -1.0;  /* ray does not intersect sphere */
}

static
void sphereImageColor (struct OBJECT *this, double hit[3], HIT_INFO *hitInfo,
		       double color[3]) {
  SPHERE_DATA *sdata = (SPHERE_DATA *) this->data;
  SPHERE_IMAGE_DATA *data = sdata->image;
  pnm_image *image = data->image;
  double scale;
  double u,v, uv[2];
  int i,j;

  /*
   * Map hit point to uv-plane.
   */
  sphereUV (sdata, &data->uvmap, hit, uv);
  u = uv[0]; v = uv[1];

  /*
   * Shift (u,v) to uv origin.
   */
  u -= data->org[0];
  v -= data->org[1];
  
  /*
   * If we are tiling the image and the hit point is outside
   * the image region we'll just return the plane's normal color.
   */
  if (!data->tile && (u < 0.0 || v < 0.0 || 
		      u > data->width || v > data->height)) {
    color[0] = sdata->color[0];
    color[1] = sdata->color[1];
    color[2] = sdata->color[2];
    return;
  }

  /*
   * Normalize to image region.
   */
  u /= data->width;
  v /= data->height;

  /*
   * Subtract of any whole terms;
   */
  u -= (double) ((int) u);
  v -= (double) ((int) v);

  /*
   * Wrap negative terms (negative floats truncate towards zero).
   */
  if (u < 0.0) u += 1.0;
  if (v < 0.0) v += 1.0;

  /*
   * Flip v axis since image vertical axis in opposite direction.
   */
  v = 1.0 - v;

  /*
   * Find corner pixel address.
   */
  u *= PNM_NC(image) - 1;
  v *= PNM_NR(image) - 1;
  j = (int) u;
  i = (int) v;

  /*
   * Compute color component normalization factor.
   */
  scale = 1.0/PNM_MAXVAL(image);

  if (data->interpolate) {  /* use bilinear interpolation */
    double w[4];

    /*
     * Determine the coordinate distances from corner pixel.
     */
    u -= j;
    v -= i;

    w[0] = (1.0 - u)*(1.0 - v); /* interpolation weights */
    w[1] = u*(1.0 - v);
    w[2] = (1.0 - u)*v;
    w[3] = u*v;
    
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
  } else {  /* just get the color of the closest pixel */
    color[0] = scale*PPM_PIXEL_R(image,i,j);
    color[1] = scale*PPM_PIXEL_G(image,i,j);
    color[2] = scale*PPM_PIXEL_B(image,i,j);
  }
}

static
void computeUVmapInfo(double pole[3], double equator[3], 
		      SPHERE_UVMAP *uvmap) {
  double scale;

  /*
   * Copy pole vector (make sure its normalized).
   */
  scale = 1.0/sqrt(pole[0]*pole[0] +
		   pole[1]*pole[1] +
		   pole[2]*pole[2]);
  uvmap->pole[0] = scale*pole[0];
  uvmap->pole[1] = scale*pole[1];
  uvmap->pole[2] = scale*pole[2];

  /*
   * Subtract off any component of equator in the pole direction
   * to make sure they are perpendicular.
   */
  scale = 
    equator[0]*uvmap->pole[0] +
    equator[1]*uvmap->pole[1] +
    equator[2]*uvmap->pole[2];
  uvmap->equator[0] = equator[0] - scale*uvmap->pole[0];
  uvmap->equator[1] = equator[1] - scale*uvmap->pole[1];
  uvmap->equator[2] = equator[2] - scale*uvmap->pole[2];

  /*
   * Make sure equator in normalized.
   */
  scale = 1.0/sqrt(uvmap->equator[0]*uvmap->equator[0] +
		   uvmap->equator[1]*uvmap->equator[1] +
		   uvmap->equator[2]*uvmap->equator[2]);
  uvmap->equator[0] *= scale;
  uvmap->equator[1] *= scale;
  uvmap->equator[2] *= scale;

  /*
   * Cache cross product pole x equator.
   */
  uvmap->pxe[0] = 
    uvmap->pole[1]*uvmap->equator[2] - 
    uvmap->pole[2]*uvmap->equator[1];
  uvmap->pxe[1] = 
    uvmap->pole[2]*uvmap->equator[0] - 
    uvmap->pole[0]*uvmap->equator[2];
  uvmap->pxe[2] = 
    uvmap->pole[0]*uvmap->equator[1] - 
    uvmap->pole[1]*uvmap->equator[0];
}

static
void sphereTextureColor(struct OBJECT *this, double hit[3], HIT_INFO *hitInfo,
			double color[3]) {
  SPHERE_DATA *sdata = (SPHERE_DATA *) this->data;
  PROCEDURAL_TEXTURE *tex = sdata->tex;
  (tex->func)(tex, hit, hitInfo, color);
}

void proceduralTextureSphere(OBJECT *object,
			     PROCEDURAL_TEXTURE *tex) {  
  SPHERE_DATA *data = (SPHERE_DATA *) object->data;

  if (data->tex == NULL)
    if ((data->tex = (PROCEDURAL_TEXTURE *) 
	 malloc(sizeof(PROCEDURAL_TEXTURE))) == NULL ||
	(tex->dataSize != 0 &&
	 (data->tex->data = malloc(tex->dataSize)) == NULL)) {
      perror("proceduralTextureSphere:malloc");
      exit(-1);
    }

  data->tex->dataSize = tex->dataSize;
  if (tex->dataSize > 0)
    memcpy(data->tex->data, tex->data, tex->dataSize);
  else
    data->tex->data = NULL;
  data->tex->func = tex->func;

  object->color = sphereTextureColor;
}

void trimSphere(OBJECT *object,              /* sphere object */
		double pole[3],              /* center to north pole */
		double equator[3],           /* center to equator pt */
		TRIM_OBJECT *trimmer) {      /* object trimmer */
  SPHERE_DATA *data = (SPHERE_DATA *) object->data;
  SPHERE_TRIM_DATA *trimdata;

  /*
   * Allocate memory for trimdata (if it isn't allocated already).
   */
  if (data->trim == NULL)
    if ((data->trim = (SPHERE_TRIM_DATA *) malloc(sizeof(SPHERE_TRIM_DATA)))
	== NULL) {
      perror("trimSphere:malloc");
      exit(-1);
    }

  trimdata = data->trim;

  /*
   * Create UV mapping data from equator and pole information.
   */
  computeUVmapInfo(pole, equator, &trimdata->uvmap);

  /*
   * Reference trimmer object.
   */
  trimdata->trimmer = trimmer;

  /*
   * Change ray intersection method.
   */
  object->rayHit = sphereTrimmedRayHit;
}

void mapImageToSphere(OBJECT *object,              /* sphere object */
		      int tile,                    /* tile image? */
		      int interpolate,             /* use bilinear interp? */
		      double pole[3],              /* center to north pole */
		      double equator[3],           /* center to equator pt */
		      double org[2],               /* image org in uv coords */
		      double width, double height, /* image uv size */
		      pnm_image *image)            /* PPM image */
{
  SPHERE_DATA *data = (SPHERE_DATA *) object->data;
  SPHERE_IMAGE_DATA *idata;
  
  /*
   * Allocate memory for image info if not allocated already.
   */
  if (data->image == NULL)
    if ((data->image = (SPHERE_IMAGE_DATA *) malloc(sizeof(SPHERE_IMAGE_DATA)))
	== NULL) {
      perror("mapImageToSphere:malloc");
      exit(-1);
    }

  idata = (SPHERE_IMAGE_DATA *) data->image;

  /*
   * Copy local parameters.
   */
  idata->tile = tile; 
  idata->interpolate = interpolate;
  idata->width = width;
  idata->height = height;
  idata->image = image;
  idata->org[0] = org[0];
  idata->org[1] = org[1];

  /*
   * Create UV mapping data from equator and pole information.
   */
  computeUVmapInfo(pole, equator, &idata->uvmap);

  /*
   * Change color method.
   */
  object->color = sphereImageColor;
}

OBJECT *createSphereObject(double center[3], double rad) {
  SPHERE_DATA *data;
  OBJECT *object;

  /*
   * Allocate memory for private sphere data and object.
   */
  if ((data = (SPHERE_DATA *) malloc(sizeof(SPHERE_DATA))) == NULL ||
      (object = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating sphere object");
    exit(-1);
  }

  /*
   * Set private sphere specific data.
   */
  data->center[0] = center[0];
  data->center[1] = center[1];
  data->center[2] = center[2];
  data->rad = rad;
  data->image = NULL;
  data->trim = NULL;

  /*
   * Assign appropriate private data and methods
   */
  object->data = data;
  object->rayHit = sphereRayHit;
  object->normal = sphereNormal;
  object->color = sphereColor;
  object->setColor = setSphereColor;

  /*
   * Assign some nice defaults for the rest.
   */
  data->color[0] = 0.2;
  data->color[1] = 0.8;
  data->color[2] = 0.2;
  object->ka = 0.15;
  object->kd = 0.40;
  object->ks = 0.45;
  object->kt = 0.55;
  object->ni = 1.52;
  object->phong = 6.0;

  return object;
}

void destroySphereObject(OBJECT *object) {
  if (((SPHERE_DATA *) (object->data))->image != NULL)
    free(((SPHERE_DATA *) (object->data))->image);
  /* XXX free up trim data and destroy trim object! */
  /* XXX free up texture data */
  free(object->data);
  free(object);
}
    
/* #define OLD_TEST_SPHERE */
#ifdef OLD_TEST_SPHERE

int main(void) {
  double center[3], eyePos[3], eyeDir[3], eyeUp[3];
  double ambient[3], background[3], color[3];
  OBJECT *spheres[9];
  LIGHT light;
  pnm_image *image;
  FILE *f;

  if ((f = fopen("spheres.ppm","wb")) == NULL) {
    perror("spheres.ppm");
    exit(-1);
  }

  /* main sphere */
  center[0] = center[1] = center[2] = 0.0;
  spheres[0] = createSphereObject(center, 1.0);
  spheres[0]->ka = 0.00;
  spheres[0]->kd = 0.00;
  spheres[0]->ks = 0.75;
  spheres[0]->kt = 0.95;
  spheres[0]->ni = 1.52;
  spheres[0]->phong = 3.0;

  /* sphere 2 */
  center[0] = -1.0; center[1] = -1.0; center[2] = -1.0;
  spheres[1] = createSphereObject(center, 0.35);
  color[0] = 0.9;  color[1] = 0.1;  color[2] = 0.1; 
  spheres[1]->setColor(spheres[1], color);

  /* sphere 3 */
  center[0] = -1.0; center[1] = -1.0; center[2] = 1.0;
  spheres[2] = createSphereObject(center, 0.35);
  spheres[2]->setColor(spheres[2], color);

  /* sphere 4 */
  center[0] = -1.0; center[1] = 1.0; center[2] = -1.0;
  spheres[3] = createSphereObject(center, 0.35);
  spheres[3]->setColor(spheres[3], color);

  /* sphere 5 */
  center[0] = -1.0; center[1] = 1.0; center[2] = 1.0;
  spheres[4] = createSphereObject(center, 0.35);
  spheres[4]->setColor(spheres[4], color);

  /* sphere 6 */
  center[0] = 1.0; center[1] = -1.0; center[2] = -1.0;
  spheres[5] = createSphereObject(center, 0.35);
  spheres[5]->setColor(spheres[5], color);

  /* sphere 7 */
  center[0] = 1.0; center[1] = -1.0; center[2] = 1.0;
  spheres[6] = createSphereObject(center, 0.35);
  spheres[6]->setColor(spheres[6], color);

  /* sphere 8 */
  center[0] = 1.0; center[1] = 1.0; center[2] = -1.0;
  spheres[7] = createSphereObject(center, 0.35);
  spheres[7]->setColor(spheres[7], color);

  /* sphere 9 */
  center[0] = 1.0; center[1] = 1.0; center[2] = 1.0;
  spheres[8] = createSphereObject(center, 0.35);
  spheres[8]->setColor(spheres[8], color);

  /* #define REFLECT_ONLY */
#ifdef REFLECT_ONLY
  {
    int i;
    for (i = 0; i < 9; i++)
      spheres[i]->kt = 0.0;
  }
#endif

  light.pos[0] = light.pos[1] = light.pos[2] = 5.0;
  light.color[0] = light.color[1] = light.color[2] = 1.0;
  light.c[0] = 1.0; light.c[1] = 0.01; light.c[2] = 0.0;

  eyePos[0] =  4.0; eyePos[1] = 0.0; eyePos[2] = 0.0; 
  eyeDir[0] = -1.0; eyeDir[1] = 0.0; eyeDir[2] = 0.0;
  eyeUp[0]  =  0.0; eyeUp[1]  = 1.0;  eyeUp[2] = 0.0;
  ambient[0] = 0.15; ambient[1] = 0.15; ambient[2] = 0.0;
  background[0] = 0.0; background[1] = 0.0; background[2] = 0.9;
  image = raytrace(eyePos, eyeDir, eyeUp,
		   1.0, 1.0, 1.0, 
		   512, 512, 2,
		   /* 256, 256, 1, */
		   /* 128, 128, 2, */
		   12, 
		   ambient, background,
		   1, &light,
		   9, spheres);

  write_pnm_image(image, f);

  return 0;
}

#endif  /* OLD_TEST_SPHERE */

#ifdef TEST_SPHERE

#include "plane.h"

int main(void) {
  double center[3], eyePos[3], eyeDir[3], eyeUp[3];
  double ambient[3], background[3];
  double point[3], normal[3];
  OBJECT *objects[2];
  LIGHT light;
  pnm_image *image, *checkerMap;
  FILE *f;

  if ((f = fopen("spheres.ppm","wb")) == NULL) {
    perror("spheres.ppm");
    exit(-1);
  }

  /* sphere */
  center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
  objects[0] = createSphereObject(center, 1.0);
  /*
  objects[0]->ka = 0.00;
  objects[0]->kd = 0.00;
  objects[0]->ks = 0.30;
  objects[0]->kt = 0.95;
  objects[0]->ni = 1.52;
  objects[0]->phong = 6.0;
  */
  objects[0]->ka = 0.05;
  objects[0]->kd = 0.05;
  objects[0]->ks = 0.50;
  objects[0]->kt = 0.80;
  objects[0]->ni = 1.52;
  objects[0]->phong = 8.0;

  /* plane */
  normal[0] = 0.0; normal[1] = 1.0; normal[2] = 0.0;
  point[0] = 0.0; point[1] = -1.00001; point[2] = 0.0;
  objects[1] = createPlaneObjectFromNormalandPoint(normal, point);
  objects[1]->ka = 0.15;
  objects[1]->kd = 0.55;
  objects[1]->ks = 0.50;
  objects[1]->kt = 0.00;
  objects[1]->ni = 1.52;
  objects[1]->phong = 6.0;

  if ((checkerMap = read_pnm_image_from_file("check.ppm")) == NULL) {
    perror("check.ppm");
    exit(-1);
  }
  
  center[0] = 0.0; center[1] = 0.0; center[2] = 0.0;
  normal[0] = -1.0; normal[1] = 0.0; normal[2] = 1.0;
  mapImageToPlane(objects[1], 1, 0, center, normal, 1.00, 1.00, checkerMap);

  light.pos[0] = light.pos[1] = light.pos[2] = 5.0;
  light.color[0] = light.color[1] = light.color[2] = 1.0;
  light.c[0] = 1.0; light.c[1] = 0.01; light.c[2] = 0.0;

  /*
  eyePos[0] =  4.0; eyePos[1] = 0.0; eyePos[2] = 0.0; 
  eyeDir[0] = -1.0; eyeDir[1] = 0.0; eyeDir[2] = 0.0;
  */
  eyePos[0] =  4.0; eyePos[1] = 0.9; eyePos[2] = 0.0; 
  eyeDir[0] = -eyePos[0]; eyeDir[1] = -eyePos[1]; eyeDir[2] = -eyePos[2];
  eyeUp[0]  =  0.0; eyeUp[1]  = 1.0;  eyeUp[2] = 0.0;
  ambient[0] = 0.15; ambient[1] = 0.15; ambient[2] = 0.0;
  background[0] = 0.0; background[1] = 0.0; background[2] = 0.9;
  image = raytrace(eyePos, eyeDir, eyeUp,
		   1.0, 1.0, 1.0, 
		   512, 512, 2,
		   /* 256, 256, 1, */
		   /* 128, 128, 2, */
		   12, 
		   ambient, background,
		   1, &light,
		   2, objects);

  write_pnm_image(image, f);

  return 0;
}

#endif  /* TEST_SPHERE */
