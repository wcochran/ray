
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef RAYTRACE_H
#define RAYTRACE_H

#include "pnmio.h"

/* 
 * Precision hack:
 *   Any t value less than this is considered to be at or behind
 *   the eye point.
 */
#define EPSILON 0.0001

/*
 * A really big number.
 */
#ifdef INFINITY  /* use my value for infinity */
#undef INFINITY
#endif
#define INFINITY 1.0e20

/*
 * Sometimes we'll use an array of doubles for a coordinate value, and
 * sometimes we'll use a structure.
 */
typedef struct {
  double x,y,z;
} POINT3;

typedef struct {
  double x,y;
} POINT2;

/*
 * The information that is computed during ray/object intersection
 * calculations can be cached for normal computations (and other
 * computations that could use precomputed intersection info). The type
 * of information needed is different for each object and the following
 * union stores the necessasy information: when an object's normal()
 * method is called it will be pased a ptr to the object's instance,
 * the hitpoint, and the information here that was computed by
 * the object's rayHit() method. The reason we do not necessarily 
 * compute the actual normal in the rayHit() method is that is
 * may be very expensive to compute and completely unnecessary
 * since another object may be closer to the eye. The reason
 * we do not cache this information in the structure itself is that
 * would create non-reintrant code. Unfortunatey this union
 * is not very opaque, but if we used something like a void *
 * pointer then we would be malloc'ing all the time which is
 * both time consuming and non-reentrant.
 */
typedef union {
  struct {
    int nada;        /* sphere has all it needs to know */
  } sphere;
  struct {
    int nada;        /* plane has all it needs to know */
  } plane;
  struct {
    double s,t;      /* parameters of intersection point */
  } bezier3;
  struct {           /* collection of bezier patches */ 
    int i;
    double s, t;
  } bezier3_mesh;
  struct {
    double h[3];     /* hit point in modelling coordinates */
  } superellipsoid;
} HIT_INFO;

/*
 * Structure used to support procedural texturing.
 * Encapsulates function to call when determining the color of an object's
 * surface at its surface.
 */
typedef struct proceduralTexture {
  size_t dataSize;               /* size (in bytes) of data field */
  void *data;                    /* texture dependant data */
  void (*func)                   /* procedural function */
       (struct proceduralTexture *this, /* ptr to instance of this struct */
	double hit[3],           /* 3D hit point to query texture value */
	HIT_INFO *hitInfo,       /* extra info about hit point */
	double color[3]);        /* return texture value (i.e. color) */
} PROCEDURAL_TEXTURE;

/*
 * All objects to be raytraced will be defined by
 * the following structure. Any private information
 * needed to model that object will be referenced with
 * the (void *data) feild below. All objects must implement
 * all the methods shown here.
 */
typedef struct OBJECT {
  double ka;            /* ambient reflection coefficient    (0 <= ka <= 1) */
  double kd;            /* diffuse reflection coefficient    (0 <= kd <= 1) */
  double ks;            /* specular reflection coefficient   (0 <= ks <= 1) */
  double kt;            /* transmission coefficient          (0 <= kt <= 1) */
  double ni;            /* material index of refraction */
  double phong;         /* specular phong exponent */
  void *data;           /* private data */
  double (*rayHit)      /* ray intersection method (return param t) */
    (struct OBJECT *this,   /* reference to current instance */
     double rayOrg[3],      /* origin of ray */
     double rayDir[3],      /* direction of ray */
     HIT_INFO *info);       /* cached info for normal() & color() methods */
  void (*normal)        /* ray intersection normal method */
    (struct OBJECT *this,   /* reference to current instance */
     double hitPoint[3],    /* intersection point */
     HIT_INFO *info,        /* cached info computed by rayHit() method */     
     double normal[3]);     /* return unit normal vector */
  void (*color)         /* get object's RGB diffuse color */
    (struct OBJECT *this,   /* reference to current instance */
     double hitPoint[3],    /* intersection point */
     HIT_INFO *info,        /* cached info computed by rayHit() method */     
     double color[3]);      /* color at intersection point */
  void (*setColor)      /* set object's color (single colored object's only) */
    (struct OBJECT *this,   /* reference to current instance */
     double color[3]);      /* RGB diffuse color */    
} OBJECT;

typedef struct LIGHT {
  double pos[3];       /* position of point light source */
  double color[3];     /* RGB color of light source */
  double c[3];         /* light attenuation coeffs for function
                          (f(d) = MIN(1/(c0 +c1*d + c2*d^2), 1)
                          d = distance light traveled, f(d) = attenuation */
} LIGHT;

#if defined(USE_PTHREADS) || defined(USE_PVM)

typedef struct {
  int threadIndex;       /* 0,1,2,...,numThreads-1 */
  int numThreads;        /* number of rendering threads being used */
  pnm_image *ppm;        /* image to render to */
  double eyePos[3];      /* position of eye */
  double eyeDir[3];      /* view direction */
  double eyeUp[3];       /* eye up vector */
  double projWidth;      /* width of projection plane */
  double projHeight;     /* height of projection plane */
  double projDist;       /* distannce from eye to projection plane */
  int imageWidth;        /* width (in pixels) of image */
  int imageHeight;       /* height of image */
  int samples;           /* sqrt of samples per pixel */
  int maxDepth;          /* maximum recursion depth (1,2,..) */
  double ambient[3];     /* RGB intensities of ambient light */
  double background[3];  /* RGB background light intensities */
  int numLights;         /* number of point light sources */
  LIGHT *lights;         /* array of point light sources */
  int numObjects;        /* number of objects in world */
  OBJECT **objects;      /* array of objects in world */
} RAYTRACE_THREAD_ARGS;

/*
 * raytraceThread:
 *   note: in pthread_create, typecast this to (void *) (*) (void *)
 *   note: the arguments are used as read only so the same structure
 *         can be used for different thread.
 */
pnm_image *raytraceThread(RAYTRACE_THREAD_ARGS *);

#else /* not USE_PTHREADS and not USE_PVM */

pnm_image * 
raytrace(double eyePos[3],      /* position of eye */
	 double eyeDir[3],      /* view direction */
	 double eyeUp[3],       /* eye up vector */
	 double projWidth,      /* width of projection plane */
	 double projHeight,     /* height of projection plane */
	 double projDist,       /* distannce from eye to projection plane */
	 int imageWidth,        /* width (in pixels) of image */
	 int imageHeight,       /* height of image */
	 int samples,           /* sqrt of samples per pixel */
	 int maxDepth,          /* maximum recursion depth (1,2,..) */
	 double ambient[3],     /* RGB intensities of ambient light */
	 double background[3],  /* RGB background light intensities */
	 int numLights,         /* number of point light sources */
	 LIGHT lights[],        /* array of point light sources */
	 int numObjects,        /* number of objects in world */
	 OBJECT *objects[]);    /* array of objects in world */

#endif /* USE_PTHREADS */

#endif

