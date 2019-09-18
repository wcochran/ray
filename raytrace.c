
/* $Author: cs548 $ $Revision: 1.2 $ $Date: 2009/10/09 20:23:17 $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "pnmio.h"
#include "raytrace.h"

/*
#define USE_HALFWAY
*/

/*
 * Factor used for each non-opaque object so even completely
 * transparent objects cast a slight shadow.
 */
#define SHADOW_ATTENUATION 0.85

static
void traceray(double rayOrg[3],      /* ray origin */ 
	      double rayDir[3],      /* direction of ray (normalized) */
	      int depth,             /* current recursion depth */
	      int maxDepth,          /* maximum recursion depth */
	      double ni,             /* index of refraction of environment */
	      double ambient[3],     /* RGB intensities of ambient light */
	      double background[3],  /* RGB background light intensities */
	      int numLights,         /* number of point light sources */
	      LIGHT lights[],        /* array of point light sources */
	      int numObjects,        /* number of objects in world */
	      OBJECT *objects[],     /* array of objects in world */
	      double color[3])       /* return RGB color intensities */
{                                    /* these values are *not* clamped */
  int i;
  double mint = 0, cosi;
  OBJECT *obj;
  double hit[3], normal[3];
  HIT_INFO hitInfo;
  double objColor[3];
  double shadow;
  int entering;

  /*
   * If we have exceeded the may raytracing depth then terminate the
   * recursion and return color with no influence.
   */
  if (depth > maxDepth) {
    color[0] = color[1] = color[2] = 0.0;
    return;
  }
  
  /*
   * Find the first object the given ray intersects.
   */
  for (i = 0, obj = NULL; i < numObjects; i++) {
    HIT_INFO info;
    double t = (*(objects[i]->rayHit)) (objects[i], rayOrg, rayDir, &info);
    if (t > EPSILON && (obj == NULL || t < mint)) {
      mint = t;
      obj = objects[i];
      memcpy(&hitInfo, &info, sizeof(HIT_INFO));
    }
  }

  /*
   * If our ray intersected no objects then return given background color.
   */
  if (obj == NULL) {
    color[0] = background[0]; 
    color[1] = background[1]; 
    color[2] = background[2];
    return;
  }

  /* 
   * Compute hit point and normal to surface.
   */
  hit[0] = rayOrg[0] + mint*rayDir[0];
  hit[1] = rayOrg[1] + mint*rayDir[1];
  hit[2] = rayOrg[2] + mint*rayDir[2];
  (*obj->normal) (obj, hit, &hitInfo, normal);

  /*
   * Fetch the object's color at the hit point.
   */
  (*obj->color) (obj, hit, &hitInfo, objColor);

  /*
   * Initialize local shading color with ambient term.
   */
  color[0] = ambient[0]*obj->ka*objColor[0];
  color[1] = ambient[1]*obj->ka*objColor[1];
  color[2] = ambient[2]*obj->ka*objColor[2];

  /*
   * Compute the cosine of the incident angle; flip the direction
   * of the normal vector if necessary. This also tells us whether
   * our ray is leaving or entering the object.
   */
  cosi = -normal[0]*rayDir[0] - normal[1]*rayDir[1] - normal[2]*rayDir[2];
  entering = 1;
  if (cosi < 0.0) {
    cosi = -cosi;
    normal[0] = -normal[0]; normal[1] = -normal[1]; normal[2] = -normal[2];
    entering = 0;
  }
    
  /*
   * Determine local shading using each light source.
   */
  for (i = 0; i < numLights; i++) {
    int j;
    double lightDist, lightDir[3];
    double atten;

    /*
     * Determine distance and direction to current light source.
     */
    lightDir[0] = lights[i].pos[0] - hit[0];
    lightDir[1] = lights[i].pos[1] - hit[1];
    lightDir[2] = lights[i].pos[2] - hit[2];
    lightDist = sqrt(lightDir[0]*lightDir[0] +
		     lightDir[1]*lightDir[1] + 
		     lightDir[2]*lightDir[2]);
    if (lightDist > 0.0) {
      double scale = 1.0/lightDist;
      lightDir[0] *= scale;
      lightDir[1] *= scale;
      lightDir[2] *= scale;
    } else {  /* we are at the light! */
      color[0] = lights[i].color[0];
      color[1] = lights[i].color[1];
      color[2] = lights[i].color[2];
      return;
    }

    /*
     * Determine the light source visability at the current
     * object intersection point.
     */
    for (shadow = 1.0, j = 0; j < numObjects; j++) {
      HIT_INFO info;
      double t = (*objects[j]->rayHit) (objects[j], hit, lightDir, &info);
      if (objects[j] == obj && t < EPSILON) 
	continue;  /* bogus self-occlusion */
      if (t > 0.0 && t < lightDist) {
	if ((shadow *= objects[j]->kt) <= 0.0)
	  break;
	shadow *= SHADOW_ATTENUATION;
      }
    }
    if (shadow <= 0.0) continue;  /* completely in shadow */

    /*
     * Compute light attenuation.
     */
    atten = 1.0/(lights[i].c[0] + 
		 lightDist*(lights[i].c[1] + lights[i].c[2]*lightDist));
    if (atten > 1.0) atten = 1.0;

    /*
     * Factor in shadow attentuation.
     */
    atten *= shadow;

    /*
     * Compute local diffuse reflection color.
     */
    if (obj->kd > 0.0) {
      double cosine = 
	normal[0]*lightDir[0] + 
	normal[1]*lightDir[1] + 
	normal[2]*lightDir[2];
      if (cosine > 0.0) {
	double prod = atten*obj->kd*cosine;
	color[0] += lights[i].color[0]*prod*objColor[0];
	color[1] += lights[i].color[1]*prod*objColor[1];
	color[2] += lights[i].color[2]*prod*objColor[2];
      }
    }

    /*
     * Compute local specular reflection color.
     */
#ifdef USE_HALFWAY
    if (obj->ks > 0.0) {
      double mag, H[3];
      H[0] = lightDir[0] - rayDir[0];  /* halfway vector */
      H[1] = lightDir[1] - rayDir[1];
      H[2] = lightDir[2] - rayDir[2];
      mag = sqrt(H[0]*H[0] + H[1]*H[1] + H[2]*H[2]);
      if (mag > 0.0) {
	double cosine, scale = 1.0/mag;
	H[0] *= scale; H[1] *= scale; H[2] *= scale;
	cosine = normal[0]*H[0] + normal[1]*H[1] + normal[2]*H[2];
	if (cosine > 0.0) {
	  scale = atten*obj->ks*pow(cosine,obj->phong);
	  color[0] += lights[i].color[0]*scale;
	  color[1] += lights[i].color[1]*scale;
	  color[2] += lights[i].color[2]*scale;
	}
      }
    }
#else  /* use reflection vector */
    if (obj->ks > 0.0) {
      double scale, cosine, R[3];
      scale = 2.0*(lightDir[0]*normal[0] +
		   lightDir[1]*normal[1] +
		   lightDir[2]*normal[2]);
      R[0] = scale*normal[0] - lightDir[0];  /* reflection vector */
      R[1] = scale*normal[1] - lightDir[1];
      R[2] = scale*normal[2] - lightDir[2];
      cosine = -(rayDir[0]*R[0] + rayDir[1]*R[1] + rayDir[2]*R[2]);
      if (cosine > 0.0) {
	scale = atten*obj->ks*pow(cosine,obj->phong);
	color[0] += lights[i].color[0]*scale;
	color[1] += lights[i].color[1]*scale;
	color[2] += lights[i].color[2]*scale;
      }
    }
#endif

  } /* end lights iteration (local shading) */

  /*
   * Global lighting computations (if we're below the max recursion level)
   */
  if (depth >= maxDepth) return;

  /*
   * Global reflection computation.
   * Compute reflection vector, recurse, and mix in result.
   */
  if (obj->ks > 0.0) {
    double scale = 2.0*cosi;
    double R[3];
    double Rcolor[3];
    R[0] = rayDir[0] + scale*normal[0];
    R[1] = rayDir[1] + scale*normal[1];
    R[2] = rayDir[2] + scale*normal[2];
    traceray(hit,R, depth+1,maxDepth, ni, ambient,background,
	     numLights, lights, numObjects, objects, Rcolor);
    color[0] += obj->ks*Rcolor[0];
    color[1] += obj->ks*Rcolor[1];
    color[2] += obj->ks*Rcolor[2];
  }

  /*
   * Global transmission computation.
   * Compute transmission vector, recurse, and mix in result.
   * Whether the current ray is entering or leaving the object
   * determines the relative index of refraction.
   */
  if (obj->kt > 0.0) {
    double nit = entering ? ni/obj->ni : obj->ni/ni; /* entering or leaving? */
    double det = 1.0 + nit*nit*(cosi*cosi - 1.0);
    if (det > 0.0) {  /* else total internal reflection */
      double T[3];
      double Tcolor[3];
      double scale = nit*cosi - sqrt(det);
      T[0] = nit*rayDir[0] + scale*normal[0];
      T[1] = nit*rayDir[1] + scale*normal[1];
      T[2] = nit*rayDir[2] + scale*normal[2];
      traceray(hit,T, depth+1,maxDepth, ni, ambient,background,
	       numLights, lights, numObjects, objects, Tcolor);
      color[0] += obj->kt*Tcolor[0];
      color[1] += obj->kt*Tcolor[1];
      color[2] += obj->kt*Tcolor[2];
    }
  }

}

#if defined(USE_PTHREADS) || defined(USE_PVM)
static
#endif
pnm_image * 
raytrace(
#if defined(USE_PTHREADS) || defined(USE_PVM)
	 int threadIndex,       /* 0,1,2,...,numThreads-1 */
	 int numThreads,        /* number of rendering threads being used */
	 pnm_image *ppm,        /* image to render to (this is returned) */
#endif
	 double eyePos[3],      /* position of eye */
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
	 OBJECT *objects[])     /* array of objects in world */
	 
{
  double mag;
  double eyeRight[3];
  double X,Y, DX,DY, dx,dy;
  double colorScale;
  int i,j;
#if !defined(USE_PTHREADS) && !defined(USE_PVM)
  pnm_image *ppm;
#endif
  
  /*
   * Normalize eye direction vector.
   */
  mag = sqrt(eyeDir[0]*eyeDir[0] + eyeDir[1]*eyeDir[1] + eyeDir[2]*eyeDir[2]);
  if (mag <= 0.0) return NULL;
  mag = 1.0/mag;
  eyeDir[0] *= mag;  eyeDir[1] *= mag;  eyeDir[2] *= mag;

  /*
   * Force up vector to be perpendicular to eye direction vector
   * and normalize it.
   */
  mag = eyeUp[0]*eyeDir[0] + eyeUp[1]*eyeDir[1] + eyeUp[2]*eyeDir[2];
  eyeUp[0] -= mag*eyeDir[0];  
  eyeUp[1] -= mag*eyeDir[1];
  eyeUp[2] -= mag*eyeDir[2];
  mag = sqrt(eyeUp[0]*eyeUp[0] + eyeUp[1]*eyeUp[1] + eyeUp[2]*eyeUp[2]);
  if (mag <= 0.0) return NULL;
  mag = 1.0/mag;
  eyeUp[0] *= mag;  eyeUp[1] *= mag;  eyeUp[2] *= mag;

  /*
   * Compute eyeRight vector as the cross product eyeDir x eyeUp (we end up
   * with a left handed coordinate system if we think of eyeDir being
   * the z-axis and eyeRight and eyeUp being the x and y axes respectively).
   */
  eyeRight[0] = eyeDir[1]*eyeUp[2] - eyeDir[2]*eyeUp[1];
  eyeRight[1] = eyeDir[2]*eyeUp[0] - eyeDir[0]*eyeUp[2];
  eyeRight[2] = eyeDir[0]*eyeUp[1] - eyeDir[1]*eyeUp[0];

#if !defined(USE_PTHREADS) && !defined(USE_PVM)
  /*
   * Allocate PPM image.
   */
  if ((ppm = allocate_ppm_image(imageHeight, imageWidth)) == NULL)
    return NULL;
#endif

  /*
   * Compute projection grid information:
   *   (DX,DY) = per pixel grid spacing;
   *   (dx,dy) = per subpixel grid spacing.
   */
  DX = projWidth/imageWidth;
  DY = projHeight/imageHeight;
  dx = DX/samples;
  dy = DY/samples;

  /*
   * Loop from left to right and from top to bottom through the image
   * and determine the color of each associated pixel by casting
   * the appropriate rays from the eye point through the projection
   * plane. The variable X and Y are used for scanning each sampled
   * region on the projection plane corresponding to a pixel. 
   * The variables x and y are for subsampling each such region
   * according to the number of samples used for each pixel.
   * The variables X,Y, DX,DY, dx,dy are all used relative to the
   * coordinate system attached to the center of the projection plane.
   */
  colorScale = 1.0/(samples*samples);
  for (Y = 0.5*projHeight, i = 0; i < imageHeight; i++, Y -= DY)
    for (X = -0.5*projWidth, j = 0; j < imageWidth; j++, X += DX) {
      int m,n;
      double x,y;
      double sumColor[3];

      if (i == 89 && j == 36)
	x = 0.0;

#if defined(USE_PTHREADS) || defined(USE_PVM)
      /*
       * Only work on pixels assigned to this thread.
       */
      if ((i + j) % numThreads != threadIndex)
	continue;
#endif

      sumColor[0] = sumColor[1] = sumColor[2] = 0.0;

      for (y = Y - 0.5*dy, n = samples; --n >= 0; y -= dy)
	for (x = X + 0.5*dx, m = samples; --m >= 0; x += dx) {
	  double color[3];
	  double rayDir[3];
	  rayDir[0] = projDist*eyeDir[0] + x*eyeRight[0] + y*eyeUp[0];
	  rayDir[1] = projDist*eyeDir[1] + x*eyeRight[1] + y*eyeUp[1];
	  rayDir[2] = projDist*eyeDir[2] + x*eyeRight[2] + y*eyeUp[2];
	  mag = sqrt(rayDir[0]*rayDir[0] +
		     rayDir[1]*rayDir[1] +
		     rayDir[2]*rayDir[2]);
	  mag = 1.0/mag;
	  rayDir[0] *= mag; rayDir[1] *= mag; rayDir[2] *= mag;
	  traceray(eyePos, rayDir, 1, maxDepth, 
		   1.0, ambient, background,
		   numLights, lights, numObjects, objects, color);
	  sumColor[0] += color[0];
	  sumColor[1] += color[1];
	  sumColor[2] += color[2];
	}
      if ((sumColor[0] *= colorScale) > 1.0) sumColor[0] = 1.0;
      if ((sumColor[1] *= colorScale) > 1.0) sumColor[1] = 1.0;
      if ((sumColor[2] *= colorScale) > 1.0) sumColor[2] = 1.0;
      PPM_PIXEL_R(ppm,i,j) = (int) (sumColor[0]*PNM_MAXVAL(ppm) + 0.5);
      PPM_PIXEL_G(ppm,i,j) = (int) (sumColor[1]*PNM_MAXVAL(ppm) + 0.5);
      PPM_PIXEL_B(ppm,i,j) = (int) (sumColor[2]*PNM_MAXVAL(ppm) + 0.5);
    }

  return ppm;
}

#if defined(USE_PTHREADS) || defined(USE_PVM)

pnm_image *raytraceThread(RAYTRACE_THREAD_ARGS *args) {
  return raytrace(args->threadIndex, args->numThreads,
		  args->ppm,
		  args->eyePos, args->eyeDir, args->eyeUp,
		  args->projWidth, args->projHeight, args->projDist,
		  args->imageWidth, args->imageHeight,
		  args->samples,
		  args->maxDepth,
		  args->ambient, args->background,
		  args->numLights, args->lights,
		  args->numObjects, args->objects);
}

#endif
