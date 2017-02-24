#ifndef VIEWPICK_H
#define VIEWPICK_H

#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#endif

typedef struct {
  int winWidth, winHeight;             /* params to glViewport() */
  GLdouble eye[3], lookat[3], up[3];   /* params to gluLookAt() */
  GLdouble fovy, aspect, zNear, zFar;  /* params to gluPerspective() */
  
  int dirty;                           /* recompute screen-to-ray info? */
  double netAspect;                    /* (w/h)*(H/W) */
  double halfWidth, halfHeight;        /* H', W' */
  double unitRight[3], unitUp[3];      /* right, up */
  double forward[3];                   /* z*forward */
} ViewPick;

void setViewport(ViewPick *vp, int w, int h);
void setLookAt(ViewPick *vp, 
               GLdouble eye[3], GLdouble lookat[3], GLdouble up[3]);
void setPerspective(ViewPick *vp,
                    GLdouble fovy, GLdouble aspect,
                    GLdouble zNear, GLdouble zFar);

void screenToRay(ViewPick *vp, int x, int y, double unitDir[3]);

double dist2PointToRay(double org[3], double unitDir[3], double p[3]);
double closestPointOnRay(double srcOrg[3], 
                         double srcDir[3],
                         double dstOrg[3], 
                         double dstDir[3],
                         double *srct);

#endif /* VIEWPICK_H */
