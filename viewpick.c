#include <math.h>
#include "viewpick.h"

#define EPSILON 0.0000001

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

#define DTOR(deg) ((M_PI/180)*(deg))

#define DOT(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])
#define CROSS(U,V,UxV) ((UxV)[0] = (U)[1]*(V)[2] - (U)[2]*(V)[1],\
                        (UxV)[1] = (U)[2]*(V)[0] - (U)[0]*(V)[2],\
                        (UxV)[2] = (U)[0]*(V)[1] - (U)[1]*(V)[0])

/*
 * computeScreenToRayInfo()
 * purpose:
 *   Precomputes info for coverting a screen coordinate (x,y)
 *   into a ray passing through the scene in world coordinates.
 * description:
 *     Here we precompute the info needed to convert a 2-D screen
 *  coordinate into a 3-D ray that passes through the scene. We
 *  assume that the user is using a perspective projection (i.e.
 *  using gluPerspective()), and is setting the viewing transformation
 *  using the gluLookAt() function. We also assume that the user
 *  is mapping points from the 3-D projection plane to the screen
 *  using the glViewport() function.
 *     We define the following variables:
 *
 *        W,H = dimensions of 2-D screen
 *        w,h = dimension of 3-D projection screen
 *        aspect = w/h
 *        theta = vertical field of view angle
 *
 *  We can map from screen coordinates (xs,ys) (where 0,0 is in the
 *  upper left-hand corner with y-values increasing as you proceed
 *  down the screen) to projection screen coordinates (x,y) (where
 *  (0,0) is in the center of the screen an y is "up") as follows:
 *
 *        x = (w/W)*(xs - W/2), y = (h/H)*(H/2 - ys)
 *
 *  If we let h = H, and thus w = aspect*H we get
 *
 *        x = aspect*H/W*(xs - W/2),  y = H/2 - ys
 *
 *  If we let A = aspect*H/W, W' = W/2, H' = H/2 we get
 *
 *        x = A*(xs - W'), y = H' - ys
 *
 *  Let theta' = theta/2, We can solve for z (distance to projection screen)
 *
 *        z = (H/2)/tan(theta/2) = H'/tan(theta')
 *
 *  which is invariant wrt to (xs,ys) so we'll precompute it.
 *
 *  The values x,y,z represent the weights of the linear combination
 *  of the unit right, up, and forward directions that yield
 *  the ray eminating from the camera through the point (xs, ys)
 *  in world coordinates.
 *
 *        ray direction vector = x*right + y*up + z*forward.
 * note:
 *   Anytime the camera moves, the perspective projection matrix
 *   changes, or the screen viewport dimensions change this
 *   method needs to be called. We use a 'dirty bit' to determine
 *   when these values need to be changed.
 */
static void computeScreenToRayInfo(ViewPick *vp) {
  double scale;

  vp->forward[0] = vp->lookat[0] - vp->eye[0];     /* forward vector */
  vp->forward[1] = vp->lookat[1] - vp->eye[1];
  vp->forward[2] = vp->lookat[2] - vp->eye[2];
  scale = 1.0/sqrt(DOT(vp->forward, vp->forward));
  vp->forward[0] *= scale;
  vp->forward[1] *= scale;
  vp->forward[2] *= scale;
  
  scale = DOT(vp->up,vp->forward);                 /* make perpendicular */
  vp->unitUp[0] = vp->up[0] - scale*vp->forward[0];
  vp->unitUp[1] = vp->up[1] - scale*vp->forward[1];
  vp->unitUp[2] = vp->up[2] - scale*vp->forward[2];
  scale = 1.0/sqrt(DOT(vp->unitUp, vp->unitUp));   /* unitize */
  vp->unitUp[0] *= scale;
  vp->unitUp[1] *= scale;
  vp->unitUp[2] *= scale;

  CROSS(vp->forward, vp->unitUp, vp->unitRight);   /* unit right vector */
  
  vp->halfWidth  = 0.5*vp->winWidth;
  vp->halfHeight = 0.5*vp->winHeight;

  vp->netAspect = vp->aspect*vp->halfHeight/vp->halfWidth;

  scale = vp->halfHeight/tan(0.5*DTOR(vp->fovy));
  vp->forward[0] *= scale;
  vp->forward[1] *= scale;
  vp->forward[2] *= scale;
  
  vp->dirty = 0;
}

void setViewport(ViewPick *vp, int w, int h) {
  vp->winWidth = w; 
  vp->winHeight = h; 
  vp->dirty = 1;
}

void setLookAt(ViewPick *vp, GLdouble e[3], GLdouble l[3], GLdouble u[3]) {
  vp->eye[0] = e[0]; vp->eye[1] = e[1]; vp->eye[2] = e[2];
  vp->lookat[0] = l[0]; vp->lookat[1] = l[1]; vp->lookat[2] = l[2];
  vp->up[0] = u[0];  vp->up[1] = u[1]; vp->up[2] = u[2];
  vp->dirty = 1;
}

void setPerspective(ViewPick *vp, GLdouble fovy, GLdouble aspect,
                    GLdouble zNear, GLdouble zFar) {
  vp->fovy = fovy;
  vp->aspect = aspect;
  vp->zNear = zNear;
  vp->zFar = zFar;
  vp->dirty = 1;
}

/*
 * screenToRay()
 * purpose:
 *  Converts a screen viewport coodinate (x,y) into a ray cast
 *  through the scene.
 * description:
 *  See comment block for computeScreenToRayInfo() method.
 * note:
 *  Invokes computeScreenToRayInfo() if dirty bit set.
 * input:
 *   x,y: integer screen coordinate (e.g. location of mouse click).
 * output:
 *   unitDir: unit direction vector of ray eminating from the eye point.
 */
void screenToRay(ViewPick *vp, int x, int y, double unitDir[3]) {
  double scale, xx, yy;
  if (vp->dirty) computeScreenToRayInfo(vp);
  xx = vp->netAspect*(x - vp->halfWidth);
  yy = vp->halfHeight - y;
  unitDir[0] = xx*vp->unitRight[0] + yy*vp->unitUp[0] + vp->forward[0];
  unitDir[1] = xx*vp->unitRight[1] + yy*vp->unitUp[1] + vp->forward[1];
  unitDir[2] = xx*vp->unitRight[2] + yy*vp->unitUp[2] + vp->forward[2];
  scale = 1.0/sqrt(DOT(unitDir,unitDir));
  unitDir[0] *= scale;
  unitDir[1] *= scale;
  unitDir[2] *= scale;
}

/*
 * dist2PointToRay()
 * purpose: Returns the distance (actually the square of the Euclidean
 *   distance) of a given point from a ray.
 * description:
 *   The distance D is calculated as
 *       D = || v - (v.u)u ||
 *   where
 *       v = p - org
 *       u = ray direction
 * input:
 *   org, unitDir: origin and unit direction vector of ray;
 *   p : point that we are using for our distance measurement.
 * returns:
 *  square of distance from point to ray.
 */
double dist2PointToRay(double org[3], double unitDir[3], double p[3]) {
  double v[3], scale;
  
  v[0] = p[0] - org[0];
  v[1] = p[1] - org[1];
  v[2] = p[2] - org[2];
  
  scale = DOT(v, unitDir);
  v[0] -= scale*unitDir[0];
  v[1] -= scale*unitDir[1];
  v[2] -= scale*unitDir[2];
  
  return DOT(v,v);
}

/*
 * closestPointOnRay()
 * purpose:
 *   Give a "source ray" and a "destination ray" this function computes
 *   the t-value for the point along the source ray corresponding to
 *   the closet point on that ray to the other ray.
 * decription:
 *   Let the source ray be r(t) = r + td (r = origin, d = unit direction).
 *   Let the dest ray be p(t) = p + tv. (p = origin, v = unit direction).
 *   The distance of a point r(t) from the the ray defined by p and v is
 *            D(t) = ||u(t) - (u(t).v)v||
 *   where
 *            u(t) = r + td - p = r-p + td
 *   so
 *            D(t) = || r-p + td - [((r-p).v) + t(d.v)v] ||
 *                 = || r-p - ((r-p).v)v + t(d - (d.v)v) ||
 *                 = || Q + tR ||
 *   where
 *            Q = r-p - ((r-p).v)v
 *            R = (d - (d.v)v)
 *
 *   so
 *            D(t) = (Qx + tRx)^2 + (Qy + tRy)^2 + (Qz + tRz)^2
 *
 *   Taking the derivative wrt to t we get
 *            D'(t) = 2[Rx(Qx + tRx) + Ry(Qy + tRy) + Rz(Qz + tRz)]
 *                  = 2[R.Q + t(R.R)]
 *   setting D'(t) equal to zero we get
 *                 -R.Q
 *            t = -------
 *                  R.R
 * input:
 *   srcOrg, srcDir: source ray (origin and unit direction vector);
 *   dstOrg, dstDir: destination ray (origin and unit direction vector);
 * output:
 *   srct: t-value along source ray closest to destination ray;
 * returns:
 *   Square of distance between rays.
 */
double closestPointOnRay(double srcOrg[3], 
                         double srcDir[3],
                         double dstOrg[3], 
                         double dstDir[3],
                         double *srct) {
  double pr[3], Q[3], R[3], V[3], scale, t;
  
  pr[0] = srcOrg[0] - dstOrg[0];         /* pr <-- r-p */
  pr[1] = srcOrg[1] - dstOrg[1];
  pr[2] = srcOrg[2] - dstOrg[2];
  
  scale = DOT(pr,dstDir);                /* Q <-- r-p - ((r-p).v)v */
  Q[0] = pr[0] - scale*dstDir[0];
  Q[1] = pr[1] - scale*dstDir[1];
  Q[2] = pr[2] - scale*dstDir[2];
  
  scale = DOT(srcDir,dstDir);            /* R <-- d - (d.v)v */
  R[0] = srcDir[0] - scale*dstDir[0];
  R[1] = srcDir[1] - scale*dstDir[1];
  R[2] = srcDir[2] - scale*dstDir[2];

  scale = DOT(R,R);                      /* t <-- R.Q/R.R */
  if (scale < EPSILON) scale = EPSILON;
  *srct = t = -DOT(R,Q)/scale;

  V[0] = Q[0] + t*R[0];                  /* dist2 <-- || Q + t*R || */
  V[1] = Q[1] + t*R[1];
  V[2] = Q[2] + t*R[2];
  return DOT(V,V);
}
