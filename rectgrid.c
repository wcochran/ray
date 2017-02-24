/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "raytrace.h"
#include "rectgrid.h"
#include "bbox.h"
#include "pnmio.h"

typedef struct {
  short tile;            /* tile image? */
  short interp;          /* use bilinear interpolation? */
  double org[2];         /* origin of image in xy-plane */
  double x[2], y[2];     /* image vertical and horizontal image */
  double width, height;  /* dimensions of mapped image */
  pnm_image *image;      /* mapped image */
} GRID_IMAGE;

typedef struct {
  int W, H;              /* width and height of grid */
  int rowStride;         /* linear distance between neighbor row values */
  double *z;             /* grid buffer stored in row-major order */
  BBOX bbox;             /* bounding box */
  double color[3];       /* rgb color of surface */
  GRID_IMAGE *image;     /* image mapped to surface (if != NULL) */
} GRID;

/*
 * color()
 * Method for retrieving solid color of grid object at hit point.
 */
static
void color(OBJECT *this,double hitPoint[3], HIT_INFO *info, double color[3]) {
  GRID *grid = (GRID *) this->data;
  color[0] = grid->color[0];
  color[1] = grid->color[1];
  color[2] = grid->color[2];
}

/*
 * setColor()
 * Method for setting solid color of grid object.
 */
static
void setColor(OBJECT *this, double color[3]) {
  GRID *grid = (GRID *) this->data;
  grid->color[0] = color[0];
  grid->color[1] = color[1];
  grid->color[2] = color[2];
}

static
double rayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
	      HIT_INFO *hitInfo) {
  GRID *grid = (GRID *) this->data;
  double trange[2];
  int r,c, W,H;
  double xorg,yorg, dx,dy;
  int zpitch;
  double *z;
  double x,y;
  int rin,cin, rout,cout;

  /*
   * If ray misses surface bounding box then return "miss".
   */
  if (!rayHitsBoundingBox(&grid->bbox, rayOrg, rayDir, trange))
    return -1.0;

  /*
   * Cache widely used grid data.
   */
  W = grid->W, H = grid->H;
  xorg = grid->bbox.min.x;
  yorg = grid->bbox.min.y;
  dx = (grid->bbox.max.x - xorg)/(W-1);
  dy = (grid->bbox.max.y - yorg)/(H-1);
  zpitch = grid->rowStride;
  z = grid->z;

#define Z(row,col) z[(row)*zpitch + (col)]
  
  /*
   * Find cell indices corresponding to where the ray
   * enters and exits the surfaces bounding box.
   *   (rin,cin) : grid row and column where ray enters bounding box;
   *   (rout,cout) : grid row and column where ray exits bounding box;
   */
  x = rayOrg[0] + trange[0]*rayDir[0];
  y = rayOrg[1] + trange[0]*rayDir[1];
  rin = (int) ((y - yorg)/dy);
  cin = (int) ((x - xorg)/dx);
  if (rin < 0) rin = 0;
  else if (rin >= H-1) rin = H-2;
  if (cin < 0) cin = 0;
  else if (cin >= W-1) cin = W-2;

  x = rayOrg[0] + trange[1]*rayDir[0];
  y = rayOrg[1] + trange[1]*rayDir[1];
  rout = (int) ((y - yorg)/dy);
  cout = (int) ((x - xorg)/dx);
  if (rout < 0) rout = 0;
  else if (rout >= H-1) rout = H-2;
  if (cout < 0) cout = 0;
  else if (cout >= W-1) cout = W-2;

  /*
   * Starting with cell (rin,cin), we march through all the cells
   * that are potentially penetrated by the ray. The first such cell
   * whose surface is penetrated by the ray yields the closest intersection.
   * We determine the next cell to try based on the direction of the 
   * ray and the relationship between the appropriate cell corner 
   * point and the ray (see below).
   *
   * To compute the actually ray/surface intersection we determine
   * the coefficients of the bilinear function f(x,y) that interpolate
   * the corner points, and then plug in the appropriate terms from
   * the parametric equation for the ray and solve the resulting
   * quadratic for t. We adjust the input ray's xy-origin and the cell's
   * xy-origin to be (0,0), so the cell support look like this:
   *
   *   dy *------------*
   *      |            |   f(x,y) = c00 + c01*x + c10*y + c11*x*y
   *      |            |
   *      |            |     z00 = f(0,0)
   *      |            |     z01 = f(dx,0)
   *      |            |     z10 = f(0,dy)
   *    0 *------------*     z11 = f(dx,dy)
   *      0            dx
   *
   *  We can then express the bilinear surface as the 
   *  function f(x,y) (see above). If we label the four corner
   *  z-values as above, we can solve for the four coefficients:
   *
   *           c00 = z00
   *           c01 = (z01 - z00)/dx
   *           c10 = (z10 - z00)/dy
   *           c11 = ((z11 - z01) - (z10 - z00))/(dx*dy)
   *
   *  Let (Ox, Oy, Oz) = (rayOrg[0] - x[0], rayOrg[1] - y[0], rayOrg[2]),
   *      (Dx, Dy, Dx) = (rayDir[0], rayDir[1], rayDir[1]),
   *  which gives us our "shifted" parametric ray equation:
   *
   *      R(t) = (Dx, Dy, Dz)*t + (Ox, Oy, Oz).
   *
   *  Let's plug each component of R(t) into f(x,y) which yields:
   *     
   *      Dz*t + Oz = c00 + c01*(Dx*t + Ox) + 
   *                    c10*(Dy*t + Oy) + c11*(Dx*t + Ox)*(Dy*t + Oy)
   *
   *  We rewrite this as a quadratic in t: A*t^2 + B*t + C = 0, where
   *
   *        A = c11*Dx*Dy,
   *        B = c01*Dx + c10*Dy + c11*(Ox*Dy + Oy*Dx) - Dz,
   *        C = c00 + c01*Ox + c10*Oy + c11*Ox*Oy - Oz.
   *  
   *  If we have complex roots, then the ray does not intersect
   *  the surface anyway. If we have 1 or 2 real roots, then we have
   *  to see if the corresponding hit points are over f's support.
   *  If both of the roots correspond to hit points over f's support,
   *  we use the smallest (non-negative) root.
   */
  for (r = rin, c = cin; ; ) {
    double x0,y0, x1,y1;
    double z00, z01, z10, z11;
    double c00, c01, c10, c11;
    double Ox, Oy;
    double A, B, C;
    int roots;
    double desc, t0,t1, xhit,yhit;

    x0 = xorg + c*dx;  /* origin of cell */
    y0 = yorg + r*dy;

    x1 = x0 + dx;      /* opposite corner of cell */
    y1 = y0 + dy;

    z00 = Z(r,c);      /* z-values at corners */
    z01 = Z(r,c+1);
    z10 = Z(r+1,c);
    z11 = Z(r+1,c+1);
    
    c00 = z00;                    /* coefficients for bilinear f(x,y) */
    c01 = (z01 - z00)/dx;
    c10 = (z10 - z00)/dy;
    c11 = ((z11 + z00) - (z10 + z01))/(dx*dy);

    Ox = rayOrg[0] - x0;          /* shifted ray xy-origin */
    Oy = rayOrg[1] - y0;

    A = c11*rayDir[0]*rayDir[1];  /* coefficients of quadratic */
    B = c01*rayDir[0] + c10*rayDir[1] + 
      c11*(Ox*rayDir[1] + Oy*rayDir[0]) - rayDir[2];
    C = c00 + c01*Ox + c10*Oy + c11*Ox*Oy - rayOrg[2];

    /*
     * Most people, including me, don't know the right way to
     * solve the quadratic equation. Section 5.5 in Numeric Recipes
     * set me straight. The "high school" solution is
     *                    
     *           -B +/- sqrt(B*B - 4*A*C)
     *       t = ------------------------,
     *                    2*A
     *
     * but if A*C << B*B then one of the roots will involve the
     * subtraction of B from a nearly equal quantity -- you'll just
     * get noise for the root! Here is the approved computation:
     *
     *       Let Q = (-1/2)(B + sign(B)*sqrt(B*B - 4*A*C)),
     *    
     *       t0 = Q/A,  t1 = C/Q.
     */

    roots = 0;

    if (A == 0.0) {
      if (B != 0.0) {
	t0 = -C/B;      
	roots = 1;
      }
    } else if ((desc = B*B - 4*A*C) == 0.0) {
      t0 = -B/(2*A);
      roots = 1;
    } else if (desc > 0.0) {
#define STABLE
#ifdef STABLE
      double Q;
      desc = sqrt(desc);
      if (B < 0) desc = -desc;
      Q = -0.5*(B + desc);
      t0 = Q/A;
      t1 = C/Q;
#else
      desc = sqrt(desc);
      t0 = (-B + desc)/(2*A);
      t1 = (-B - desc)/(2*A);
#endif
      if (t0 > t1) {              /* make sure t's are in order */
	double tmp = t0;
	t0 = t1;
	t1 = tmp;
      }
      roots = 2;
    }

    if (roots > 0) {
      if (t0 >= EPSILON &&
	  (xhit = Ox + t0*rayDir[0]) >= 0.0 && xhit <= dx &&
	  (yhit = Oy + t0*rayDir[1]) >= 0.0 && yhit <= dy)
	return t0;                  /* closest root works! */      
      if (--roots > 0) {            	
	if (t1 >= EPSILON &&
	    (xhit = Ox + t1*rayDir[0]) >= 0.0 && xhit <= dx &&
	    (yhit = Oy + t1*rayDir[1]) >= 0.0 && yhit <= dy)
	  return t1;                  /* further root works! */
      }
    }
    
    /*
     * If that was our last cell to check then the ray
     * does not intersect the surface.
     */
    if (r == rout && c == cout)
      break;

    /*
     * The previous cell yielded no intersection, so let's
     * try the next cell. We are at cell (r,c) and we need to
     * make a decision on which is the next cell. There are
     * four cases based on the direction of the ray as shown
     * below; Each of these cases leaves two choices left.
     *
     *     rayDir[0] < 0         rayDir[0] > 0 
     *     rayDir[1] > 0         rayDir[1] > 0
     *             *-----*       *-----*
     *    (x0,y1)  |     |       |     |  (x1,y1)
     *           \ |r+1,c|       |r+1,c| /
     *       *-----*-----*       *-----*-----*
     *       |     |     |       |     |     |
     *       |r,c-1| r,c |       | r,c |r,c+1|
     *       *-----*-----*       *-----*-----*
     *
     *
     *     rayDir[0] < 0         rayDir[0] > 0
     *     rayDir[1] < 0         rayDir[1] < 0
     *       *-----*-----*       *-----*-----*
     *       |     |     |       |     |     |
     *       |r,c-1| r,c |       | r,c |r,c+1|
     *       *-----*-----*       *-----*-----*
     *           / |     |       |     | \
     *    (x0,y0)  |r-1,c|       |r-1,c|  (x1,y0)
     *             *-----*       *-----*
     *
     * In each subcase we need to choose between making a horizontal
     * or vertical "move." Note that in each of the four figures
     * above we label the corner that is going to help us make this
     * decision. Using the relative position of this point with the ray,
     * we determine the final move; We can determine if a point
     * (x,y) is "to the left" of the ray (dx,dy)*t + (ox, oy)
     * by examining the sign of the z-component of the 
     * cross prosspduct U x V, where U = (dx,dy) and V = (x-ox, y-oy):
     *
     *      (U x V) = dx*(y - oy) - dy*(x - ox) > 0
     *
     * For example, if rayDir[0] < 0 and rayDir[1] > 0 and
     * (x0, y1) is to the left of the ray, then we must move
     * to (r+1,c). If the cross product iz zero then we make
     * a diaganol move.
     */

    if (rayDir[0] > 0 && rayDir[1] > 0) {
      double cross = rayDir[0]*(y1 - rayOrg[1]) - rayDir[1]*(x1 - rayOrg[0]);
      if (cross > 0)
	c++;
      else if (cross < 0)
	r++;
      else
	c++, r++;
    } else if (rayDir[0] < 0 && rayDir[1] > 0) {
      double cross = rayDir[0]*(y1 - rayOrg[1]) - rayDir[1]*(x0 - rayOrg[0]);
      if (cross > 0)
	r++;
      else if (cross < 0)
	c--;
      else
	r++, c--;
    } else if (rayDir[0] < 0 && rayDir[1] < 0) {
      double cross = rayDir[0]*(y0 - rayOrg[1]) - rayDir[1]*(x0 - rayOrg[0]);
      if (cross > 0)
	c--;
      else if (cross < 0)
	r--;
      else
	c--, r--;
    } else {
      double cross = rayDir[0]*(y0 - rayOrg[1]) - rayDir[1]*(x1 - rayOrg[0]);
      if (cross > 0)
	r--;
      else if (cross < 0)
	c++;
      else
	r--, c++;
    }

    /*
     * Safety check: Let's not step off the end of the world!
     */
    if (r < 0 || r >= H-1 ||
	c < 0 || c >= W-1)
      break;    
  }

  return -1.0;  /* ray missed surface */
}
	
/*  z10                    (x,y,z): hit point
 *    *---------*          x0,x1: min/max x-values of cell
 *    |  (x,y,z)|          y0,y1: min/max y-values of cell
 *    *-----*   |          z00, z01, z10: z-values at corners of cell
 *    |     |   |          
 *    |     |   |          U = (x - x0, 0, z - lerp(z00, z10, (y-y0)/(y1-y0)))
 *    *-----*---*          V = (0, y - y0, z - lerp(z00, z01, (x-x0)/(x1-x0)))
 *  z00        z01         N = U x V
 */
static
void normal(OBJECT *this, double hitPoint[3], HIT_INFO *info, 
	    double normal[3]) {
  GRID *grid = (GRID *) this->data;
  double x0,y0, dx,dy;
  int r,c;
  double U[3], V[3], scale;
  double z[2][2];

  x0 = grid->bbox.min.x;
  y0 = grid->bbox.min.y;
  dx = (grid->bbox.max.x - x0)/(grid->W-1);
  dy = (grid->bbox.max.y - y0)/(grid->H-1);

  c = (int) ((hitPoint[0] - x0)/dx);
  if (c < 0) c = 0;
  else if (c >= grid->W-1) c = grid->W-2;
  r = (int) ((hitPoint[1] - y0)/dy);
  if (r < 0) r = 0;
  else if (r >= grid->H-1) r = grid->H-2;

  z[0][0] = grid->z[r*grid->rowStride + c];
  z[1][0] = grid->z[(r+1)*grid->rowStride + c];
  z[0][1] = grid->z[r*grid->rowStride + c+1];

  x0 += c*dx;
  y0 += r*dy;

  U[0] = hitPoint[0] - x0;
  U[1] = 0.0;
  U[2] = hitPoint[2] - (z[0][0] + (z[1][0] - z[0][0])*(hitPoint[1] - y0)/dy);

  V[0] = 0.0;
  V[1] = hitPoint[1] - y0;
  V[2] = hitPoint[2] - (z[0][0] + (z[0][1] - z[0][0])*(hitPoint[0] - x0)/dx);
  
  normal[0] = -U[2]*V[1];
  normal[1] = -U[0]*V[2];
  normal[2] =  U[0]*V[1];

  scale = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
  if (scale > 0.0) {
    scale = 1.0/sqrt(scale);
    normal[0] *= scale;
    normal[1] *= scale;
    normal[2] *= scale;
  } else {
    normal[0] = normal[1] = 0.0;  /* rare case: hit point on seam */
    normal[1] = 1.0;
  }
}

/*
 * Creates a surface z = f(x,y) that interpolates WxN knots (xi, yj, zij)
 * where zij = f(xi, yj), 0 <= i <= W-1, 0 <= j <= N-1. Bilinear
 * interpolation is used between four neighboring "pixels."
 * Note: We don't copy the given 2D z-array, but merely reference it,
 * so don't trash it -- if you modify the grid, creates a new object
 * (really a new bounding box needs to be computed).
 */
OBJECT *createSurfaceSupportedByRectangularGrid(int W, int H, 
						double orgx, double orgy,
						double width, double height,
						double *z, int rowStride) {
  GRID *grid;
  OBJECT *obj;
  int r,c;

  if ((grid = (GRID *) malloc(sizeof(GRID))) == NULL ||
      (obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("allocating rectaangular grid object");
    exit(-1);
  }

  grid->W = W, grid->H = H;
  grid->rowStride = rowStride;
  grid->z = z;
  grid->bbox.min.x = orgx;
  grid->bbox.min.y = orgy;
  grid->bbox.max.x = orgx + width;
  grid->bbox.max.y = orgy + height;

  grid->bbox.min.z = grid->bbox.max.z = z[0];

  for (r = 0; r < H; r++)
    for (c = 0; c < W; c++) {
      double zz = z[r*rowStride + c];
      if (zz < grid->bbox.min.z)
	grid->bbox.min.z = zz;
      else if (zz > grid->bbox.max.z)
	grid->bbox.max.z = zz;
    }

  grid->color[0] = 0.2;   /* nice blue surface> */
  grid->color[1] = 0.2;
  grid->color[2] = 0.6;

  grid->image = NULL;

  obj->data = grid;       /* hook grid data to ray trace object */
  
  obj->ka = 0.15;         /* some nice material coefficients */
  obj->kd = 0.40;
  obj->ks = 0.45;
  obj->kt = 0.55;
  obj->ni = 1.52;
  obj->phong = 6.0;

  obj->rayHit = rayHit;   /* set object's polymorphic method operators */
  obj->normal = normal;
  obj->color = color;
  obj->setColor = setColor;

  return obj;
}

static
void imageColor(OBJECT *this,double hit[3], HIT_INFO *info, double color[3]) {
  GRID *grid = (GRID *) this->data;
  GRID_IMAGE *image = grid->image;
  pnm_image *img = image->image;
  double p[2], x,y, scale;
  int i,j;
  
  /*
   * Adjust hit point to the image's origin.
   * We ignore the z-coordinate.
   */
  p[0] = hit[0] - image->org[0];
  p[1] = hit[1] - image->org[1];

  /*
   * Find the hit point in the image's coordinate system via projection.
   */
  x = p[0]*image->x[0] + p[1]*image->x[1];
  y = p[0]*image->y[0] + p[1]*image->y[1];

  /*
   * If we are tiling the image and the hit point is outside
   * the image region we'll just return the plane's normal color.
   */
  if (!image->tile && (x < 0.0 || y < 0.0 || 
		       x > image->width || y > image->height)) {
    color[0] = grid->color[0];
    color[1] = grid->color[1];
    color[2] = grid->color[2];
    return;
  }
  
  /*
   * Scale x,y to unit box.
   */
  x /= image->width;
  y /= image->height;

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
  /* XXX
  y = 1.0 - y;
  */
  
  /*
   * Find corner pixel address.
   */
  x *= PNM_NC(img) - 1;
  y *= PNM_NR(img) - 1;
  j = (int) x;
  i = (int) y;

  /*
   * Compute component normalization factor.
   */
  scale = 1.0/PNM_MAXVAL(img);
  
  if (image->interp) { /* use bilinear interpolation */
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

    color[0] = scale*(w[0]*PPM_PIXEL_R(img,i,j) +
		      w[1]*PPM_PIXEL_R(img,i,j+1) +
		      w[2]*PPM_PIXEL_R(img,i+1,j) +
		      w[3]*PPM_PIXEL_R(img,i+1,j+1));
    color[1] = scale*(w[0]*PPM_PIXEL_G(img,i,j) +
		      w[1]*PPM_PIXEL_G(img,i,j+1) +
		      w[2]*PPM_PIXEL_G(img,i+1,j) +
		      w[3]*PPM_PIXEL_G(img,i+1,j+1));
    color[2] = scale*(w[0]*PPM_PIXEL_B(img,i,j) +
		      w[1]*PPM_PIXEL_B(img,i,j+1) +
		      w[2]*PPM_PIXEL_B(img,i+1,j) +
		      w[3]*PPM_PIXEL_B(img,i+1,j+1));
  } else {  /* use closest pixel color */
    color[0] = scale*(PPM_PIXEL_R(img,i,j));
    color[1] = scale*(PPM_PIXEL_G(img,i,j));
    color[2] = scale*(PPM_PIXEL_B(img,i,j));
  }
}
  
void mapImageToSurfaceSupportedByRectangularGrid(OBJECT *obj,
						 int tile,
						 int interpolate,
						 double xyorg[2],
						 double yaxis[2],
						 double width, double height,
						 pnm_image *image) {
  GRID *grid = (GRID *) obj->data;
  GRID_IMAGE *gridImage;
  double mag;
  
  /*
   * Allocate memory for image info if not allocated already.
   */
  if (grid->image == NULL &&
      (grid->image = (GRID_IMAGE *) malloc(sizeof(GRID_IMAGE))) == NULL) {
    perror("mapImageToSurfaceSupportedByRectangularGrid:malloc()");
    exit(-1);
  }

  gridImage = grid->image;

  gridImage->tile = tile;
  gridImage->interp = interpolate;
  gridImage->org[0] = xyorg[0];
  gridImage->org[1] = xyorg[1];
  mag = sqrt(yaxis[0]*yaxis[0] + yaxis[1]*yaxis[1]);
  if (mag > 0.0)
    gridImage->y[0] = yaxis[0]/mag, gridImage->y[1] = yaxis[1]/mag;
  else
    gridImage->y[0] = 1.0, gridImage->y[1] = 0.0; /* caller a dope! */
  gridImage->x[0] = gridImage->y[1];
  gridImage->x[1] = -gridImage->y[0];
  gridImage->width = width;
  gridImage->height = height;
  gridImage->image = image;

  obj->color = imageColor;
}
