
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

/*
 * The algorithm presented here was extracted from
 *       Nishita, Sederberg, Kakimoto
 *       "Ray Tracing Trimmed Rational Surface Patches"
 *       Computer Graphics (SigGraph Proceedings) August 1990
 *       Volume 24, Number 4, pp 337-345
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "bezier3.h"

/* 
 * Maximum ray intersection with surface; there is no known
 * theoretical bound.
 */
#define MAX_ZEROS 10

/*
 * If a bezier clip fails to reduce the parameter interval
 * by 20% then we split the patch in half.
 */
#define SPLIT_TOLERANCE 0.8

/*
 * Given four points that represent the upper or lower boundary
 * of a 2D polygon that is initially evenly spaced in x going from
 * left to right, we remove either of the middle points (or both)
 * that cause the polygon to not be convex.
 */
static
int forceConvexHull(POINT2 hull[4], int upper) {
  int n;
  double cross;

  /*
   * Test vertex hull[1] to see it causes a concavity.
   * The cross product simplifies since we know the points
   * are evenly space in x.
   */
  n = 4;
  cross = (hull[1].y - hull[2].y) + (hull[1].y - hull[0].y);
  if (upper) cross = -cross;
  if (cross >= 0) {
    n--;
    hull[1] = hull[2];
    hull[2] = hull[3];
  }

  /*
   * Test vertex hull[n-2] to see it causes a concavity.
   * If it does, then we need to go back and see if hull[1] now
   * causes a concavity if we didn't remove it already.
   */
  cross = (hull[n-3].x - hull[n-2].x)*(hull[n-1].y - hull[n-2].y) -
          (hull[n-1].x - hull[n-2].x)*(hull[n-3].y - hull[n-2].y);
  if (upper) cross = -cross;
  if (cross >= 0) {
    n--;
    hull[n-1] = hull[n];
    if (n > 2) {
      cross = (hull[0].x - hull[1].x)*(hull[2].y - hull[1].y) -
	      (hull[2].x - hull[1].x)*(hull[0].y - hull[1].y);
      if (upper) cross = -cross;
      if (cross >= 0) {
	n = 2;
	hull[1] = hull[2];
      }
    }
  }

  return n;
}

static
void splitInHalf(POINT2 P[4][4],     /* 2D Bezier Patch to split */
		 int split_s,        /* split along s or t? */
		 POINT2 L[4][4],     /* left half */
		 POINT2 R[4][4]) {   /* right half */
  int i;
  POINT2 P01, P12, P23;
  POINT2 P012, P123;
  POINT2 P0123;

  if (split_s) {
    for (i = 0; i <= 3; i++) {
      P01.x = 0.5*(P[0][i].x + P[1][i].x);
      P01.y = 0.5*(P[0][i].y + P[1][i].y);
      P12.x = 0.5*(P[1][i].x + P[2][i].x);
      P12.y = 0.5*(P[1][i].y + P[2][i].y);
      P23.x = 0.5*(P[2][i].x + P[3][i].x);
      P23.y = 0.5*(P[2][i].y + P[3][i].y);

      P012.x = 0.5*(P01.x + P12.x);
      P012.y = 0.5*(P01.y + P12.y);
      P123.x = 0.5*(P12.x + P23.x);
      P123.y = 0.5*(P12.y + P23.y);
      P0123.x = 0.5*(P012.x + P123.x);
      P0123.y = 0.5*(P012.y + P123.y);

      L[0][i].x = P[0][i].x;
      L[0][i].y = P[0][i].y;
      L[1][i].x = P01.x;
      L[1][i].y = P01.y;
      L[2][i].x = P012.x;
      L[2][i].y = P012.y;
      L[3][i].x = R[0][i].x = P0123.x;
      L[3][i].y = R[0][i].y = P0123.y;
      R[1][i].x = P123.x;
      R[1][i].y = P123.y;
      R[2][i].x = P23.x;
      R[2][i].y = P23.y;
      R[3][i].x = P[3][i].x;
      R[3][i].y = P[3][i].y;      
    }
  } else {
    for (i = 0; i <= 3; i++) {
      P01.x = 0.5*(P[i][0].x + P[i][1].x);
      P01.y = 0.5*(P[i][0].y + P[i][1].y);
      P12.x = 0.5*(P[i][1].x + P[i][2].x);
      P12.y = 0.5*(P[i][1].y + P[i][2].y);
      P23.x = 0.5*(P[i][2].x + P[i][3].x);
      P23.y = 0.5*(P[i][2].y + P[i][3].y);

      P012.x = 0.5*(P01.x + P12.x);
      P012.y = 0.5*(P01.y + P12.y);
      P123.x = 0.5*(P12.x + P23.x);
      P123.y = 0.5*(P12.y + P23.y);
      P0123.x = 0.5*(P012.x + P123.x);
      P0123.y = 0.5*(P012.y + P123.y);

      L[i][0].x = P[i][0].x;
      L[i][0].y = P[i][0].y;
      L[i][1].x = P01.x;
      L[i][1].y = P01.y;
      L[i][2].x = P012.x;
      L[i][2].y = P012.y;
      L[i][3].x = R[i][0].x = P0123.x;
      L[i][3].y = R[i][0].y = P0123.y;
      R[i][1].x = P123.x;
      R[i][1].y = P123.y;
      R[i][2].x = P23.x;
      R[i][2].y = P23.y;
      R[i][3].x = P[i][3].x;
      R[i][3].y = P[i][3].y;      
    }
  }
}

static
void split(POINT2 P[4][4],     /* 2D Bezier Patch to split */
	   int split_s,        /* split along s or t? */
	   double alpha,       /* split value for t or s */
	   POINT2 L[4][4],     /* left half */
	   POINT2 R[4][4]) {   /* right half */
  int i;
  double beta = 1.0 - alpha;
  POINT2 P01, P12, P23;
  POINT2 P012, P123;
  POINT2 P0123;

  if (split_s) {
    for (i = 0; i <= 3; i++) {
      P01.x = beta*P[0][i].x + alpha*P[1][i].x;
      P01.y = beta*P[0][i].y + alpha*P[1][i].y;
      P12.x = beta*P[1][i].x + alpha*P[2][i].x;
      P12.y = beta*P[1][i].y + alpha*P[2][i].y;
      P23.x = beta*P[2][i].x + alpha*P[3][i].x;
      P23.y = beta*P[2][i].y + alpha*P[3][i].y;

      P012.x = beta*P01.x + alpha*P12.x;
      P012.y = beta*P01.y + alpha*P12.y;
      P123.x = beta*P12.x + alpha*P23.x;
      P123.y = beta*P12.y + alpha*P23.y;
      P0123.x = beta*P012.x + alpha*P123.x;
      P0123.y = beta*P012.y + alpha*P123.y;

      L[0][i].x = P[0][i].x;
      L[0][i].y = P[0][i].y;
      L[1][i].x = P01.x;
      L[1][i].y = P01.y;
      L[2][i].x = P012.x;
      L[2][i].y = P012.y;
      L[3][i].x = R[0][i].x = P0123.x;
      L[3][i].y = R[0][i].y = P0123.y;
      R[1][i].x = P123.x;
      R[1][i].y = P123.y;
      R[2][i].x = P23.x;
      R[2][i].y = P23.y;
      R[3][i].x = P[3][i].x;
      R[3][i].y = P[3][i].y;      
    }
  } else {
    for (i = 0; i <= 3; i++) {
      P01.x = beta*P[i][0].x + alpha*P[i][1].x;
      P01.y = beta*P[i][0].y + alpha*P[i][1].y;
      P12.x = beta*P[i][1].x + alpha*P[i][2].x;
      P12.y = beta*P[i][1].y + alpha*P[i][2].y;
      P23.x = beta*P[i][2].x + alpha*P[i][3].x;
      P23.y = beta*P[i][2].y + alpha*P[i][3].y;

      P012.x = beta*P01.x + alpha*P12.x;
      P012.y = beta*P01.y + alpha*P12.y;
      P123.x = beta*P12.x + alpha*P23.x;
      P123.y = beta*P12.y + alpha*P23.y;
      P0123.x = beta*P012.x + alpha*P123.x;
      P0123.y = beta*P012.y + alpha*P123.y;

      L[i][0].x = P[i][0].x;
      L[i][0].y = P[i][0].y;
      L[i][1].x = P01.x;
      L[i][1].y = P01.y;
      L[i][2].x = P012.x;
      L[i][2].y = P012.y;
      L[i][3].x = R[i][0].x = P0123.x;
      L[i][3].y = R[i][0].y = P0123.y;
      R[i][1].x = P123.x;
      R[i][1].y = P123.y;
      R[i][2].x = P23.x;
      R[i][2].y = P23.y;
      R[i][3].x = P[i][3].x;
      R[i][3].y = P[i][3].y;      
    }
  }
}

static
void clip(POINT2 P[4][4], int clip_s, double min, double max) {
  POINT2 L[4][4], R[4][4];
  int i,j;

  /*
   * Split P along current parameterization using min 
   * (if min > 0) into left and right patches. Discard
   * the left patch.
   */
  if (min > 0.0) {
    split(P, clip_s, min, L, R);
  } else
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
	R[i][j] = P[i][j];
  
  /*
   * Now split R using max. Since we have already adjusted
   * The size of the interval we need to compute what the
   * new max should be for the already trimmed interval.
   */
  if (max < 1.0)
    split(R, clip_s, (max - min)/(1.0 - min), P,L);
  else
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
	P[i][j] = R[i][j];
}
  

/*
 * P[i][j]
 *   0 <= i <= 3 : s parameter control points;
 *   0 <= j <= 3 : t parameter control points;
 */
static
int findZeros(POINT2 P[4][4],       /* 2D projected control points */
	      double s[2],          /* current range on s parameter */
	      double t[2],          /* current range on t parameter */
	      int clip_s,           /* are we clipping in s or t */
	      double tolerance,     /* tells us when to quit subdividing */
	      double zeros[][2])  { /* return array of s,t intersections */
  POINT2 line;
  POINT2 L[4][4], R[4][4];
  double scale;
  double D[4][4];
  int upperCount, lowerCount;
  POINT2 upperHull[4], lowerHull[4];
  static double param[4] = {0.0, 1.0/3.0, 2.0/3.0, 1.0};
  int i,j;
  double minx, maxx;
  double smin,smax, tmin,tmax;
  int numSplits;
  int extraClips = 0;

  do {
    /*
     * Create a line L that is parallel to V0+V1 and passes through
     * the origin where V0 and V1 are the vectors P[0][3] - P[0][0]
     * and P[3][3] - P[3][0] if we are clipping along the s axis
     * otherwise V0 and V1 are the vectors P[3][0] - P[0][0] and
     * P[3][3] - P[0][3] when we are clipping along the t axis.
     * The line is implicitly defined as
     *      D(x,y) = L.x*x + L.y*y
     * where L.x^2 + L.y^2 = 1 so that D(x,y) represents the distance
     * of each control point to the line.
     */
    if (clip_s) {
      line.y = (P[0][3].x - P[0][0].x) + (P[3][3].x - P[3][0].x);
      line.x = (P[0][0].y - P[0][3].y) + (P[3][0].y - P[3][3].y);
    } else {
      line.y = (P[3][0].x - P[0][0].x) + (P[3][3].x - P[0][3].x);
      line.x = (P[0][0].y - P[3][0].y) + (P[0][3].y - P[3][3].y);
    }
    scale = line.x*line.x + line.y*line.y;
    if (scale == 0.0) goto split_in_half;  /* extremely rare case */
    scale = 1.0/sqrt(scale);
    line.x *= scale;
    line.y *= scale;
    
    /*
     * Create 4x4 array of D[i][j] control points where the z-value
     * each control point (s,t,D[i][j]) is explicilty represented
     * in the array and represents the distance to the line L
     * and s = i/3, t = j/3.
     */
    for (i = 0; i <= 3; i++)
      for (j = 0; j <= 3; j++)
	D[i][j] = line.x*P[i][j].x + line.y*P[i][j].y;
    
    /*
     * Now we need to compute the convex hull of the D(s,t) patch
     * as it projected onto either the t=0 or s=0 plane.
     * We will do this by first determining the minimum and maximum
     * of each s value (t=0 projection) or each t value (s=0 projection)
     * which will give us initial upper and lower hull values. 
     */
    if (clip_s)  /* t=0 projection */
      for (i = 0; i <= 3; i++) {  /* loop through each s value */
	double min, max; 
	min = max = D[i][0];
	for (j = 1; j <= 3; j++)
	  if (D[i][j] < min) 
	    min = D[i][j];
	  else if (D[i][j] > max)
	    max = D[i][j];
	upperHull[i].x = lowerHull[i].x = param[i];
	upperHull[i].y = max;
	lowerHull[i].y = min;
      }
    else   /* s=0 projection */
      for (j = 0; j <= 3; j++) {  /* loop through each t value */
	double min, max; 
	min = max = D[0][j];
	for (i = 1; i <= 3; i++)
	  if (D[i][j] < min) 
	    min = D[i][j];
	  else if (D[i][j] > max)
	    max = D[i][j];
	upperHull[j].x = lowerHull[j].x = param[j];
	upperHull[j].y = max;
	lowerHull[j].y = min;
      }      
    
    /*
     * Now we scan the boundary edges of the convex hull to get 
     * a range where the projected convex hull intersects the x-axis.
     * The convex hull will intersect the x-axis either 0 times
     * or exactly 2 times. We'll first check the left and right
     * boundaries of the convex hull which are vertical edges.
     * If the left vertical boundary intersects the x-axis we'll set the
     * minx value to zero once and for all; If it doesn't we'll
     * set it to a value > 1 (we'll set it to 1). If the right
     * vertical boundary intersects the x-axis then we'll set maxx to 1.
     * (otherwise we'll set it to a value less than 0). 
     * If both the left and right vertical boundaries intersect 
     * the x-axis then we are done.
     */
    numSplits = 0;
    minx = 1.1; maxx = -0.1;
    if (upperHull[0].y > 0.0 && lowerHull[0].y < 0.0) {
      numSplits++;
      minx = 0.0;
    }
    if (upperHull[3].y > 0.0 && lowerHull[3].y < 0.0) {
      numSplits++;
      maxx = 1.0;
    }
    
    /*
     *   Now we scan the boundary edges of the convex hull to get 
     * a range where the projected convex hull intersects the x-axis.
     * We'll first scan the upper hull edges assuming we have not
     * found both splits yet.
     */
    if (numSplits < 2) {
      
      /*
       *   Now we need to remove the portions of the convex hull that
       * result in concavities to make sure its convex. We didn't
       * do this earlier since we knew the first and last points
       * of both the upper and lower hull arrays were 
       * part of the convex hull. 
       */
      upperCount = forceConvexHull(upperHull, 1);
      
      for (i = 1; i < upperCount; i++) {
	double y0 = upperHull[i-1].y, y1 = upperHull[i].y;
	double x0, x1;
	double x, dx,dy;

	/*
	 * Check to see if the edge is completely on one side
	 * of the x-axis.
	 */
	if ((y0 < 0.0 && y1 < 0.0) || (y0 > 0.0 && y1 > 0.0))
	  continue;
	
	/*
         * Check for horizontal edge on the x-axis.
	 * If this is the case then we have found our intersections.
	 */
	if ((dy = y1 - y0) == 0.0) {
	  minx = upperHull[i-1].x;
	  maxx = upperHull[i].x;
	  numSplits = 2;
	  break;
	}
	
	/*
	 * Find x value of edge y=0 intersection.
	 */
	x0 = upperHull[i-1].x; x1 = upperHull[i].x;
	dx = x1 - x0;
	x = x0 - dx*y0/dy;
	
	/*
	 * Update boundary values.
	 */
	if (x < minx) minx = x;
	if (x > maxx) maxx = x;
	if (++numSplits == 2) break;  /* could use goto here (nah!) */
      }
      
      /*
       * If we still have not found 2 intersections
       * then we'll scan the lower hull.
       */
      if (numSplits < 2) {
	
	/*
	 * Make sure the lower hull values are convex.
	 */
	lowerCount = forceConvexHull(lowerHull, 0);
	
	for (i = 1; i < lowerCount; i++) {
	  double y0 = lowerHull[i-1].y, y1 = lowerHull[i].y;
	  double x0, x1;
	  double x, dx,dy;
	  
	  /*
	   * Check to see if the edge is completely on one side
	   * of the x-axis.
	   */
	  if ((y0 < 0.0 && y1 < 0.0) || (y0 > 0.0 && y1 > 0.0))
	    continue;
	  
	  /*
	   * Check for horizontal edge on the x-axis.
	   * If this is the case then we have found our intersections.
	   */
	  if ((dy = y1 - y0) == 0.0) {
	    minx = lowerHull[i-1].x;
	    maxx = lowerHull[i].x;
	    numSplits = 2;
	    break;
	  }

	  /*
	   * Find x value of edge y=0 intersection.
	   */
	  x0 = lowerHull[i-1].x; x1 = lowerHull[i].x;
	  dx = x1 - x0;
	  x = x0 - dx*y0/dy;
	  
	  /*
	   * Update boundary values.
	   */
	  if (x < minx) minx = x;
	  if (x > maxx) maxx = x;
	  if (++numSplits == 2) break;
	}
      }
    }
    
    /*
     * No intersection: We are outta here!
     */
    if (numSplits == 0)
      return 0;
    
    /*
     * Possible multiple intersections.
     *   Split the patch in half along the current parameterization.
     *   Find the number of zeros in each half and merge the result.
     */
    if (maxx - minx > SPLIT_TOLERANCE) {
      int numLeft, numRight;
      double leftZeros[MAX_ZEROS][2], rightZeros[MAX_ZEROS][2];
      double srange[2], trange[2];
      
    split_in_half:
      splitInHalf(P, clip_s, L,R);
      
      if (clip_s) {
	double mids = 0.5*(s[0] + s[1]);
	
	srange[0] = s[0]; srange[1] = mids;
	trange[0] = t[0]; trange[1] = t[1];
	numLeft = findZeros(L,srange,trange,0,tolerance,leftZeros);
	
	srange[0] = mids; srange[1] = s[1];
	trange[0] = t[0]; trange[1] = t[1];
	numRight = findZeros(R,srange,trange,0,tolerance,rightZeros);
      } else {
	double midt = 0.5*(t[0] + t[1]);
	double trange[2];
	
	srange[0] = s[0]; srange[1] = s[1];
	trange[0] = t[0]; trange[1] = midt;
	numLeft = findZeros(L,srange,trange,1,tolerance,leftZeros);
	
	srange[0] = s[0]; srange[1] = s[1];
	trange[0] = midt; trange[1] = t[1];
	numRight = findZeros(R,srange,trange,1,tolerance,rightZeros);
      }
      
      for (i = 0; i < numLeft; i++) {
	zeros[i][0] = leftZeros[i][0];
	zeros[i][1] = leftZeros[i][1];
      }
      
      for (i = 0; i < numRight; i++) {
	zeros[numLeft+i][0] = rightZeros[i][0];
	zeros[numLeft+i][1] = rightZeros[i][1];
      }
      
      return numLeft + numRight;
    }
    
    /*
     * Compute the current global bounds on s & t.
     */
    if (clip_s) {
      double ds = s[1] - s[0];
      smin = s[0] + minx*ds;
      smax = s[0] + maxx*ds;
      tmin = t[0];
      tmax = t[1];
    } else {
      double dt = t[1] - t[0];
      smin = s[0];
      smax = s[1];
      tmin = t[0] + minx*dt;
      tmax = t[0] + maxx*dt;
    }
    
    /*
     * If both parameter range differences are within tolerance
     * and we are not on the boundary of the *original* patch then we have
     * found a zero point! If we are on the boundary and we have
     * performed at least two extra clips then we'll assume
     * we have made contact with the surface.
     */
    if ((smax - smin) < tolerance && (tmax - tmin) < tolerance) {
      if ((smin != 0.0 && smax != 1.0 && tmin != 0.0 && tmax != 1.0)
	  || extraClips > 1) {

	/*
	 * Estimate s. We can do a better job of estimating s and t
	 * by assuming that the patch is roughly planar and rectangular
	 * after all the clipping.
	 */
	line.y = (P[0][3].x - P[0][0].x) + (P[3][3].x - P[3][0].x);
	line.x = (P[0][0].y - P[0][3].y) + (P[3][0].y - P[3][3].y);
	scale = line.x*line.x + line.y*line.y;
	if (scale > 0.0) {
	  scale = 1.0/sqrt(scale);
	  line.x *= scale;
	  line.y *= scale;
	  smin = line.x*P[0][0].x + line.y*P[0][0].y;
	  smax = line.x*P[3][0].x + line.y*P[3][0].y;
	  zeros[0][0] = s[0] - smin*(s[1] - s[0])/(smax - smin);
	} else {
	  zeros[0][0] = 0.5*(s[0] + s[1]);
	}
	  
	/*
	 * Estimate t.
	 */
	line.y = (P[3][0].x - P[0][0].x) + (P[3][3].x - P[0][3].x);
	line.x = (P[0][0].y - P[3][0].y) + (P[0][3].y - P[3][3].y);
	scale = line.x*line.x + line.y*line.y;
	if (scale > 0.0) {
	  scale = 1.0/sqrt(scale);
	  line.x *= scale;
	  line.y *= scale;
	  tmin = line.x*P[0][0].x + line.y*P[0][0].y;
	  tmax = line.x*P[0][3].x + line.y*P[0][3].y;
	  zeros[0][1] = t[0] - tmin*(t[1] - t[0])/(tmax - tmin);
	} else {
	  zeros[0][1] = 0.5*(t[0] + t[1]);
	}

	return 1;
      }
      extraClips++;
    }

    /*
     * Here we clip away the portion of the surface that we know
     * can *not* contain a zero.
     * We make the interval slightly wider to avoid numerical
     * roundoff problems and we adjust the global bounds
     * accordingly. Then we tail recurse -- or loop back to the top.
     */
    minx *= 0.99; 
    maxx = maxx*0.99 + 0.01;
    if (clip_s) {
      double ds = s[1] - s[0];
      s[1] = s[0] + maxx*ds;
      s[0] += minx*ds;
    } else {
      double dt = t[1] - t[0];
      t[1] = t[0] + maxx*dt;
      t[0] += minx*dt;
    }
    clip(P, clip_s, minx,maxx);
    clip_s = !clip_s;
  } while(1);
}

static
int findAllZeros(double rayOrg[3], double rayDir[3],
		 POINT3 P[4][4], double st_tolerance,
		 double zeros[][2]) {
  POINT2 P2[4][4];
  double N1[3], e1;
  double N2[3], e2;
  double scale;
  int i,j;
  double s[2], t[2];

  /*
   * Contstruct a plane that the given ray lies in.
   * Find a normal vector, unitize it, and find the plane's
   * distance from the origin. Here we assume that rayDir
   * is normalized.
   */
  if (fabs(rayDir[0]) > 0.5 || fabs(rayDir[1]) > 0.5) {
    N1[0] = rayDir[1]; N1[1] = -rayDir[0]; N1[2] = 0.0;
    scale = 1.0/sqrt(N1[0]*N1[0] + N1[1]*N1[1]);
    N1[0] *= scale;
    N1[1] *= scale;
    e1 = -(N1[0]*rayOrg[0] + N1[1]*rayOrg[1]);
  } else {
    N1[0] = rayDir[2]; N1[1] = 0.0; N1[2] = -rayDir[0];
    scale = 1.0/sqrt(N1[0]*N1[0] + N1[2]*N1[2]);
    N1[0] *= scale;
    N1[2] *= scale;
    e1 = -(N1[0]*rayOrg[0] + N1[2]*rayOrg[2]);
  }

  /*
   * Construct another such plane that is perpendicular to 
   * the one above.
   */
  N2[0] = rayDir[1]*N1[2] - rayDir[2]*N1[1];
  N2[1] = rayDir[2]*N1[0] - rayDir[0]*N1[2];
  N2[2] = rayDir[0]*N1[1] - rayDir[1]*N1[0];
  e2 = -(rayOrg[0]*N2[0] + rayOrg[1]*N2[1] + rayOrg[2]*N2[2]);

  /*
   * Now use these two planes to form the projection of the
   * input Bezier patch to 2D. The projected control point x
   * values are the distance to the first plane; the y values
   * are the distance to the second plane.
   */
  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++) {
      P2[i][j].x = P[i][j].x*N1[0] + P[i][j].y*N1[1] + P[i][j].z*N1[2] + e1;
      P2[i][j].y = P[i][j].x*N2[0] + P[i][j].y*N2[1] + P[i][j].z*N2[2] + e2;
    }

  /*
   * Now we just need to find the (s,t) parameter values that
   * give us zero in the projected mesh.
   */
  s[0] = t[0] = 0.0;
  s[1] = t[1] = 1.0;
  return findZeros(P2, s,t, 1, st_tolerance, zeros);
}

#ifdef COMMENT_OUT

static
void pointOnSurface(POINT3 P[4][4], double s, double t, POINT3 *C) {
  double s_ = 1.0 - s;
  double ss = s*s;
  double sss = s*ss;
  double s_s_ = s_*s_;
  double s_s_s_ = s_*s_s_;
  double s_s_s = 3.0*s_s_*s;
  double s_ss = 3.0*s_*ss;

  double t_ = 1.0 - t;
  double tt = t*t;
  double ttt = t*tt;
  double t_t_ = t_*t_;
  double t_t_t_ = t_*t_t_;
  double t_t_t = 3.0*t_t_*t;
  double t_tt = 3.0*t_*tt;

  int i;
  POINT3 Q[4];

  for (i = 0; i <= 3; i++) {
    Q[i].x = 
      t_t_t_*P[i][0].x + t_t_t*P[i][1].x + t_tt*P[i][2].x + ttt*P[i][3].x;
    Q[i].y = 
      t_t_t_*P[i][0].y + t_t_t*P[i][1].y + t_tt*P[i][2].y + ttt*P[i][3].y;
    Q[i].z = 
      t_t_t_*P[i][0].z + t_t_t*P[i][1].z + t_tt*P[i][2].z + ttt*P[i][3].z;
  }

  C->x = s_s_s_*Q[0].x + s_s_s*Q[1].x + s_ss*Q[2].x + sss*Q[3].x;
  C->y = s_s_s_*Q[0].y + s_s_s*Q[1].y + s_ss*Q[2].y + sss*Q[3].y;
  C->z = s_s_s_*Q[0].z + s_s_s*Q[1].z + s_ss*Q[2].z + sss*Q[3].z;
}

#endif /* COMMENT_OUT */

void pointOnSurface(POINT3 P[4][4], double s, double t, POINT3 *C) {
  int i;
  POINT3 Q[4];
  double s_ = 1.0 - s, t_ = 1.0 - t;
  POINT3 P01, P12, P23;

  for (i = 0; i <= 3; i++) {
    P01.x = t_*P[i][0].x + t*P[i][1].x;
    P01.y = t_*P[i][0].y + t*P[i][1].y;
    P01.z = t_*P[i][0].z + t*P[i][1].z;
    
    P12.x = t_*P[i][1].x + t*P[i][2].x;
    P12.y = t_*P[i][1].y + t*P[i][2].y;
    P12.z = t_*P[i][1].z + t*P[i][2].z;

    P23.x = t_*P[i][2].x + t*P[i][3].x;
    P23.y = t_*P[i][2].y + t*P[i][3].y;
    P23.z = t_*P[i][2].z + t*P[i][3].z;

    P01.x = t_*P01.x + t*P12.x;
    P01.y = t_*P01.y + t*P12.y;
    P01.z = t_*P01.z + t*P12.z;

    P12.x = t_*P12.x + t*P23.x;
    P12.y = t_*P12.y + t*P23.y;
    P12.z = t_*P12.z + t*P23.z;

    Q[i].x = t_*P01.x + t*P12.x;
    Q[i].y = t_*P01.y + t*P12.y;
    Q[i].z = t_*P01.z + t*P12.z;
  }

  Q[0].x = s_*Q[0].x + s*Q[1].x;
  Q[0].y = s_*Q[0].y + s*Q[1].y;
  Q[0].z = s_*Q[0].z + s*Q[1].z;
  
  Q[1].x = s_*Q[1].x + s*Q[2].x;
  Q[1].y = s_*Q[1].y + s*Q[2].y;
  Q[1].z = s_*Q[1].z + s*Q[2].z;

  Q[2].x = s_*Q[2].x + s*Q[3].x;
  Q[2].y = s_*Q[2].y + s*Q[3].y;
  Q[2].z = s_*Q[2].z + s*Q[3].z;

  Q[0].x = s_*Q[0].x + s*Q[1].x;
  Q[0].y = s_*Q[0].y + s*Q[1].y;
  Q[0].z = s_*Q[0].z + s*Q[1].z;
  
  Q[1].x = s_*Q[1].x + s*Q[2].x;
  Q[1].y = s_*Q[1].y + s*Q[2].y;
  Q[1].z = s_*Q[1].z + s*Q[2].z;

  C->x = s_*Q[0].x + s*Q[1].x;
  C->y = s_*Q[0].y + s*Q[1].y;
  C->z = s_*Q[0].z + s*Q[1].z;
}


static
void derivative(POINT3 P[4][4], double s, double t, int s_dir, POINT3 *D) {
  POINT3 Q[4], dQ[4];
  POINT3 P01, P12, P23;
  POINT3 P012, P123;
  double u,u_;
  int i;

  /*
   * If we want the derivative in the s_dir then we perform
   * deCasteljau iteration along each row of control points and
   * build a set of control points for a curve in the s
   * direction.
   */
  if (s_dir) {    
    double t_ = 1.0 - t;
    u = s; u_ = 1.0 - s;
    for (i = 0; i <= 3; i++) {
      P01.x = t_*P[i][0].x + t*P[i][1].x;
      P01.y = t_*P[i][0].y + t*P[i][1].y;
      P01.z = t_*P[i][0].z + t*P[i][1].z;
    
      P12.x = t_*P[i][1].x + t*P[i][2].x;
      P12.y = t_*P[i][1].y + t*P[i][2].y;
      P12.z = t_*P[i][1].z + t*P[i][2].z;
  
      P23.x = t_*P[i][2].x + t*P[i][3].x;
      P23.y = t_*P[i][2].y + t*P[i][3].y;
      P23.z = t_*P[i][2].z + t*P[i][3].z;

      P012.x = t_*P01.x + t*P12.x;
      P012.y = t_*P01.y + t*P12.y;
      P012.z = t_*P01.z + t*P12.z;
      
      P123.x = t_*P12.x + t*P23.x;
      P123.y = t_*P12.y + t*P23.y;
      P123.z = t_*P12.z + t*P23.z;

      Q[i].x = t_*P012.x + t*P123.x;
      Q[i].y = t_*P012.y + t*P123.y;
      Q[i].z = t_*P012.z + t*P123.z;
    }
  } else {
    double s_ = 1.0 - s;
    u = t; u_ = 1.0 - t;
    for (i = 0; i <= 3; i++) {
      P01.x = s_*P[0][i].x + s*P[1][i].x;
      P01.y = s_*P[0][i].y + s*P[1][i].y;
      P01.z = s_*P[0][i].z + s*P[1][i].z;
    
      P12.x = s_*P[1][i].x + s*P[2][i].x;
      P12.y = s_*P[1][i].y + s*P[2][i].y;
      P12.z = s_*P[1][i].z + s*P[2][i].z;
  
      P23.x = s_*P[2][i].x + s*P[3][i].x;
      P23.y = s_*P[2][i].y + s*P[3][i].y;
      P23.z = s_*P[2][i].z + s*P[3][i].z;

      P012.x = s_*P01.x + s*P12.x;
      P012.y = s_*P01.y + s*P12.y;
      P012.z = s_*P01.z + s*P12.z;
      
      P123.x = s_*P12.x + s*P23.x;
      P123.y = s_*P12.y + s*P23.y;
      P123.z = s_*P12.z + s*P23.z;

      Q[i].x = s_*P012.x + s*P123.x;
      Q[i].y = s_*P012.y + s*P123.y;
      Q[i].z = s_*P012.z + s*P123.z;
    }
  }

  /*
   * To compute the derivative we build a quadratic difference
   * curve.
   */
  for (i = 0; i < 3; i++) {
    dQ[i].x = Q[i+1].x - Q[i].x;
    dQ[i].y = Q[i+1].y - Q[i].y;
    dQ[i].z = Q[i+1].z - Q[i].z;
  }

  /*
   * Now we evaluate this quadratic difference curve at the
   * appropriate point.
   */
  P01.x = u_*dQ[0].x + u*dQ[1].x;
  P01.y = u_*dQ[0].y + u*dQ[1].y;
  P01.z = u_*dQ[0].z + u*dQ[1].z;

  P12.x = u_*dQ[1].x + u*dQ[2].x;
  P12.y = u_*dQ[1].y + u*dQ[2].y;
  P12.z = u_*dQ[1].z + u*dQ[2].z;

  /*
   * The actual derivative would need to be scaled by 3
   * but we ultimately only want the direction.
   */
  D->x = u_*P01.x + u*P12.x;
  D->y = u_*P01.y + u*P12.y;
  D->z = u_*P01.z + u*P12.z;
}

static
void hitNormal(OBJECT *this, double hit[3], HIT_INFO *hitInfo, 
	       double normal[3]) {
  BEZIER3_DATA *data = (BEZIER3_DATA *) this->data;
  POINT3 Ds, Dt;
  double mag;

  /*
   * Compute the derivative along both the s and t
   * parameterization.
   */
  derivative(data->P, hitInfo->bezier3.s, hitInfo->bezier3.t, 1, &Ds);
  derivative(data->P, hitInfo->bezier3.s, hitInfo->bezier3.t, 0, &Dt);
  
  /*
   * Compute the cross product of the two derivatives to
   * get the normal.
   */
  normal[0] = Ds.y*Dt.z - Ds.z*Dt.y;
  normal[1] = Ds.z*Dt.x - Ds.x*Dt.z;
  normal[2] = Ds.x*Dt.y - Ds.y*Dt.x;

  /*
   * Normalize the normal vector (if its non-zero).
   */
  mag = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
  if (mag > 0.0) {
    mag = 1.0/sqrt(mag);
    normal[0] *= mag;
    normal[1] *= mag;
    normal[2] *= mag;
  }
}

static
double rayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
	      HIT_INFO *hitInfo) {
  BEZIER3_DATA *data = (BEZIER3_DATA *) this->data;
  int i, numZeros;
  double zeros[MAX_ZEROS][2];
  double mint;

  /*
   * Find all the (s,t) parameterizations that intersect the patch.
   */
  numZeros = findAllZeros(rayOrg, rayDir,
			  data->P, data->st_tolerance, zeros);
  if (numZeros <= 0) return -1.0;  /* no intersection */

  /*
   * Find the minimum positive distance along the ray
   * to the point. The distance has to be greater than
   * EPSILON to avoid self-intersection stuff.
   * If there is an intersection in front of the eye we want
   * save the u,v coords of the intersection for 
   * future normal computations.
   */
  mint = INFINITY;
  for (i = 0; i < numZeros; i++) {
    POINT3 C,V;
    double t;
    pointOnSurface(data->P, zeros[i][0], zeros[i][1], &C);
    V.x = C.x - rayOrg[0];
    V.y = C.y - rayOrg[1];
    V.z = C.z - rayOrg[2];
    t = V.x*rayDir[0] + V.y*rayDir[1] + V.z*rayDir[2];
    if (t > EPSILON && t < mint) {
      mint = t;
      hitInfo->bezier3.s = zeros[i][0];
      hitInfo->bezier3.t = zeros[i][1];
    }
  }

  return (mint < INFINITY) ? mint : -1.0;
}

static
double trimRayHit(OBJECT *this, double rayOrg[3], double rayDir[3],
		  HIT_INFO *hitInfo) {
  double t;
  BEZIER3_DATA *data;
  BEZIER3_TRIM_DATA *trimData;
  double u,v, u_,v_;

  /*
   * Call super class rayHit method, and if we don't get a hit
   * then we are done.
   */
  if ((t = rayHit(this, rayOrg, rayDir, hitInfo)) < EPSILON)
    return -1.0;

  /*
   * Get object specific data which includes trimming data.
   */
  data = (BEZIER3_DATA *) this->data;
  trimData = data->trim;

  /*
   * Here we take our s and t parameter values and map them 
   * to our u,v values for trimming. First we adjust to the
   * trimming origin. Note that we map s -> u, t -> v.
   */
  u_ = hitInfo->bezier3.s - trimData->uvmap.org[0];
  v_ = hitInfo->bezier3.t - trimData->uvmap.org[1];

  /*
   * Next we map u,v to align with our u,v axis.
   */
  u = u_*trimData->uvmap.u[0] + v_*trimData->uvmap.u[1];
  v = u_*trimData->uvmap.v[0] + v_*trimData->uvmap.v[1];

  /*
   * If we are not tiling the trimming region and the hit point
   * is outside the region we'll return the normal hit information.
   */
  if (!trimData->uvmap.tile && 
      (u < 0.0 || v < 0.0 || 
       u > trimData->uvmap.width || v > trimData->uvmap.height)) 
    return t;

  /*
   * Scale (u,v) to unit box.
   */
  u /= trimData->uvmap.width;
  v /= trimData->uvmap.height;

  /*
   * Subtract off any whole terms.
   */
  u -= (double) ((int) u);
  v -= (double) ((int) v);

  /*
   * Wrap negative terms (negative floats truncate towards zero).
   */
  if (u < 0.0) u += 1.0;
  if (v < 0.0) v += 1.0;

  return ((trimData->trimmer->trim)(trimData->trimmer, u,v)) ? t : -1.0;
}
  
static
void solidColor(OBJECT *this, double hit[3], HIT_INFO *hitInfo, 
		double color[3]) {
  BEZIER3_DATA *data = (BEZIER3_DATA *) this->data;
  color[0] = data->color[0];
  color[1] = data->color[1];
  color[2] = data->color[2];
}

static
void setSolidColor(OBJECT *this, double color[3]) {
  BEZIER3_DATA *data = (BEZIER3_DATA *) this->data;
  data->color[0] = color[0];
  data->color[1] = color[1];
  data->color[2] = color[2];
}

OBJECT *createBezier3Object(POINT3 P[4][4]) {
  OBJECT *obj;
  BEZIER3_DATA *data;
  int i,j;

  /*
   * Allocate memory for object.
   */
  if ((obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL ||
      (data = (BEZIER3_DATA *) malloc(sizeof(BEZIER3_DATA))) == NULL) {
    perror("createBezier3Object:malloc");
    exit(-1);
  }
  
  /*
   * This should be computed by finding a bound on the largest
   * first derivative of screen space x and y with respect to the
   * parameter space s and t. For now we'll just set it to a small
   * value.
   */
  data->st_tolerance = ST_TOLERANCE;

  /*
   * Initially no trimming information.
   */
  data->trim = NULL;

  /*
   * Copy control points.
   */
  for (i = 0; i <= 3; i++)
    for (j = 0; j <= 3; j++)
      data->P[i][j] = P[i][j];

  /*
   * Pick some nice defaults.
   */
  data->color[0] = 0.8;
  data->color[1] = 0.1;
  data->color[2] = 0.1;

  obj->ka = 0.15;
  obj->kd = 0.60;
  obj->ks = 0.25;
  obj->kt = 0.0;
  obj->ni = 0.0;
  obj->phong = 10.0;

  obj->data = data;
  obj->rayHit = rayHit;
  obj->normal = hitNormal;
  obj->color = solidColor;
  obj->setColor = setSolidColor;

  return obj;
}

void destroyBezierObject(OBJECT *obj) {
  free(obj->data);
  free(obj);
  /* XXX need to free trim data if there is any */
}

void trimBezier3(OBJECT *object,             /* bezier3 object */
		 TRIM_OBJECT *trimObject,    /* trimming object */
		 int tile,                   /* tile trimming? */
		 double origin[2],           /* u,v origin on st-plane */
		 double uaxis[2],            /* u axis on st-plane */
		 double vaxis[2],            /* v axis on st-plane */
		 double width,               /* width of trim */
		 double height) {            /* height of trim */
  BEZIER3_DATA *data = (BEZIER3_DATA *) object->data;
  BEZIER3_TRIM_DATA *trimData;
  double scale;
  double u[2],v[2];

  /*
   * Allocate memory for trimdata (if its not there already)
   */
  if (data->trim == NULL)
    if ((data->trim = (BEZIER3_TRIM_DATA *) 
	 malloc(sizeof(BEZIER3_TRIM_DATA))) == NULL) {
      perror("trimBezier3:malloc");
      exit(-1);
    }

  trimData = data->trim;

  /*
   * Unitize u,v coordinate axes and make sure they are orthoganol.
   */
  u[0] = uaxis[0];  u[1] = uaxis[1];
  v[0] = vaxis[0];  v[1] = vaxis[1];
  scale = 1.0/sqrt(u[0]*u[0] + u[1]*u[1]);
  u[0] *= scale;
  u[1] *= scale;
  scale = v[0]*u[0] + v[1]*u[1];
  v[0] -= scale*u[0];
  v[1] -= scale*u[1];
  scale = 1.0/sqrt(v[0]*v[0] + v[1]*v[1]);
  v[0] *= scale;
  v[1] *= scale;

  /*
   * Copy trimming data.
   */
  trimData->uvmap.tile = tile;
  trimData->uvmap.org[0] = origin[0];
  trimData->uvmap.org[1] = origin[1];
  trimData->uvmap.u[0] = u[0];
  trimData->uvmap.u[1] = u[1];
  trimData->uvmap.v[0] = v[0];
  trimData->uvmap.v[1] = v[1];
  trimData->uvmap.width = width;
  trimData->uvmap.height = height;

  /*
   * Here we would like to make a clone of the object, but for
   * now we will just reference the given one (and hope it
   * isn't destroyed).
   */
  trimData->trimmer = trimObject;

  /*
   * Update the rayHit method.
   */
  object->rayHit = trimRayHit;
}
