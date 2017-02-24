/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef RECTGRID_H
#define RECTGRID_H

#include "raytrace.h"

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
						double *z, int rowStride);

void mapImageToSurfaceSupportedByRectangularGrid(OBJECT *obj,
						 int tile,
						 int interpolate,
						 double xyorg[2],
						 double yaxis[2],
						 double width, double height,
						 pnm_image *image);
#endif /* RECTGRID_H */
