#include <stdlib.h>
#include <math.h>
#include "checker.h"
#include "noise.h"

static void color(struct OBJECT *this,
		  double hitPoint[3],
		  HIT_INFO *info,
		  double color[3]) {
  const double n = noise3(hitPoint[0],hitPoint[1],hitPoint[2]);
  const double t = turbulence3(hitPoint[0],hitPoint[1],hitPoint[2],8);
  // const int i = (int) (6*n*sin(hitPoint[0]));
  // const int j = (int) (0.5*hitPoint[1]);
  // const int k = (int) (3*hitPoint[2]);
  static double color0[3] = {0, 0, 1};
  static double color1[3] = {0.9, 0.8, 0};
  static double color2[3] = {0, 0.4, 1};
  static double color3[3] = {1, 0.3, 0.5};
  const double f = sin(1.5*n*hitPoint[0] + 1.6*t);
  const double u = (f + 1)/2;
  const double uu = u*u;
  const double v = 1-u;
  const double vv = v*v;
  color[0] = vv*v*color0[0] + 3*vv*u*color1[0] + 3*v*uu*color2[0] + uu*u*color3[0];
  color[1] = vv*v*color0[1] + 3*vv*u*color1[1] + 3*v*uu*color2[1] + uu*u*color3[1];
  color[2] = vv*v*color0[2] + 3*vv*u*color1[2] + 3*v*uu*color2[2] + uu*u*color3[2];
}

#ifdef COMMENT_OUT_V1

static void color(struct OBJECT *this,
		  double hitPoint[3],
		  HIT_INFO *info,
		  double color[3]) {
  const double n = noise3(hitPoint[0],hitPoint[1],hitPoint[2]);
  /*
    const double n = 
    turbulence3(hitPoint[0],hitPoint[1],hitPoint[2], 16);
  const int i = (int) (6*n*sin(hitPoint[0]));
  const int j = (int) (0.5*hitPoint[1]);
  const int k = (int) (3*hitPoint[2]);
  */
  const int i = (int) (3.5*n*hitPoint[0]);
  const int j = (int) hitPoint[1];
  const int k = (int) hitPoint[2];
  static double color0[3] = {1, 0, 0};
  static double color1[3] = {1, 1, 0};
  if (i  % 2 == 0) {
    color[0] = color0[0];
    color[1] = color0[1];
    color[2] = color0[2];
  } else {
    color[0] = color1[0];
    color[1] = color1[1];
    color[2] = color1[2];
  }
}

#endif

#ifdef COMMENT_OUT_XXXXXXXXXX

static void color(struct OBJECT *this,
		  double hitPoint[3],
		  HIT_INFO *info,
		  double color[3]) {

  /* marble_color[sin((x,y,z)*dir + turbulence(x,y,z))] */

  const double x = hitPoint[0];
  const double y = hitPoint[1];
  const double z = hitPoint[2];
  const double f = sin(3*x + y + turbulence3(x,y,z,16));
  const double w = (f + 1)/2;
  
  static double color0[3] = {0, 0, 1};
  static double color1[3] = {1, 1, 0};

  color[0] = (1-w)*color0[0] + w*color1[0];
  color[1] = (1-w)*color0[1] + w*color1[1];
  color[2] = (1-w)*color0[2] + w*color1[2];
}

static void color(struct OBJECT *this,
		  double hitPoint[3],
		  HIT_INFO *info,
		  double color[3]) {

  /* marble_color[sin((x,y,z)*dir + turbulence(x,y,z))] */

  const double x = hitPoint[0];
  const double y = hitPoint[1];
  const double z = hitPoint[2];
  const double f = sin(3*x + y + turbulence3(x,y,z,16));
  const double w = (f + 1)/2;
  
  static double color0[3] = {0, 0, 1};
  static double color1[3] = {1, 1, 0};

  color[0] = (1-w)*color0[0] + w*color1[0];
  color[1] = (1-w)*color0[1] + w*color1[1];
  color[2] = (1-w)*color0[2] + w*color1[2];
}


#endif

OBJECT *checkerObject(OBJECT *old) {
  OBJECT *obj;
  
  if ((obj = (OBJECT *) malloc(sizeof(OBJECT))) == NULL) {
    perror("checkerObject():malloc()");
    exit(-1);
  }

  *obj = *old;
  obj->color = color;

  return obj;
}
