#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "unimodalroot.h"

/*
 * g(t) = s*(x - xc)^2 + yc
 */
typedef struct {
  double s, xc, yc;
} PARABOLA_DATA;

double parabola(void *data, double x) {
  PARABOLA_DATA *p = (PARABOLA_DATA *) data;
  double dx = x - p->xc;
  return p->s*dx*dx + p->yc;
}

int main(void) {
  PARABOLA_DATA pdata;
  double a, b, ga, gb;
  double r, gr;

  /*
   * Test case 1.
   * g(t) = t^2 - 1.
   * tmin = 0;
   * [a,b] = -10, 20.
   */
  pdata.xc = 0;
  pdata.yc = -1;
  pdata.s = 1;
  a = -10;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  /*
   * Test case 2 (no roots)
   * g(t) = t^2 + 1 
   * tmin = 0;
   * [a,b] = -10, 20.
   */
  pdata.xc = 0;
  pdata.yc = 1;
  pdata.s = 1;
  a = -10;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f (should be -1)\n", r);

  /*
   * Test case 3.
   * g(t) = -t^2 + 1.
   * tmin = 0;
   * [a,b] = -10, 20.
   */
  pdata.xc = 0;
  pdata.yc = +1;
  pdata.s = -1;
  a = -10;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  /*
   * Test case 4.
   * g(t) = 10*(t - 5)^2 - 7.
   * tmin = 0;
   * [a,b] = -10, 20.
   */
  pdata.xc = 5;
  pdata.yc = -7;
  pdata.s = 10;
  a = -10;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  /*
   * Test case 5.
   * g(t) = 10*(t - 5)^2 - 7.
   * tmin = 0;
   * [a,b] = 0.0, 20.
   */
  pdata.xc = 5;
  pdata.yc = -7;
  pdata.s = 10;
  a = 0.0;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  /*
   * Test case 6.
   * g(t) = t^2 - 1.
   * tmin = 0;
   * [a,b] = -1.1, 20.
   */
  pdata.xc = 0;
  pdata.yc = -1;
  pdata.s = 1;
  a = -1.1;
  b = 20;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  /*
   * Test case 7.
   * g(t) = t^2 - 1.
   * tmin = 0;
   * [a,b] = -20, 1.1.
   */
  pdata.xc = 0;
  pdata.yc = -1;
  pdata.s = 1;
  a = -20;
  b = 1.1;
  ga = parabola(&pdata, a);
  gb = parabola(&pdata, b);
  r = unimodalRoot(&pdata, parabola, a,b, ga,gb, 0);
  gr = parabola(&pdata, r);
  printf("r=%0.10f, g(r)=%0.10f\n", r, gr);

  return 0;
}
