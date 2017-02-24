#include "unimodalroot.h"
#include <math.h>

#define TOL 1e-8
#define GOLD_TOL 1e-5

/*
 * Root of g(t) bracketed on interval [a,b].
 * Employs Bisection Method to isolate root.
 * If b < minRoot at any point we short-circuit
 * the process since we're convering on a root
 * we are not interested in.
 */
static
double bisection(void *data, double (*g)(void *data, double t),
		 double a, double b, double ga, double gb, 
		 double minRoot) {
  double c, gc;
  do {
    c = 0.5*(a + b);
    gc = g(data, c);
    if (ga*gc < 0) {
      b = c;
      gb = gc;
      if (b < minRoot)  /* short-circuit : not interested in this root */
	return b;
    } else {
      a = c;
      ga = gc;
    }
  } while (b - a > TOL);
  return c;
}

/*
 * Given a bracketed minimum (a,b,c) we initiate a search
 * for the minimum via the Golden Section Search method.
 * If we encounted a t such that g(t) < 0 along the way,
 * we immediately return that. 
 * gb initially holds g(b) on input.
 * gb holds the minimum or first encountered negative value
 * upon return.
 * The minimal t-value is returned.
 * The Golden Section Search converges when (x3 - x0) < sqrt(TOL)
 * -- any more function evaluations is a waste of time.
 */
static
double goldenSearch(void *data, double (*g)(void *data, double t),
		    double a, double b, double c, double *gb) {
#define W 0.381966
  double x0, x1, x2, x3;
  double g1, g2;

  x0 = a;
  x3 = c;
  if (c - b > b - a) {
    x1 = b;
    g1 = *gb;
    x2 = b + W*(c - b);
    g2 = g(data, x2);
    if (g2 < 0) {
      *gb = g2;
      return x2;
    }
  } else {
    x1 = b - W*(b - a);
    g1 = g(data, x1);
    if (g1 < 0) {
      *gb = g1;
      return x1;
    }
    x2 = b;
    g2 = *gb; 
  }

  while (x3 - x0 > GOLD_TOL)
    if (g1 < g2) {
      x3 = x2;
      x2 = x1;
      g2 = g1;
      x1 = (1-W)*x1 + W*x0;
      g1 = g(data, x1);
      if (g1 < 0) {
	*gb = g1;
	return x1;
      }
    } else {
      x0 = x1;
      x1 = x2;
      g1 = g2;
      x2 = (1-W)*x2 + W*x3;
      g2 = g(data, x2);
      if (g2 < 0) {
	*gb = g2;
	return x2;
      }
    }

  if (g1 < g2) {
    *gb = g1;
    return x1;
  }

  *gb = g2;
  return x2;
}

/*
 * Define state datatype and function
 * for h(t) = -g(t).
 */
typedef struct {
  double (*g)(void *, double);  /* original function g */
  double *gdata;                /* original state info for g */
} H_DATA;

double h(void *data, double t) {  /* h(t) = -g(t) */
  H_DATA *hdata = (H_DATA*) data;
  return -hdata->g(hdata->gdata, t);
}

/*
 * input:
 *   data : state information to pass to g when evaluating g(t).
 *   g(data,t) : function for which we are seeking roots.
 *   a, b : a < b input interval
 *   ga, gb : g(a) and g(b) -- we assume neither are 0.
 *   tmin : minimal root we are seeking
 * returns:
 *   smallest root >= tmin or -1 if no root found.
 */
double unimodalRoot(void *data, double (*g)(void *data, double t),
		    double a, double b, double ga, double gb, double tmin) {
  H_DATA hdata;
  double r, c, gc;

  /*
   * If root bracketed, use bisection method
   * to find root.
   */
  if (ga*gb < 0) {
    r = bisection(data, g, a, b, ga, gb, tmin);
    return (r >= tmin) ? r : -1;
  }

  /*
   * If ga < 0 and gb < 0 we switch problem to
   * finding roots of h(t) = -g(t)
   */
  if (ga < 0 && gb < 0) {
    hdata.g = g;
    hdata.gdata = data;
    data = &hdata;
    g = h;
    ga = -ga;
    gb = -gb;
  }

  /*
   * At this point we ga > 0 and gb > 0.
   * We'll begin a search for the minimum of g(t) on [a,b]
   * and see if there exists a value c on the interval
   * where g(c) < 0. If so we will have bracketed two
   * roots with intervals [a,c] and [c,b].
   */
  do {
    c = 0.5*(a + b);
    gc = g(data, c);

    if (gc < 0) {
      r = bisection(data, g, a, c, ga, gc, tmin);
      if (r >= tmin) return r;
      r = bisection(data, g, c, b, gc, gb, tmin);
      return (r >= tmin) ? r : -1;
    }

    if (gc < ga) {

      if (gc < gb) {
	/*
	 * gc < ga and gc < gb.
	 * We have a minimum bracketed so we begin the Golden Section search.
	 * If search yields a c such thatn g(c) < 0, then we
	 * have bracketed both roots and we bisect to find them.
	 * Otherwise, there are no roots.
	 */
	c = goldenSearch(data, g, a,c,b, &gc);
	if (gc < 0.0) {
	  r = bisection(data, g, a, c, ga, gc, tmin);
	  if (r >= tmin) return r;
	  r = bisection(data, g, c, b, gc, gb, tmin);
	  if (r >= tmin) return r;
	} else
	  return -1;
      }

      /*
       * Minimum (and our roots) will be in [c,b] if they exist.
       */
      a = c;
      ga = gc;

    } else if (gc < gb) {

      /*
       * Minimum (and our roots) will be in [a,c] if they exist.
       */
      b = c;
      gb = gc;

    } else {

      /*
       * Function increased or stayed the same.
       * No roots possible.
       */
      return -1;
    }

  } while (b - a > TOL);

  return -1;
}
