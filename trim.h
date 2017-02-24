
/* $Author: cs548 $ $Revision: 1.1 $ $Date: 2009/10/09 19:02:19 $ */

#ifndef TRIM_H
#define TRIM_H

typedef struct TRIM_OBJECT {
  void *data;                     /* private data */
  int (*trim)                     /* trim method (0 => out, 1 => in) */
    (struct TRIM_OBJECT *this,    /* pointer to this instance */
     double u, double v);         /* (u,v) point to test (0 <= u,v <= 1) */
} TRIM_OBJECT;

#endif /* TRIM_H */
