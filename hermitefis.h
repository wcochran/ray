#ifndef HERMITEFIS_H
#define HERMITEFIS_H

#include "rectfis.h"

RectFis *createHermiteFis(int M, int N, double *X, double *Y);
void setHermiteFis(RectFis *fis, double **alpha, 
		   double **Z, double **Zx, double **Zy, double **Zxy);

#endif /* HERMITEFIS_H */
