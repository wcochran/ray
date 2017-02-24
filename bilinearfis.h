#ifndef BILINEARFIS_H
#define BILINEARFIS_H

#include "bilinear.h"
#include "rectfis.h"

RectFis *createBilinearFis(int M, int N, double *X, double *Y);
void setBilinearFis(RectFis *fis, double **Z, double **alpha);

#endif /* BILINEARFIS_H */
