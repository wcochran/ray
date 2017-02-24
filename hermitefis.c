#include <stdio.h>
#include <stdlib.h>
#include "rectfis.h"
#include "bilinear.h"
#include "hermite.h"
#include "hermitefis.h"

RectFis *createHermiteFis(int M, int N, double *X, double *Y) {
  RectFis *fis = rectFisCreate(M, N, X, Y);
  int i,j;

  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++) {
      RectFisMap *map = &fis->maps[j][i];
      map->alpha = createBilinearFunc2(fis->X[i], fis->X[i+1],
				       fis->Y[j], fis->Y[j+1]);
      map->q = createGeneralHermiteFunc2(fis->X[i], fis->X[i+1],
					 fis->Y[j], fis->Y[j+1]);
    }

  return fis;
}

#define EVEN(n) (((n) & 1) == 0)

/*
 * ui = (x - bi)/ai
 * vj = (y - dj)/cj
 * f(x,y) = qij(x,y) + alphaij(x,y)*f(ui, vj)
 *
 * d/dx f(x,y) = d/dx qij(x,y) + alphaij(x,y)/ai * d/du f(ui, vj)
 *                             + d/dx alphaij(x,y) * f(ui, vj)
 *
 * d/dy f(x,y) = d/dy qij(x,y) + alphaij(x,y)/cj * d/dv f(ui, vj)
 *                             + d/dy alphaij(x,y) * f(ui, vj)
 *
 * d2/dxdy f(x,y) = d2/dxdy qij(x,y) + alphaij(x,y)/(ai*cj) * d2/dudv f(ui, vj)
 *                                   + d/dy alphaij(x,y)/ai * d/du f(ui, vj)
 *                                   + d/dx alpjaij(x,y)/cj * d/dv f(ui, vj)
 *                                   + d2/dxdy alphaij(x,y) * f(ui, vj)
 */
void setHermiteFis(RectFis *fis, double **alpha, 
		   double **Z, double **Zx, double **Zy, double **Zxy) {
  int M = fis->M, N = fis->N;
  int i,j;

  for (j = 0; j <= N; j++)
    for (i = 0; i <= M; i++)
      fis->Z[j][i] = Z[j][i];

  for (j = 0; j < N; j++) {
    int jj[2];

    if (EVEN(j)) {
      jj[0] = 0;
      jj[1] = N;
    } else {
      jj[0] = N;
      jj[1] = 0;
    }
    
    for (i = 0; i < M; i++) {
      RectFisMap *map = &fis->maps[j][i];
      double inva = 1.0/map->a, invc = 1.0/map->c;
      int ii[2];
      double z[2][2], zx[2][2], zy[2][2], zxy[2][2];
      int r, c;

      if (EVEN(i)) {
	ii[0] = 0;
	ii[1] = M;
      } else {
	ii[0] = M;
	ii[1] = 0;
      }
    
      z[0][0] = alpha[j][i];
      z[0][1] = alpha[j][i+1];
      z[1][0] = alpha[j+1][i];
      z[1][1] = alpha[j+1][i+1];
      setBilinearFunc2(map->alpha, z);
  
      for (r = 0; r < 2; r++)
        for (c = 0; c < 2; c++) {
          double F = Z[j+r][i+c];
          double Fx = Zx[j+r][i+c];
          double Fy = Zy[j+r][i+c];
          double Fxy = Zxy[j+r][i+c];
          double FF = Z[jj[r]][ii[c]];
          double FFx = Zx[jj[r]][ii[c]];
          double FFy = Zy[jj[r]][ii[c]];
          double FFxy = Zxy[jj[r]][ii[c]];
          double A = alpha[j+r][i+c];
          double Ax = (*map->alpha->dx)(map->alpha, fis->X[i+c],fis->Y[j+r]);
          double Ay = (*map->alpha->dy)(map->alpha, fis->X[i+c],fis->Y[j+r]);
          double Axy = (*map->alpha->dxdy)(map->alpha, fis->X[i+c],fis->Y[j+r]);

          z[r][c] = F - A*FF;
          zx[r][c] = Fx - (A*inva*FFx + Ax*FF);
          zy[r][c] = Fy - (A*invc*FFy + Ay*FF);
          zxy[r][c] = Fxy - 
            (A*inva*invc*FFxy + Ay*inva*FFx + Ax*invc*FFy + Axy*FF);
        }

      setGeneralHermiteFunc2(map->q, z, zx, zy, zxy);
    }
  }
}
