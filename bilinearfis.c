#include <stdio.h>
#include <stdlib.h>
#include "bilinearfis.h"

#define EVEN(n) (((n) & 1) == 0)

RectFis *createBilinearFis(int M, int N, double *X, double *Y) {
  RectFis *fis = rectFisCreate(M, N, X, Y);
  int i,j;

  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++) {
      RectFisMap *map = &fis->maps[j][i];
      map->alpha = createBilinearFunc2(fis->X[i], fis->X[i+1],
				       fis->Y[j], fis->Y[j+1]);
      map->q = createBilinearFunc2(fis->X[i], fis->X[i+1],
				   fis->Y[j], fis->Y[j+1]);
    }

  return fis;
}

void setBilinearFis(RectFis *fis, double **Z, double **alpha) {
  int M = fis->M, N = fis->N;
  int i,j;

  for (j = 0; j <= N; j++)
    for (i = 0; i <= M; i++)
      fis->Z[j][i] = Z[j][i];

  for (j = 0; j < N; j++) {
    int J0, J1;

    if (EVEN(j)) {
      J0 = 0;
      J1 = N;
    } else {
      J0 = N;
      J1 = 0;
    }
    
    for (i = 0; i < M; i++) {
      RectFisMap *map = &fis->maps[j][i];
      int I0, I1;
      double z[2][2];

      if (EVEN(i)) {
	I0 = 0;
	I1 = M;
      } else {
	I0 = M;
	I1 = 0;
      }
    
      z[0][0] = alpha[j][i];
      z[0][1] = alpha[j][i+1];
      z[1][0] = alpha[j+1][i];
      z[1][1] = alpha[j+1][i+1];
      setBilinearFunc2(map->alpha, z);

      z[0][0] = Z[j][i] - alpha[j][i]*Z[J0][I0];
      z[0][1] = Z[j][i+1] - alpha[j][i+1]*Z[J0][I1];
      z[1][0] = Z[j+1][i] - alpha[j+1][i]*Z[J1][I0];
      z[1][1] = Z[j+1][i+1] - alpha[j+1][i+1]*Z[J1][I1];
      setBilinearFunc2(map->q, z);
    }
  }
}
