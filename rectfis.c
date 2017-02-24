#include <stdio.h>
#include <stdlib.h>
#include "func2.h"
#include "rectfis.h"

static void *myalloc(size_t sz) {
  void *p = malloc(sz);
  if (p == NULL) {
    perror("malloc()");
    exit(-1);
  }
  return p;
}

double rectFisTransform(RectFisMap *map, double src[3], double dst[3]) {
  double x_ = dst[0] = map->a*src[0] + map->b;
  double y_ = dst[1] = map->c*src[1] + map->d;
  dst[2] = (*map->q->eval)(map->q, x_, y_) + 
    (*map->alpha->eval)(map->alpha, x_, y_)*src[2];
  return dst[2];
}

#define EVEN(n) (((n) & 1) == 0)

RectFis *rectFisCreate(int M, int N, double *X, double *Y) {
  RectFis *fis = (RectFis *) myalloc(sizeof(RectFis));
  int i,j;

  fis->M = M;
  fis->N = N;

  fis->X = (double *) myalloc((M+1)*sizeof(double));
  fis->Y = (double *) myalloc((N+1)*sizeof(double));
  fis->Z = (double **) myalloc((N+1)*sizeof(double *));
  for (j = 0; j <= N; j++)
    fis->Z[j] = (double *) myalloc((M+1)*sizeof(double));
  for (i = 0; i <= M; i++)
    fis->X[i] = X[i];
  for (j = 0; j <= N; j++) {
    fis->Y[j] = Y[j];
    for (i = 0; i <= M; i++)
      fis->Z[j][i] = 0.0;
  }

  fis->maps = (RectFisMap **) myalloc(N*sizeof(RectFisMap *));
  for (j = 0; j < N; j++)
    fis->maps[j] = (RectFisMap *) myalloc(M*sizeof(RectFisMap));

  for (j = 0; j < N; j++) {
    double c, d;

    if (EVEN(j)) {
      c = (Y[j+1] - Y[j])/(Y[N] - Y[0]);
      d = Y[j] - c*Y[0];
    } else {
      c = (Y[j] - Y[j+1])/(Y[N] - Y[0]);
      d = Y[j] - c*Y[N];
    }

    for (i = 0; i < M; i++) {
      RectFisMap *map = &fis->maps[j][i];
      double a, b;

      if (EVEN(i)) {
	a = (X[i+1] - X[i])/(X[M] - X[0]);
	b = X[i] - a*X[0];
      } else {
	a = (X[i] - X[i+1])/(X[M] - X[0]);
	b = X[i] - a*X[M];
      }

      map->a = a;
      map->b = b;
      map->c = c;
      map->d = d;
      map->q = map->alpha = NULL;  /* pure virtual functions */
    }
  }

  return fis;
}

void rectFisMapSequenceX(RectFis *fis, double x, int len, int seq[]) {
  double a = 1.0, b = 0.0;
  int M = fis->M;
  double *X = fis->X;
  int k;
  
  for (k = 0; k < len; k++) {
    int i;

    if (a >= 0.0)
      for (i = 0; i < M-1 && x > a*X[i+1] + b; i++)
        ;
    else
      for (i = M-1; i > 0 && x > a*X[i] + b; i--)
        ;

    b += a*fis->maps[0][i].b;
    a *= fis->maps[0][i].a;

    seq[len-k-1] = i;
  }  
}

void rectFisMapSequenceY(RectFis *fis, double y, int len, int seq[]) {
  double c = 1.0, d = 0.0;
  int N = fis->N;
  double *Y = fis->Y;
  int k;
  
  for (k = 0; k < len; k++) {
    int j;

    if (c >= 0.0)
      for (j = 0; j < N-1 && y > c*Y[j+1] + d; j++)
        ;
    else
      for (j = N-1; j > 0 && y > c*Y[j] + d; j--)
        ;

    d += c*fis->maps[j][0].d;
    c *= fis->maps[j][0].c;

    seq[len-k-1] = j;
  }
}

double rectFisEvalAux(RectFis *fis, double x0, double y0, double z0,
		      int len, int seqX[], int seqY[]) {
  RectFisMap **maps = fis->maps;
  int k;
  double p[3];

  p[0] = x0; p[1] = y0; p[2] = z0;
  for (k = 0; k < len; k++)
    rectFisTransform(&maps[seqY[k]][seqX[k]], p, p);
    
  return p[2];
}

#define SEQLEN 15

double rectFisEval(RectFis *fis, double x, double y) {
  int seqX[SEQLEN], seqY[SEQLEN];
  int i = fis->M/2, j = fis->N/2;
  rectFisMapSequenceX(fis, x, SEQLEN, seqX);
  rectFisMapSequenceY(fis, y, SEQLEN, seqY);
  return rectFisEvalAux(fis, fis->X[i], fis->Y[j], fis->Z[j][i],
			SEQLEN, seqX, seqY);
}

static double eval(Func2 *f, double x, double y) {
  return rectFisEval((RectFis *) f->data, x, y);
}

Func2 *rectFisToFunc2(RectFis *fis) {
  Func2 *f = (Func2 *) myalloc(sizeof(Func2));
  f->eval = eval;
  f->dx = f->dy = f->dxdy = NULL; /* XXXX -- later */
  f->data = fis;
  return f;
}
