#include <stdlib.h>
#include <stdio.h>
#include "func2.h"
#include "bilinear.h"
#ifdef USE_BEZIER
#include "bezfunc.h"
#else
#include "hermite.h"
#endif
#include "ris.h"

/*
 * Allocates and initializes bilinear RIS.
 * M, N : grid size
 * Sx, Sy: size of domain cell
 * x0, y0 : origin
 * W, H : size of support
 * alpha : "roughness" function (note: we only store alpha reference
 *   so don't go destroyin' it while its being used here).
 * All z[j][i] are initialized to 0.
 * zx = zy = zxy = NULL since we ignore them.
 */
Ris *createBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
                       double W, double H, Func2 *alpha) {
  Ris *ris;
  int i,j;
  double x,y, dx,dy;

  if ((ris = (Ris *) malloc(sizeof(Ris))) == NULL ||
      (ris->x = (double *) malloc((M+1)*sizeof(double))) == NULL ||
      (ris->y = (double *) malloc((N+1)*sizeof(double))) == NULL ||
      (ris->z = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (ris->maps = (RisMap **) malloc(N*sizeof(RisMap))) == NULL) {
    perror("createBilinearRis():malloc()");
    exit(-1);
  }

  ris->M = M;
  ris->N = N;

  dx = W/M;
  for (i = 0, x = x0; i <= M; i++, x += dx)
    ris->x[i] = x;

  dy = H/N;
  for (j = 0, y = y0; j <= N; j++, y += dy)
    ris->y[j] = y;

  for (j = 0; j <= N; j++) {
    if ((ris->z[j] = (double *) malloc((M+1)*sizeof(double))) == NULL) {
      perror("createBilinearRis():malloc()");
      exit(-1);
    }
    for (i = 0; i <= M; i++)
      ris->z[j][i] = 0.0;
  }

  ris->zx = ris->zy = ris->zxy = NULL;

  ris->Sx = Sx;
  ris->Sy = Sy;

  for (j = 0; j < N; j++) {
    int B;
    double c, d;

    B = j*Sy % N;
    c = 1.0/Sy; 
    d = ris->y[j] - c*ris->y[B];

    if ((ris->maps[j] = (RisMap *) malloc(M*sizeof(RisMap))) == NULL) {
      perror("createBilinearRis():malloc()");
      exit(-1);
    }

    for (i = 0; i < M; i++) {
      RisMap *map = &ris->maps[j][i];
      int L;
      double a, b;

      L = i*Sx % M;
      a = 1.0/Sy;
      b = ris->x[i] - a*ris->x[L];
      
      map->a = a;
      map->b = b;
      map->c = c;
      map->d = d;
      map->alpha = alpha;
      map->q = createBilinearFunc2(ris->x[i], ris->x[i+1], 
                                   ris->y[j], ris->y[j+1]);
    }
  }

  return ris;
}

void destroyBilinearRis(Ris *ris) {
  int i,j;

  for (j = 0; j < ris->N; j++) {
    for (i = 0; i < ris->M; i++)
      destroyBilinearFunc2(ris->maps[j][i].q);
    free(ris->maps[j]);
  }

  free(ris->maps);

  for (j = 0; j <= ris->N; j++)
    free(ris->z[j]);

  free(ris->z);
  free(ris->y);
  free(ris->x);
  free(ris);
}

/*
 * Anytime and of the z[j][i] values change or alpha() is modified
 * externally, this needs to be called.
 */
void computeBilinearRisMaps(Ris *ris) { 
  int M = ris->M, N = ris->N;
  int Sx = ris->Sx, Sy = ris->Sy;
  int i,j;

  for (j = 0; j < N; j++) {
    int B = j*Sy % N;
    for  (i = 0; i < M; i++) {
      RisMap *map = &ris->maps[j][i];
      int L = i*Sx % M;
      double q[2][2];
      int m,n, mm,nn;

      for (n = 0, nn = B; n < 2; n++, nn += Sy)
        for (m = 0, mm = L; m < 2; m++, mm += Sx)
          q[n][m] = ris->z[j+n][i+m] -
            (*map->alpha->eval)(map->alpha,ris->x[i+m],ris->y[j+n])*
            ris->z[nn][mm];
      
      setBilinearFunc2(map->q, q);
    }
  }
}

Ris *createHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                      double W, double H, Func2 *alpha) {
  Ris *ris;
  int i,j;
  double x,y, dx,dy;

  if ((ris = (Ris *) malloc(sizeof(Ris))) == NULL ||
      (ris->x = (double *) malloc((M+1)*sizeof(double))) == NULL ||
      (ris->y = (double *) malloc((N+1)*sizeof(double))) == NULL ||
      (ris->z = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (ris->zx = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (ris->zy = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (ris->zxy = (double **) malloc((N+1)*sizeof(double *))) == NULL ||
      (ris->maps = (RisMap **) malloc(N*sizeof(RisMap))) == NULL) {
    perror("createHermiteRis():malloc()");
    exit(-1);
  }

  ris->M = M;
  ris->N = N;

  dx = W/M;
  for (i = 0, x = x0; i <= M; i++, x += dx)
    ris->x[i] = x;

  dy = H/N;
  for (j = 0, y = y0; j <= N; j++, y += dy)
    ris->y[j] = y;

  for (j = 0; j <= N; j++) {
    if ((ris->z[j] = (double *) malloc((M+1)*sizeof(double))) == NULL ||
        (ris->zx[j] = (double *) malloc((M+1)*sizeof(double))) == NULL ||
        (ris->zy[j] = (double *) malloc((M+1)*sizeof(double))) == NULL ||
        (ris->zxy[j] = (double *) malloc((M+1)*sizeof(double))) == NULL) {
      perror("createHermiteRis():malloc()");
      exit(-1);
    }
    for (i = 0; i <= M; i++)
      ris->z[j][i] = ris->zx[j][i] = ris->zy[j][i] = ris->zxy[j][i] = 0.0;
  }
  
  ris->Sx = Sx;
  ris->Sy = Sy;

  for (j = 0; j < N; j++) {
    int B;
    double c, d;

    B = j*Sy % N;
    c = 1.0/Sy; 
    d = ris->y[j] - c*ris->y[B];

    if ((ris->maps[j] = (RisMap *) malloc(M*sizeof(RisMap))) == NULL) {
      perror("createHermiteRis():malloc()");
      exit(-1);
    }

    for (i = 0; i < M; i++) {
      RisMap *map = &ris->maps[j][i];
      int L;
      double a, b;

      L = i*Sx % M;
      a = 1.0/Sx;
      b = ris->x[i] - a*ris->x[L];
      
      map->a = a;
      map->b = b;
      map->c = c;
      map->d = d;
      map->alpha = alpha;
#ifdef USE_BEZIER
      map->q = createBezFunc2(ris->x[i], ris->x[i+1],
                              ris->y[j], ris->y[j+1]);
#else
      map->q = createGeneralHermiteFunc2(ris->x[i], ris->x[i+1], 
                                         ris->y[j], ris->y[j+1]);
#endif
    }
  }

  return ris;
}

void destroyHermiteRis(Ris *ris) {
  int i,j;

  for (j = 0; j < ris->N; j++) {
    for (i = 0; i < ris->M; i++)
#ifdef USE_BEZIER
      destroyBezFunc2(ris->maps[j][i].q);
#else
      destroyGeneralHermiteFunc2(ris->maps[j][i].q);
#endif
    free(ris->maps[j]);
  }

  free(ris->maps);

  for (j = 0; j <= ris->N; j++) {
    free(ris->z[j]);
    free(ris->zx[j]);
    free(ris->zy[j]);
    free(ris->zxy[j]);
  }

  free(ris->z);
  free(ris->zx);
  free(ris->zy);
  free(ris->zxy);
  free(ris->y);
  free(ris->x);
  free(ris);  
}

void computeHermiteRisMaps(Ris *ris) { 
  int M = ris->M, N = ris->N;
  int Sx = ris->Sx, Sy = ris->Sy;
  int i,j;

  for (j = 0; j < N; j++) {
    int B = j*Sy % N;
    for  (i = 0; i < M; i++) {
      RisMap *map = &ris->maps[j][i];
      int L = i*Sx % M;
      double q[2][2], qx[2][2], qy[2][2], qxy[2][2];
      int m,n, mm,nn;

      for (n = 0, nn = B; n < 2; n++, nn += Sy)
        for (m = 0, mm = L; m < 2; m++, mm += Sx) {
          double x_ = ris->x[i+m];
          double y_ = ris->y[j+n];
          double z_ = ris->z[j+n][i+m];
          double zx_ = ris->zx[j+n][i+m];
          double zy_ = ris->zy[j+n][i+m];
          double z = ris->z[nn][mm];
          double zx = ris->zx[nn][mm];
          double zy = ris->zy[nn][mm];
          double alpha = (*map->alpha->eval)(map->alpha, x_, y_);
          double alpha_x = (*map->alpha->dx)(map->alpha, x_, y_);
          double alpha_y = (*map->alpha->dy)(map->alpha, x_, y_);
          q[n][m] = z_ - alpha*z;
          qx[n][m] = zx_ - (Sx*alpha*zx + alpha_x*z);
          qy[n][m] = zy_ - (Sy*alpha*zy + alpha_y*z);
          qxy[n][m] = ris->zxy[j+n][i+m] - 
            (Sx*Sy*alpha*ris->zxy[nn][mm] + 
             Sx*alpha_y*zx + Sy*alpha_x*zy +
             (*map->alpha->dxdy)(map->alpha, x_, y_)*z);
        }

#ifdef USE_BEZIER
      setBezFunc2FromDerivs(map->q, q, qx, qy, qxy); 
#else
      setGeneralHermiteFunc2(map->q, q, qx, qy, qxy);
#endif
    }
  }
}

void risMapSeqX(Ris *ris, double x, int n, int seq[]) {
  double invdx = 1.0/(ris->x[1] - ris->x[0]);
  int k;

  for (k = 0; k < n; k++) {
    int L, i = (int) ((x - ris->x[0])*invdx);
    if (i < 0) i = 0;  /* avoid numerical error */
    else if (i >= ris->M) i = ris->M-1;
    L = i*ris->Sx % ris->M;
    seq[n-k-1] = i;
    x = (x - ris->x[i])*ris->Sx + ris->x[L];
  }
}

void risMapSeqY(Ris *ris, double y, int n, int seq[]) {
  double invdy = 1.0/(ris->y[1] - ris->y[0]);
  int k;

  for (k = 0; k < n; k++) {
    int B, j = (int) ((y - ris->y[0])*invdy);
    if (j < 0) j = 0;  /* avoid numerical error */
    else if (j >= ris->N) j = ris->N-1;
    B = j*ris->Sy % ris->N;
    seq[n-k-1] = j;
    y = (y - ris->y[j])*ris->Sy + ris->y[B];
  }
}

double risMapTransform(RisMap *map, double src[3], double dst[3]) {
  dst[0] = map->a*src[0] + map->b;
  dst[1] = map->c*src[1] + map->d;
  dst[2] = (*map->q->eval)(map->q, dst[0], dst[1]) +
    (*map->alpha->eval)(map->alpha, dst[0], dst[1])*src[2];
  return dst[2];
}

/*
 * Transform src = (x, y, z, zx, zy, zxy)
 * into dst = (x', y', z', zx', zy', zxy')
 * where
 * x' = a*x + b
 * y' = c*y + d
 * z' = q(x',y') + alpha(x',y')*z
 * zx' = qx(x',y') + Sx*alpha(x',y')*zx + alpha_x(x',y')*z
 * zy' = qy(x',y') + Sy*alpha(x',y')*zy + alpha_y(x',y')*z
 * zxy' = qxy(x',y') + Sx*Sy*alpha(x',y')*zxy +
 *           Sx*alpha_y(x',y')*zx + Sy*alpha_x(x',y')*zy
 *           alpha_xy(x',y')*z
 */
void risMapTransformWithDerivs(RisMap *map, int Sx, int Sy,
                               double src[6], double dst[6]) {
  double x_ = map->a*src[0] + map->b;
  double y_ = map->c*src[1] + map->d;
  double alpha = (*map->alpha->eval)(map->alpha, x_, y_);
  double alpha_x = (*map->alpha->dx)(map->alpha, x_, y_);
  double alpha_y = (*map->alpha->dy)(map->alpha, x_, y_);
  double z_ = (*map->q->eval)(map->q, x_, y_) + alpha*src[2];
  double zx_ = (*map->q->dx)(map->q, x_, y_) + 
    Sx*alpha*src[3] + alpha_x*src[2];
  double zy_ = (*map->q->dy)(map->q, x_, y_) + 
    Sy*alpha*src[4] + alpha_y*src[2];
  double zxy_ = (*map->q->dxdy)(map->q, x_, y_) + Sx*Sy*alpha*dst[5] +
    Sx*alpha_y*src[3] + Sy*alpha_x*src[4] +
    (*map->alpha->dxdy)(map->alpha, x_, y_)*src[2];
  dst[0] = x_;
  dst[1] = y_;
  dst[2] = z_;
  dst[3] = zx_;
  dst[4] = zy_;
  dst[5] = zxy_;
}

double risEvalSeq(Ris *ris, int n, int seqx[], int seqy[]) {
  int I = seqx[0]*ris->Sx % ris->M + ris->Sx/2;
  int J = seqy[0]*ris->Sy % ris->N + ris->Sy/2;
  int k;
  double p[3];

  p[0] = ris->x[I];
  p[1] = ris->y[J];
  p[2] = ris->z[J][I];

  for (k = 0; k < n; k++)
    risMapTransform(&ris->maps[seqy[k]][seqx[k]], p, p);
    
  return p[2];
}

/*
 * Given map sequence, computes
 * z = f(x,y)
 * zx = fx(x,y)
 * zy = fy(x,y)
 * zxy = fxy(x,y)
 */
void risEvalSeqWithDerivs(Ris *ris, int n, int seqx[], int seqy[],
                          double *z, double *zx, double *zy, double *zxy) {
  int I = seqx[0]*ris->Sx % ris->M + ris->Sx/2;
  int J = seqy[0]*ris->Sy % ris->N + ris->Sy/2;
  int k;
  double p[6];

  p[0] = ris->x[I];
  p[1] = ris->y[J];
  p[2] = ris->z[J][I];
  p[3] = (ris->zx == NULL) ? 0.0 : ris->zx[J][I];
  p[4] = (ris->zy == NULL) ? 0.0 : ris->zy[J][I];
  p[5] = (ris->zxy == NULL) ? 0.0 : ris->zxy[J][I];

  for (k = 0; k < n; k++)
    risMapTransformWithDerivs(&ris->maps[seqy[k]][seqx[k]], 
			      ris->Sx, ris->Sy, p, p);

  *z = p[2];
  *zx = p[3];
  *zy = p[4];
  *zxy = p[5];
}

#define SEQLEN 12

double risEval(Ris *ris, double x, double y) {
  int seqx[SEQLEN], seqy[SEQLEN];
  risMapSeqX(ris, x, SEQLEN, seqx);
  risMapSeqY(ris, y, SEQLEN, seqy);
  return risEvalSeq(ris, SEQLEN, seqx, seqy);
}

void resEvalWithDerivs(Ris *ris, double x, double y, 
                       double *z, double *zx, double *zy, double *zxy) {
  int seqx[SEQLEN], seqy[SEQLEN];
  risMapSeqX(ris, x, SEQLEN, seqx);
  risMapSeqY(ris, y, SEQLEN, seqy);
  risEvalSeqWithDerivs(ris, SEQLEN, seqx, seqy, z, zx, zy, zxy);
}

static double eval(Func2 *this, double x, double y) {
  Ris *ris = (Ris *) this->data;
  return risEval(ris, x, y);
}

/*
 * Derivative computation is very redundant.
 * Perhaps cash (x,y) and (z,zx,zy,zx) data if we
 * ever need a speed up. Try to preserve reentrancy.
 */ 
static double dx(Func2 *this, double x, double y) {
  Ris *ris = (Ris *) this->data;
  double z, zx, zy, zxy;
  resEvalWithDerivs(ris, x, y, &z, &zx, &zy, &zxy);
  return zx;
}

static double dy(Func2 *this, double x, double y) {
  Ris *ris = (Ris *) this->data;
  double z, zx, zy, zxy;
  resEvalWithDerivs(ris, x, y, &z, &zx, &zy, &zxy);
  return zy;
}

static double dxdy(Func2 *this, double x, double y) {
  Ris *ris = (Ris *) this->data;
  double z, zx, zy, zxy;
  resEvalWithDerivs(ris, x, y, &z, &zx, &zy, &zxy);
  return zxy;
}

Func2 *risToFunc2(Ris *ris) {
  Func2 *f;

  if ((f = (Func2 *) malloc(sizeof(Func2))) == NULL) {
    perror("risToFunc2():malloc()");
    exit(-1);
  }

  f->eval = eval;
  f->dx = dx;
  f->dy = dy;
  f->dxdy = dxdy;
  f->data = ris;
  
  return f;
}
