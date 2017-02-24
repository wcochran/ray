#ifdef WIN32
#include <windows.h>
#endif
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "func2.h"
#include "ris.h"
#include "viewpick.h"

/* #define FOR_HARDCOPY */  /* white background, gray-level image */

#define GRID_SIZE   6    /* default grid size */
#define DOMAIN_SIZE 2    /* default domain cell size */
#define ALPHA     -0.25    /* default alpha value */

#define SAMPLES      50  /* number of samples in low resolution surface */
#define BIG_SAMPLES 200  /* number of samples in high resolution surface */

#ifdef COMMENT_OUT_XXXX
#define SIZE 6.0         /* width of function support */
#define ORG (-SIZE/2)    /* x0 = y0 = ORG = origin of function support */
#endif
#define SIZE 1.0         /* width of function support */
#define ORG (-SIZE/2)    /* x0 = y0 = ORG = origin of function support */

#define EPSILON 0.0000001

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

int useHermiteRis = GL_FALSE;  /* Using bilinear or hermite RIS? */
Func2 *alpha;           /* "roughness" function. */
Ris *ris;               /* Recurrent Interpolation Surface (RIS) */

/*
 * Surface material properties
 */
GLfloat material_ambient[4]  = {0.2, 0.2, 0.2, 1.0};
#ifdef FOR_HARDCOPY
GLfloat material_diffuse[4]  = {0.4, 0.4, 0.4, 1.0};
#else
GLfloat material_diffuse[4]  = {0.1, 0.5, 0.8, 1.0};
#endif
GLfloat material_specular[4] = {0.1, 0.1, 0.1, 1.0};
GLfloat material_phong = 30.0;

/*
 * Lighting
 */
GLfloat lmodel_ambient[4]    = {0.1, 0.1, 0.1, 1.0};
GLfloat light0_position[4]   = {1.0, 1.0, 1.0, 0.0};
GLfloat light0_ambient[4]    = {0.3, 0.3, 0.3, 1.0};
GLfloat light0_diffuse[4]    = {1.0, 1.0, 1.0, 1.0};
GLfloat light0_specular[4]   = {1.0, 1.0, 1.0, 1.0};  

/*
 * Grids used for rendering the surface.
 * The 'samples' grid is the smaller grid which can (hopefully) be
 * rendered at interactive frame rates. The 'bigSamples' grid
 * contains more detail.
 */
int useBigGrid = GL_TRUE;                  /* render w/ "hi res" grid? */
int dirtyGrid = GL_TRUE;                   /* recompute "lo res" grid? */
double grid[SAMPLES*SAMPLES];              /* lo-res samples grid */
int dirtyBigGrid = GL_TRUE;                /* recompute "hi res" grid? */
double bigGrid[BIG_SAMPLES*BIG_SAMPLES];   /* lo-res samples grid */

int drawControls = GL_TRUE;                /* draw surface controls? */
GLfloat controlColor[4] = {0.7, 0.7, 0.7, 1.0}; /* control color */

/*
 * Camera setup:
 * Eye stored in spherical coords with the lookat point at the center;
 * Eye converted to cartestion coords for lookat transformation.
 */ 
GLdouble eyeRadius, eyePhi, eyeTheta;
GLdouble eye[3], lookat[3], up[3];

/*
 * Utility for mouse selection and dragging of controls based
 * on glViewport(), gluLookAt(), and gluPerspective() info.
 */
ViewPick viewpick;

/*
 * For controlling displayed text message
 */
GLfloat messagePos[3];
char *message = NULL;
GLfloat messageColor[4] = {0.7, 0.7, 0.7, 1.0};

/*
 * Init RIS with a sinusoidal surface.
 */
void initRis(void) {
  int M = ris->M, N = ris->N;
  double x,y, dx,dy, DX,DY;
  int i,j;

  DX = ris->x[M] - ris->x[0];
  DY = ris->y[N] - ris->y[0];

  dx = DX/M;
  dy = DY/N;

#define SQRT_AMP SIZE/6
#define CYCLES   0.5
  
  for (j = 0, y = ris->y[0]; j < N; j++, y += dy) {
    double fy = SQRT_AMP*sin(CYCLES*2.0*M_PI*(y - ris->x[0])/DY);
    for (i = 0, x = ris->x[0]; i < M; i++, x += dx) {
      double fx = SQRT_AMP*sin(CYCLES*2.0*M_PI*(x - ris->y[0])/DX);
      ris->z[j][i] = fx + fy;
    }
  }

  for (j = 0; j <= N; j++)   /* match boundary data */
    ris->z[j][M] = ris->z[j][0];

  for (i = 0; i <= M; i++)
    ris->z[N][i] = ris->z[0][i];

  computeBilinearRisMaps(ris);
}

/*
 * This gives random values to the RIS controls.
 */
void randomizeRis(void) {
  int M = ris->M, N = ris->N;
  int i,j;

#define ZSCALE SIZE/6
  
  for (j = 0; j < N; j++)
    for (i = 0; i < M; i++)
      ris->z[j][i] = (2*drand48() - 1)*ZSCALE;

  for (j = 0; j <= N; j++)   /* match boundary data */
    ris->z[j][M] = ris->z[j][0];

  for (i = 0; i <= M; i++)
    ris->z[N][i] = ris->z[0][i];
    
#define DZSCALE SIZE/12
  
  if (ris->zx != NULL) {
    for (j = 0; j < N; j++) {
      for (i = 0; i < M; i++) {
        ris->zx[j][i] = (2*drand48() - 1)*DZSCALE;
        ris->zy[j][i] = (2*drand48() - 1)*DZSCALE;
        for (j = 0; j <= N; j++) {  /* match boundary data */
          ris->zx[j][M] = ris->zx[j][0];
          ris->zy[j][M] = ris->zy[j][0];
          ris->zxy[j][M] = ris->zxy[j][0];
        }
        for (i = 0; i <= M; i++) {
          ris->zx[j][M] = ris->zx[j][0];
          ris->zy[j][M] = ris->zy[j][0];
          ris->zxy[j][M] = ris->zxy[j][0];
        }
      }
    }
    computeHermiteRisMaps(ris);
  } else {
    computeBilinearRisMaps(ris);
  }  
}

#define SEQLEN 10      /* default sequence length for surface evaluation */

/*
 * Converts Bilinear to Hermite ris.
 * Global ris is altered; old ris is destroyed.
 */
void bilinearToHermiteRis(void) {
  int M = ris->M, N = ris->N;
  Ris *hris;
  int i,j;

  hris = createHermiteRis(M, N, ris->Sx, ris->Sy,
                          ris->x[0], ris->y[0],
                          ris->x[M] - ris->x[0], ris->y[N] - ris->y[0],
                          ris->maps[0][0].alpha);

  for (j = 0; j < N; j++) {
    int prevj = (j == 0) ? N : j-1;
    int nextj = j+1;
    double dy = ris->y[nextj] - ris->y[prevj];
    for (i = 0; i < M; i++) {
      int previ = (i == 0) ? M : i-1;
      int nexti = i+1;
      double dx = ris->x[nexti] - ris->x[previ];
      hris->z[j][i] = ris->z[j][i];
      hris->zx[j][i] = (ris->z[j][nexti] - ris->z[j][previ])/dx;
      hris->zy[j][i] = (ris->z[nextj][i] - ris->z[prevj][i])/dy;
      hris->zxy[j][i] = 0.0;
    }
  }
      
  for (i = 0; i < M; i++) {        /* match top edge with bottom */
    hris->z[N][i] = hris->z[0][i];
    hris->zx[N][i] = hris->zx[0][i];
    hris->zy[N][i] = hris->zy[0][i];
    hris->zxy[N][i] = hris->zxy[0][i];
  }

  for (j = 0; j <= N; j++) {        /* match right edge with left */
    hris->z[j][M] = hris->z[j][0];
    hris->zx[j][M] = hris->zx[j][0];
    hris->zy[j][M] = hris->zy[j][0];
    hris->zxy[j][M] = hris->zxy[j][0];
  }

  computeHermiteRisMaps(hris);
  
  destroyBilinearRis(ris);
  useHermiteRis = GL_TRUE;
  ris = hris;
}

/*
 * Convert Hermite RIS to Bilinear RIS
 * Global ris is altered; old ris is destroyed.
 */ 
void hermiteToBilinearRis(void) {
  int M = ris->M, N = ris->N;
  Ris *bris;
  int i,j;

  bris = createBilinearRis(M, N, ris->Sx, ris->Sy,
                           ris->x[0], ris->y[0],
                           ris->x[M] - ris->x[0], ris->y[N] - ris->y[0],
                           ris->maps[0][0].alpha);

  for (j = 0; j <= N; j++)
    for (i = 0; i <= M; i++)
      bris->z[j][i] = ris->z[j][i];

  computeBilinearRisMaps(bris);

  destroyHermiteRis(ris);
  useHermiteRis = GL_FALSE;
  ris = bris;
}

void getSurfaceSamples(int samples, double *grid) {
  int M = ris->M, N = ris->N;
  double x,y, dx,dy;
  int i,j;

  dx = (ris->x[M] - ris->x[0])/(samples-1);
  dy = (ris->y[N] - ris->y[0])/(samples-1);

  for (j = 0, y = ris->y[0]; j < samples; j++, y += dy) {
    int seqy[SEQLEN];
    risMapSeqY(ris, y, SEQLEN, seqy);
    for (i = 0, x = ris->x[0]; i < samples; i++, x += dx) {
      int seqx[SEQLEN];
      risMapSeqX(ris, x, SEQLEN, seqx);
      grid[samples*j + i] = risEvalSeq(ris, SEQLEN, seqx, seqy);
    }
  }    
}

#define DOT(A,B) ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2])
#define SCALE(k,V) ((V)[0] *= k, (V)[1] *= k, (V)[2] *= k)

/*
 *   Renders surface as triangle strips row by row as follows:
 *
 *     i    i+1
 *     *----*----*----*----*----*----*----*  j
 *     |   /|   /|   /|   /|   /|   /|   /|
 *     | /  | /  | /  | /  | /  | /  | /  |
 *     *----*----*----*----*----*----*----*  j+1
 *
 *   Triangles are rendered using flat shading (i.e. all 3 vertices
 *   share the same surface normal).
 */
void drawGrid(int samples, double *grid) {
  int M = ris->M, N = ris->N;
  double x,y, dx,dy, dxdy;
  int i,j;

  dx = (ris->x[M] - ris->x[0])/(samples-1);
  dy = (ris->y[N] - ris->y[0])/(samples-1);
  dxdy = dx*dy;

  for (j = 0, y = ris->y[0]; j < samples-1; j++, y += dy) {
    glBegin(GL_TRIANGLES);
    for (i = 0, x = ris->x[0]; i < samples-1; i++, x += dx) {
      double s;
      double n[3];
      double z[2][2];

#define Z(i,j) grid[(j)*samples + (i)]
      z[0][0] = Z(i,j);
      z[0][1] = Z(i+1,j);
      z[1][0] = Z(i,j+1);
      z[1][1] = Z(i+1,j+1);
#undef Z

      n[0] = -dy*(z[0][1] - z[0][0]);
      n[1] = -dx*(z[1][0] - z[0][0]);
      n[2] = dxdy;
      s = 1.0/sqrt(DOT(n,n));
      SCALE(s,n);
      glNormal3dv(n);
      glVertex3d(x, y, z[0][0]);
      glVertex3d(x, y+dy, z[1][0]);
      glVertex3d(x+dx, y, z[0][1]);

      n[0] = dy*(z[1][0] - z[1][1]);
      n[1] = dx*(z[0][1] - z[1][1]);
      n[2] = dxdy;
      s = 1.0/sqrt(DOT(n,n));
      SCALE(s,n);
      glNormal3dv(n);
      glVertex3d(x+dx, y, z[0][1]);
      glVertex3d(x, y+dy, z[1][0]);
      glVertex3d(x+dx, y+dy, z[1][1]);
    }
    glEnd();
  }    
}


int drawq = GL_FALSE;

/*
 * Render q surface's only -- for debugging purposes.
 */
void qSurfaces(int samples) {
  int i,j;
  double x,y, dx,dy;

  dx = (ris->x[1] - ris->x[0])/(samples-1);
  dy = (ris->y[1] - ris->y[0])/(samples-1);
  
  for (j = 0; j < ris->N; j++)
    for (i = 0; i < ris->M; i++) {
      Func2 *q = ris->maps[j][i].q;
      int ii, jj;
      for (jj = 1, y = ris->y[j]; jj < samples; jj++, y += dy) {
        if (jj == samples-1)
          y = ris->y[j+1] - dy;
        
        glBegin(GL_QUAD_STRIP);
        for (ii = 1, x = ris->x[i]; ii <= samples; ii++, x += dx) {
          double norm[3], s;

          if (ii == samples)
            x = ris->x[i+1];
          
          norm[0] = -(*q->dx)(q, x, y);
          norm[1] = -(*q->dy)(q, x, y);
          norm[2] = 1.0;
          s = 1.0/sqrt(DOT(norm, norm));
          SCALE(s, norm);
          glNormal3dv(norm);
          glVertex3d(x, y, (*q->eval)(q, x, y));

          norm[0] = -(*q->dx)(q, x, y+dy);
          norm[1] = -(*q->dy)(q, x, y+dy);
          norm[2] = 1.0;
          s = 1.0/sqrt(DOT(norm, norm));
          SCALE(s, norm);
          glNormal3dv(norm);
          glVertex3d(x, y+dy, (*q->eval)(q, x, y+dy));
        }
        glEnd();
        
      }      
    }
}

#define Q_SAMPLES 240
int dirty_qgrid = GL_TRUE;
double qgrid[Q_SAMPLES*Q_SAMPLES];

void get_qsamples(void) {
  int i,j;
  double x,y, dx,dy;

  dx = (ris->x[ris->M] - ris->x[0])/(Q_SAMPLES-1);
  dy = (ris->y[ris->N] - ris->y[0])/(Q_SAMPLES-1);

  for (j = 0, y = ris->y[0]; j < Q_SAMPLES; j++, y += dy) {
    int mj = (int) ((y - ris->y[0])/(ris->y[1] - ris->y[0]));
    if (mj >= ris->N) mj = ris->N-1;
    for (i = 0, x = ris->x[0]; i < Q_SAMPLES; i++, x += dx) {
      int mi = (int) ((x - ris->x[0])/(ris->x[1] - ris->x[0]));
      RisMap *map;
      if (mi >= ris->M) mi = ris->M-1;
      map = &ris->maps[mj][mi];
      qgrid[j*Q_SAMPLES + i] = (*map->q->eval)(map->q, x, y);
    }
  }
}

/*
 * Render RIS surface with random surface points generated
 * by the (modified) chaos game.
 */
void chaos(int iterations) {
  double p[6];
  double DX, DY;
  int domx, domy;

  /*
   * Choose initial random point on surface.
   * The point is enhanced to hold derivative info as well.
   */
  p[0] = ris->x[0];
  p[1] = ris->y[0];
  p[2] = ris->z[0][0];
  if (ris->zx != NULL) {    /* hermite */
    p[3] = ris->zx[0][0];
    p[4] = ris->zy[0][0];
    p[5] = ris->zxy[0][0];
  } else {                  /* bilinear */
    p[3] = p[4] = p[5] = 0.0;
  }

  DX = ris->x[ris->Sx] - ris->x[0];   /* size of domain cell */
  DY = ris->y[ris->Sy] - ris->y[0];
    
  /*
   * domx, domy : number of unique domain cells in x & y;
   * This is also the vertical and horizontal spread between
   * range cells that correspond to the same domain cell.
   */
  domx = ris->M/ris->Sx;
  domy = ris->N/ris->Sy;

  /*
   * Play the chaos game and plot points with normals.
   */
  glBegin(GL_POINTS);
  while (--iterations >= 0) {
    int i,j, I,J, rx, ry;
    double s, n[3];
    
    /*
     * Find (lower-left) range cell index.
     */
    I = (int) ((p[0] - ris->x[0])/DX);
    if (I >= domx) I = 0;
    J = (int) ((p[1] - ris->y[0])/DY);
    if (J >= domy) J = 0;

    /*
     * Choose of the corresponding range cells (there are Sx*Sy of them).
     */
#define RANDOM(n) ((int) ((rand()/(double)RAND_MAX)*(n)))
    rx = RANDOM(ris->Sx);
    ry = RANDOM(ris->Sy);
    i = I + rx*domx;
    j = J + ry*domy;

    /*
     * Transform point with map corresponding to selected range cell.
     */
    risMapTransformWithDerivs(&ris->maps[j][i], ris->Sx, ris->Sy, p, p);

    /*
     * Compute point normal based on derivative info,
     * and send the normal and point down the OpenGL pipeline.
     */
    n[0] = -p[3];   /* normal = (-fx, -fy, 1) */
    n[1] = -p[4];
    n[2] = 1.0;
    s = 1.0/sqrt(DOT(n,n));
    SCALE(s,n);
    glNormal3dv(n);
    glVertex3dv(p);
  }
  glEnd();
}

int useChaos = GL_FALSE;
int chaosIters = 2e6;

/*
 * display callback
 */
void display(void) {

  /*
   * Clear frame buffer and z-buffer
   */
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  /*
   * Render either the lo-res or hi-res surface sample grid
   * using z-buffer and flat shading.
   */
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_LIGHTING);
  glShadeModel(GL_FLAT);

  if (drawq) {
    /* XXX
    qSurfaces(40);
    */
    if (dirty_qgrid) {
      get_qsamples();
      dirty_qgrid = GL_FALSE;
    }
    drawGrid(Q_SAMPLES, qgrid);
    glutSwapBuffers();
    return;
  }

  if (useChaos) {
    glPointSize(1.0);
    chaos(chaosIters);
  } else if (useBigGrid) {
    if (dirtyBigGrid) {
      getSurfaceSamples(BIG_SAMPLES, bigGrid);
      dirtyBigGrid = GL_FALSE;
    }
    drawGrid(BIG_SAMPLES, bigGrid);
  } else {
    if (dirtyGrid) {
      getSurfaceSamples(SAMPLES, grid);
      dirtyGrid = GL_FALSE;
    }
    drawGrid(SAMPLES, grid);
  }

  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  /*
   * Unless disabled, draw surface interpolation points.
   */
  if (drawControls) {
    int M = ris->M, N = ris->N;
    int i,j;

    glColor4fv(controlColor);
    glPointSize(5.0);
    glBegin(GL_POINTS);
    for (j = 0; j <= N; j++)
      for (i = 0; i <= M; i++)
        glVertex3d(ris->x[i], ris->y[j], ris->z[j][i]);
    glEnd();

    if (useHermiteRis) {  /* draw surface derivative controls */
#define DELTA (0.05*SIZE)
      glBegin(GL_LINES);
      for (j = 0; j <= N; j++)
        for (i = 0; i <= M; i++) {
          glVertex3d(ris->x[i], ris->y[j], ris->z[j][i]);
          glVertex3d(ris->x[i] + DELTA, ris->y[j],
                     ris->z[j][i] + ris->zx[j][i]*DELTA);
          glVertex3d(ris->x[i], ris->y[j], ris->z[j][i]);
          glVertex3d(ris->x[i], ris->y[j] + DELTA,
                     ris->z[j][i] + ris->zy[j][i]*DELTA);
        }
      glEnd();
    }
  }
  
  /*
   * Display text message if there is one
   */
  if (message != NULL) {
    char *ch;
    /* glColor4fv(controlColor); XXX */
    glColor3f(1.0, 1.0, 1.0);
    glRasterPos3fv(messagePos);
    for (ch = message; *ch != '\0'; ch++)
      glutBitmapCharacter(GLUT_BITMAP_9_BY_15, *ch);
  }
  
  glutSwapBuffers();
}

/*
 *Convert point in sperical coordinates to cartesion coordinate
 */
void sphericalToCartesian(double r, double theta, double phi,
                          double *x, double *y, double *z) {
  double sin_phi = sin(phi);
  *x = r*cos(theta)*sin_phi;
  *y = r*sin(theta)*sin_phi;
  *z = r*cos(phi);
}

double fovy, aspect, zNear, zFar;   /* perspective projecion parameters */

/*
 * Initialize rendering setup.
 */
void initGraphics() {
  int M = ris->M, N = ris->N;
  double diag,dx,dy;

  /*
   * Determine reasonable initial lookat transformation
   */
  dx = ris->x[M] - ris->x[0];
  dy = ris->y[N] - ris->y[0];
  diag = sqrt(dx*dx + dy*dy);

  lookat[0] = 0.5*(ris->x[M] + ris->x[0]);
  lookat[1] = 0.5*(ris->y[N] + ris->y[0]);
  lookat[2] = 0.0;

  /* eyeRadius = 2*diag; XXX */
  eyeRadius = 1.5*diag;
  eyeTheta = M_PI/4.0;
  eyePhi = M_PI/4.0;
  sphericalToCartesian(eyeRadius, eyeTheta, eyePhi, 
                       &eye[0], &eye[1], &eye[2]);
  eye[0] += lookat[0];
  eye[1] += lookat[1];
  eye[2] += lookat[2];

  up[0] = up[1] = 0.0;
  up[2] = 1.0;
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye[0], eye[1], eye[2],
            lookat[0], lookat[1], lookat[2],
            up[0], up[1], up[2]);
  setLookAt(&viewpick, eye, lookat, up);

  /*
   * Choose nice perspective projection
   */
  fovy = 40.0;
  aspect = 1.0;
  zNear = 0.01*diag;
  zFar = 10*diag;
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(fovy, aspect, zNear, zFar);
  setPerspective(&viewpick, fovy, aspect, zNear, zFar);

  /*
   * Flat shade surface
   */
  glFrontFace(GL_CCW);
  glPolygonMode(GL_FRONT, GL_FILL);
  glShadeModel(GL_FLAT); /* actually set is display() */

  /*
   * Set up surface material and lighting parameters
   */
  glMaterialfv(GL_FRONT, GL_AMBIENT, material_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, &material_phong);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glEnable(GL_LIGHT0);

#ifdef FOR_HARDCOPY
  glClearColor(1.0, 1.0, 1.0, 1.0);
#endif
}

/*
 * Reshape callback.
 * Proper initialization done on first call.
 * From then on we just change the viewport.
 */
void reshape(int w, int h) {
  static int first = GL_TRUE;

  if (first) {
    initGraphics();
    first = GL_FALSE;
  }

  setViewport(&viewpick, w, h);
  glViewport(0,0, w,h);
}

int mouseRotate = GL_FALSE;  /* currently rotating camera view? */
int mousex, mousey;          /* used to track mouse position */
int knoti, knotj;            /* current knot being "dragged" */
enum {
  DRAG_NONE, DRAG_KNOT, DRAG_DX, DRAG_DY, DRAG_ALPHA
} drag = DRAG_NONE;          /* how knot is being altered */

#define ALPHA_BAR 200        /* floating vertical slider length */
int alphaBarMin;             /* top of slider => alpha = 1 */

/*
 * mouse callback
 * middle-button rotate camera about center of view
 * left-button drags controls
 * shift-left-button drags alpha-val
 */
void mouse(int button, int state, int x, int y) {
  if (state == GLUT_DOWN) {
    if (button == GLUT_MIDDLE_BUTTON) {  /* rotate camera */
      mouseRotate = GL_TRUE;
      mousex = x;
      mousey = y;
      useBigGrid = GL_FALSE;
    } else if (button == GLUT_LEFT_BUTTON) { 
      if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) { /* drag alpha val */
        int yy = viewpick.winHeight - y;
        alphaBarMin = yy - (int) ((*alpha->eval)(alpha, 0,0)*ALPHA_BAR);
        useBigGrid = GL_FALSE;
        drag = DRAG_ALPHA;
      } else { /* see if knot selected */
#define KNOT_RADIUS 0.001*SIZE /* 0.05 XXXX */
        int i,j, M = ris->M, N = ris->N;
        double dir[3];
        screenToRay(&viewpick, x, y, dir);
        
        for (j = 0; j <= N; j++) {
          for (i = 0; i <= M; i++) {
            double dist, knot[3];
            
            knot[0] = ris->x[i];        /* interpolation point */
            knot[1] = ris->y[j];
            knot[2] = ris->z[j][i];
            dist = dist2PointToRay(eye, dir, knot);
            if (dist < KNOT_RADIUS) { 
              knoti = i;
              knotj = j;
              useBigGrid = GL_FALSE;
              drag = DRAG_KNOT;
              return;
            }

            if (useHermiteRis) {
              knot[0] = ris->x[i] + DELTA;   /* df/dx control */
              knot[1] = ris->y[j];
              knot[2] = ris->z[j][i] + DELTA*ris->zx[j][i];
              dist = dist2PointToRay(eye, dir, knot);
              if (dist < KNOT_RADIUS) {
                knoti = i;
                knotj = j;
                useBigGrid = GL_FALSE;
                drag = DRAG_DX;
                return;
              }
              
              knot[0] = ris->x[i];          /* df/dy control */
              knot[1] = ris->y[j] + DELTA;
              knot[2] = ris->z[j][i] + DELTA*ris->zy[j][i];
              dist = dist2PointToRay(eye, dir, knot);
              if (dist < KNOT_RADIUS) {
                knoti = i;
                knotj = j;
                useBigGrid = GL_FALSE;
                drag = DRAG_DY;
                return;
              }
            }
            
          }
        }
      }
    } 
  } else if (state == GLUT_UP) {
    if (mouseRotate) {
      mouseRotate = GL_FALSE;
      useBigGrid = GL_TRUE;
      glutPostRedisplay();
    } else if (drag != DRAG_NONE) {
      drag = DRAG_NONE;
      message = NULL;
      useBigGrid = GL_TRUE;
      dirtyBigGrid = GL_TRUE;
      dirty_qgrid = GL_TRUE;
      glutPostRedisplay();      
    }
  }
}

/*
 * Given possible RIS edge index (i,j), this makes
 * sure that the proper data is duplicated.
 * This got ugly -- oh well.
 */
void matchEdgeVal(int i, int j) {
  if (j == 0) {
    ris->z[ris->N][i] = ris->z[j][i];
    if (useHermiteRis) {
      ris->zx[ris->N][i] = ris->zx[j][i];
      ris->zy[ris->N][i] = ris->zy[j][i];
      ris->zxy[ris->N][i] = ris->zxy[j][i];
    }
    if (i == 0) {
      ris->z[ris->N][ris->M] = ris->z[j][i];
      if (useHermiteRis) {
        ris->zx[ris->N][ris->M] = ris->zx[j][i];
        ris->zy[ris->N][ris->M] = ris->zy[j][i];
        ris->zxy[ris->N][ris->M] = ris->zxy[j][i];
      }
    } else if (i == ris->M) {
      ris->z[ris->N][0] = ris->z[j][i];
      if (useHermiteRis) {
        ris->zx[ris->N][0] = ris->zx[j][i];
        ris->zy[ris->N][0] = ris->zy[j][i];
        ris->zxy[ris->N][0] = ris->zxy[j][i];
      }
    }
  } else if (j == ris->N) {
    ris->z[0][i] = ris->z[j][i];
    if (useHermiteRis) {
      ris->zx[0][i] = ris->zx[j][i];
      ris->zy[0][i] = ris->zy[j][i];
      ris->zxy[0][i] = ris->zxy[j][i];
    }
    if (i == 0) {
      ris->z[0][ris->M] = ris->z[j][i];
      if (useHermiteRis) {
        ris->zx[0][ris->M] = ris->zx[j][i];
        ris->zy[0][ris->M] = ris->zy[j][i];
        ris->zxy[0][ris->M] = ris->zxy[j][i];
      }
    } else if (i == ris->M) {
      ris->z[0][0] = ris->z[j][i];
      if (useHermiteRis) {
        ris->zx[0][0] = ris->zx[j][i];
        ris->zy[0][0] = ris->zy[j][i];
        ris->zxy[0][0] = ris->zxy[j][i];
      }
    }
  }
  if (i == 0) {
    ris->z[j][ris->M] = ris->z[j][i];
    if (useHermiteRis) {
      ris->zx[j][ris->M] = ris->zx[j][i];
      ris->zy[j][ris->M] = ris->zy[j][i];
      ris->zxy[j][ris->M] = ris->zxy[j][i];
    }
  } else if (i == ris->M) {
    ris->z[j][0] = ris->z[j][i];
    if (useHermiteRis) {
      ris->zx[j][0] = ris->zx[j][i];
      ris->zy[j][0] = ris->zy[j][i];
      ris->zxy[j][0] = ris->zxy[j][i];
    }
  }
}    

/*
 * mouse motion callback
 * rotating view
 * dragging knot
 */
void mouseMotion(int x, int y) {
  if (mouseRotate) {
#define RADIANS_PER_PIXEL (M_PI/(2*90.0))
    int dx = x - mousex, dy = y - mousey;
    eyeTheta -= dx*RADIANS_PER_PIXEL;
    eyePhi -= dy*RADIANS_PER_PIXEL;
    if (eyePhi >= M_PI)
      eyePhi = M_PI - EPSILON;
    else if (eyePhi <= 0.0)
      eyePhi = EPSILON;
    sphericalToCartesian(eyeRadius, eyeTheta, eyePhi,
                         &eye[0], &eye[1], &eye[2]);
    eye[0] += lookat[0];
    eye[1] += lookat[1];
    eye[2] += lookat[2];
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              lookat[0], lookat[1], lookat[2],
              up[0], up[1], up[2]);
    setLookAt(&viewpick, eye, lookat, up);
    mousex = x;
    mousey = y;
    glutPostRedisplay();
  } else if (drag == DRAG_ALPHA) {
    int yy = viewpick.winHeight - y;
    double alphaVal = ((double) yy - alphaBarMin)/ALPHA_BAR;
    static char msg[80];
    sprintf(msg, "alpha = %0.4f", alphaVal);
    messagePos[0] = 0.0;
    messagePos[1] = 0.0;
    messagePos[2] = 0.0;
    message = msg;
    setConstantFunc2(alpha, alphaVal);
    if (useHermiteRis)
      computeHermiteRisMaps(ris);
    else
      computeBilinearRisMaps(ris);
    dirtyGrid = GL_TRUE;
    dirty_qgrid = GL_TRUE;
    glutPostRedisplay();
  } else if (drag == DRAG_KNOT) {
    double dir[3], knot[3], t;
    screenToRay(&viewpick, x, y, dir);
    knot[0] = ris->x[knoti];
    knot[1] = ris->y[knotj];
    knot[2] = ris->z[knotj][knoti];
    closestPointOnRay(knot, up, eye, dir, &t);
    ris->z[knotj][knoti] += t;
    matchEdgeVal(knoti, knotj);
    if (useHermiteRis)
      computeHermiteRisMaps(ris);
    else
      computeBilinearRisMaps(ris);
    dirtyGrid = GL_TRUE;
    dirty_qgrid = GL_TRUE;
    glutPostRedisplay();
  } else if (drag == DRAG_DX) {
    double dir[3], knot[3], t, dz;
    screenToRay(&viewpick, x, y, dir);
    knot[0] = ris->x[knoti] + DELTA;
    knot[1] = ris->y[knotj];
    knot[2] = ris->z[knotj][knoti] + (dz = DELTA*ris->zx[knotj][knoti]);
    closestPointOnRay(knot, up, eye, dir, &t);
    ris->zx[knotj][knoti] = (dz + t)/DELTA;
    matchEdgeVal(knoti, knotj);
    computeHermiteRisMaps(ris);
    dirtyGrid = GL_TRUE;
    dirty_qgrid = GL_TRUE;
    glutPostRedisplay();
  } else if (drag == DRAG_DY) {
    double dir[3], knot[3], t, dz;
    screenToRay(&viewpick, x, y, dir);
    knot[0] = ris->x[knoti];
    knot[1] = ris->y[knotj] + DELTA;
    knot[2] = ris->z[knotj][knoti] + (dz = DELTA*ris->zy[knotj][knoti]);
    closestPointOnRay(knot, up, eye, dir, &t);
    ris->zy[knotj][knoti] = (dz + t)/DELTA;
    matchEdgeVal(knoti, knotj);
    computeHermiteRisMaps(ris);
    dirtyGrid = GL_TRUE;
    dirty_qgrid = GL_TRUE;
    glutPostRedisplay();
  }
}

void printLists(FILE *f, int M, int N, double **a) {
  int i, j;
  fprintf(f, "\t{");
  for (j = 0; j <= N; j++) {
    fprintf(f, "{");
    for (i = 0; i <= M; i++)
      fprintf(f, " %g%c", a[j][i], (i < M) ? ',' : '}');
    if (j < N)
      fprintf(f, ",\n\t ");
  }
  fprintf(f, "}");
}

#define DUMP_FILE "/tmp/xxx.in"

void dumpRayTraceFile(void) {
  FILE *f;

  if ((f = fopen(DUMP_FILE, "w")) == NULL) {
    perror(DUMP_FILE);
    exit(-1);
  }

  fprintf(f, "# generated by edris.c()\n\n");
  fprintf(f, "LOOKAT = (%g, %g, %g), (%g, %g, %g), (%g, %g, %g)\n",
          eye[0], eye[1], eye[2],
          lookat[0], lookat[1], lookat[2],
          up[0], up[1], up[2]);
  fprintf(f, "PROJECTION = %g, %g\n", fovy, aspect);
  fprintf(f, "IMAGE = \"xxx.ppm\", %d, %d, 2\n", (int) (512*aspect + 0.5), 512);
  fprintf(f, "RECURSEDEPTH = 5\n");
  fprintf(f, "AMBIENT = (%g, %g, %g)\n",
          lmodel_ambient[0], lmodel_ambient[1], lmodel_ambient[2]);
  fprintf(f, "BACKGROUND = (0, 0.1, 0.6)\n");
  fprintf(f, "LIGHT = (%g, %g, %g), (1, 1, 1), 1, 0, 0\n",
          light0_position[0], light0_position[1], light0_position[2]);
  fprintf(f, "KA = %g\n",
          (material_ambient[0] + material_ambient[1] + material_ambient[2])/3);
  fprintf(f, "KD = %g\n",
          (material_diffuse[0] + material_diffuse[1] + material_diffuse[2])/3);
  fprintf(f, "KS = %g\n",
          (material_specular[0] + material_specular[1] + material_specular[2])/3);
  fprintf(f, "KT = 0\n");
  fprintf(f, "NI = 1.52\n");
  fprintf(f, "PHONG = %g\n", material_phong);
  fprintf(f, "COLOR = (%g, %g, %g)\n",
          material_diffuse[0], material_diffuse[1], material_diffuse[2]);
  if (useHermiteRis) {
    fprintf(f, "HERMITERIS = %d, %d, %d, %d, %g, %g, %g, %g,\n",
            ris->M, ris->N, ris->Sx, ris->Sy, ris->x[0], ris->y[0],
            ris->x[ris->M] - ris->x[0], ris->y[ris->N] - ris->y[0]);
    fprintf(f, "\t1000, 1000,\n");
    fprintf(f, "\t%g,\n", (alpha->eval)(alpha, 0, 0)); /* xxx constant alpha */
    printLists(f, ris->M, ris->N, ris->z);
    fprintf(f, ",\n");
    printLists(f, ris->M, ris->N, ris->zx);
    fprintf(f, ",\n");
    printLists(f, ris->M, ris->N, ris->zy);
    fprintf(f, ",\n");
    printLists(f, ris->M, ris->N, ris->zxy);
    fprintf(f, "\n");
  } else {
    fprintf(f, "BILINEARRIS = %d, %d, %d, %d, %g, %g, %g, %g,\n",
            ris->M, ris->N, ris->Sx, ris->Sy, ris->x[0], ris->y[0],
            ris->x[ris->M] - ris->x[0], ris->y[ris->N] - ris->y[0]);
    fprintf(f, "\t1000, 1000,\n");
    fprintf(f, "\t%g,\n", (alpha->eval)(alpha, 0, 0)); /* xxx constant alpha */    	printLists(f, ris->M, ris->N, ris->z);
    fprintf(f, "\n");
  }

  fclose(f);
}

enum {
  MENU_DRAW_CONTROLS = 1,    /* enable/disable display of grid controls */
  MENU_BILINEAR_TO_HERMITE,  /* convert bilinear RIS to hermite RIS */
  MENU_HERMITE_TO_BILINEAR,  /* convert hermite RIS to bilinear RIS */
  MENU_DRAWQ,                /* draw q surfaces only */
  MENU_RANDOMIZE,            /* give random values to surface controls */
  MENU_DUMP_RAYTRACE,
  MENU_QUIT = 666            /* terminate application */
};

/*
 * menu call back
 */
void menu(int value) {
  switch(value) {
    case MENU_DRAW_CONTROLS:
      drawControls = !drawControls;
      glutPostRedisplay();
      break;
    case MENU_BILINEAR_TO_HERMITE:
      if (!useHermiteRis) {
        bilinearToHermiteRis();
        dirtyBigGrid = dirtyGrid = GL_TRUE;
        dirty_qgrid = GL_TRUE;
        glutPostRedisplay();
      }
      break;
    case MENU_HERMITE_TO_BILINEAR:
      if (useHermiteRis) {
        hermiteToBilinearRis();
        dirtyBigGrid = dirtyGrid = GL_TRUE;
        dirty_qgrid = GL_TRUE;
        glutPostRedisplay();
      }
      break;
    case MENU_DRAWQ:
      drawq = !drawq;
      glutPostRedisplay();
      break;
    case MENU_RANDOMIZE:
      randomizeRis();
      dirtyBigGrid = dirtyGrid = GL_TRUE;
      dirty_qgrid = GL_TRUE;
      glutPostRedisplay();
      break;
    case MENU_DUMP_RAYTRACE:
      dumpRayTraceFile();
      break;
    case MENU_QUIT:
      exit(0);
  }
}

int chaosItersList[] = {
  0, 0, 1e4, 5e4, 1e5, 5e5, 1e6, 2e6, 5e6, 1e7   /* first 2 entries not used */
};

#define MENU_CHAOS_OFF 1

void chaosMenu(int val) {
  if (val == MENU_CHAOS_OFF) {
    useChaos = GL_FALSE;
  } else {
    chaosIters = chaosItersList[val];
    useChaos = GL_TRUE;
  }
  glutPostRedisplay();
}

/*
 * keyboard callback
 * ESC : kill application
 * z : zoom in
 * Z : zoom out
 */
void keyboard(unsigned char key, int x, int y) {
#define ESC 27
  int eyeChanged = GL_FALSE;
  
  switch(key) {
    case ESC: 
      exit(0);
    case 'z':
      eyeRadius *= 0.97;
      eyeChanged = GL_TRUE;
      break;
    case 'Z':
      eyeRadius *= 1.03;
      eyeChanged = GL_TRUE;
      break;
  }

  if (eyeChanged) {
    sphericalToCartesian(eyeRadius, eyeTheta, eyePhi,
                         &eye[0], &eye[1], &eye[2]);
    eye[0] += lookat[0];
    eye[1] += lookat[1];
    eye[2] += lookat[2];
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye[0], eye[1], eye[2],
              lookat[0], lookat[1], lookat[2],
              up[0], up[1], up[2]);
    setLookAt(&viewpick, eye, lookat, up);
    glutPostRedisplay();
  }
}

int main(int argc, char *argv[]) {
  int i;
  int menuId, chaosMenuId;
  int M = GRID_SIZE, S = DOMAIN_SIZE;

#ifdef USE_COMMAND_LINE_ARGS
  
  /*
   * Check to see if user has specified any parameters.
   */
  if (argc > 1) {
    M = atoi(argv[1]);
    if (M < 4) {
      fprintf(stderr, "Grid size sucks eggs!\n");
      exit(-1);
    }
    if (argc > 2) {
      S = atoi(argv[2]);
      if (S < 2 || M % S != 0) {
        fprintf(stderr, "Domain size not compatible with grid size!\n");
        exit(-1);
      }
    }
  }

#endif

  /*
   * Create and initialize RIS.
   */
  alpha = createConstantFunc2(ALPHA);
  ris = createBilinearRis(M, M, S, S, ORG, ORG, SIZE, SIZE, alpha);
  initRis();

  randomizeRis(); /* XXXX */

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(700,700);
  glutInitWindowPosition(10,10);
  glutCreateWindow("Recurrent Interpolation Surface");
  glutReshapeFunc(reshape);
  glutDisplayFunc(display);
  glutMouseFunc(mouse);
  glutMotionFunc(mouseMotion);
  menuId = glutCreateMenu(menu);
  glutAddMenuEntry("Draw Controls", MENU_DRAW_CONTROLS);
  chaosMenuId = glutCreateMenu(chaosMenu);
  glutAddMenuEntry("Off", MENU_CHAOS_OFF);
  for (i = 2; i < sizeof(chaosItersList)/sizeof(int); i++) {
    char ibuf[30];
    sprintf(ibuf, "%d iters", chaosItersList[i]);
    glutAddMenuEntry(ibuf, i);
  }
  glutSetMenu(menuId);
  glutAddSubMenu("Chaos Game", chaosMenuId);
  glutAddMenuEntry("Bilinear to Hermite RIS", MENU_BILINEAR_TO_HERMITE);
  glutAddMenuEntry("Hermite to Bilinear RIS", MENU_HERMITE_TO_BILINEAR);
  glutAddMenuEntry("Draw q's only", MENU_DRAWQ);
  glutAddMenuEntry("Random Controls", MENU_RANDOMIZE);
  glutAddMenuEntry("Dump Raytrace Input File", MENU_DUMP_RAYTRACE);
  glutAddMenuEntry("Quit", MENU_QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
  glutKeyboardFunc(keyboard);
  glutMainLoop();        

  return 0;
}

