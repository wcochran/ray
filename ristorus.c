/*
 * Recurrent Torus.
 *
 * Torus is defined by its "big" and "little" radii R and r.
 * Given a surface point parameterized by (u,v), we convert
 * to cartesian coordinates as follows:
 *   x(u,v) = (R + r*cos(u))*cos(v),   0 <= u,v <= 2*pi
 *   y(u,v) = (R + r*cos(u))*sin(v)
 *   z(u,v) = r*sin(u)
 *
 * To create our fractal torus, we allow r to vary as follows
 *    r = r' + d(u,v)
 * where r' is the nominal "little" radius and 
 * d(u,v) is a Recurrent Interpolation Surface (RIS).
 */

#ifdef WIN32
#include <windows.h>
#endif
#if defined(__APPLE__) || defined(MACOSX)
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ris.h"
#include "func2.h"
#include "noise.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ALPHA -0.3

double bigRadius = 10.0;
double littleRadius = 4.0;

Func2 *alpha;
Ris *ris;
Func2* risFunc;

int useHermiteRis = GL_FALSE;

typedef struct {      /* torus surface mesh */
  int w, h;
  double (**verts)[3];
} Mesh;

int useBigMesh = GL_TRUE;
int dirtyMesh = GL_TRUE;
int dirtyBigMesh = GL_TRUE;
Mesh *mesh;
Mesh *bigMesh;

void *myalloc(size_t sz) {
  void *p = malloc(sz);
  if (p == NULL) {
    perror("malloc()");
    exit(-1);
  }
  return p;
}

Mesh *createMesh(int w, int h) {
  Mesh *mesh = (Mesh *) myalloc(sizeof(Mesh));
  int j;
  mesh->w = w;
  mesh->h = h;
  mesh->verts = (double (**)[3]) myalloc(h*sizeof(void *));
  for (j = 0; j < h; j++)
    mesh->verts[j] = (double (*)[3]) myalloc(3*w*sizeof(double));
  return mesh;
}

void torusToMesh(double R, double r, Func2* delta, Mesh *mesh) {
  int i,j;
  double u,v;
  double du = 2*M_PI/(mesh->w);
  double dv = 2*M_PI/(mesh->h);
  for (j = 0, v = 0.0; j < mesh->h; j++, v += dv) {
    double sinv = sin(v);
    double cosv = cos(v);
    for (i = 0, u = 0.0; i < mesh->w; i++, u += du) {
      double *p = mesh->verts[j][i];
      double tweek = (*delta->eval)(delta, u,v);
      double A = (R + (r + tweek)*cos(u));
      p[0] = A*cosv;
      p[1] = A*sinv;
      p[2] = (r + tweek)*sin(u);
    }
  }
}

void init(void) {
  int i,j;
  double u,v, du,dv;
  alpha = createConstantFunc2(ALPHA);
  ris = createBilinearRis(24, 24, 2, 2, 0.0, 0.0, 2*M_PI, 2*M_PI, alpha);
  /*
  du = 10*2*M_PI/ris->M;
  dv = 16*2*M_PI/ris->N;
  */
  du = 2*2*M_PI/ris->M;
  dv = 3*2*M_PI/ris->N;
  for (j = 0, v = 0.0; j <= ris->N; j++, v += dv) {
    /* double fv = sin(v); */
    for (i = 0, u = 0.0; i <= ris->M; i++, u += du) {
      /* double fu = sin(u); */
      /* ris->z[j][i] = fv + fu; */
      /* ris->z[j][i] = 2.5*noise3(27*sin(u), 31*sin(v), 0.0); */
      ris->z[j][i] = 3.5*noise3(1*sin(u), 3*sin(v), 0.0);
    }
  }
  computeBilinearRisMaps(ris);
  risFunc = risToFunc2(ris);
  mesh = createMesh(40, 40);
  torusToMesh(bigRadius, littleRadius, risFunc, mesh);
  dirtyMesh = GL_FALSE;
  bigMesh = createMesh(450, 450);
  torusToMesh(bigRadius, littleRadius, risFunc, bigMesh);
  dirtyBigMesh = GL_FALSE;
  useBigMesh = GL_TRUE;
}

void toHermite(void) {
  int M = ris->M, N = ris->N;
  int i,j;
  
  Ris *hris = createHermiteRis(M, N, ris->Sx, ris->Sy,
                               0.0, 0.0, 2*M_PI, 2*M_PI, alpha);
  for (j = 0; j <= N; j++) {
    int prevj = (j == 0) ? N : j-1;
    int nextj = (j+1) % N;
    double dy = ris->y[nextj] - ris->y[prevj];
    for (i = 0; i <= M; i++) {
      int previ = (i == 0) ? M : i-1;
      int nexti = (i+1) % M;
      double dx = ris->x[nexti] - ris->x[previ];
      hris->z[j][i] = ris->z[j][i];
      hris->zx[j][i] = (ris->z[j][nexti] - ris->z[j][previ])/dx;
      hris->zy[j][i] = (ris->z[nextj][i] - ris->z[prevj][i])/dy;
      hris->zxy[j][i] = 0.0;
    }
  }

  computeHermiteRisMaps(hris);
  risFunc = risToFunc2(hris);  /* XXX memory leak */
  destroyBilinearRis(ris);
  useHermiteRis = GL_TRUE;
  ris = hris;  
}

#define DIFF(u,v,d) (d[0] = u[0] - v[0],\
                     d[1] = u[1] - v[1],\
                     d[2] = u[2] - v[2])
#define DOT(u,v) (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])
#define CROSS(u,v,uxv) (uxv[0] = u[1]*v[2] - u[2]*v[1],\
                        uxv[1] = u[2]*v[0] - u[0]*v[2],\
                        uxv[2] = u[0]*v[1] - u[1]*v[0])
#define SCALE(k,u) (u[0] *= k, u[1] *= k, u[2] *= k)

void drawMeshAsQuads(Mesh *mesh) {
  int i,j;
  for (j = 0; j < mesh->h; j++) {
    int nextj = (j+1) % mesh->h;
    glBegin(GL_QUAD_STRIP);
    for (i = 0; i <= mesh->w; i++) {
      int thisi = i % mesh->w;
      int nexti = (i+1) % mesh->w;
      double *p00 = mesh->verts[j][thisi];
      double *p01 = mesh->verts[j][nexti];
      double *p10 = mesh->verts[nextj][thisi];
      double u[3], v[3], n[3], s;
      DIFF(p10, p00, u);
      DIFF(p01, p00, v);
      CROSS(u, v, n);
      s = 1.0/sqrt(DOT(n,n));
      SCALE(s, n);
      glNormal3dv(n);
      glVertex3dv(p00);
      glVertex3dv(p10);
    }
    glEnd();
  }
}

void drawMeshAsTriangles(Mesh *mesh) {
  int i,j;
  for (j = 0; j < mesh->h; j++) {
    int nextj = (j+1) % mesh->h;
    glBegin(GL_TRIANGLE_STRIP);
    for (i = 0; i <= mesh->w; i++) {
      int thisi = i % mesh->w;
      int nexti = (i+1) % mesh->w;
      double *p00 = mesh->verts[j][thisi];
      double *p01 = mesh->verts[j][nexti];
      double *p10 = mesh->verts[nextj][thisi];
      double *p11 = mesh->verts[nextj][nexti];
      double u[3], v[3], n[3], s;
      DIFF(p10, p00, u);
      DIFF(p01, p00, v);
      CROSS(u, v, n);
      s = 1.0/sqrt(DOT(n,n));
      SCALE(s, n);
      glNormal3dv(n);
      glVertex3dv(p00);
      DIFF(p01, p11, u);
      DIFF(p10, p11, v);
      CROSS(u, v, n);
      s = 1.0/sqrt(DOT(n,n));
      SCALE(s, n);
      glNormal3dv(n);
      glVertex3dv(p10);
    }
    glEnd();
  }
}
    
void display(void) {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  if (useBigMesh) {
    if (dirtyBigMesh) {
      torusToMesh(bigRadius, littleRadius, risFunc, bigMesh);
      dirtyBigMesh = GL_FALSE;
    }
    drawMeshAsTriangles(bigMesh);
  } else {
    if (dirtyMesh) {
      torusToMesh(bigRadius, littleRadius, risFunc, mesh);
      dirtyMesh = GL_FALSE;
    }
    drawMeshAsQuads(mesh);
  }
  glutSwapBuffers();
}

double eyeRadius, eyeTheta, eyePhi;
GLdouble eye[3], lookat[3], up[3];

void sphericalToCartesian(double r, double theta, double phi,
				 double *x, double *y, double *z) {
  double sin_phi = sin(phi);
  *x = r*cos(theta)*sin_phi;
  *y = r*sin(theta)*sin_phi;
  *z = r*cos(phi);
}

GLfloat lmodel_ambient[4]    = {0.2, 0.2, 0.2, 2.0};

GLfloat light0_position[4]   = {30, 10, 80, 1.0};
GLfloat light0_ambient[4]    = {0.2, 0.2, 0.2, 1.0};
GLfloat light0_diffuse[4]    = {1.0, 1.0, 1.0, 1.0};
GLfloat light0_specular[4]   = {1.0, 1.0, 1.0, 1.0};

GLfloat light1_position[4]   = {1.0, 1.0, 1.0, 0.0}; /* infinite */
GLfloat light1_ambient[4]    = {0.2, 0.2, 0.2, 1.0};
GLfloat light1_diffuse[4]    = {0.4, 0.4, 0.4, 1.0};
GLfloat light1_specular[4]   = {0.6, 0.6, 0.6, 1.0};

GLfloat material_ambient[4]  = {0.24725, 0.1995, 0.0745, 1.0};
GLfloat material_diffuse[4]  = {0.75164, 0.60648, 0.22648, 1.0};
GLfloat material_specular[4] = {0.628281, 0.555802, 0.366025, 1.0};
GLfloat material_shininess   = 51.2;

void graphicsInit(void) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(40.0, 1.0, 0.01*bigRadius, 10*bigRadius);
  
  eyeRadius = 4.5*bigRadius;
  eyeTheta = 0.0;
  eyePhi = M_PI/3;
  lookat[0] = lookat[1] = lookat[2] = 0.0;
  up[0] = up[1] = 0.0; up[2] = 1.0;

  sphericalToCartesian(eyeRadius, eyeTheta, eyePhi, &eye[0], &eye[1], &eye[2]);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(eye[0], eye[1], eye[2], 
	    lookat[0], lookat[1], lookat[2],
	    up[0], up[1], up[2]);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);  
  glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light1_specular);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);  

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_LIGHT1);

  glMaterialfv(GL_FRONT, GL_AMBIENT, material_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, material_diffuse);
  glMaterialfv(GL_FRONT, GL_SPECULAR, material_specular);
  glMaterialfv(GL_FRONT, GL_SHININESS, &material_shininess);

  glShadeModel(/* GL_SMOOTH */ GL_FLAT);

  glFrontFace(GL_CCW);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  
  glEnable(GL_DEPTH_TEST);
}

enum {
  MENU_BILINEAR_HERMITE_TOGGLE = 1,
  MENU_QUIT = 666
};

void menu(int value) {
  switch(value) {
    case MENU_BILINEAR_HERMITE_TOGGLE:
      if (!useHermiteRis) {
        toHermite();
        useHermiteRis = GL_TRUE;
        dirtyMesh = dirtyBigMesh = GL_TRUE;
        glutPostRedisplay();
      }
      break;
    case MENU_QUIT:
      exit(0);
  }
}      

#define ESC 27

void keyboard(unsigned char key, int x, int y) {
  if (key == ESC) exit(0);
}

int main(int argc, char *argv[]) {
  init();
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
  glutInitWindowSize(700,700);
  glutInitWindowPosition(10,10);
  glutCreateWindow("Donut");
  glutDisplayFunc(display);
  /*
  glutMouseFunc(mouse);  
  glutMotionFunc(mouseMotion);
  */
  glutCreateMenu(menu);
  glutAddMenuEntry("Bilinear/Hermite", MENU_BILINEAR_HERMITE_TOGGLE);
  glutAddMenuEntry("Quit", MENU_QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);  
  glutKeyboardFunc(keyboard);
  graphicsInit();
  glutMainLoop();
  return 0;
}
