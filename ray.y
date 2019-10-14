%{
#include <stdlib.h>
#include <math.h>

extern int yylex(void);
void yyerror(char *msg);
void *myalloc(size_t sz);

typedef struct List {  /* base type for all lists */
  struct List *next;
} List;

int length(List *list);
List *insert(List *list, List *node);  /* insert onto tail */

typedef struct VecList {
  struct VecList *next;
  double v[3];
} VecList;

typedef struct VecListList {
  struct VecListList *next;
  VecList *list;
} VecListList;

typedef struct NumList {
  struct NumList *next;
  double num;
} NumList;

typedef struct Matrix {
  struct Matrix *next;
  NumList *row;
} Matrix;

VecList *newVecList(double v[3]);
VecListList *newVecListList(VecList *);
NumList *newNumList(double num);
Matrix *newMatrixRow(NumList *row);

void go(void);
void setLookAt(double eye[3], double lookat[3], double up[3]);
void setProjection(double fovy, double aspect);
void setImage(char *fname, int w, int h, int samples);
void setRecursionDepth(int depth);
void setAmbient(double ambient[3]);
void setBackground(double background[3]);
void setLight(double pos[3], double color[3], 
              double c0, double c1, double c2);
void setVariable(char *name, double val);
double lookupVariable(char *name);
void setMaterialAmbient(double);
void setMaterialDiffuse(double);
void setMaterialSpecular(double);
void setMaterialTransmission(double);
void setMaterialRefraction(double);
void setMaterialPhong(double);
void setMaterialColor(double color[3]);
void setMarbleTexture(VecList *colors, double veinDir[3], int octaves);
void setCheckerTexture(void);
void makeSphere(double center[3], double rad);
void makePlane(double normal[3], double point[3]);
void makeImagePlane(double normal[3], double point[3],
		    double org[3], double yaxis[3],
		    double width, double height,
		    char *imagename);
void makeBezier3(VecListList *mesh);
void makeTeapot(double org[3], double up[3], double spout[3], double scale);
void makeElevationMap(double orgxy[], double widthheight[], 
                      double zminmax[], char *fname);
void makeHermiteFunc(double x[], double y[], 
                     int W, int H,
                     double *f, double *fx, double *fy, double *fxy);
void makeBilinearFis(double x[], double y[],
		     int W, int H,
		     double **alpha, double **z);
void makeHermiteFis(double x[], double y[],
		    int W, int H,
		    double **alpha, 
		    double **z, double **zx, double **zy, double **zxy);
void makeBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
	             double width, double height, int W, int H,
		     double alpha, Matrix *z);
void makeHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                    double width, double height, int W, int H, double alpha, 
                    Matrix *z, Matrix *zx, Matrix *zy, Matrix *zxy);
void makeSuperellipsoid(double n, double m, 
			double size[3], double z[3], double x[3], double center[3]);
void makePointCloud(const char *fname, double maxLevel, double sigma, double halfThickness);
#ifdef YYDEBUG
int yydebug=1;
#endif
%}
%union {
  char *s;
  int i;
  double d;
  double v[3];
  double *a;
  double **m;
  struct VecList *veclist;
  struct VecListList *vecllist;
  struct NumList *numlist;
  struct Matrix *matrix;
}
%token <s> STR
%token <s> IDENT
%token <i> INT
%token <d> REAL
%token LOOKAT PROJ IMAGE RECDEPTH AMBIENT BACKGROUND LIGHTSRC
%token SIN COS TAN SQRT
%token KA KD KS KT NI PHONG COLOR MARBLE CHECKER
%token SPHERE PLANE IMAGEPLANE BEZIER3 TEAPOT ELEVMAP HERMITEFUNC 
%token BILINEARFIS HERMITEFIS
%token BILINEARRIS HERMITERIS
%token SUPERELLIPSOID POINTCLOUD
%left '+' '-'
%left '*' '/'
%right '^'
%right UMINUS
%type <d> expr num
%type <v> vec pair
%type <a> a4 mat2x2
%type <m> mat4x4
%type <veclist> vec_list
%type <vecllist> vec_llist
%type <numlist> num_list
%type <matrix> matrix
%%
args       : params           {go();}
           ;

params     : param
           | param params
           ;

param      : lookat
           | projection
           | image
           | rec_depth
           | ambient
           | background
           | light
           | assign
           | material
           | object
           ;

lookat     : LOOKAT '=' vec ',' vec ',' vec    {setLookAt($3, $5, $7);}
           ;

projection : PROJ '=' expr ',' expr            {setProjection($3, $5);}
           ;

image      : IMAGE '=' STR ',' expr ',' expr ',' expr 
                             {setImage($3, (int) $5, (int) $7, (int) $9); 
			      free($3);}
           | IMAGE '=' STR ',' expr ',' expr
                             {setImage($3, (int) $5, (int) $7, 1); 
			      free($3);}
           | IMAGE '=' STR ',' expr
                             {setImage($3, (int) $5, -1, 1); 
			      free($3);}
           ;

rec_depth  : RECDEPTH '=' expr          {setRecursionDepth((int) $3);}
           ;

ambient    : AMBIENT '=' vec            {setAmbient($3);}
           ;

background : BACKGROUND '=' vec         {setBackground($3);}
           ;

light      : LIGHTSRC '=' vec ',' vec ',' expr ',' expr ',' expr
                                        {setLight($3, $5, $7, $9, $11);}
           | LIGHTSRC '=' vec ',' vec   {setLight($3, $5, 1, 0, 0);}
           ;

vec        : '(' expr ',' expr ',' expr ')'
                                        {$$[0] = $2; $$[1] = $4; $$[2] = $6;}
           ;

vec_list   : vec                        {$$ = newVecList($1);}
	   | vec_list ',' vec           {$$ = (VecList *) 
                                            insert((List *) $1, 
	                                           (List *) newVecList($3));}
	   ;

vec_llist  : '{' vec_list '}'                 {$$ = newVecListList($2);}
           |  vec_llist ',' '{' vec_list '}'  {$$ = (VecListList *)
                                                  insert((List *) $1,
                                                         (List *)
                                                          newVecListList($4));}
           ;

pair       : '(' expr ',' expr ')'      {$$[0] = $2; $$[1] = $4;}
           ;

mat2x2     : '(' pair ',' pair ')'      {$$ = myalloc(4*sizeof(double));
                                         $$[0] = $2[0]; $$[1] = $2[1];
                                         $$[2] = $4[0]; $$[3] = $4[1];}
           ;

a4         : '(' expr ',' expr ',' expr ',' expr ')'
                                        {$$ = myalloc(4*sizeof(double));
                                         $$[0] = $2; $$[1] = $4;
                                         $$[2] = $6; $$[3] = $6;}
           ;

mat4x4     : '(' a4 ',' a4 ',' a4 ',' a4 ')'
                                        {$$ = myalloc(4*sizeof(double *));
                                         $$[0] = $2; $$[1] = $4;
                                         $$[2] = $6; $$[3] = $6;}
           ;

assign     : IDENT '=' expr             {setVariable($1, $3); free($1);}
           ;

expr       : IDENT                      {$$ = lookupVariable($1); free($1);}
           | num
           | '(' expr ')'               {$$ = $2;}
           | '-' expr     %prec UMINUS  {$$ = -$2;}
           | '+' expr     %prec UMINUS  {$$ = $2;}
           | expr '+' expr              {$$ = $1 + $3;}
           | expr '-' expr              {$$ = $1 - $3;}
           | expr '*' expr              {$$ = $1 * $3;}
           | expr '/' expr              {$$ = $1 / $3;}
           | expr '^' expr              {$$ = pow($1,$3);}
           | SIN '(' expr ')'           {$$ = sin($3);}
           | COS '(' expr ')'           {$$ = cos($3);}
           | TAN '(' expr ')'           {$$ = tan($3);}
           | SQRT '(' expr ')'          {$$ = sqrt($3);}
           ;

num        : INT                        {$$ = (double) $1;}
           | REAL
           ;

num_list   : expr                       {$$ = newNumList($1);}
           | num_list ',' expr          {$$ = (NumList *)
                                          insert((List *) $1,
                                                 (List *) newNumList($3));}
           ;

matrix     : '{' num_list '}'             {$$ = newMatrixRow($2);}
           | matrix ',' '{' num_list '}'  {$$ = (Matrix *)
                                             insert((List *) $1,
                                                  (List *) newMatrixRow($4));}
           ;

material   : KA '=' expr                {setMaterialAmbient($3);}
           | KD '=' expr                {setMaterialDiffuse($3);}
           | KS '=' expr                {setMaterialSpecular($3);}
           | KT '=' expr                {setMaterialTransmission($3);}
           | NI '=' expr                {setMaterialRefraction($3);}
           | PHONG '=' expr             {setMaterialPhong($3);}
           | COLOR '=' vec              {setMaterialColor($3);}
           | MARBLE '=' '(' vec_list ')' ',' vec ',' expr
                                        {setMarbleTexture($4, $7, $9);}
           | CHECKER                    {setCheckerTexture();}
           ;

object     : SPHERE '=' vec ',' expr    {makeSphere($3, $5);}
           | PLANE '=' vec ',' vec      {makePlane($3, $5);}
           | IMAGEPLANE '=' vec ',' vec ','
                vec ',' vec ',' expr ',' expr ',' STR
                                        {makeImagePlane($3, $5, 
                                                $7, $9, $11, $13, $15);}
           | BEZIER3 '=' vec_llist      {makeBezier3($3);}
           | TEAPOT '=' vec ',' vec ',' vec ',' expr
                                        {makeTeapot($3, $5, $7, $9);}
           | ELEVMAP '=' pair ',' pair ',' pair ',' STR
                                        {makeElevationMap($3, $5, 
                                                          $7, $9);} 
           | HERMITEFUNC '=' pair ',' pair ',' expr ',' expr ','
               mat2x2 ',' mat2x2 ',' mat2x2 ',' mat2x2
                                        {makeHermiteFunc($3, $5, 
                                                         (int) $7, (int) $9, 
                                                         $11, $13, $15, $17);}
           | BILINEARFIS '=' pair ',' pair ',' expr ',' expr ',' 
               mat4x4 ',' mat4x4
                                        {makeBilinearFis($3, $5, $7, $9,
							 $11, $13);}
           | HERMITEFIS '=' pair ',' pair ',' expr ',' expr ',' 
               mat4x4 ',' mat4x4 ',' mat4x4 ',' mat4x4 ',' mat4x4
                                        {makeHermiteFis($3, $5, $7, $9,
                                                        $11, 
                                                        $13, $15, $17, $19);}
           | BILINEARRIS '=' expr ',' expr ',' expr ',' expr ','
               expr ',' expr ',' expr ',' expr ',' expr ',' expr ','
               expr ',' '{' matrix '}'
	                                {makeBilinearRis($3, $5, $7, $9,
                                                         $11, $13, $15, $17,
                                                         $19, $21,
                                                         $23, 
							 $26);}
           | HERMITERIS '=' expr ',' expr ',' expr ',' expr ','
               expr ',' expr ',' expr ',' expr ',' expr ',' expr ','
               expr ',' '{' matrix '}' ',' 
               '{' matrix '}' ',' '{' matrix '}' ',' '{' matrix '}'
	                                {makeHermiteRis($3, $5, $7, $9,
                                                        $11, $13, $15, $17,
                                                        $19, $21,
							$23,
							$26, $30, $34, $38);}
           /* n, m, size[3], z[3], x[3], center[3] */
           | SUPERELLIPSOID '=' expr ',' expr ',' vec ',' vec ',' vec ',' vec
                                        {makeSuperellipsoid($3,$5,$7,$9,$11,$13);}

           | POINTCLOUD '=' STR ',' expr ',' expr ',' expr
                   {makePointCloud($3, $5, $7, $9);}
           ;
%%

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/errno.h>
#include <sys/param.h>
#include <assert.h>
#include <unistd.h>
#include "pnmio.h"
#include "raytrace.h"
#include "marble.h"
#include "checker.h"
#include "sphere.h"
#include "plane.h"
#include "bezier3.h"
#include "teapot.h"
#include "func2.h"
#include "rectgrid.h"
#include "hermite.h"
#include "bilinear.h"
#include "rectfis.h"
#include "bilinearfis.h"
#include "hermitefis.h"
#include "ris.h"
#include "superellipsoid.h"
#include "pointcloud.h"

/*
 * Bitmaps to track which parameters were set at least once.
 */
enum {
  SET_LOOKAT        = 0x00001,
  SET_PROJECTION    = 0x00002,
  SET_IMAGE         = 0x00004,
  SET_RECDEPTH      = 0x00008,
  SET_AMBIENT       = 0x00010,
  SET_BACKGROUND    = 0x00020,
  SET_MAT_AMBIENT   = 0x00080,
  SET_MAT_DIFFUSE   = 0x00100,
  SET_MAT_SPECULAR  = 0x00200,
  SET_MAT_TRANS     = 0x00400,
  SET_MAT_REFRACT   = 0x00800,
  SET_MAT_PHONG     = 0x01000,  
  SET_MAT_COLOR     = 0x02000  
};

static int setParams = 0;

/*
 * Arguments to raytrace().
 */
typedef struct {
  double eyePos[3];      /* position of eye */
  double eyeDir[3];      /* view direction */
  double eyeUp[3];       /* eye up vector */
  double projWidth;      /* width of projection plane */
  double projHeight;     /* height of projection plane */
  double projDist;       /* distannce from eye to projection plane */
  int imageWidth;        /* width (in pixels) of image */
  int imageHeight;       /* height of image */
  int samples;           /* sqrt of samples per pixel */
  int maxDepth;          /* maximum recursion depth (1,2,..) */
  double ambient[3];     /* RGB intensities of ambient light */
  double background[3];  /* RGB background light intensities */
  int numLights;         /* number of point light sources */
  LIGHT *lights;         /* array of point light sources */
  int numObjects;        /* number of objects in world */
  OBJECT **objects;      /* array of objects in world */
} RAYTRACE_ARGS;

static RAYTRACE_ARGS args;

#define MAX_LIGHTS 100
static int numLights = 0;
static LIGHT lights[MAX_LIGHTS];

#define MAX_OBJECTS 200
static int numObjects = 0;
static OBJECT *objects[MAX_OBJECTS];

static char *filename;
static char *imagename;
extern int lineno;      /* lineno declared in ray.l */

void yyerror(char *msg) {
  fprintf(stderr, "%s [%d] : %s\n", filename, lineno, msg);
  exit(-1);
}

void fatalError(char *msg) {
  fprintf(stderr, "%s : %s\n", filename, msg);
  exit(-1);
}

void warning(char *msg) {
  fprintf(stderr, "%s : %s\n", filename, msg);
}

static void checkParams(void) {
  if ((setParams & SET_LOOKAT) == 0)
    fatalError("Look-at transformation not specified!");
  if ((setParams & SET_PROJECTION) == 0)
    fatalError("Perspective projection not specified!");
  if ((setParams & SET_IMAGE) == 0)
    fatalError("Output image not specified!");
  if (numLights <= 0)
    fatalError("No lights specified!");
  if (numObjects <= 0)
    fatalError("No objects specified!");

  if ((setParams & SET_RECDEPTH) == 0)
    warning("Ray tracing recursion depth not specified!");
  if ((setParams & SET_AMBIENT) == 0)
    warning("Ambient light in scene not specified!");
  if ((setParams & SET_BACKGROUND) == 0)
    warning("Background color not specified!");
  if ((setParams & SET_MAT_AMBIENT) == 0)
    warning("Ambient material coefficient never specified!");
  if ((setParams & SET_MAT_DIFFUSE) == 0)
    warning("Diffuse material coefficient never specified!");
  if ((setParams & SET_MAT_SPECULAR) == 0)
    warning("Specular material coefficient never specified!");
  if ((setParams & SET_MAT_TRANS) == 0)
    warning("Transmission material coefficient never specified!");
  if ((setParams & SET_MAT_REFRACT) == 0)
    warning("Index of refraction material coefficient never specified!");
  if ((setParams & SET_MAT_PHONG) == 0)
    warning("Phong material exponent never specified!");
}

void go(void) {
  FILE *f;
  pnm_image *image;

#ifdef VERBOSE
  printf("go()\n");
#endif

  checkParams();

  if ((f = fopen(imagename, "wb")) == NULL) {
    perror(imagename);
    exit(-1);
  }

  if (args.imageHeight < 1) {  /* match aspect ratio of proj screen */
    double aspect = args.projWidth/args.projHeight;
    args.imageHeight = (int) (args.imageWidth/aspect + 0.5);
  }

  args.numLights = numLights;
  args.lights = lights;
  args.numObjects = numObjects;
  args.objects = objects;

  image = raytrace(args.eyePos, args.eyeDir, args.eyeUp,
		   args.projWidth, args.projHeight, args.projDist,
		   args.imageWidth, args.imageHeight, args.samples,
		   args.maxDepth,
		   args.ambient, args.background,
		   numLights, lights,
		   numObjects, objects);

  write_pnm_image(image, f);

  fclose(f);
}

void setLookAt(double eye[3], double lookat[3], double up[3]) {
#ifdef VERBOSE
  printf("setLookAt((%f,%f,%f), (%f,%f,%f), (%f,%f,%f))\n",
	 eye[0], eye[1], eye[2],
	 lookat[0], lookat[1], lookat[2],
	 up[0], up[1], up[2]);
#endif

  args.eyePos[0] = eye[0];
  args.eyePos[1] = eye[1];
  args.eyePos[2] = eye[2];
  args.eyeDir[0] = lookat[0] - eye[0];
  args.eyeDir[1] = lookat[1] - eye[1];
  args.eyeDir[2] = lookat[2] - eye[2];
  args.eyeUp[0] = up[0];
  args.eyeUp[1] = up[1];
  args.eyeUp[2] = up[2];

  if (fabs(args.eyeDir[0]) <= EPSILON  &&
      fabs(args.eyeDir[1]) <= EPSILON  &&
      fabs(args.eyeDir[2]) <= EPSILON)
    fatalError("Eye position and lookat coincide!");

  if (fabs(args.eyeDir[1]*up[2] - args.eyeDir[2]*up[1]) <= EPSILON &&
      fabs(args.eyeDir[2]*up[0] - args.eyeDir[0]*up[2]) <= EPSILON &&
      fabs(args.eyeDir[0]*up[1] - args.eyeDir[1]*up[0]) <= EPSILON)
    fatalError("Eye direction and up coincide!");

  if ((setParams & SET_LOOKAT) != 0)
    warning("Lookat specified more than once!");
  
  setParams |= SET_LOOKAT;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define D2R(d) ((d)*M_PI/180.0)

void setProjection(double fovy, double aspect) {

#ifdef VERBOSE
  printf("setProjection(%f,%f)\n", fovy, aspect);
#endif

  if (fovy <= EPSILON)
    fatalError("Bogus field of view!");

  if (aspect <= EPSILON)
    fatalError("Bogus projection aspect ratio!");

  args.projHeight = 2*tan(D2R(fovy)/2);
  args.projWidth = aspect*args.projHeight;
  args.projDist = 1.0;

  if ((setParams & SET_PROJECTION) != 0)
    warning("Projection specified more than once!");
  
  setParams |= SET_PROJECTION;
}

void setImage(char *fname, int w, int h, int samples) {
  struct stat sb;
  int i = 0;

#ifdef VERBOSE
  printf("setImage(%s,%d,%d,%d)\n", fname, w, h, samples);
#endif

  if (setParams & SET_IMAGE)
    fatalError("Output image specified more than once!");

  imagename = strdup(fname);

  while (stat(imagename, &sb) >= 0 && i < 999) {
    static char buf[MAXPATHLEN+1];
    sprintf(buf, "%s-%03d", fname, ++i);
    imagename = buf;
  }

  if (i >= 999)
    fatalError("Too many output images with the same name!");

  if (errno != ENOENT) {
    perror(imagename);
    exit(-1);
  }

  if (w < 1)
    fatalError("Bogus image dimensions!");
  /* h < 1 is ok
   * We'll adjust according to aspect ratio of
   * projetion screen later on.
   */

  if (samples < 0 || samples > 4)
    fatalError("Bogus pixel sampling!");

  args.imageWidth = w;
  args.imageHeight = h;
  args.samples = samples;

  setParams |= SET_IMAGE;
}

void setRecursionDepth(int depth) {

#ifdef VERBOSE
  printf("setRecursionDepth(%d)\n", depth);
#endif

  if (depth < 0 || depth > 12)
    fatalError("Bogus recusions depth!");

  args.maxDepth = depth;

  if (setParams & SET_RECDEPTH)
    warning("Recursion depth set more than once!");

  setParams |= SET_RECDEPTH;
}

void setAmbient(double ambient[3]) {
  args.ambient[0] = ambient[0];
  args.ambient[1] = ambient[1];
  args.ambient[2] = ambient[2];

  if (setParams & SET_AMBIENT)
    warning("Ambient color set more than once!");

  setParams |= SET_AMBIENT;
}

void setBackground(double background[3]) {

#ifdef VERBOSE
  printf("setBackground((%f,%f,%f)\n", 
	 background[0], background[1], background[2]);
#endif

  args.background[0] = background[0];
  args.background[1] = background[1];
  args.background[2] = background[2];

  if (setParams & SET_BACKGROUND)
    warning("Background color set more than once!");

  setParams |= SET_BACKGROUND;
}

int length(List *list) {
  int n;
  List *node;
  for (n = 0, node = list; node != NULL; n++, node = node->next)
    ;
  return n;
}

List *insert(List *list, List *node) {
  List *tail;
  if (list == NULL)
    return node;
  for (tail = list; tail->next != NULL; tail = tail->next)
    ;
  tail->next = node;
  node->next = NULL;
  return list;
}

VecList *newVecList(double v[3]) {
  VecList *node = (VecList *) myalloc(sizeof(VecList));
  node->v[0] = v[0];
  node->v[1] = v[1];
  node->v[2] = v[2];
  node->next = NULL;
  return node;
}

VecListList *newVecListList(VecList *vlist) {
  VecListList *node = (VecListList *) myalloc(sizeof(VecListList));
  node->list = vlist;
  node->next = NULL;
  return node;
}

NumList *newNumList(double num) {
  NumList *list = (NumList *) myalloc(sizeof(NumList));
  list->num = num;
  list->next = NULL;
  return list;
}

Matrix *newMatrixRow(NumList *row) {
  Matrix *m = (Matrix *) myalloc(sizeof(Matrix));
  m->row = row;
  m->next = NULL;
  return m;
}

void setLight(double pos[3], double color[3], 
              double c0, double c1, double c2) {

#ifdef VERBOSE
  printf("setLight((%f,%f,%f), (%f,%f,%f), %f, %f, %f)\n", 
	 pos[0], pos[1], pos[2],
	 color[0], color[1], color[2],
	 c0, c1, c2);
#endif

  if (numLights >= MAX_LIGHTS)
    fatalError("Too many friggin' lights!");

  lights[numLights].pos[0] = pos[0];
  lights[numLights].pos[1] = pos[1];
  lights[numLights].pos[2] = pos[2];

  lights[numLights].color[0] = color[0];
  lights[numLights].color[1] = color[1];
  lights[numLights].color[2] = color[2];
  
  lights[numLights].c[0] = c0;
  lights[numLights].c[1] = c1;
  lights[numLights].c[2] = c2;

  numLights++;
}

void *myalloc(size_t sz) {
  void *p = malloc(sz);
  if (p == NULL) {
    perror("malloc()");
    exit(-1);
  }
  return p;
}

typedef struct Var {
  char *name;
  double val;
  struct Var *next;
} Var;

Var *vars = NULL;

void setVariable(char *name, double val) {
  Var *v;

#ifdef VERBOSE
  printf("setVariable(%s, %f)\n", name, val);
#endif

  for (v = vars; v != NULL; v = v->next)
    if (strcmp(name, v->name) == 0) {
      v->val = val;
      return;
    }

  v = (Var *) myalloc(sizeof(Var));

  v->name = strdup(name);

  v->val = val;
  v->next = vars;
  vars = v;
}

double lookupVariable(char *name) {
  Var *v;

#ifdef VERBOSE
  printf("lookupVariable(%s)\n", name);
#endif

  for (v = vars; v != NULL; v = v->next)
    if (strcmp(name, v->name) == 0)
      return v->val;

  yyerror("Use of unknown variable!");

  return 0.0;
}

typedef struct {
  double ka;
  double kd;
  double ks;
  double kt;
  double ni;
  double phong;
  double color[3];
} Material;

Material material = {  /* current material settings */
  0.1, 0.7, 0.2, 0.55, 1.52, 6.0, {0.2, 0.6, 1.0}
};

void setMaterialAmbient(double ka) {
#ifdef VERBOSE
  printf("setMaterialAmbient(%f)\n", ka);
#endif

  material.ka = ka;
  setParams |= SET_MAT_AMBIENT;
}

void setMaterialDiffuse(double kd) {
#ifdef VERBOSE
  printf("setMaterialDiffuse(%f)\n", kd);
#endif

  material.kd = kd;
  setParams |= SET_MAT_DIFFUSE;
}

void setMaterialSpecular(double ks) {
#ifdef VERBOSE
  printf("setMaterialSpecular(%f)\n", ks);
#endif

  material.ks = ks;
  setParams |= SET_MAT_SPECULAR;
}

void setMaterialTransmission(double kt) {
#ifdef VERBOSE
  printf("setMaterialTransmission(%f)\n", kt);
#endif

  material.kt = kt;
  setParams |= SET_MAT_TRANS;
}

void setMaterialRefraction(double ni) {
#ifdef VERBOSE
  printf("setMaterialRefraction(%f)\n", ni);
#endif

  material.ni = ni;
  setParams |= SET_MAT_REFRACT;
}

void setMaterialPhong(double phong) {
#ifdef VERBOSE
  printf("setMaterialPhong(%f)\n", phong);
#endif

  material.phong = phong;
  setParams |= SET_MAT_PHONG;
}

enum {
  USE_SOLID_COLOR,
  USE_MARBLE,
  USE_CHECKER
} colorMode = USE_SOLID_COLOR;

void setMaterialColor(double color[3]) {
#ifdef VERBOSE
  printf("setMaterialColor((%f,%f,%f))\n", color[0], color[1], color[2]);
#endif

  material.color[0] = color[0];
  material.color[1] = color[1];
  material.color[2] = color[2];
  setParams |= SET_MAT_COLOR;

  colorMode = USE_SOLID_COLOR;
}

struct {
  int n;              /* number of colors */
  double (*rgb)[3];   /* n rgb colors */
  double veinDir[3];  /* direction of marble "grain" */
  int octaves;        /* number of octaves of noise */
} marble = {0, NULL, {0,0,0}, 0};

void setMarbleTexture(VecList *colors, double veinDir[3], int octaves) {
  int i, n = length((List *) colors);
  double (*rgb)[3];
  VecList *list;
  
  rgb = (double (*)[3]) myalloc(3*n*sizeof(double));
  for (i = 0, list = colors; i < n; i++, list = list->next) {
    rgb[i][0] = list->v[0];
    rgb[i][1] = list->v[1];
    rgb[i][2] = list->v[2];
  }

  marble.n = n;
  marble.rgb = rgb;
  marble.veinDir[0] = veinDir[0];
  marble.veinDir[1] = veinDir[1];
  marble.veinDir[2] = veinDir[2];
  marble.octaves = octaves;

  colorMode = USE_MARBLE;
}

void setCheckerTexture(void) {
  colorMode = USE_CHECKER;
}

OBJECT *setMaterial(OBJECT *object) {
  object->ka = material.ka;
  object->kd = material.kd;
  object->ks = material.ks;
  object->kt = material.kt;
  object->ni = material.ni;
  object->setColor(object, material.color);

  if (colorMode == USE_MARBLE) {
    object = marbleObject(object, marble.n, marble.rgb, 
			  marble.veinDir, marble.octaves);
  } else if (colorMode == USE_CHECKER) {
    object = checkerObject(object);
  }

  return object;
}

void makeSphere(double center[3], double rad) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeSphere((%f,%f,%f), %f)\n", 
	 center[0], center[1], center[2], rad);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");

  object = createSphereObject(center, rad);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeSuperellipsoid(double n, double m, 
			double size[3], double z[3], double x[3], double center[3]) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeSuperellipsoid(%f,%f, (%f,%f,%f), (%f,%f,%f), (%f,%f,%f), (%f,%f,%f))\n",
	 n, m, size[0],size[1],size[2], z[0],z[1],z[2], x[0],x[1],x[2],
	 center[0],center[1],center[2]);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");

  object = createSuperellipsoidObject(n,m, size, z, x, center);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makePointCloud(const char *fname, double maxLevel, double sigma, double halfThickness) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makePointCloud(%s,%f,%f,%f)\n", 
	 fname, maxLevel, sigma, halfThickness);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");

  object = createPointCloudObject(fname, (int) maxLevel, (float) sigma, (float) halfThickness);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makePlane(double normal[3], double point[3]) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makePlane((%f,%f,%f), (%f,%f,%f))\n",
	 normal[0], normal[1], normal[2], point[0], point[1], point[2]);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");
    
  object = createPlaneObjectFromNormalandPoint(normal, point);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeImagePlane(double normal[3], double point[3],
		    double org[3], double yaxis[3],
		    double width, double height,
		    char *imagename) {
  OBJECT *object;
  pnm_image *image;

#ifdef VERBOSE
  printf("makePlane((%f,%f,%f), (%f,%f,%f))\n",
	 normal[0], normal[1], normal[2], point[0], point[1], point[2]);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");
    
  object = createPlaneObjectFromNormalandPoint(normal, point);
  object = setMaterial(object);

  image = read_pnm_image_from_file(imagename);
  mapImageToPlane(object, 1, 1, org, yaxis, width, height, image);

  objects[numObjects++] = object;
}

void makeBezier3(VecListList *mesh) {
  OBJECT *object;
  VecListList *vllist;
  POINT3 P[4][4];
  int r,c;

  if (length((List *) mesh) != 4)
    fatalError("Cubic Bezier mesh must have exactly 4 rows!");
  for (vllist = mesh; vllist != NULL; vllist = vllist->next)
    if (length((List *) vllist->list) != 4)
      fatalError("Cubic Bezier mesh must have exactly 4 columns per row!");

  for (r = 0, vllist = mesh; r < 4; r++, vllist = vllist->next) {
    VecList *vlist;
    for (c = 0, vlist = vllist->list; c < 4; c++, vlist = vlist->next) {
      P[r][c].x = vlist->v[0];
      P[r][c].y = vlist->v[1];
      P[r][c].z = vlist->v[2];
    }
  }

  object = (OBJECT *) createBezier3Object(P);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeTeapot(double org[3], double up[3], double spout[3], double scale) {
  OBJECT *object;

#ifdef VERBOSE
  printf("makeTeapot((%f,%f,%f), (%f,%f,%f), (%f,%f,%f), %f)\n",
	 org[0], org[1], org[2], 
	 up[0], up[1], up[2], 
	 spout[0], spout[1], spout[2], scale);
#endif

  if (numObjects >= MAX_OBJECTS)
    fatalError("Too friggin many objects!");
    
  object = createTeapotObject(org, up, spout, scale);
  object = setMaterial(object);

  objects[numObjects++] = object;
}

void makeElevationMap(double orgxy[], double widthheight[], 
                      double zminmax[], char *fname) {
  OBJECT *object;
  pnm_image *image = read_pnm_image_from_file(fname);
  double zscale, zshift;
  double *z;
  int r,c;

  if (!IS_PGM(image))
    fatalError("Elevation map is not a PGM image!");

  zscale = (zminmax[1] - zminmax[0])/PNM_MAXVAL(image);
  zshift = zminmax[0];
  z = (double *) myalloc(PNM_NR(image)*PNM_NC(image)*sizeof(double));
  for (r = 0; r < PNM_NR(image); r++)
    for (c = 0; c < PNM_NC(image); c++)
      z[r*PNM_NC(image) + c] = PGM_PIXEL(image, r,c)*zscale + zshift;

  object = 
    createSurfaceSupportedByRectangularGrid(PNM_NC(image), PNM_NR(image),
					    orgxy[0], orgxy[1],
					    widthheight[0], widthheight[1],
					    z, PNM_NC(image));
  object = setMaterial(object);
  objects[numObjects++] = object;
}

void makeHermiteFunc(double x[], double y[], 
                     int W, int H,
                     double *f, double *fx, double *fy, double *fxy) {  
  OBJECT *object;
  double *z = (double *) myalloc(W*H*sizeof(double));
  Func2 *func = createGeneralHermiteFunc2(x[0], x[1], y[0], y[1]);
  setGeneralHermiteFunc2(func, (double (*)[2]) f,
			 (double (*)[2]) fx, (double (*)[2]) fy,
			 (double (*)[2]) fxy);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], z, W);
  object =  createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						    x[1] - x[0], y[1] - y[0],
						    z, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;
}

void makeBilinearFis(double x[], double y[],
		     int W, int H,
		     double **alpha, double **z) {
  double X[4], Y[4];
  RectFis *fis;
  double *buf = (double *) myalloc(W*H*sizeof(double));
  Func2 *func;
  double dx = (x[1] - x[0])/3, dy = (y[1] - y[0])/3;
  OBJECT *object;
  X[0] = x[0];
  X[1] = x[0] + dx;
  X[2] = x[0] + 2*dx;
  X[3] = x[1];
  Y[0] = y[0];
  Y[1] = y[0] + dy;
  Y[2] = y[0] + 2*dy;
  Y[3] = y[1];
  fis = createBilinearFis(3, 3, X, Y);
  setBilinearFis(fis, z, alpha);
  func = rectFisToFunc2(fis);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						   x[1] - x[0], y[1] - y[0],
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}


void makeHermiteFis(double x[], double y[],
		    int W, int H,
		    double **alpha, 
		    double **z, double **zx, double **zy, double **zxy) {
  double X[4], Y[4];
  RectFis *fis;
  double *buf = (double *) myalloc(W*H*sizeof(double));
  Func2 *func;
  double dx = (x[1] - x[0])/3, dy = (y[1] - y[0])/3;
  OBJECT *object;
  X[0] = x[0];
  X[1] = x[0] + dx;
  X[2] = x[0] + 2*dx;
  X[3] = x[1];
  Y[0] = y[0];
  Y[1] = y[0] + dy;
  Y[2] = y[0] + 2*dy;
  Y[3] = y[1];
  fis = createHermiteFis(3, 3, X, Y);
  setHermiteFis(fis, alpha, z, zx, zy, zxy);
  func = rectFisToFunc2(fis);
  getFunc2Samples(func, W, H, x[0], y[0], x[1] - x[0], y[1] - y[0], buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x[0], y[0],
						   x[1] - x[0], y[1] - y[0],
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}

void makeBilinearRis(int M, int N, int Sx, int Sy, double x0, double y0,
	             double width, double height, int W, int H,
		     double alpha, Matrix *z) {
  Func2 *alphaFunc = createConstantFunc2(alpha);
  Ris *ris;
  Func2 *func;
  double *buf;
  int i,j;
  OBJECT *object;

  assert(W > 10 && H > 10);
  buf = (double *) myalloc(W*H*sizeof(double));

  assert(M >= 2 && N >= 2 &&
	 Sx >= 2 && (M % Sx) == 0 && 
	 Sy >= 2 && (N % Sy) == 0);

  ris = createBilinearRis(M, N, Sx, Sy, x0, y0, width, height, alphaFunc);
  for (j = 0; z != NULL; z = z->next, j++) {
    NumList *n = z->row;
    assert(j <= N);
    for (i = 0; n != NULL; n = n->next, i++) {
      assert(i <= M);
      ris->z[j][i] = n->num;
    }
  }

  assert(i > M && j > N);
  computeBilinearRisMaps(ris);
  func = risToFunc2(ris);

  getFunc2Samples(func, W, H, x0, y0, width, height, buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x0, y0,
						   width, height, 
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;  
}

void makeHermiteRis(int M, int N, int Sx, int Sy, double x0, double y0,
                    double width, double height, int W, int H, double alpha, 
                    Matrix *z, Matrix *zx, Matrix *zy, Matrix *zxy) {
  
  Func2 *alphaFunc = createConstantFunc2(alpha);
  Ris *ris;
  Func2 *func;
  double *buf;
  int i,j;
  OBJECT *object;

  assert(W > 10 && H > 10);
  buf = (double *) myalloc(W*H*sizeof(double));

  assert(M >= 2 && N >= 2 &&
	 Sx >= 2 && (M % Sx) == 0 && 
	 Sy >= 2 && (N % Sy) == 0);

  ris = createHermiteRis(M, N, Sx, Sy, x0, y0, width, height, alphaFunc);
  for (j = 0; z != NULL; j++,
       z = z->next, zx = zx->next, zy = zy->next, zxy = zxy->next) {
    NumList *n = z->row, *nx, *ny, *nxy;
    assert(j <= N && zx != NULL && zy != NULL && zxy != NULL);
    nx = zx->row; ny = zy->row; nxy = zxy->row;
    for (i = 0; n != NULL; i++, 
	   n = n->next, nx = nx->next, ny = ny->next, nxy = nxy->next) {
      assert(i <= M && nx != NULL && ny != NULL && nxy != NULL);
      ris->z[j][i] = n->num;
      ris->zx[j][i] = nx->num;
      ris->zy[j][i] = ny->num;
      ris->zxy[j][i] = nxy->num;
    }
  }

  assert(i > M && j > N);
  computeHermiteRisMaps(ris);
  func = risToFunc2(ris);

  getFunc2Samples(func, W, H, x0, y0, width, height, buf, W);
  object = createSurfaceSupportedByRectangularGrid(W, H, x0, y0,
						   width, height, 
						   buf, W);  
  object = setMaterial(object);
  objects[numObjects++] = object;    
}

extern FILE *yyin;  /* from ray.l */
extern int yylex(void);
extern int yyparse(void);

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "usage: %s <scene>\n", argv[0]);
    exit(-1);
  }

  if ((yyin = fopen(filename = argv[1], "r")) == NULL) {
    perror(filename);
    exit(-1);
  }

  setVariable("pi", M_PI);

  yyparse();

  fclose(yyin);

  return 0;
}
