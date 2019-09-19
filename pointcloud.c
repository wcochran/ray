#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "bbox.h"
#include "pointcloud.h"

typedef struct GridData {
    struct GridData *next;   // link to next bucket in hash tbl
    uint8_t i,j,k;           // grid indix
    uint8_t color[3];        // RGB (8-bit fixed pt)
    float f;                 // signed distance function val
    float normal[3];         // surface normal
} GridData;

typedef struct {  // holds data indexed on (i,j,k)
    int numElems;
    GridData **table;
} GridHashTable;

static
GridHashTable *createGridHashTable(int size) {
    GridHashTable *hashTable = (GridHashTable*)malloc(sizeof(GridHashTable));
    hashTable->numElems = size;
    hashTable->table = (GridData**) malloc(size*sizeof(GridData*));
    for (int i = 0; i < size; i++)
        hashTable->table[i] = NULL;
    return hashTable;
}

//https://en.wikipedia.org/wiki/Pairing_function 
#define HASH2(i,j) (((i) + (j))*((i) + (j) + 1)/2 + (j))
#define HASH3(i,j,k) HASH2(HASH2(i,j),k)

static
GridData *gridDataLookup(GridHashTable *hashTable, int i, int j, int k) {
    const int index = HASH3(i,j,k) % hashTable->numElems;
    GridData *n = hashTable->table[index];
    while (n != NULL && (n->i != i || n->j != j || n->k != k))
        n = n->next;
    return n;
}

// Note: reuses gridData node!
static
GridData *gridDataInsert(GridHashTable *hashTable, GridData *gridData) { 
    const int i = gridData->i, j = gridData->j, k = gridData->k;
    const int index = HASH3(i,j,k) % hashTable->numElems;
    GridData *n = hashTable->table[index];
    gridData->next = n;
    hashTable->table[index] = gridData;
    return gridData;
}

typedef struct {
    int i,j,k;
} Index;

typedef struct {
    Index min;
    Index max;
} IndexBox;

//
// Three case for octree node
//   node = NULL => empty leaf
//   node != NULL and node->child == NULL => filled leaf
//   node != NULL and node->child != NULL => internal node with 
//              8 children (some of which may be NULL)
//  The actual data for a filled leaf is stored in a GridDataHashTable
//  which is indexed at its eight corners.
//
typedef struct OctreeNode {
    struct OctreeNode **child;  // NULL => filled leaf
} OctreeNode;                   // not NULL => internal node

typedef struct {
    int maxLevel;
    BBOX bbox;
    IndexBox indexBox;
    GridHashTable *gridHashTable;
    OctreeNode *root;
} Octree;

typedef struct {
    float x,y,z;
} Vec3;

typedef struct {
    float r, g, b;
} Color;

typedef struct {
    Vec3 pos;
    Color color;
    Vec3 normal;
} Vert;

typedef struct {
    int capacity;
    int num;
    Vert *verts;
} VertArray;

static
VertArray *createVertArray(int initialSize, int initialCapacity) {
    assert(initialSize <= initialCapacity);
    VertArray *vertArray = (VertArray*) malloc(sizeof(VertArray));
    vertArray->capacity = initialCapacity;
    vertArray->num = initialSize;
    vertArray->verts = (Vert *) calloc(initialCapacity,sizeof(Vert));
    return vertArray;
}

static
void destroyVertArray(VertArray *vertArray) {
    free(vertArray->verts);
    free(vertArray);
}

static
void addVert(VertArray *vertArray, Vert *vert) { // copies vertex
    if (vertArray->num >= vertArray->capacity) {
        vertArray->capacity *= 2;
        vertArray->verts = (Vert *) realloc(vertArray->verts, vertArray->capacity * sizeof(Vert));
        assert(vertArray->verts != NULL);
    }
    vertArray->verts[vertArray->num++] = *vert;
}

// Returns NULL if not vert's in bbox.
static
VertArray *vertsInsideExpandedBBox(VertArray *inputVertArray, BBOX *bbox, double r) {
    VertArray *vertArray = NULL;
    for (int i = 0; i < inputVertArray->num; i++) {
        Vert *vert = &inputVertArray->verts[i];
        if (bbox->min.x - r <= vert->pos.x && vert->pos.x < bbox->max.x + r &&
            bbox->min.y - r <= vert->pos.y && vert->pos.y < bbox->max.y + r &&
            bbox->min.z - r <= vert->pos.z && vert->pos.z < bbox->max.z + r) {
            if (vertArray == NULL) {
                vertArray = createVertArray(0, 10);
            }
            addVert(vertArray, vert);
        }
    }
    return vertArray;
}

static
double signedDistanceFunction(VertArray *verts, POINT3 *P, POINT3 *N, double color[3], 
        double radius, double sigma) {
    const double r2 = radius*radius;
    const double b = -1/(sigma*sigma);
    float weight = 0;
    POINT3 weightP = {0,0,0};
    POINT3 weightN = {0,0,0};
    double weightColor[3] = {0,0,0};
    int n = 0;
    for (int i = 0; i < verts->num; i++) {
        Vert *vert = &verts->verts[i];
        const double dx = vert->pos.x - P->x;
        const double dy = vert->pos.y - P->y;
        const double dz = vert->pos.z - P->z;
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < r2) {
            const double w = exp(b*d2);
            weight += w;
            weightP.x += w*vert->pos.x;
            weightP.y += w*vert->pos.y;
            weightP.z += w*vert->pos.z;
            weightN.x += w*vert->normal.x;
            weightN.y += w*vert->normal.y;
            weightN.z += w*vert->normal.z;
            weightColor[0] += w*vert->color.r;
            weightColor[1] += w*vert->color.g;
            weightColor[2] += w*vert->color.b;
            n++;
        }
    }
    const double scale = weight > 0 ? 1/weight : 1;
    POINT3 Pave = {scale*weightP.x, scale*weightP.y, scale*weightP.z};
    N->x = scale*weightN.x;
    N->y = scale*weightN.y;
    N->z = scale*weightN.z;
    const double mag2 = N->x*N->x + N->y*N->y + N->z*N->z; // renormalize N
    const double s = mag2 > 0 ? 1/sqrt(mag2) : 1;
    N->x *= s;
    N->y *= s;
    N->z *= s;
    color[0] = scale*weightColor[0];
    color[1] = scale*weightColor[1];
    color[2] = scale*weightColor[2];
    POINT3 D = {P->x - Pave.x, P->y - Pave.y, P->z - Pave.z};
    const double f = D.x*N->x + D.y*N->y + D.z*N->z;
    return f;
}

static
void createFilledLeaf(VertArray *verts, 
        BBOX *bbox, IndexBox *indexBox, 
        float radius, float sigma, 
        GridHashTable *hashTable) {
    for (int k = 0; k < 2; k++) {
        const double z = (k == 0) ? bbox->min.z : bbox->max.z;
        const int K = (k == 0) ? indexBox->min.k : indexBox->max.k;
        for (int j = 0; j < 2; j++) {
            const double y = (j == 0) ? bbox->min.y : bbox->max.y;
            const int J = (j == 0) ? indexBox->min.j : indexBox->max.j;
            for (int i = 0; i < 2; i++) {
                const double x = (i == 0) ? bbox->min.x : bbox->max.x;
                const int I = (i == 0) ? indexBox->min.i : indexBox->max.i;
                GridData *gridData = gridDataLookup(hashTable, I,J,K);
                if (gridData == NULL) {
                    gridData = (GridData*) malloc(sizeof(GridData));
                    POINT3 P = {x, y, z};
                    POINT3 N;
                    double color[3];
                    const double f = signedDistanceFunction(verts, &P, &N, color, 
                                        radius, sigma);
                    gridData->f = (float)f;
                    for (int ch = 0; ch < 3; ch++)
                        gridData->color[ch] = (uint8_t) (255*color[ch] + 0.5);
                    gridData->normal[0] = (float) N.x;
                    gridData->normal[1] = (float) N.y;
                    gridData->normal[2] = (float) N.z;
                    gridData->i = (uint8_t) I;
                    gridData->j = (uint8_t) J;
                    gridData->k = (uint8_t) K;
                    gridDataInsert(hashTable, gridData);
                }
            }
        }
    }
}

Vec3 bboxCenter(BBOX *bbox) {
    Vec3 C;
    C.x = (bbox->min.x + bbox->max.x)/2;
    C.y = (bbox->min.y + bbox->max.y)/2;
    C.z = (bbox->min.z + bbox->max.z)/2;
    return C;
}

Index IndexBoxCenter(IndexBox *iBox) {
    Index centerIndex;
    centerIndex.i = (iBox->min.i + iBox->max.i)/2;
    centerIndex.j = (iBox->min.j + iBox->max.j)/2;
    centerIndex.k = (iBox->min.k + iBox->max.k)/2;
    return centerIndex;
}

BBOX createBoundingBox(double xmin, double ymin, double zmin,
        double xmax, double ymax, double zmax) {
    BBOX bbox = {{xmin, ymin, zmin}, {xmax, ymax, zmax}};
    return bbox;
}

static
void subdivide(OctreeNode *node, int level, int maxLevel,
        BBOX *bbox, IndexBox *iBox, GridHashTable *hashTable,
        VertArray *verts, float sigma, float radius) {
    Vec3 C = bboxCenter(bbox);
    double X[3] = {bbox->min.x, C.x, bbox->max.x};
    double Y[3] = {bbox->min.y, C.y, bbox->max.y};
    double Z[3] = {bbox->min.z, C.z, bbox->max.z};

    Index centerIndex = IndexBoxCenter(iBox);
    int I[3] = {iBox->min.i, centerIndex.i, iBox->max.i};
    int J[3] = {iBox->min.j, centerIndex.j, iBox->max.j};
    int K[3] = {iBox->min.k, centerIndex.k, iBox->max.k};
    
//    // XXX
//    printf("bbox (%0.4f, %0.4f, %0.4f) - (%0.4f, %0.4f, %0.4f)\n",
//           bbox->min.x, bbox->min.y, bbox->min.z,
//           bbox->max.x, bbox->max.y, bbox->max.z);

    int n = 0;

    node->child = (OctreeNode **) malloc(8*sizeof(OctreeNode*));

    BBOX childBBox;
    IndexBox childIndexBox;
    
    for (int k = 0; k < 2; k++) {
        childBBox.min.z = Z[k];
        childBBox.max.z = Z[k+1];
        childIndexBox.min.k = K[k];
        childIndexBox.max.k = K[k+1];
        for (int j = 0; j < 2; j++) {
            childBBox.min.y = Y[j];
            childBBox.max.y = Y[j+1];
            childIndexBox.min.j = J[j];
            childIndexBox.max.j = J[j+1];
            for (int i = 0; i < 2; i++) {
                childBBox.min.x = X[i];
                childBBox.max.x = X[i+1];
                childIndexBox.min.i = I[i];
                childIndexBox.max.i = I[i+1];
                
//                printf("   child bbox %d: (%0.4f, %0.4f, %0.4f) - (%0.4f, %0.4f, %0.4f)\n",
//                       n,
//                       childBBox.min.x, childBBox.min.y, childBBox.min.z,
//                       childBBox.max.x, childBBox.max.y, childBBox.max.z);

                VertArray *childVerts = 
                    vertsInsideExpandedBBox(verts, &childBBox, radius);
                if (childVerts == NULL) {
                    node->child[n] = NULL; // empty 
                } else if (level >= maxLevel) {
                    node->child[n] = (OctreeNode*) malloc(sizeof(OctreeNode));
                    node->child[n]->child = NULL; // filled leaf
                    createFilledLeaf(childVerts, &childBBox, &childIndexBox, 
                            radius, sigma, hashTable);
                } else { // internal node;
                    node->child[n] = (OctreeNode*) malloc(sizeof(OctreeNode));
                    subdivide(node->child[n], level + 1, maxLevel,
                        &childBBox, &childIndexBox, hashTable, childVerts,
                        radius, sigma);
                }

                if (childVerts != NULL)
                    destroyVertArray(childVerts);

                n++;
            }
        }
    }
    
}

static
void setIndexBox(IndexBox *ibox, int maxLevel) {
    ibox->min.i = ibox->min.j = ibox->min.k = 0;
    ibox->max.i = ibox->max.j = ibox->max.k = 1 << maxLevel;
}

void getBoundingBox(BBOX *bbox, VertArray *pointCloud) {
    assert(pointCloud->num > 0);
    bbox->min.x = bbox->max.x = pointCloud->verts[0].pos.x;
    bbox->min.y = bbox->max.y = pointCloud->verts[0].pos.y;
    bbox->min.z = bbox->max.z = pointCloud->verts[0].pos.z;
    for (int i = 1; i < pointCloud->num; i++) {
        Vec3 *p = &pointCloud->verts[i].pos;
        if (p->x < bbox->min.x) bbox->min.x = p->x;
        if (p->y < bbox->min.y) bbox->min.y = p->y;
        if (p->z < bbox->min.z) bbox->min.z = p->z;
        if (p->x > bbox->max.x) bbox->max.x = p->x;
        if (p->y > bbox->max.y) bbox->max.y = p->y;
        if (p->z > bbox->max.z) bbox->max.z = p->z;
    }
}

void inflateBoundingBox(BBOX *bbox, double fraction) {
    const double dx = fraction*(bbox->max.x - bbox->min.x)/2;
    const double dy = fraction*(bbox->max.y - bbox->min.y)/2;
    const double dz = fraction*(bbox->max.z - bbox->min.z)/2;
    bbox->min.x -= dx;
    bbox->min.y -= dy;
    bbox->min.z -= dz;
    bbox->max.x += dx;
    bbox->max.y += dy;
    bbox->max.z += dz;
}

void cubeBoundingBox(BBOX *bbox) {
    const double w = bbox->max.x - bbox->min.x;
    const double h = bbox->max.y - bbox->min.y;
    const double d = bbox->max.z - bbox->min.z;
    float size = w;
    if (h > size) size = h;
    if (d > size) size = d;
    const double dx = (size - w)/2;
    const double dy = (size - h)/2;
    const double dz = (size - d)/2;
    bbox->min.x -= dx;
    bbox->min.y -= dy;
    bbox->min.z -= dz;
    bbox->max.x += dx;
    bbox->max.y += dy;
    bbox->max.z += dz;
}

static
Octree *createOctree(VertArray *pointCloud, 
        int maxLevel,
        float sigma, float radius) {
    Octree *octree = (Octree*) malloc(sizeof(Octree));
    octree->maxLevel = maxLevel;
    getBoundingBox(&octree->bbox, pointCloud);
    inflateBoundingBox(&octree->bbox, 0.01);
    cubeBoundingBox(&octree->bbox);
    setIndexBox(&octree->indexBox, maxLevel);
    octree->gridHashTable = createGridHashTable(10003);
    octree->root = (OctreeNode*) malloc(sizeof(OctreeNode));
    octree->root->child = NULL;
    subdivide(octree->root, 1, maxLevel,
            &octree->bbox, &octree->indexBox,
            octree->gridHashTable, pointCloud, 
            sigma, radius);
    return octree;
}

//
// Convert ray from "world coordinates" to 
// "unit bounding box coordinates."
//
// Let (x0,y0,z0) = bbox->min, (x1,y1,z1) = bbox->max
// Given a point P = (x,y,z), the point U relative to the box is
// U = ((x - x0)/(x1 - x0), (y - y0)/(y1 - y0), (z - z0)/(z1 - z0)).
//
// Let R(t) = a + bt be a point on the ray at a distance t
// where a = rayOrg, b = rayDir. The corresponding point U(t) is
//    U(t) = ((ax + bx*t - x0)/(x1 - x0), ...)
//         = ((ax - x0)/(x1 - x0) + bx/(x1 - x0)*t, ....)
// So the normalized ray is U(t) = A + B*t where
//    Ax = (ax - x0)/(x1 - x0), Bx = bx/(x1 - x0)
// and similarly for Ay,By and Az,z
//
// Note that unitized ray direction no longer has length 1.
//
static
void unitRay(BBOX *bbox, 
        double rayOrg[3], double rayDir[3],
        double unitRayOrg[3], double unitRayDir[3]) {
    const double sx = 1/(bbox->max.x - bbox->min.x);
    const double sy = 1/(bbox->max.y - bbox->min.y);
    const double sz = 1/(bbox->max.z - bbox->min.z);
    // note sx = sy = sz if bbox is a cube
    unitRayOrg[0] = (bbox->max.x - rayOrg[0])*sx;
    unitRayOrg[1] = (bbox->max.y - rayOrg[1])*sy;
    unitRayOrg[2] = (bbox->max.z - rayOrg[2])*sz;
    unitRayDir[0] = rayDir[0]*sx;
    unitRayDir[1] = rayDir[1]*sy;
    unitRayDir[2] = rayDir[2]*sz;
}

//
// Given the ray direction and bounding box for a filled
// octree leaaf, we find the coefficeints for a cubic polynomial
//    f(t) = A*t^3 + B*t^2 + C*t + D
// whose zeros occur at the t-values where the ray intersects
// the underlying trilinear implicit surface.
// See Section 2 of "Fast and Accurate Ray-Voxel Intersection Techniques
// fir Iso Surface Ray Tracing" by Marmitt, et al 
// ftp://ftp.mpi-sb.mpg.de/pub/conferences/vmv04/submissions/153/isoisec.pdf
// https://www.cs.utah.edu/~shirley/papers/iso/iso.pdf
//
//
static
void rayCubicTrilinearInterpolater(double rayOrg[3], double rayDir[3],
        BBOX *bbox, 
        IndexBox *ibox, GridHashTable *hashTable, // grid information
        double coeff[4]) {
    assert(ibox->max.i - ibox->min.i == 1); // leafs all at the same bottom level
    assert(ibox->max.j - ibox->min.j == 1);
    assert(ibox->max.k - ibox->min.k == 1);
    
    //
    // Convert input ray R(t) to "unit ray" U(t) that yields points
    // relative to "unit cube."
    //
    double unitRayOrg[3], unitRayDir[3];
    unitRay(bbox, rayOrg, rayDir, unitRayOrg, unitRayDir);
    
    //
    // We have two rays who components are used as weights
    // for trilinear interpolation of 8 corner values:
    //    U0(t) = 1 - (org + dir*t) = (1 - org) + (-dir)*t
    //    U1(t) = org + dir*t
    //
    const double uorg[2] = {1 - unitRayOrg[0], unitRayOrg[0]};
    const double vorg[2] = {1 - unitRayOrg[1], unitRayOrg[1]};
    const double worg[2] = {1 - unitRayOrg[2], unitRayOrg[2]};
    
    const double udir[2] = {-unitRayDir[0], unitRayDir[0]};
    const double vdir[2] = {-unitRayDir[1], unitRayDir[1]};
    const double wdir[2] = {-unitRayDir[2], unitRayDir[2]};

    //
    // Compute coefficients for cubic polynomal that yields
    // the appropriate interpolatied signed distance value
    // at each point along the ray.
    //
    double A = 0, B = 0, C = 0, D = 0;
    for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
                GridData *g = gridDataLookup(hashTable, // assume filled leaf
                    ibox->min.i+i, ibox->min.j+j, ibox->min.k+k);
                assert(g != NULL);
                double f = g->f;
                A += udir[i]*vdir[j]*wdir[k]*f;
                B += (uorg[i]*vdir[j]*wdir[k] + 
                      udir[i]*vorg[j]*wdir[k] +
                      udir[i]*vdir[j]*worg[k])*f;
                C += (udir[i]*vorg[j]*worg[k] + 
                      uorg[i]*vdir[j]*worg[k] +
                      uorg[i]*vorg[j]*wdir[k])*f;
                D += uorg[i]*vorg[j]*worg[k]*f;
            }

    coeff[0] = A;
    coeff[1] = B;
    coeff[2] = C;
    coeff[3] = D;
}

static
int quadratic(double B, double C, double roots[2]) {
  double det = B*B - 4*C;
  if (det < 0.0) return 0;  /* imaginary */
  det = sqrt(det);
  roots[0] = 0.5*(-B + det);
  roots[1] = 0.5*(-B - det);
  return 1;
}

//
// Eval cubic polynomial
//    F(t) = A*t^3 + B*t^2 + C*t + D
//
static
inline double F(double t, double A, double B, double C, double D) {
    return ((A*t + B)*t + C)*t + D;
}

//
// See Section 4 (Algorithm 3) of "Fast and Accurate Ray-Voxel Intersection Techniques
// for Iso Surface Ray Tracing" by Marmitt, et al 
// ftp://ftp.mpi-sb.mpg.de/pub/conferences/vmv04/submissions/153/isoisec.pdf
//
static
double cubicRootFinder(double C[4], double tin, double tout) {

    //
    // Handle constant, linear, or quadratic case.
    //
    if (fabs(C[0]) <= 0) {
        if (fabs(C[1]) <= 0) {
            if (fabs(C[2]) <= 0)  
                return -1;  // constant -- no root
            double t = C[3]/C[2];
            if (tin <= t && t <= tout) return t; // linear -- one root
            return -1.0;
        }
        double roots[2];
        if (quadratic(C[2]/C[1], C[3]/C[1], roots)) { // quadratic -- two roots
            double t;
            if (roots[0] < tin && roots[1] < tin) return -1;
            if (roots[0] < tin) t = roots[1];
            else if (roots[1] < tin) t = roots[0];
            else t = roots[0] < roots[1] ? roots[0] : roots[1];
            if (t <= tout) return t;
        }
        return -1.0;
    }

    //
    // Follow Algorithm 3 in paper.
    //

    double t0 = tin;
    double t1 = tout;
    double f0 = F(t0, C[0],C[1],C[2],C[3]);
    double f1 = F(t1, C[0],C[1],C[2],C[3]);

    double droots[2];
    if (quadratic(2*C[1]/(3*C[0]), C[2]/(3*C[0]), droots)) { // has 2 real roots
        if (droots[0] > droots[1]) { // sort roots in ascending order
            double tmp = droots[0];
            droots[0] = droots[1];
            droots[1] = tmp;
        }
        for (int i = 0; i < 2; i++) {
            double t = droots[i];
            if (t0 <= t && t <= t1) {
                double f = F(t, C[0],C[1],C[2],C[3]);
                if (signbit(f) == signbit(f0)) {
                    t0 = t;
                    f0 = f;
                } else {
                    t1 = t;
                    f1 = f;
                }
            }
        }
    }

    if (signbit(f0) == signbit(f1)) return -1.0;

    const int numIters = 3;
    for (int i = 0; i < numIters; i++) {
        assert(f0 != f1);
        double t = t0 + (t1 - t0)*(-f0/(f1 - f0));
        double f = F(t, C[0],C[1],C[2],C[3]);
        if (signbit(f) == signbit(f0)) {
            t0 = t;
            f0 = f;
        } else {
            t1 = t;
            f1 = f;
        }
    }

    assert(f0 != f1);
    double thit = t0 + (t1 - t0)*(-f0/(f1 - f0));
    return thit;
}

typedef struct {    // Child ray intersection info
    double t[2];    // entry / exit ray parameter
    int n;          // octree index
    BBOX bbox;      // bounding box (interior node)
    IndexBox ibox;  // grid indices
} ChildInfo;

//
// Sorting octree children info on t[0] = smallest ray
// entry point.
// https://stackoverflow.com/questions/3886446/problem-trying-to-use-the-c-qsort-function
//
static
int childInfoPtrCompare(const ChildInfo **A, const ChildInfo **B) {
    const double fa = (*A)->t[0];
    const double fb = (*B)->t[0];
    return (fa > fb) - (fa < fb);
}

static
double rayIntersection(double rayOrg[3], double rayDir[3],
        OctreeNode *node, BBOX *bbox, IndexBox *ibox,
        GridHashTable *hashTable,
        double tinterval[2], Index *hitIndex) {
    double thit = -1.0;

    if (node->child == NULL) {  // filled leaf
        thit = tinterval[0]; // XXXX
//        double C[4];
//        rayCubicTrilinearInterpolater(rayOrg, rayDir,
//            bbox, ibox, hashTable, C);
//        thit = cubicRootFinder(C, tinterval[0], tinterval[1]);
        if (thit >= EPSILON) {
            *hitIndex = ibox->min;
        }
    } else { // internal node
        Vec3 C = bboxCenter(bbox);
        double X[3] = {bbox->min.x, C.x, bbox->max.x};
        double Y[3] = {bbox->min.y, C.y, bbox->max.y};
        double Z[3] = {bbox->min.z, C.z, bbox->max.z};

        Index centerIndex = IndexBoxCenter(ibox);
        int I[3] = {ibox->min.i, centerIndex.i, ibox->max.i};
        int J[3] = {ibox->min.j, centerIndex.j, ibox->max.j};
        int K[3] = {ibox->min.k, centerIndex.k, ibox->max.k};

        int n = 0;

        int numChildren = 0;
        ChildInfo childInfo[8];

        int childIntersections = 0;
        
        for (int k = 0; k < 2; k++) {
            for (int j = 0; j < 2; j++) {
                for (int i = 0; i < 2; i++) {

                    if (node->child[n] != NULL) {
                        BBOX cbox = {{X[i], Y[j], Z[k]}, {X[i+1], Y[j+1], Z[k+1]}};
                        double tchild[2];
                        if (rayHitsBoundingBox(&cbox, rayOrg, rayDir, tchild)) {
                            childIntersections++;
                            assert(tchild[0] >= 0.0);
                            if (tchild[0] <= tinterval[1] && tchild[1] >= tinterval[0]) {
                                ChildInfo cinfo = {
                                    {tchild[0], tchild[1]},
                                    n,
                                    cbox,
                                    {{I[i], J[j], K[k]}, {I[i+1], J[j+1], K[k+1]}}
                                };
                                childInfo[numChildren++] = cinfo;
                            }
                        }
                    }

                    n++;
                }
            }
        }
        
        // XXX assert(childIntersections > 0);
        // XXX assert(numChildren > 0);
        if (childIntersections <= 0) {
            printf("rayOrg = (%f, %f, %f)\n", rayOrg[0], rayOrg[1], rayOrg[2]);
            printf("rayDir = (%f, %f, %f)\n", rayOrg[0], rayOrg[1], rayOrg[2]);
            double hitEntry[3], hitExit[3];
            for (int i = 0; i < 3; i++) {
                hitEntry[i] = rayOrg[i] + tinterval[0]*rayDir[i];
                hitExit[i] = rayOrg[i] + tinterval[1]*rayDir[i];
            }
            printf("hitEntry = (%f, %f, %f) t=%f\n",
                   hitEntry[0], hitEntry[1], hitEntry[2], tinterval[0]);
            printf("hitExit = (%f, %f, %f) t=%f\n",
                   hitExit[0], hitExit[1], hitExit[2], tinterval[1]);
            printf("bbox.min = (%f, %f, %f)\n", bbox->min.x, bbox->min.y, bbox->min.z);
            printf("bbox.max = (%f, %f, %f)\n", bbox->max.x, bbox->max.y, bbox->max.z);
            // XXX assert(childIntersections > 0);
        }
        
        if (numChildren > 0) {
            ChildInfo *childInfoPtr[8];
            for (int i = 0; i < numChildren; i++)
                childInfoPtr[i] = &childInfo[i];
            qsort(childInfoPtr, numChildren, sizeof(ChildInfo*), 
                (int (*)(const void *, const void *)) childInfoPtrCompare);
            for (int i = 0; i < numChildren; i++) {
                double t0 = childInfoPtr[i]->t[0];
                if (thit < 0 || t0 < thit) {
                    int n = childInfoPtr[i]->n;
                    BBOX *cbox = &childInfoPtr[i]->bbox;
                    IndexBox *ibox = &childInfoPtr[i]->ibox;
                    double t1 = childInfoPtr[i]->t[1];
                    double tint[2] = {t0, t1};
                    Index closestHitIndex;
                    double t = rayIntersection(rayOrg, rayDir, node->child[n], 
                                    cbox, ibox, hashTable, tint, &closestHitIndex);
                    if (t > 0 && (thit < 0 || t < thit)) {
                        thit = t;
                        *hitIndex = closestHitIndex;
                    }
                }
            }
        }

    }

    return thit;
}

VertArray *readPLY(const char *fname) {
    FILE *f = fopen(fname, "rb'");
    if (f == NULL) {perror(fname); return NULL;}
    char buf[200];
    if (fgets(buf, sizeof(buf), f) == NULL || strncmp("ply", buf, 3) != 0) {
        fprintf(stderr, "%s: bogus ply file!\n", fname);
        fclose(f);
        return NULL;
    }
    while (fgets(buf, sizeof(buf), f) != NULL && 
        strncmp("element vertex", buf, 14) != 0)
        ;
    int numVerts = 0;
    if (sscanf(buf + 14, "%d", &numVerts) != 1 || numVerts < 1) {
        fprintf(stderr, "%s: bogus number of vertices!\n", fname);
        fclose(f);
        return NULL;
    }
    int numProperties = 0;
    while (fgets(buf, sizeof(buf), f) != NULL &&
            strncmp("property", buf, 8) == 0)
        numProperties++;
    if (numProperties != 10) {
        fprintf(stderr, 
            "%s: should be exactly 10 properties (verts:3 + colors:4 + normals:3)!\n", 
            fname);
        fclose(f);
        return NULL;
    }
    
    while (fgets(buf, sizeof(buf), f) != NULL &&
            strncmp("end_header", buf, 8) != 0)
        ;

    VertArray *pointCloud = createVertArray(0, numVerts);

    for (int n = 0; n < numVerts; n++) {
        float pos[3];
        uint8_t color[4];
        float normal[3];
        if (fread(pos, sizeof(float), 3, f) != 3 ||
            fread(color, sizeof(uint8_t), 4, f) != 4 ||
            fread(normal, sizeof(float), 3, f) != 3) {
            fprintf(stderr, 
                "%s: Error reading vertex %d info!\n", 
                fname, n);
            fclose(f);
            return NULL;
        }
        Vert vert = {
            {pos[0], pos[1], pos[2]},
            {color[0]/255.0, color[1]/255.0, color[2]/255.0}, 
            {normal[0], normal[1], normal[2]}
        };
        addVert(pointCloud, &vert);
    }
    fclose(f);

    return pointCloud;
}

typedef struct {
    VertArray *pointCloud;
    Octree *octree;
} POINTCLOUD_DATA;

static
double rayHit(struct OBJECT *this, 
            double rayOrg[3],
            double rayDir[3],
            HIT_INFO *info) {
    POINTCLOUD_DATA *data = (POINTCLOUD_DATA*) this->data;
    Octree *octree = data->octree;
    double tbox[2];
    if (rayHitsBoundingBox(&octree->bbox, rayOrg, rayDir, tbox)) {
        Index hitIndex;
        double t = rayIntersection(rayOrg, rayDir,
            octree->root, 
            &octree->bbox, 
            &octree->indexBox, octree->gridHashTable, 
            tbox, &hitIndex);
        info->pointcloud.i = hitIndex.i;
        info->pointcloud.j = hitIndex.j;
        info->pointcloud.k = hitIndex.k;
        return t;
    }
    return -1.0;
}

static
void getTrilinearInterpolationWeights(Octree *octree,
        Index *index, double hitPoint[3],
        double U[2], double V[2], double W[2]) {
    const int M = 1 << octree->maxLevel;
    const double dx = (octree->bbox.max.x - octree->bbox.min.x)/M;
    const double dy = (octree->bbox.max.y - octree->bbox.min.y)/M;
    const double dz = (octree->bbox.max.z - octree->bbox.min.z)/M;
    // if bbox a cube then dx = dy = dz
    const double x0 = index->i*dx + octree->bbox.min.x;
    const double y0 = index->j*dy + octree->bbox.min.y;
    const double z0 = index->k*dz + octree->bbox.min.z;
    const double u = (hitPoint[0] - x0)/dx;
    const double v = (hitPoint[1] - y0)/dy;
    const double w = (hitPoint[2] - z0)/dz;
    U[0] = 1-u; U[1] = u;
    V[0] = 1-v; V[1] = v;
    W[0] = 1-w; W[1] = w;
}

static
void getNormal(struct OBJECT *this, 
            double hitPoint[3],
            HIT_INFO *info,  
            double normal[3]) {
    POINTCLOUD_DATA *data = (POINTCLOUD_DATA*) this->data;
    Octree *octree = data->octree;
    GridHashTable *hashTable = octree->gridHashTable;
    Index index = {
        info->pointcloud.i,
        info->pointcloud.j,
        info->pointcloud.k
    };

    //
    // Determine hitPoint relative to the grid cell that
    // the surface intersection occurred in.
    // If the bounding box is a cube, then dx = dy = dz.
    //
    double U[2], V[2], W[2];
    getTrilinearInterpolationWeights(octree, &index, hitPoint,
        U, V, W);

    //
    // Determine resulting normal via trilinear interpolation
    // of the 8 corner normals of the grid cell.
    // The grid normals are stored in our grid hash table.
    //
    normal[0] = normal[1] = normal[2] = 0;
    for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
                GridData *g = gridDataLookup(hashTable,
                                index.i,index.j,index.k);
                assert(g != NULL);
                const double w = U[i]*V[j]*W[k];
                for (int d = 0; d < 3; d++)
                    normal[d] += w*g->normal[d];
            }

    //
    // Normalize the normal so that it has unit length.
    //
    double m2 = normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2];
    double s = m2 > 0 ? 1/sqrt(m2) : 1;
    normal[0] *= s;
    normal[1] *= s;
    normal[2] *= s;
}

static
void getColor(struct OBJECT *this, 
            double hitPoint[3],
            HIT_INFO *info,  
            double color[3]) {
   POINTCLOUD_DATA *data = (POINTCLOUD_DATA*) this->data;
    Octree *octree = data->octree;
    GridHashTable *hashTable = octree->gridHashTable;
    Index index = {
        info->pointcloud.i,
        info->pointcloud.j,
        info->pointcloud.k
    };

    //
    // Determine hitPoint relative to the grid cell that
    // the surface intersection occurred in.
    // If the bounding box is a cube, then dx = dy = dz.
    //
    double U[2], V[2], W[2];
    getTrilinearInterpolationWeights(octree, &index, hitPoint,
        U, V, W);

    //
    // Determine resulting normal via trilinear interpolation
    // of the 8 corner normals of the grid cell.
    // The grid normals are stored in our grid hash table.
    //
    color[0] = color[1] = color[2] = 0;
    for (int k = 0; k < 2; k++)
        for (int j = 0; j < 2; j++)
            for (int i = 0; i < 2; i++) {
                GridData *g = gridDataLookup(hashTable,
                                index.i,index.j,index.k);
                assert(g != NULL);
                const double w = U[i]*V[j]*W[k];
                for (int d = 0; d < 3; d++)
                    color[d] += w*g->color[d]/255.0; // convert from 8-bit
            }
}

static
void setColor(OBJECT *this, double color[3]) {
 // XXX    POINTCLOUD_DATA *data = (POINTCLOUD_DATA *) this->data;
}

//typedef struct {
//    Vec3 point;
//    double distance;
//} PointDist;
//
//typedef struct {
//    int capacity;
//    int num;
//    PointDist *points;
//} PointDistArray;
//
//static
//PointDistArray *createPointDistArray(int maxCapacity) {
//    PointDistArray *array = (PointDistArray*) malloc(sizeof(PointDistArray));
//    array->capacity = maxCapacity;
//    array->num = 0;
//    array->points = (PointDist*) calloc(maxCapacity, sizeof(PointDist));
//    return array;
//}
//
//static
//void destroyPointDistanceArray(PointDistArray *array) {
//    free(array->points);
//    free(array);
//}
//
//static
//void insertPointDistance(PointDistArray *array, Vec3 *point, double distance) {
//    const int N = array->num;
//    const int C = array->capacity;
//    int i = 0;
//    while (i < N && array->points[i].distance <= distance)
//        i++;
//    if (i < C) {
//        int L = N;
//        if (L > C-1)
//            L = C-1;
//        for (int j = L; j > i; --j)
//            array->points[j] = array->points[j-1];
//        array->points[i].point = *point;
//        array->points[i].distance = distance;
//        if (N < C)
//            array->num++;
//    }
//}
//
//static
//void nearestNeighborsAux(OctreeNode *node, IndexBox *ibox, int level, int maxLevel) {
//    
//}

//static
//double meanDistanceBetweenPoints(VertArray *pointCloud, double *stdev) {
//    const int numPoints = pointCloud->num;
//    double dist = 0;
//    double dist2 = 0;
//    for (int i = 0; i < n; i++) {
//        const Vert *vert = pointCloud->verts[i];
//
//    }
//}

OBJECT *createPointCloudObject(const char *plyfile, int maxLevel, float sigma) {
    VertArray *pointcloud = readPLY(plyfile);
    if (pointcloud == NULL) return NULL;
    float radius = 2.5*sigma;
    Octree *octree = createOctree(pointcloud, maxLevel, sigma, radius);

    OBJECT *object = (OBJECT *) malloc(sizeof(OBJECT));

    POINTCLOUD_DATA *data = (POINTCLOUD_DATA*) malloc(sizeof(POINTCLOUD_DATA));
    data->pointCloud = pointcloud;
    data->octree = octree;

    object->data = data;
    object->rayHit = rayHit;
    object->normal = getNormal;
    object->color = getColor;
    object->setColor = setColor;
    
    object->ka = 0.1;
    object->kd = 0.80;
    object->ks = 0.1;
    object->kt = 0.0;
    object->ni = 1.52;
    object->phong = 6.0;

    return object;
}

// destroyPointCloudObject
// We'll leak memory for now
