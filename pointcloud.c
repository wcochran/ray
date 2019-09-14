#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
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
    hashTable->table = (GridData*) malloc(size*sizeof(GridData*));
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
//   node = NULL => empty child
//   node != NULL and node->child == NULL => filled leaf
//   node != NULL and node->child != NULL => internal node with 
//              8 children (some of which may be NULL)
//  The actually data for a filled leaf is stored in a GridDataHashTable
//  which is indexed at its eight corners.
//
typedef struct OctreeNode {
    struct OctreeNode **child;  // NULL => filled leaf
} OctreeNode;                   // not NULL => internal node

typedef struct {
    int maxLevel;
    BBOX bbox;
    IndexBox indexBox;
    GridHashTable *GridHashTable;
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
            bbox->min.z - r <= vert->pos.z && vert->pos.z < bbox->max.y + r) {
            if (vertArray == NULL) {
                vertArray = createVertArray(0, 10);
            }
            addVert(vertArray, &inputVertArray->verts[i]);
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
            weightColor[0] *= w*vert->color.r;
            weightColor[1] *= w*vert->color.g;
            weightColor[2] *= w*vert->color.b;
            n++;
        }
    }
    const double scale = weight > 0 ? 1/weight : 1;
    POINT3 Pave = {scale*weightP.x, scale*weightP.y, scale*weightP.z};
    N->x = scale*weightN.x;
    N->y = scale*weightN.y;
    N->z = scale*weightN.z;
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

static
void subdivide(OctreeNode *node, int level, int maxLevel,
        BBOX *bbox, IndexBox *iBox, GridHashTable *hashTable,
        VertArray *verts, float radius, float sigma) {
    Vec3 C;
    C.x = (bbox->min.x + bbox->max.x)/2;
    C.y = (bbox->min.y + bbox->max.y)/2;
    C.z = (bbox->min.z + bbox->max.z)/2;
    double X[3] = {bbox->min.x, C.x, bbox->max.x};
    double Y[3] = {bbox->min.y, C.y, bbox->max.y};
    double Z[3] = {bbox->min.z, C.z, bbox->max.z};

    Index centerIndex;
    centerIndex.i = (iBox->min.i + iBox->max.i)/2;
    centerIndex.j = (iBox->min.j + iBox->max.j)/2;
    centerIndex.k = (iBox->min.k + iBox->max.k)/2;
    int I[3] = {iBox->min.i, centerIndex.i, iBox->max.i};
    int J[3] = {iBox->min.j, centerIndex.j, iBox->max.j};
    int K[3] = {iBox->min.k, centerIndex.k, iBox->max.k};

    int n = 0;

    node->isFilled = true; // don't care
    node->child = (OctreeNode **) malloc(8*sizeof(OctreeNode*));

    BBOX childBBox;
    IndexBox childIndexBox;
    
    for (int k = 0; k < 2; k++) {
        childBBox.min.z = Z[k];
        childBBox.max.z = Z[k+1];
        childIndexBox.min.k = K[k];
        childIndexBox.max.k = K[k+1];
        for (int j = 0; j < 2; j++) {
            childBBox.min.y = Y[k];
            childBBox.max.y = Y[k+1];
            childIndexBox.min.j = J[j];
            childIndexBox.max.j = J[j+1];
            for (int i = 0; i < 2; i++) {
                childBBox.min.x = X[i];
                childBBox.max.x = X[i+1];
                childIndexBox.min.i = I[i];
                childIndexBox.max.i = I[i+1];

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