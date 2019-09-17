#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include "raytrace.h"

OBJECT *createPointCloudObject(const char *plyfile, int maxLevel, float sigma);

#endif /* POINTCLOUD_H */