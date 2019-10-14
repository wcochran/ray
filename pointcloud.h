#ifndef POINTCLOUD_H
#define POINTCLOUD_H

#include "raytrace.h"

#define THICKNESS_SIGNED_DISTANCE_FUNCTION

OBJECT *createPointCloudObject(const char *plyfile, int maxLevel, 
 #ifdef THICKNESS_SIGNED_DISTANCE_FUNCTION
        float halfThickness,
#endif   
        float sigma);

#endif /* POINTCLOUD_H */