#ifndef TEAPOT_H
#define TEAPOT_H

#include "raytrace.h"

OBJECT *createTeapotObject(double bottom[3], 
			   double lidDir[3], double spoutDir[3],
			   double scale);

#endif
