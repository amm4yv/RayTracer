#include <math.h>
#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
#include "rayScene.h"
#include "rayBox.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayBox::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
    Point3D pos = ray.position;
    Point3D dir = ray.direction;
//    double tmin = (min.x - pos.p[0]) / dir.p[0];
//    double tmax = (max.x - pos.p[0]) / dir.p[0];
//    if (tmin > tmax){
//        double temp = tmin;
//        tmin = tmax;
//        tmax = temp;
//    }
//    
//    double tymin = (min.y - pos.p[1]) / dir.p[1];
//    double tymax = (max.y - pos.p[1]) / dir.p[1];
//    if (tymin > tymax){
//        double temp = tymin;
//        tymin = tymax;
//        tymax = temp;
//    }
//    if ((tmin > tymax) || (tymin > tmax))
//        return -1;
//    if (tymin > tmin) tmin = tymin;
//    if (tymax < tmax) tmax = tymax;
//    
//    double tzmin = (min.z - pos.p[2]) / dir.p[1];
//    double tzmax = (max.z - pos.p[2]) / dir.p[1];
//    if (tzmin > tzmax){
//        double temp = tzmin;
//        tzmin = tzmax;
//        tzmax = temp;
//    }
//    if ((tmin > tzmax) || (tzmin > tmax))
//        return -1;
//    
//    if (tzmin > tmin)
//        tmin = tzmin;
//    if (tzmax < tmax)
//        tmax = tzmax;
//    if ((tmin > r.tmax) || (tmax < r.tmin)) return -1;
//    if (r.tmin < tmin) r.tmin = tmin;
//    if (r.tmax > tmax) r.tmax = tmax;
//    return true;
//    
	return -1;
}
BoundingBox3D RayBox::setBoundingBox(void){
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayBox::drawOpenGL(int materialIndex){
	return -1;
}
