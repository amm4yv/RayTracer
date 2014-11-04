#include <math.h>
#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
#include "rayScene.h"
#include "rayCone.h"


////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayCone::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
    
    double a, b, c;
    Point3D P0 = ray.position;
    Point3D V = ray.direction;
    /** The center of the cylinder */
    Point3D O = center;
	/** The height of the cylinder */
    double h = height;
	/** The radius of the cylinder */
    double r = radius;
    // a=xD2+yD2-zD2, b=2xExD+2yEyD-2zEzD, and c=xE2+yE2-zE2.
    
    a = (V.p[0] * V.p[0]) + (V.p[1] * V.p[1]) - (V.p[2] * V.p[2]);
    b = 2 * (P0.p[0] * V.p[0] + P0.p[1] * V.p[1] - P0.p[2] * V.p[2]);
    c = (P0.p[0] * P0.p[0]) + (P0.p[1] * P0.p[1]) - (P0.p[2] * P0.p[2]);
    
    double dis = b * b - 4 * a * c;
    //printf("Cylinder success\n");

    if (dis < 0) return -1;

    double t;
    if (-b - sqrt(dis) <= -b + sqrt(dis)) t = (-b - sqrt(dis)) / 2;
    else t = (-b + sqrt(dis)) / 2;
    Point3D P = P0 + V * t;

    if (t <= 0) return -1;

    iInfo.iCoordinate = P;
    iInfo.normal = (P - O).unit();
    iInfo.material = material;
    return t;
    
}

BoundingBox3D RayCone::setBoundingBox(void){
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayCone::drawOpenGL(int materialIndex){
	return -1;
}
