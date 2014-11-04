#include <math.h>
#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
#include "rayScene.h"
#include "rayCylinder.h"


////////////////////////
//  Ray-tracing stuff //
////////////////////////
double RayCylinder::intersect(Ray3D ray,RayIntersectionInfo& iInfo,double mx){
    
    double a, b, c;
    
    //The ends lie at y=cy-h/2 and y=cy+h/2.
    
    Point3D P0 = ray.position;
    Point3D V = ray.direction;
    
    /* The center of the cylinder */
    Point3D O = center;
	/** The hieght of the cylinder */

    double h = height;
	/** The radius of the cylinder */

    double r = radius;
        
    /*a=xD2+yD2, b=2xExD+2yEyD, and c=xE2+yE2-1.*/

//    double a = (P0.p[0] * P0.p[0]) + (P0.p[1] * P0.p[1]);
//    double b = 2 * V.p[0] * P0.p[0] + 2 * V.p[1] * P0.p[1];
//    double c = (V.p[0] * V.p[0]) + (V.p[1] * V.p[1]) - (r * r);
        
    
    //P is the start of the ray and V is the direction. Substituting this formua for X in the implicit equation yields:
    //Here, X is a point on the sphere or circle, C is the center, and r is the radius. The equation for the ray is:
    //(V.V)*t^2 + 2*(P.V - C.V)*t + P.P + C.C - 2*P.C - r^2 = 0
    
    P0.p[1] = O.p[1];
    V.p[1] = O.p[1];
    
    a = (V.dot(V));
    b = 2 * (P0.dot(V) - O.dot(V));
    c = P0.dot(P0) + O.dot(O) - 2 * (P0.dot(O)) - (r * r);
       
    double dis = b * b - 4 * a * c;
    printf("Cylinder success\n");

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
        
    /* (p - pa + vt - (va
. (p - pa + vt))va)2 - r2 = 0
     */
    
    
	return -1;
}

BoundingBox3D RayCylinder::setBoundingBox(void){
	return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////
int RayCylinder::drawOpenGL(int materialIndex){
	return -1;
}
