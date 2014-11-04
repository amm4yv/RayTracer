#include <math.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "rayScene.h"
#include "raySphere.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

double RaySphere::intersect(Ray3D ray, RayIntersectionInfo& iInfo, double mx) {

    Point3D P0 = ray.position;
    Point3D V = ray.direction.unit();
    Point3D O = center;

    double b = (V * 2).dot(P0 - O);
    double c = ((P0 - O).length() * (P0 - O).length()) - (radius * radius);
    double dis = b * b - 4 * c;

    if (dis < 0) return -1;

    double t;
    if (-b - sqrt(dis) <= -b + sqrt(dis)) t = (-b - sqrt(dis)) / 2;
    else t = (-b + sqrt(dis)) / 2;
    Point3D P = P0 + V * t;

    if (t <= 0) return -1;
    
    // TEXTURE  
    Point3D d = (O - P).unit();  
    double u = 0.5 + (atan2(d.p[2], d.p[0]) / 2*PI);
    double v = 0.5 - (asin(d.p[1]) / PI);   
    iInfo.texCoordinate = Point2D(u, v);   

    iInfo.iCoordinate = P;
    iInfo.normal = (P - O).unit();
    iInfo.material = material;
    return t;
}

BoundingBox3D RaySphere::setBoundingBox(void) {
    
    Point3D O = center;
    double r = radius;
    
//    bBox.p[0] = O + Point3D(-r, -r, -r);
//    bBox.p[1] = O + Point3D(r, -r, -r);
//    bBox.p[2] = O + Point3D(-r, -r, r);
//    bBox.p[3] = O + Point3D(r, -r, r);
//    bBox.p[4] = O + Point3D(-r, r, -r);
//    bBox.p[5] = O + Point3D(r, r, -r);
//    bBox.p[6] = O + Point3D(-r, r, r);
//    bBox.p[7] = O + Point3D(r, r, r);
    
//    bBox.p[0] = O + Point3D(-r, -r, -r);
//    bBox.p[1] = O + Point3D(r, r, r);
//       
//    return bBox;
    
    bBox = BoundingBox3D(O + Point3D(-r, -r, -r),O + Point3D(r, r, r));
    return bBox;
}

//////////////////
// OpenGL stuff //
//////////////////

int RaySphere::drawOpenGL(int materialIndex) {
    return -1;
}
