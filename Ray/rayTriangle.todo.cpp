#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rayTriangle.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

void RayTriangle::initialize(void) {
    v1 = v[1]->position - v[2]->position;
    v2 = v[0]->position - v[2]->position;
    Point3D n = v1.crossProduct(v2).unit();
    plane.normal = v[0]->normal;

    double d = -1 * ((n).dot(v[2]->position));

    plane.distance = d;

}

double RayTriangle::intersect(Ray3D ray, RayIntersectionInfo& iInfo, double mx) {

    initialize();

    Point3D N = plane.normal.unit();
    Point3D t0 = v[0]->position;
    Point3D t1 = v[1]->position;
    Point3D t2 = v[2]->position;
    double d = plane.distance;
    Point3D V = ray.direction.unit();
    Point3D P0 = ray.position;

    if (V.dot(N) == 0) return -1;

    double t = -(P0.dot(N) + d) / (V.dot(N));

    Point3D P = P0 + V * t;

    Point3D w0 = t0 - P0;
    Point3D w1 = t1 - P0;
    Point3D w2 = t2 - P0;

    Point3D N1 = w1.crossProduct(w0);
    if ((P - P0).dot(N1) < 0) return -1;

    Point3D N2 = w0.crossProduct(w2);
    if ((P - P0).dot(N2) < 0) return -1;

    Point3D N3 = w2.crossProduct(w1);
    if ((P - P0).dot(N3) < 0) return -1;

    if (t <= 0) return -1;
    
    //TEXTURE
    double areaT = ((t1 - t0).crossProduct(t2 - t0)).length() / 2;
    double area1 = ((P - t0).crossProduct(t2 - t0)).length() / 2;
    double area2 = ((P - t0).crossProduct(t1 - t0)).length() / 2;
    //double area3 = ((t1 - t0).crossProduct(t2 - t0)).length() / 2;  
    double u = area1 / areaT;
    double v = area2 / areaT;    
    double w = 1 - u - v;    
    iInfo.texCoordinate = Point2D(u, v);  

    iInfo.iCoordinate = P;
    iInfo.material = material;
    iInfo.normal = N;
    return t;

}

BoundingBox3D RayTriangle::setBoundingBox(void) {

//    bBox.p[0] = v[0];
//    bBox.p[1] = v[0];
//
//    for (int i = 0; i < 3; ++i) {
//        if (v[i] < bBox.p[0]) bBox.p[0] = v[i];
//        if (v[i] > bBox.p[1]) bBox.p[1] = v[i];
//    }
//        return bBox;

    Point3D pList[] = {v[0]->position, v[1]->position, v[2]->position};
    
    //pList->p[0] = v[0];
    
    return BoundingBox3D(pList, 3);

}

//////////////////
// OpenGL stuff //
//////////////////

int RayTriangle::drawOpenGL(int materialIndex) {
    return -1;
}
