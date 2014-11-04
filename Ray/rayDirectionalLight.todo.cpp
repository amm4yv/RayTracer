#include <math.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "rayDirectionalLight.h"
#include "rayScene.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

Point3D RayDirectionalLight::getDiffuse(Point3D cameraPosition, RayIntersectionInfo& iInfo) {

    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_d = m->diffuse;

    //Normal vector to surface
    Point3D N = iInfo.normal.unit();

    //Color of the light source
    Point3D I0 = color;

    /** The direction the outgoing light rays */
    Point3D L = direction.negate().unit();
    
    if (N.dot(L) < 0) return Point3D(0, 0, 0);

    double rI_d = fmax(0, 1 * (fmin(1, K_d.p[0] * N.dot(L) * I0.p[0])));
    double gI_d = fmax(0, 1 * (fmin(1, K_d.p[1] * N.dot(L) * I0.p[1])));
    double bI_d = fmax(0, 1 * (fmin(1, K_d.p[2] * N.dot(L) * I0.p[2])));
    

    return Point3D(rI_d, gI_d, bI_d);

}

Point3D RayDirectionalLight::getSpecular(Point3D cameraPosition, RayIntersectionInfo& iInfo) {
    
    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_s = m->specular.unit();   
    
    //Normal vector to surface
    Point3D N = iInfo.normal.unit();
    //Shine
    double n = m->specularFallOff;
    
    //Color of the light source
    Point3D I0 = color;
    
    /** The direction the outgoing light rays */
    Point3D L = direction.negate().unit();

    //Reflected vector
    Point3D R = ((N * 2) * (N.dot(L)) - L).unit();
    //View vector
    Point3D V = (cameraPosition - iInfo.iCoordinate).unit();
    
    if (V.dot(R) < 0) return Point3D(0, 0, 0);

    double rI_s = fmax(0, (fmin(1, K_s.p[0] * (pow(V.dot(R), n)) * I0.p[0])));
    double gI_s = fmax(0, (fmin(1, K_s.p[1] * (pow(V.dot(R), n)) * I0.p[1])));
    double bI_s = fmax(0, (fmin(1, K_s.p[2] * (pow(V.dot(R), n)) * I0.p[2])));
    
    return Point3D(rI_s, gI_s, bI_s);
}

int RayDirectionalLight::isInShadow(RayIntersectionInfo& iInfo, RayShape* shape, int& isectCount) {   
    Point3D P0 = iInfo.iCoordinate;
    Point3D V = (direction).negate().unit();
    RayIntersectionInfo temp = iInfo;
    if (shape->intersect(Ray3D(P0, V), temp, -1) <= 0){
        //printf("%f\n", shape->intersect(Ray3D(P0, V), temp, -1));
        return 1;    
    }
    return 0;
}

/*This virtual method tests if the intersection point represented by iInfo 
 * is in partial shadow from the light source. A ray is cast from the hit 
 * location to the light source, and the transparency values are accumulated. 
 * If the transparency value falls below cLimit, the testing terminates. 
 * In computing the transparency value, isectCount is incremented according 
 * to the number of intersection tests performed and rayCount is incremented 
 * according to the number of rays that need to be cast. The returned value 
 * is a 3D point whose coefficients should all be in the range [0,1].*/
Point3D RayDirectionalLight::transparency(RayIntersectionInfo& iInfo, RayShape* shape, Point3D cLimit) {
    
    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_s = m->specular.unit();   
    
    //Normal vector to surface
    Point3D N = iInfo.normal.unit();
    //Shine
    double n = m->specularFallOff;
    
    //Color of the light source
    Point3D I0 = color;
    
    /** The direction to the light rays */
    Point3D L = direction.negate().unit();
    
    
    
    return Point3D(1, 1, 1);
}

//////////////////
// OpenGL stuff //
//////////////////

void RayDirectionalLight::drawOpenGL(int index) {
}
