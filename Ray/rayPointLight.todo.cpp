#include <math.h>
#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
#include "rayPointLight.h"
#include "rayScene.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

/*This virtual method returns the diffuse contribution of the light source 
 * to the specified hit location. It computes the amount of diffuse light 
 * reaching the hit location and scales that using the material properties 
 * of the hit location. The returned value is a 3D point whose coefficients 
 * should all be in the range [0,1].*/
Point3D RayPointLight::getDiffuse(Point3D cameraPosition,RayIntersectionInfo& iInfo){
    
    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_d = m->diffuse;
        
    //Normal vector to surface
    Point3D N = iInfo.normal;
   
    //Color of the light source
    Point3D I0 = color;
    
    /** The position of the spot-light */
    //Point3D locat = location;
    
    //Vector between light and surface
    Point3D L = location - iInfo.iCoordinate;
    
    /** The constant term of the attenuation equation */
    double kc = constAtten;
    /** The linear term of the attenuation equation */
    double kl = linearAtten;
    /** The quadratic term of the attenuation equation */
    double kq = quadAtten;
    
    double d = L.length();

    double rI_l = I0.p[0] / (kc + kl*d + kq * d * d);
    double gI_l = I0.p[1] / (kc + kl*d + kq * d * d);   
    double bI_l = I0.p[2] / (kc + kl*d + kq * d * d);    
    
    double rI_d = fmax(0, -1 * (fmin(1, K_d.p[0] * N.dot(L) * rI_l)));
    double gI_d = fmax(0, -1 * (fmin(1, K_d.p[1] * N.dot(L) * gI_l)));
    double bI_d = fmax(0, -1 * (fmin(1, K_d.p[2] * N.dot(L) * bI_l)));
    
    
    return Point3D(rI_d, gI_d, bI_d);
}

/*This virtual method returns the specular contribution of the light source 
 * to the specified hit location. It computes the amount of diffuse light 
 * reaching the hit location, using the normal direction at the hit location 
 * (in iInfo), and scales that using the material properties of the hit location. 
 * The returned value is a 3D point whose coefficients should all be in the range [0,1].*/
Point3D RayPointLight::getSpecular(Point3D cameraPosition,RayIntersectionInfo& iInfo){
      
    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_s = m->specular;
        
    //Normal vector to surface
    Point3D N = iInfo.normal;
    //Shine
    double n = m->specularFallOff;
        
    //Color of the light source
    Point3D I0 = color;
    
    /** The position of the spot-light */
    //Point3D locat = location;
    
    //Vector between light and surface
    Point3D L = location - iInfo.iCoordinate;
    //Reflected vector
    Point3D R = ((N * 2) * (N.dot(L)) - L).unit() ;
    //View vector
    Point3D V = (cameraPosition - iInfo.iCoordinate).unit() ;
    
    /** The constant term of the attenuation equation */
    double kc = constAtten;
    /** The linear term of the attenuation equation */
    double kl = linearAtten;
    /** The quadratic term of the attenuation equation */
    double kq = quadAtten;
    
    double d = L.length();    
    
    double rI_l = I0.p[0] / (kc + kl*d + kq * d * d);
    double gI_l = I0.p[1] / (kc + kl*d + kq * d * d);   
    double bI_l = I0.p[2] / (kc + kl*d + kq * d * d);    
    
    double rI_s = fmax(0, -1 * (fmin(1, K_s.p[0] * (pow(V.dot(R), n)) * rI_l)));
    double gI_s = fmax(0, -1 * (fmin(1, K_s.p[1] * (pow(V.dot(R), n)) * gI_l)));
    double bI_s = fmax(0, -1 * (fmin(1, K_s.p[2] * (pow(V.dot(R), n)) * bI_l)));
      
    return Point3D(rI_s, gI_s, bI_s);
    
 }

/*This virtual method tests if the intersection point represented by iInfo is 
 * in shadow from the light source. The value of isectCount is incremented 
 * according to the number of intersection tests performed in order to determine 
 * if the point of intersection is in shadow. The returned value is either 0 if 
 * the the intersection point is not in shadow or 1 if it is.*/
int RayPointLight::isInShadow(RayIntersectionInfo& iInfo,RayShape* shape,int& isectCount){
    
    Point3D P0 = iInfo.iCoordinate;
    Point3D P1 = location; 
    
    if (shape->intersect(Ray3D(P0, (P1 - P0).unit()), iInfo, 0) != -1) return 0;
    
    return 1;
}

/*This virtual method tests if the intersection point represented by iInfo is 
 * in partial shadow from the light source. A ray is cast from the hit location 
 * to the light source, and the transparency values are accumulated. If the 
 * transparency value falls below cLimit, the testing terminates. In computing 
 * the transparency value, isectCount is incremented according to the number of 
 * intersection tests performed and rayCount is incremented according to the 
 * number of rays that need to be cast. The returned value is a 3D point whose 
 * coefficients should all be in the range [0,1].*/
Point3D RayPointLight::transparency(RayIntersectionInfo& iInfo,RayShape* shape,Point3D cLimit){

    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_d = m->diffuse;
        
    //Normal vector to surface
    Point3D N = iInfo.normal;
   
    //Color of the light source
    Point3D I0 = color;
    
    //Vector between light and surface
    Point3D L = location - iInfo.iCoordinate;
    
    return Point3D(1,1,1);
}


//////////////////
// OpenGL stuff //
//////////////////
void RayPointLight::drawOpenGL(int index){
}
