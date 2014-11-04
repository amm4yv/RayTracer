#include <math.h>
#ifdef __APPLE__
	#include <GLUT/glut.h>
#else
	#include <GL/glut.h>
#endif
#include "rayScene.h"
#include "raySpotLight.h"


////////////////////////
//  Ray-tracing stuff //
////////////////////////
Point3D RaySpotLight::getDiffuse(Point3D cameraPosition,RayIntersectionInfo& iInfo){
    
    double rI_l, gI_l, bI_l;
    
    RayMaterial * m = iInfo.material;
    //Surface property
    Point3D K_d = m->diffuse;
        
    //Normal vector to surface
    Point3D N = iInfo.normal;
   
    //Color of the light source
    Point3D I0 = color;
    
    /** The position of the spot-light */
    Point3D L = location;
    double d = L.length();

    /** The preferred direction the outgoing light rays */
    Point3D D = direction;
    double Dmag = sqrt(D.dot(D));
    double Lmag = sqrt(L.dot(L));

    /** The constant term of the attenuation equation */
    double kc = constAtten;
    /** The linear term of the attenuation equation */
    double kl = linearAtten;
    /** The quadratic term of the attenuation equation */
    double kq = quadAtten;

    /** The cut-off angle for the spot light (should be in the range [0,Pi/2]) */
    double gamma = cutOffAngle;
    /** The rate at which the intensity falls off as light travels in the non-preferred direction (should be in the range [0,128]) */
    double alpha = dropOffRate;
    
    double DL = acos((D.dot(L)) / (Dmag * Lmag));
    
    if (DL <= cos(gamma)){
        rI_l = 0;
        gI_l = 0;
        bI_l = 0;
    }
    
    else{        
        rI_l = (I0.p[0] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
        gI_l = (I0.p[1] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
        bI_l = (I0.p[2] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
    }

    double rI_d = fmax(0, -1 * (fmin(1, K_d.p[0] * N.dot(L) * rI_l)));
    double gI_d = fmax(0, -1 * (fmin(1, K_d.p[1] * N.dot(L) * gI_l)));
    double bI_d = fmax(0, -1 * (fmin(1, K_d.p[2] * N.dot(L) * bI_l)));
    
    return Point3D(rI_d, gI_d, bI_d);
}


Point3D RaySpotLight::getSpecular(Point3D cameraPosition,RayIntersectionInfo& iInfo){
    
    double rI_l, gI_l, bI_l;
      
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
    Point3D locat = location;
    double d = locat.length();

    //Vector between light and surface
    Point3D L = locat - iInfo.iCoordinate;
    //Reflected vector
    Point3D R = ((N * 2) * (N.dot(L)) - L).unit() ;
    //View vector
    Point3D V = (cameraPosition - iInfo.iCoordinate).unit();
   
    /** The preferred direction the outgoing light rays */
    Point3D D = direction;
    double Dmag = sqrt(D.dot(D));
    double Lmag = sqrt(L.dot(L));

    /** The constant term of the attenuation equation */
    double kc = constAtten;
    /** The linear term of the attenuation equation */
    double kl = linearAtten;
    /** The quadratic term of the attenuation equation */
    double kq = quadAtten;

    /** The cut-off angle for the spot light (should be in the range [0,Pi/2]) */
    double gamma = cutOffAngle;
    /** The rate at which the intensity falls off as light travels in the 
     * non-preferred direction (should be in the range [0,128]) */
    double alpha = dropOffRate;
    
    double DL = acos((D.dot(L)) / (Dmag * Lmag));
    
    if (DL <= cos(gamma)){
        rI_l = 0;
        gI_l = 0;
        bI_l = 0;
    }
    
    else{        
        rI_l = (I0.p[0] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
        gI_l = (I0.p[1] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
        bI_l = (I0.p[2] * (pow(DL, alpha))) / (kc + kl*d + kq*d*d);
    }
       
    double rI_s = fmax(0, -1 * (fmin(1, K_s.p[0] * (pow(V.dot(R), n)) * rI_l)));
    double gI_s = fmax(0, -1 * (fmin(1, K_s.p[1] * (pow(V.dot(R), n)) * gI_l)));
    double bI_s = fmax(0, -1 * (fmin(1, K_s.p[2] * (pow(V.dot(R), n)) * bI_l)));
    
    return Point3D(rI_s, gI_s, bI_s);
}
int RaySpotLight::isInShadow(RayIntersectionInfo& iInfo,RayShape* shape,int& isectCount){
    
    Point3D P0 = iInfo.iCoordinate;
    Point3D P1 = location;
    
    if (shape->intersect(Ray3D(P0, (P1 - P0).unit()), iInfo, 0) != -1) return 0;
    
    return 1;
}


Point3D RaySpotLight::transparency(RayIntersectionInfo& iInfo,RayShape* shape,Point3D cLimit){
       
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
    Point3D locat = location;
    double d = locat.length();
    
    return Point3D(1,1,1);
}

//////////////////
// OpenGL stuff //
//////////////////
void RaySpotLight::drawOpenGL(int index){
}
