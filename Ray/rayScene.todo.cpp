#include "rayScene.h"
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include <math.h>


///////////////////////
// Ray-tracing stuff //
///////////////////////

Point3D RayScene::Reflect(Point3D v, Point3D n) {
    return v - ((n) * (v.dot(n)) * 2);
}

/*if computing the refraction direction would require evaluating the arcsin 
 * of a number with magnitude larger than 1 so that the angle of the ray is 
 * greater than the critical angle
 */
int RayScene::Refract(Point3D v, Point3D n, double ir, Point3D& refract) {
    //    double d = v.dot(n);    
    //    double cos = 1 - ir * (1 - d * d);
    //    
    //    if (cos < 0) return 0;
    //    
    //    refract = (v * ir) + n * ((v.dot(n)) * ir - cos);
    //    return ir;
    
//    float rindex = prim->GetMaterial()->GetRefrIndex();
//	float n = a_RIndex / rindex;
//	vector3 N = prim->GetNormal( pi ) * (float)result;
//	float cosI = -DOT( N, a_Ray.GetDirection() );
//	float cosT2 = 1.0f - n * n * (1.0f - cosI * cosI);
//	if (cosT2 > 0.0f)
//	{
//		vector3 T = (n * a_Ray.GetDirection()) + (n * cosI - sqrtf( cosT2 )) * N;

    double cosI = -1 * n.dot(v);
    double cosT2 = 1.0 - (ir * ir * (1.0 - (cosI * cosI)));
    if (cosT2 > 0){
        refract = (v * ir) + (n * (ir * cosI - sqrt(cosT2)));
        return ir;
    }
    
    return 0;
    
//    double c1 = -1 * v.dot(n);
//    double cs2 = 1 - ir * ir * (1 - c1 * c1);
//    if (cs2 < 0)
//        return 0;
//    refract = (v * ir) + (n * (ir * c1 - sqrt(cs2)));
//    return ir;
}

/*Generate rays from the camera's position through (i,j)-th pixel of a width x height view plane.
 */
Ray3D RayScene::GetRay(RayCamera* camera, int i, int j, int width, int height) {

    Point3D p = camera->position;
    Point3D V = camera->direction;
    Point3D up = camera->up;
    Point3D right = camera->right;

    double t = (height / 2) / tan(camera->heightAngle / 2);
    Point3D p1 = p + (V * t) + up * (j - (height / 2)) + right * (i - (width / 2));

    return Ray3D(p, p1.unit());
}

Point3D RayScene::GetColor(Ray3D ray, int rDepth, Point3D cLimit) {

    if (rDepth == 0 || cLimit.length() > 1) return Point3D(0, 0, 0);

    Point3D color = Point3D(0, 0, 0);

    RayIntersectionInfo info;
    RayIntersectionInfo info2;
    int isectCount = 0;

    double mx = group->intersect(ray, info, -1);

    //printf("%d ", group->sNum);
    
    if (mx != -1) {

        color = info.material->emissive + ambient * info.material->ambient;

        for (int i = 0; i < lightNum; i++) {

            int shadow = 1;
            for (int j = 0; j < group->sNum; j++) {
                info2 = info;
                shadow = lights[i]->isInShadow(info2, group->shapes[j], isectCount);
                if (shadow == 0) break;
            }
            color += (lights[i]->getDiffuse(camera->position, info)) * shadow;
            color += (lights[i]->getSpecular(camera->position, info)) * shadow;

        }

        if ((info.material->specular.length() != 0)) {
            // compute reflection
            Point3D reflect = (Reflect(ray.direction.unit(), info.normal.unit())).unit();
            // recurse
            Point3D reflectionColor;
            reflectionColor += (GetColor(Ray3D(info.iCoordinate, reflect), rDepth - 1, cLimit * info.material->specular.length())) / info.material->specular.length();

            color += reflectionColor;
        }

        double t = info.material->transparent.length();

        int ir = info.material->refind;
                if (t != 0){
                    Point3D refract;
                    
                    Point3D N;
                    if (mx > 0) N = info.normal.unit();
                    else N = info.normal.unit().negate();
                    
                    int ir = Refract(ray.direction.unit(), N, info.material->refind, refract);
                    
                    if (ir != 0){
                        Point3D refractColor;
                        info.material->refind = 1/ir;
                        t = info.material->transparent.length();
                        refractColor += GetColor(Ray3D(info.iCoordinate, refract), rDepth - 1, cLimit / t) / info.material->specular.length();
                        color += refractColor;
                    }
               
                }


        return color;

    }
    //reflect not on background
    if (rDepth < 5)
        return Point3D(0, 0, 0);
    
    return background;
}

//////////////////
// OpenGL stuff //
//////////////////

void RayMaterial::drawOpenGL(void) {
}

void RayTexture::setUpOpenGL(void) {
}
