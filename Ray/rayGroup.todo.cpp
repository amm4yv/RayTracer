#include <stdlib.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "rayGroup.h"

////////////////////////
//  Ray-tracing stuff //
////////////////////////

double RayGroup::intersect(Ray3D ray, RayIntersectionInfo& iInfo, double mx) {

    Matrix4D M = getMatrix();
    Matrix4D N = getNormalMatrix();
    Matrix4D I = getInverseMatrix();

    Point3D c = ray.position;
    RayIntersectionInfo temp = iInfo;

    Ray3D tray = I.mult(ray);
    double dt, qL, mxW;

    //printf("%d ", sNum);

    
    for (int i = 0; i < sNum; i++) {

        
        /*Modify this method to only test for an intersection or a ray 
         * with a child RayShape if the ray intersects the bounding box 
         * of the child. Keep in mind that if the value of the mx 
         * parameter is greater than 0 and the distance to the intersection 
         * of the ray with a child's bounding box is great than mx, you do 
         * not need to test for intersection with that child.*/

        //dt = shapes[i]->bBox.intersect(ray);

        //&& (mx > 0 && dt > mx)
        //if (dt != -1) {

            //printf("Intersect w/ bbox\n");
            
            //dt = shapes[i]->intersect(ray, temp, mx);

            //TRANSFORMATION      
            dt = shapes[i]->intersect(tray, temp, mx);
            temp.iCoordinate = M.multPosition(temp.iCoordinate);
            temp.normal = N.multDirection(temp.normal).unit();
            if (dt != -1) dt = (temp.iCoordinate - c).length();

            if ((mx < 0 || dt < mx) && (dt > 0)) {
                mx = dt;
                iInfo.iCoordinate = temp.iCoordinate;
                iInfo.material = temp.material;
                iInfo.normal = temp.normal;
                return dt;
            }
        //}
    }

    return -1;
}

/*To implement the RayShape::getBoundingBox method for For a RayGroup, 
 * accumulate the bounding boxes of all the child RayShapes, compute 
 * the bounding box of the transformed accumulation of bounding boxes, 
 * store and return it.
 * (Note: When the parser is done reading the .ray file it automatically 
 * calls the RayShape::setBoundingBox method for the root node, so if you 
 * have implemented this method for all of the subclasses of RayShape, 
 * the bounding boxes are already in place to be used for intersection 
 * queries, and you do not have to reset them.)*/

BoundingBox3D RayGroup::setBoundingBox(void) {



    for (int i = 0; i < sNum; i++) {

        bBox += shapes[i]->bBox;

    }

    bBox.transform(getMatrix());

    return bBox;
}

int StaticRayGroup::set(void) {

    inverseTransform = getMatrix().invert();
    normalTransform = (getMatrix().transpose()).invert();

    return 1;
}
//////////////////
// OpenGL stuff //
//////////////////

int RayGroup::getOpenGLCallList(void) {
    return 0;
}

int RayGroup::drawOpenGL(int materialIndex) {
    return -1;
}

/////////////////////
// Animation Stuff //
/////////////////////

Matrix4D ParametrizedEulerAnglesAndTranslation::getMatrix(void) {
    return Matrix4D::IdentityMatrix();
}

Matrix4D ParametrizedClosestRotationAndTranslation::getMatrix(void) {
    return Matrix4D::IdentityMatrix();
}

Matrix4D ParametrizedRotationLogarithmAndTranslation::getMatrix(void) {
    return Matrix4D::IdentityMatrix();
}

Matrix4D ParametrizedQuaternionAndTranslation::getMatrix(void) {
    return Matrix4D::IdentityMatrix();
}
