#include <stdlib.h>
#include <math.h>

#include <SVD/SVDFit.h>
#include <SVD/MatrixMNTC.h>

#include "geometry.h"


///////////////////////
// Ray-tracing stuff //
///////////////////////

/*Return the distance along the ray to the nearest point of 
 * intersection with the interior of the bounding box. 
 * (If the ray starts off inside the bounding box, then 
 * a value of 0 should be returned.) If the ray does not 
 * intersect the bounding box, a negative value should be returned.*/

double BoundingBox3D::intersect(const Ray3D& ray) const {

    double tnear, tfar, temp;


    Point3D O = ray.position;
    Point3D V = ray.direction;
    Point3D p0 = p[0];
    Point3D p1 = p[1];

    if ((O.p[2] >= p1.p[2] && O.p[1] <= p1.p[1] && O.p[0] <= p1.p[0]) ||
            (O.p[2] <= p0.p[2] && O.p[1] >= p0.p[1] && O.p[0] >= p0.p[0])) {

        printf("Ray inside\n");

        return 0;
    }

    double t0x = (p0.p[0] - O.p[0]) / V.p[0];
    double t1x = (p1.p[0] - O.p[0]) / V.p[0];
    double t0y = (p0.p[1] - O.p[1]) / V.p[1];
    double t1y = (p1.p[1] - O.p[1]) / V.p[1];
    double t0z = (p0.p[2] - O.p[2]) / V.p[2];
    double t1z = (p1.p[2] - O.p[2]) / V.p[2];

    if (t0x > t1x) {
        temp = t0x;
        t0x = t1x;
        t1x = temp;
    }
    if (t0x > tnear) {
        tnear = t1x;
    }
    tfar = t1x;
    if (tnear > tfar) {
        return -1;
    }
    if (tfar < 0) return -1;

    if (t0y > t1y) {
        temp = t0y;
        t0y = t1y;
        t1y = temp;
    }
    if (t0y > tnear) {
        tnear = t1y;
    }
    if (t1y < tfar) {
        tfar = t1y;
    }
    if (tnear > tfar) {
        return -1;
    }
    if (tfar < 0) return -1;

    if (t0z > t1z) {
        temp = t0z;
        t0z = t1z;
        t1z = temp;
    }
    if (t0z > tnear) {
        tnear = t1z;
    }
    if (t1z < tfar) {
        tfar = t1z;
    }
    if (tnear > tfar) {
        return -1;
    }
    if (tfar < 0) return -1;


    //printf("Passed tests: %f\n", tnear);
    
    return tnear;


    //    
    //    if (t0x > t0y) {
    //        tmin = t0x;
    //    } else tmin = t0y;
    //
    //    if (t1x < t1y) {
    //        tmax = t1x;
    //    } else tmax = t1y;
    //
    //    if (t0x > t1y || t0y > t1x){
    //        //printf("cond1\n");
    //        return -1;
    //    }
    //    
    //    
    //    if (tmin > t1z || t0z > tmax){
    //        printf("cond2\n");
    //        return -1;
    //    }
    //
    //    if (t0z > tmin) tmin = t0z;
    //
    //    if (t1z < tmax) tmax = t1z;
    //    
    //    //if (tmin > ray.tmax || tmax < ray.tmin) return false
    //    
    //    double t = tmax - tmin;
    //
    //    return t;
}

/////////////////////
// Animation stuff //
/////////////////////

Matrix3D::Matrix3D(const Point3D& e) {
    (*this) = Matrix3D();
}

Matrix3D::Matrix3D(const Quaternion& q) {
    (*this) = Matrix3D();
}

Matrix3D Matrix3D::closestRotation(void) const {
    return (*this);
}

/* While these Exp and Log implementations are the direct implementations of the Taylor series, the Log
 * function tends to run into convergence issues so we use the other ones:*/
Matrix3D Matrix3D::Exp(const Matrix3D& m, int iter) {
    return m;
}
