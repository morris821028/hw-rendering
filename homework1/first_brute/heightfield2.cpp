
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield2.cpp*
#include "stdafx.h"
#include "shapes/heightfield2.h"
#include "shapes/trianglemesh.h"
#include "paramset.h"

// Heightfield2 Method Definitions
Heightfield2::Heightfield2(const Transform *o2w, const Transform *w2o,
        bool ro, int x, int y, const float *zs)
    : Shape(o2w, w2o, ro) {
    nx = x;
    ny = y;
    z = new float[nx*ny];
    memcpy(z, zs, nx*ny*sizeof(float));
}


Heightfield2::~Heightfield2() {
    delete[] z;
}


BBox Heightfield2::ObjectBound() const {
    float minz = z[0], maxz = z[0];
    for (int i = 1; i < nx*ny; ++i) {
        if (z[i] < minz) minz = z[i];
        if (z[i] > maxz) maxz = z[i];
    }
    return BBox(Point(0,0,minz), Point(1,1,maxz));
}


bool Heightfield2::CanIntersect() const {
	// modify by morris821028
    // return false;
	return true;
}

// not work if Heightfield2::CanIntersect() return false;
void Heightfield2::Refine(vector<Reference<Shape> > &refined) const {
    int ntris = 2*(nx-1)*(ny-1);
    refined.reserve(ntris);
    int *verts = new int[3*ntris];
    Point *P = new Point[nx*ny];
    float *uvs = new float[2*nx*ny];
    int nverts = nx*ny;
    int x, y;
    // Compute heightfield vertex positions
    int pos = 0;
    for (y = 0; y < ny; ++y) {
        for (x = 0; x < nx; ++x) {
            P[pos].x = uvs[2*pos]   = (float)x / (float)(nx-1);
            P[pos].y = uvs[2*pos+1] = (float)y / (float)(ny-1);
            P[pos].z = z[pos];
            ++pos;
        }
    }

    // Fill in heightfield vertex offset array
    int *vp = verts;
    for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
#define VERT(x,y) ((x)+(y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y);
            *vp++ = VERT(x+1, y+1);
    
            *vp++ = VERT(x, y);
            *vp++ = VERT(x+1, y+1);
            *vp++ = VERT(x, y+1);
        }
#undef VERT
    }
    ParamSet paramSet;
    paramSet.AddInt("indices", verts, 3*ntris);
    paramSet.AddFloat("uv", uvs, 2 * nverts);
    paramSet.AddPoint("P", P, nverts);
    refined.push_back(CreateTriangleMeshShape(ObjectToWorld, WorldToObject, ReverseOrientation, paramSet));
    delete[] P;
    delete[] uvs;
    delete[] verts;
}

/*
	grid[nx][ny] = height
 */
Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params) {
    int nu = params.FindOneInt("nu", -1);
    int nv = params.FindOneInt("nv", -1);
    int nitems;
    const float *Pz = params.FindFloat("Pz", &nitems);
    Assert(nitems == nu*nv);
    Assert(nu != -1 && nv != -1 && Pz != NULL);
    return new Heightfield2(o2w, w2o, reverseOrientation, nu, nv, Pz);
}

// custom implement R04922067 begin
bool Heightfield2::Intersect(const Ray &r, float *tHit, float *rayEpsilon,
        DifferentialGeometry *dg) const {
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

	// simplify 3D DDA to 2D DDA
	
	int x, y;
	bool has = false;
	for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
			float tx[2], ty[2];
			float tz[4] = {z[y*nx+x], z[y*nx+(x+1)], z[(y+1)*nx+x], z[(y+1)*nx+x+1]};
            tx[0] = (float)x / (float) (nx-1);
            ty[0] = (float)y / (float) (ny-1);
			tx[1] = (float)(x+1) / (float) (nx-1);
            ty[1] = (float)(y+1) / (float) (ny-1);
			Point p[4];
			p[0].x = tx[0], p[0].y = ty[0], p[0].z = tz[0]; // (0, 0)
			p[1].x = tx[1], p[1].y = ty[0], p[1].z = tz[1]; // (1, 0)
			p[2].x = tx[0], p[2].y = ty[1], p[2].z = tz[2]; // (0, 1)
			p[3].x = tx[1], p[3].y = ty[1], p[3].z = tz[3]; // (1, 1)
			{
				Point triangle[3] = {p[0], p[1], p[3]};
				if (triangleIntersect(r, tHit, rayEpsilon, dg, triangle)) {
					has = true;
				}
			}
			{
				Point triangle[3] = {p[0], p[3], p[2]};
				if (triangleIntersect(r, tHit, rayEpsilon, dg, triangle)) {
					has = true;
				}
			}
        }
    }
	return has;
}
// copy shapes/tianglemesh.cpp -> bool Triangle::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const
bool Heightfield2::triangleIntersect (const Ray &ray, float *tHit, float *rayEpsilon, 
									  DifferentialGeometry *dg, Point triangle[]) const {
	const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
	Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
	
	// Compute triangle partial derivatives
    Vector dpdu, dpdv;
	
	// Compute deltas for triangle partial derivatives
    float du1 = triangle[0].x - triangle[2].x;
    float du2 = triangle[1].x - triangle[2].x;
    float dv1 = triangle[0].y - triangle[2].y;
    float dv2 = triangle[1].y - triangle[2].y;
    Vector dp1 = p1 - p3, dp2 = p2 - p3;
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f) {
        // Handle zero determinant for triangle partial derivative matrix
        CoordinateSystem(Normalize(Cross(e2, e1)), &dpdu, &dpdv);
    }
    else {
        float invdet = 1.f / determinant;
        dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet;
        dpdv = (-du2 * dp1 + du1 * dp2) * invdet;
    }

	// Interpolate $(u,v)$ triangle parametric coordinates
    float b0 = 1 - b1 - b2;
    float tu = b0*triangle[0].x + b1*triangle[1].x + b2*triangle[2].x;
    float tv = b0*triangle[0].y + b1*triangle[1].y + b2*triangle[2].y;

    // Fill in _DifferentialGeometry_ from triangle hit
    *dg = DifferentialGeometry(ray(t), dpdu, dpdv,
                               Normal(0,0,0), Normal(0,0,0),
                               tu, tv, this);
	
	*tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
	return true;
}
// copy shapes/tianglemesh.cpp -> bool Triangle::IntersectP(const Ray &ray) const
bool Heightfield2::triangleIntersectP(const Ray &ray, Point triangle[]) const {
	const Point p1 = (*ObjectToWorld)(triangle[0]);
    const Point p2 = (*ObjectToWorld)(triangle[1]);
    const Point p3 = (*ObjectToWorld)(triangle[2]);
	Vector e1 = p2 - p1;
    Vector e2 = p3 - p1;
    Vector s1 = Cross(ray.d, e2);
    float divisor = Dot(s1, e1);
    
    if (divisor == 0.)
        return false;
    float invDivisor = 1.f / divisor;
	// Compute first barycentric coordinate
    Vector s = ray.o - p1;
    float b1 = Dot(s, s1) * invDivisor;
    if (b1 < 0. || b1 > 1.)
        return false;
    // Compute second barycentric coordinate
    Vector s2 = Cross(s, e1);
    float b2 = Dot(ray.d, s2) * invDivisor;
    if (b2 < 0. || b1 + b2 > 1.)
        return false;

    // Compute _t_ to intersection point
    float t = Dot(e2, s2) * invDivisor;
    if (t < ray.mint || t > ray.maxt)
        return false;
	return true;
}

bool Heightfield2::IntersectP(const Ray &r) const {
    // Transform _Ray_ to object space
    Ray ray;
    (*WorldToObject)(r, &ray);

	int x, y;
	bool has = false;
	for (y = 0; y < ny-1; ++y) {
        for (x = 0; x < nx-1; ++x) {
			float tx[2], ty[2];
			float tz[4] = {z[y*nx+x], z[y*nx+(x+1)], z[(y+1)*nx+x], z[(y+1)*nx+x+1]};
            tx[0] = (float)x / (float) (nx-1);
            ty[0] = (float)y / (float) (ny-1);
			tx[1] = (float)(x+1) / (float) (nx-1);
            ty[1] = (float)(y+1) / (float) (ny-1);
			Point p[4];
			p[0].x = tx[0], p[0].y = ty[0], p[0].z = tz[0]; // (0, 0)
			p[1].x = tx[1], p[1].y = ty[0], p[1].z = tz[1]; // (1, 0)
			p[2].x = tx[0], p[2].y = ty[1], p[2].z = tz[2]; // (0, 1)
			p[3].x = tx[1], p[3].y = ty[1], p[3].z = tz[3]; // (1, 1)
			{
				Point triangle[3] = {p[0], p[1], p[3]};
				if (triangleIntersectP(r, triangle)) {
					has = true;
				}
				if (has)	return true;
			}
			{
				Point triangle[3] = {p[0], p[3], p[2]};
				if (triangleIntersectP(r, triangle)) {
					has = true;
				}
				if (has)	return true;
			}
        }
    }
	return false;
}