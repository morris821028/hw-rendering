
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
	// Custom R04922067 begin
	bounds = ObjectBound();
	nVoxels[0] = nx-1;
	nVoxels[1] = ny-1;
	width = Vector(1.f / (nx-1), 1.f / (ny-1), 0);
	triangles = vector<mTriangle>(nx*ny*2);
	int pos = 0;
	for (int y = 0; y < ny-1; ++y) {
        for (int x = 0; x < nx-1; ++x) {
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
			mTriangle tmp;
			tmp.dx = width[0], tmp.dy = width[1], tmp.bx = tx[0], tmp.by = ty[0];
			tmp.p[0] = p[0], tmp.p[1] = p[1], tmp.p[2] = p[3];
			triangles[pos++] = tmp;
			tmp.p[0] = p[0], tmp.p[1] = p[3], tmp.p[2] = p[2];
			triangles[pos++] = tmp;
        }
    }
	Grid2D tmp = {nVoxels, bounds, width, Vector(nx-1, ny-1, 0)};
	GRID2D = tmp;

	// Phong interpolation, computing normal for each point on the grid.
	m_normals = new Normal[nx * ny];
	initFieldNormals();
	// Custom R04922067 end
}


Heightfield2::~Heightfield2() {
    delete[] z;
	delete[] m_normals;
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

	// Start 3D DDA

	// Check ray against overall grid bounds
    float rayT;
	if (bounds.Inside(ray(ray.mint))) 
		rayT = ray.mint;
	else if (!bounds.IntersectP(ray, &rayT))
		return false;
	Point gridIntersect = ray(rayT);

	// Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = GRID2D.posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (GRID2D.voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (GRID2D.voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }
	
	// Walk grid for shadow ray
	bool has = false;
    for (;;) {
        int o = GRID2D.offset(Pos[0], Pos[1]) * 2;
        
		if (triangles[o|0].Intersect(ray, tHit, rayEpsilon, dg, this, *ObjectToWorld))
            return true;
		if (triangles[o|1].Intersect(ray, tHit, rayEpsilon, dg, this, *ObjectToWorld))
            return true;
        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
        if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
	return false;
}
// copy shapes/tianglemesh.cpp -> bool Triangle::Intersect(const Ray &ray, float *tHit, float *rayEpsilon, DifferentialGeometry *dg) const
bool mTriangle::Intersect (const Ray &ray, float *tHit, float *rayEpsilon, 
									  DifferentialGeometry *dg, const Shape *belong, const Transform &objToWorld) const {
	PBRT_RAY_TRIANGLE_INTERSECTIONP_TEST(const_cast<Ray *>(&ray), const_cast<Triangle *>(this));
	const Point p1 = p[0], p2 = p[1], p3 = p[2];
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
    float du1 = p1.x - p3.x;
    float du2 = p2.x - p3.x;
    float dv1 = p1.y - p3.y;
    float dv2 = p2.y - p3.y;
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
    float tu = b0*p1.x + b1*p2.x + b2*p3.x;
    float tv = b0*p1.y + b1*p2.y + b2*p3.y;
	Point tp = (ray(t));
	tu = tp.x, tv = tp.y; 
    *dg = DifferentialGeometry((objToWorld) (ray(t)), (objToWorld) (dpdu), (objToWorld) (dpdv),
                               (objToWorld) (Normal(0,0,0)), (objToWorld) (Normal(0,0,0)),
                               tu, tv, belong);
	*tHit = t;
    *rayEpsilon = 1e-3f * *tHit;
	 PBRT_RAY_TRIANGLE_INTERSECTIONP_HIT(const_cast<Ray *>(&ray), t);
	return true;
}
// copy shapes/tianglemesh.cpp -> bool Triangle::IntersectP(const Ray &ray) const
bool mTriangle::IntersectP(const Ray &ray, const Shape *belong) const {
	// Get triangle vertices in _p1_, _p2_, and _p3_
	const Point p1 = p[0], p2 = p[1], p3 = p[2];
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

	// Start 3D DDA

	// Check ray against overall grid bounds
    float rayT;
	if (bounds.Inside(ray(ray.mint))) 
		rayT = ray.mint;
	else if (!bounds.IntersectP(ray, &rayT))
		return false;
	Point gridIntersect = ray(rayT);
	
	// Set up 2D DDA for ray
    float NextCrossingT[2], DeltaT[2];
    int Step[2], Out[2], Pos[2];
    for (int axis = 0; axis < 2; ++axis) {
        // Compute current voxel for axis
        Pos[axis] = GRID2D.posToVoxel(gridIntersect, axis);
        if (ray.d[axis] >= 0) {
            // Handle ray with positive direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (GRID2D.voxelToPos(Pos[axis]+1, axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = width[axis] / ray.d[axis];
            Step[axis] = 1;
            Out[axis] = nVoxels[axis];
        }
        else {
            // Handle ray with negative direction for voxel stepping
            NextCrossingT[axis] = rayT +
                (GRID2D.voxelToPos(Pos[axis], axis) - gridIntersect[axis]) / ray.d[axis];
            DeltaT[axis] = -width[axis] / ray.d[axis];
            Step[axis] = -1;
            Out[axis] = -1;
        }
    }
	
	// Walk grid for shadow ray
	bool has = false;
    for (;;) {
        int o = GRID2D.offset(Pos[0], Pos[1]) * 2;
        
        if (triangles[o|0].IntersectP(ray, this))
            return true;
		if (triangles[o|1].IntersectP(ray, this))
            return true;
        // Advance to next voxel

        // Find _stepAxis_ for stepping to next voxel
        int stepAxis = NextCrossingT[0] < NextCrossingT[1] ? 0 : 1;
        if (ray.maxt < NextCrossingT[stepAxis])
            break;
        Pos[stepAxis] += Step[stepAxis];
        if (Pos[stepAxis] == Out[stepAxis])
            break;
        NextCrossingT[stepAxis] += DeltaT[stepAxis];
    }
	return false;
}

// Phong interpolation of normals across, initFieldNormals & GetShadingGeometry

Normal Heightfield2::getGridNormal(int i,int j) const {
#define VERT(x,y) ((x)+(y)*nx)
   Vector dx(1.f * width[0], 0, 0), dy(0, 1.f * width[1], 0);
   if (i == nx-1) 
		dx.z = z[VERT(i, j)] - z[VERT(i-1, j)];
   else if (i == 0) 
		dx.z = z[VERT(i+1, j)] - z[VERT(i, j)];
   else 
		dx.z = 0.5f * (z[VERT(i+1, j)] - z[VERT(i-1, j)]);
   
   if (j == ny-1) 
		dy.z = z[VERT(i, j)] - z[VERT(i, j-1)];
   else if (j == 0)
		dy.z = z[VERT(i, j+1)] - z[VERT(i, j)];
   else 
		dy.z = 0.5f * (z[VERT(i, j+1)] - z[VERT(i, j-1)]);
#undef VERT
   return Normalize(Normal(Cross(dx, dy)));
}
void Heightfield2::initFieldNormals() {
	for (int j = 0; j < ny; ++j) {
		for (int i = 0; i < nx; ++i) {
			m_normals[i + nx * j] = getGridNormal(i, j);
		}
	}
}

void Heightfield2::GetShadingGeometry(const Transform &obj2world,
        const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const {
	// *dgShading = dg;
	// return ;
	int x = Clamp(Float2Int(dg.u * (nx-1)), 0, nx - 1), 
		y = Clamp(Float2Int(dg.v * (ny-1)), 0, ny - 1);

	int type = dg.u - x * width[0] <= dg.v - y * width[1];	// upper or lower triangle

	Point p1 = Point(x * width[0], y * width[1], z[x + nx * y]),
		  p2 = Point((x + 1) * width[0], (y + 1) * width[1], z[(x + 1) + nx * (y + 1)]),
		  p3 = Point(x * width[0], (y + 1) * width[1], z[x + nx * (y + 1)]),
		  p4 = Point(dg.u, p1.y + (dg.u - p1.x), 0);	// on diagonal line

	if (type == 0)
		p3 = Point((x + 1) * width[0], y * width[1], z[(x + 1) + nx * y]);
			
	// Compute the normal at the hit point
	int q = type == 0 ? (x + 1) + nx * (y) : (x) + nx * (y + 1);
	Normal normals[3] = {m_normals[x + nx * y], m_normals[(x + 1) + nx * (y + 1)], m_normals[q]};	// normal of p1, p2, p3

	// Compute normal with Phong interpolation
	Point a = Point(p1.x, p1.y, 0), 
		  b = Point(p2.x, p2.y, 0), 
		  c = Point(p4.x, p4.y, 0);
	Normal n1 = ((a - c).Length() * normals[1] + (b - c).Length() * normals[0]),	// diagonal normal
		   n2 = type == 0 ? (fabs(dg.u - p3.x) * normals[0] + fabs(p1.x - dg.u) * normals[2]) : 
						    (fabs(dg.u - p2.x) * normals[2] + fabs(p3.x - dg.u) * normals[1]) ;
	float d1 = type == 0 ? fabs(p4.y - dg.v) : fabs(p4.y - dg.v),
		  d2 = type == 0 ? fabs(p1.y - dg.v) : fabs(p2.y - dg.v);
	Normal hitNormal = (*ObjectToWorld) (Normalize(d2 * Normalize(n1) + d1 * Normalize(n2)));

	// Use _n_ and _s_ to compute shading tangents for triangle, _ss_ and _ts_
    Normal ns = hitNormal;
    Vector ss, ts;
    ss = Normalize(dg.dpdu);
    
    ts = Cross(ss, ns);
    if (ts.LengthSquared() > 0.f) {
        ts = Normalize(ts);
        ss = Cross(ts, ns);
    }
    else
        CoordinateSystem((Vector)ns, &ss, &ts);

	// Compute deltas for triangle partial derivatives of normal
	Normal dndu, dndv;
    float du1 = p2.x - p1.x;
    float du2 = p3.x - p1.x;
    float dv1 = p2.y - p1.y;
    float dv2 = p3.y - p1.y;
	Normal dn1 = normals[1] - normals[0], dn2 = normals[2] - normals[0];
    float determinant = du1 * dv2 - dv1 * du2;
    if (determinant == 0.f)
        dndu = dndv = Normal(0,0,0);
    else {
        float invdet = 1.f / determinant;
        dndu = ( dv2 * dn1 - dv1 * dn2) * invdet;
        dndv = (-du2 * dn1 + du1 * dn2) * invdet;
    }

	*dgShading = DifferentialGeometry(dg.p, ss, ts,
				(*ObjectToWorld)(dndu), (*ObjectToWorld)(dndv), dg.u, dg.v, dg.shape);
	dgShading->dudx = dg.dudx;  dgShading->dvdx = dg.dvdx;
	dgShading->dudy = dg.dudy;  dgShading->dvdy = dg.dvdy;
	dgShading->dpdx = dg.dpdx;  dgShading->dpdy = dg.dpdy;
}