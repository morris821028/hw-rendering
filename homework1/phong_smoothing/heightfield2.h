
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_SHAPES_HEIGHTFIELD2_H
#define PBRT_SHAPES_HEIGHTFIELD2_H

// shapes/heightfield2.h*
#include "shape.h"

// Heightfield2 Forward Declarations
struct mTriangle;

// mTriangle Declarations
struct mTriangle {
	// Public Member
	Point p[3];
	float dx, dy, bx, by;
    // Public Methods
    mTriangle() { }
    bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon, 
									  DifferentialGeometry *dg, const Shape *belong, const Transform &objToWorld) const;
	bool IntersectP(const Ray &ray, const Shape *belong) const;
};

// Grid2D Declarations
struct Grid2D {
	// Public Member
	int *nVoxels;
	BBox bounds;
	Vector width, invWidth;
	// Public Methods
	int posToVoxel(const Point &P, int axis) const {
		int v = Float2Int((P[axis] - bounds.pMin[axis]) *
							invWidth[axis]);
		return Clamp(v, 0, nVoxels[axis]-1);
	}
	float voxelToPos(int p, int axis) const {
		return bounds.pMin[axis] + p * width[axis];
	}
	inline int offset(int x, int y) const {
		return y*nVoxels[0] + x;
	}
};

// Heightfield2 Declarations
class Heightfield2 : public Shape {
public:
    // Heightfield2 Public Methods
    Heightfield2(const Transform *o2, const Transform *w2o, bool ro, int nu, int nv, const float *zs);
    ~Heightfield2();
    bool CanIntersect() const;
    void Refine(vector<Reference<Shape> > &refined) const;
	// custom implement R04922067 begin
	bool Intersect(const Ray &ray, float *tHit, float *rayEpsilon,
                   DifferentialGeometry *dg) const;
    bool IntersectP(const Ray &ray) const;
	virtual void GetShadingGeometry(const Transform &obj2world,
									const DifferentialGeometry &dg, DifferentialGeometry *dgShading) const;
	// custom implement R04922067 end
    BBox ObjectBound() const;
private:
    // Heightfield Private Data
    float *z;
    int nx, ny;
	// custom implement R04922067 begin
	bool triangleIntersectP(const Ray &ray, Point * triangle) const;
	vector<mTriangle> triangles;
	Grid2D GRID2D;
	BBox bounds;
	int nVoxels[2];
	Vector width;
	// Phong interpolation of normals across
	Normal *m_normals;
	Normal Heightfield2::getGridNormal(int idxX,int idxY) const;
	void initFieldNormals();
	// custom implement R04922067 end;
};

Heightfield2 *CreateHeightfield2Shape(const Transform *o2w, const Transform *w2o,
        bool reverseOrientation, const ParamSet &params);

#endif // PBRT_SHAPES_HEIGHTFIELD_H