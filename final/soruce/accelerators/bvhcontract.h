
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

#ifndef PBRT_ACCELERATORS_BVH_CONTRACT_H
#define PBRT_ACCELERATORS_BVH_CONTRACT_H

// accelerators/bvhcontract.h*
#include "pbrt.h"
#include "primitive.h"

// BVHAccel Forward Declarations
#include "accelerators/bvh.h"
struct LinearBVHContractNode;

// BVHAccel Declarations
class BVHContractAccel : public Aggregate {
public:
    // BVHAccel Public Methods
    BVHContractAccel(const vector<Reference<Primitive> > &p, uint32_t maxPrims = 1,
             const string &sm = "sah");
    BBox WorldBound() const;
    bool CanIntersect() const { return true; }
    ~BVHContractAccel();
    bool Intersect(const Ray &ray, Intersection *isect) const;
    bool IntersectP(const Ray &ray) const;
	void tick();
private:
    // BVHAccel Private Methods
    BVHBuildNode *recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start, uint32_t end,
        uint32_t *totalNodes, vector<Reference<Primitive> > &orderedPrims);
    uint32_t flattenBVHTree(BVHBuildNode *node, uint32_t *offset, uint32_t parentOffset);
	// MORRIS ADD
	void recursiveContractSA(uint32_t uoffset);
	void recursiveContractRD(uint32_t uoffset);
	bool contractionCriterionSA(LinearBVHContractNode *node, LinearBVHContractNode *pnode);
	bool contractionCriterionRD(LinearBVHContractNode *node, LinearBVHContractNode *pnode);
	uint32_t flattenLinearBVHTree(uint32_t uoffset, LinearBVHContractNode mem[], uint32_t *offset, uint32_t parentOffset);
	void reorderChild(uint32_t uoffset, LinearBVHContractNode mem[]);
    // BVHAccel Private Data
    uint32_t maxPrimsInNode;
	uint32_t realNodes, intersectTest;
    enum SplitMethod { SPLIT_MIDDLE, SPLIT_EQUAL_COUNTS, SPLIT_SAH };
    SplitMethod splitMethod;
    vector<Reference<Primitive> > primitives;
    LinearBVHContractNode *nodes;
};


BVHContractAccel * CreateBVHContractAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps);

#endif // PBRT_ACCELERATORS_BVH_H
