
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


// accelerators/bvhcontract.cpp*
#include "stdafx.h"
#include "accelerators/bvhcontract.h"
#include "accelerators/bvh.h"
#include "probes.h"
#include "paramset.h"

// BVHAccel Local Declarations

extern struct BVHPrimitiveInfo {
    BVHPrimitiveInfo() { }
    BVHPrimitiveInfo(int pn, const BBox &b)
        : primitiveNumber(pn), bounds(b) {
        centroid = .5f * b.pMin + .5f * b.pMax;
    }
    int primitiveNumber;
    Point centroid;
    BBox bounds;
};


extern struct BVHBuildNode {
    // BVHBuildNode Public Methods
    BVHBuildNode() { children[0] = children[1] = NULL; }
    void InitLeaf(uint32_t first, uint32_t n, const BBox &b) {
        firstPrimOffset = first;
        nPrimitives = n;
        bounds = b;
    }
    void InitInterior(uint32_t axis, BVHBuildNode *c0, BVHBuildNode *c1) {
        children[0] = c0;
        children[1] = c1;
        bounds = Union(c0->bounds, c1->bounds);
        splitAxis = axis;
        nPrimitives = 0;
    }
    BBox bounds;
    BVHBuildNode *children[2];
    uint32_t splitAxis, firstPrimOffset, nPrimitives;
};

extern struct CompareToMid {
    CompareToMid(int d, float m) { dim = d; mid = m; }
    int dim;
    float mid;
    bool operator()(const BVHPrimitiveInfo &a) const {
        return a.centroid[dim] < mid;
    }
};


extern struct ComparePoints {
    ComparePoints(int d) { dim = d; }
    int dim;
    bool operator()(const BVHPrimitiveInfo &a,
                    const BVHPrimitiveInfo &b) const {
        return a.centroid[dim] < b.centroid[dim];
    }
};


extern struct CompareToBucket {
    CompareToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const BVHPrimitiveInfo &p) const;

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};

extern struct LinearBVHNode {
    BBox bounds;
    union {
        uint32_t primitivesOffset;    // leaf
        uint32_t secondChildOffset;   // interior
    };

    uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};

struct LinearBVHContractNode {
    BBox bounds;
	uint32_t primitivesOffset;    // leaf
	uint32_t parentOffset;
	uint32_t childOffsetHead, childOffsetTail;	// interior
	uint32_t siblingOffsetNext, siblingOffsetPrev;	// interior

	uint32_t numChild;
	uint32_t visitCount;	// contract record

	uint8_t nPrimitives;  // 0 -> interior node
    uint8_t axis;         // interior node: xyz
    uint8_t pad[2];       // ensure 32 byte total size
};

static inline bool IntersectP(const BBox &bounds, const Ray &ray,
        const Vector &invDir, const uint32_t dirIsNeg[3]) {
    // Check for ray intersection against $x$ and $y$ slabs
    float tmin =  (bounds[  dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tmax =  (bounds[1-dirIsNeg[0]].x - ray.o.x) * invDir.x;
    float tymin = (bounds[  dirIsNeg[1]].y - ray.o.y) * invDir.y;
    float tymax = (bounds[1-dirIsNeg[1]].y - ray.o.y) * invDir.y;
    if ((tmin > tymax) || (tymin > tmax))
        return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    // Check for ray intersection against $z$ slab
    float tzmin = (bounds[  dirIsNeg[2]].z - ray.o.z) * invDir.z;
    float tzmax = (bounds[1-dirIsNeg[2]].z - ray.o.z) * invDir.z;
    if ((tmin > tzmax) || (tzmin > tmax))
        return false;
    if (tzmin > tmin)
        tmin = tzmin;
    if (tzmax < tmax)
        tmax = tzmax;
    return (tmin < ray.maxt) && (tmax > ray.mint);
}


#define MORRISDEBUG
// #define TRAVERSAL_REC
#define TRAVERSAL_LOOP
// BVHAccel Method Definitions
BVHContractAccel::BVHContractAccel(const vector<Reference<Primitive> > &p,
                   uint32_t mp, const string &sm) {

#ifdef MORRISDEBUG
	fprintf(stderr, "sizeof(LinearBVHContractNode) = %d, sizeof(LinearBVHNode) = %d\n", sizeof(LinearBVHContractNode), sizeof(LinearBVHNode));
	fprintf(stderr, "BVHContractAccle Create\n");
#endif
    maxPrimsInNode = min(255u, mp);
    for (uint32_t i = 0; i < p.size(); ++i)
        p[i]->FullyRefine(primitives);
    if (sm == "sah")         splitMethod = SPLIT_SAH;
    else if (sm == "middle") splitMethod = SPLIT_MIDDLE;
    else if (sm == "equal")  splitMethod = SPLIT_EQUAL_COUNTS;
    else {
        Warning("BVH split method \"%s\" unknown.  Using \"sah\".",
                sm.c_str());
        splitMethod = SPLIT_SAH;
    }

    if (primitives.size() == 0) {
        nodes = NULL;
        return;
    }
    // Build BVH from _primitives_
    PBRT_BVH_STARTED_CONSTRUCTION(this, primitives.size());

    // Initialize _buildData_ array for primitives
    vector<BVHPrimitiveInfo> buildData;
    buildData.reserve(primitives.size());
    for (uint32_t i = 0; i < primitives.size(); ++i) {
        BBox bbox = primitives[i]->WorldBound();
        buildData.push_back(BVHPrimitiveInfo(i, bbox));
    }

    // Recursively build BVH tree for primitives
    MemoryArena buildArena;
    uint32_t totalNodes = 0;
    vector<Reference<Primitive> > orderedPrims;
    orderedPrims.reserve(primitives.size());
    BVHBuildNode *root = recursiveBuild(buildArena, buildData, 0,
                                        primitives.size(), &totalNodes,
                                        orderedPrims);
    primitives.swap(orderedPrims);
        Info("BVH created with %d nodes for %d primitives (%.2f MB)", totalNodes,
             (int)primitives.size(), float(totalNodes * sizeof(LinearBVHNode))/(1024.f*1024.f));

    // Compute representation of depth-first traversal of BVH tree
	fprintf(stderr, "Total Nodes = %d\n", totalNodes);
    nodes = AllocAligned<LinearBVHContractNode>(totalNodes);
    for (uint32_t i = 0; i < totalNodes; ++i)
        new (&nodes[i]) LinearBVHContractNode;
    uint32_t offset = 0;
    flattenBVHTree(root, &offset, -1);

	if (totalNodes) {
		intersectTest = 0;
		realNodes = totalNodes;
		recursiveContractSA(0);
		fprintf(stderr, "Preprocessing Contract %d / %d\n", realNodes, totalNodes);
		fprintf(stderr, "re-flattenBVH tree %d / %d\n", realNodes, totalNodes);
		LinearBVHContractNode *mem = AllocAligned<LinearBVHContractNode>(realNodes);
		for (uint32_t i = 0; i < realNodes; ++i)
			new (&mem[i]) LinearBVHContractNode;
		uint32_t toffset = 0;
		flattenLinearBVHTree(0, mem, &toffset, -1);
		Assert(toffset == realNodes);
		totalNodes = realNodes, offset = totalNodes;
		FreeAligned(nodes);
		nodes = mem;
	}

    Assert(offset == totalNodes);
    PBRT_BVH_FINISHED_CONSTRUCTION(this);
}


BBox BVHContractAccel::WorldBound() const {
    return nodes ? nodes[0].bounds : BBox();
}


BVHBuildNode *BVHContractAccel::recursiveBuild(MemoryArena &buildArena,
        vector<BVHPrimitiveInfo> &buildData, uint32_t start,
        uint32_t end, uint32_t *totalNodes,
        vector<Reference<Primitive> > &orderedPrims) {
    Assert(start != end);
    (*totalNodes)++;
    BVHBuildNode *node = buildArena.Alloc<BVHBuildNode>();
    // Compute bounds of all primitives in BVH node
    BBox bbox;
    for (uint32_t i = start; i < end; ++i)
        bbox = Union(bbox, buildData[i].bounds);
    uint32_t nPrimitives = end - start;
    if (nPrimitives == 1) {
        // Create leaf _BVHBuildNode_
        uint32_t firstPrimOffset = orderedPrims.size();
        for (uint32_t i = start; i < end; ++i) {
            uint32_t primNum = buildData[i].primitiveNumber;
            orderedPrims.push_back(primitives[primNum]);
        }
        node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
    }
    else {
        // Compute bound of primitive centroids, choose split dimension _dim_
        BBox centroidBounds;
        for (uint32_t i = start; i < end; ++i)
            centroidBounds = Union(centroidBounds, buildData[i].centroid);
        int dim = centroidBounds.MaximumExtent();

        // Partition primitives into two sets and build children
        uint32_t mid = (start + end) / 2;
        if (centroidBounds.pMax[dim] == centroidBounds.pMin[dim]) {
            // If nPrimitives is no greater than maxPrimsInNode,
            // then all the nodes can be stored in a compact bvh node.
            if (nPrimitives <= maxPrimsInNode) {
                // Create leaf _BVHBuildNode_
                uint32_t firstPrimOffset = orderedPrims.size();
                for (uint32_t i = start; i < end; ++i) {
                    uint32_t primNum = buildData[i].primitiveNumber;
                    orderedPrims.push_back(primitives[primNum]);
                }
                node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                return node;
            }
            else {
                // else if nPrimitives is greater than maxPrimsInNode, we
                // need to split it further to guarantee each node contains
                // no more than maxPrimsInNode primitives.
                node->InitInterior(dim,
                                   recursiveBuild(buildArena, buildData, start, mid,
                                                  totalNodes, orderedPrims),
                                   recursiveBuild(buildArena, buildData, mid, end,
                                                  totalNodes, orderedPrims));
                return node;
            }
        }

        // Partition primitives based on _splitMethod_
        switch (splitMethod) {
        case SPLIT_MIDDLE: {
            // Partition primitives through node's midpoint
            float pmid = .5f * (centroidBounds.pMin[dim] + centroidBounds.pMax[dim]);
            BVHPrimitiveInfo *midPtr = std::partition(&buildData[start],
                                                      &buildData[end-1]+1,
                                                      CompareToMid(dim, pmid));
            mid = midPtr - &buildData[0];
            if (mid != start && mid != end)
                // for lots of prims with large overlapping bounding boxes, this
                // may fail to partition; in that case don't break and fall through
                // to SPLIT_EQUAL_COUNTS
                break;
        }
        case SPLIT_EQUAL_COUNTS: {
            // Partition primitives into equally-sized subsets
            mid = (start + end) / 2;
            std::nth_element(&buildData[start], &buildData[mid],
                             &buildData[end-1]+1, ComparePoints(dim));
            break;
        }
        case SPLIT_SAH: default: {
            // Partition primitives using approximate SAH
            if (nPrimitives <= 4) {
                // Partition primitives into equally-sized subsets
                mid = (start + end) / 2;
                std::nth_element(&buildData[start], &buildData[mid],
                                 &buildData[end-1]+1, ComparePoints(dim));
            }
            else {
                // Allocate _BucketInfo_ for SAH partition buckets
                const int nBuckets = 12;
                struct BucketInfo {
                    BucketInfo() { count = 0; }
                    int count;
                    BBox bounds;
                };
                BucketInfo buckets[nBuckets];

                // Initialize _BucketInfo_ for SAH partition buckets
                for (uint32_t i = start; i < end; ++i) {
                    int b = nBuckets *
                        ((buildData[i].centroid[dim] - centroidBounds.pMin[dim]) /
                         (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
                    if (b == nBuckets) b = nBuckets-1;
                    Assert(b >= 0 && b < nBuckets);
                    buckets[b].count++;
                    buckets[b].bounds = Union(buckets[b].bounds, buildData[i].bounds);
                }

                // Compute costs for splitting after each bucket
                float cost[nBuckets-1];
                for (int i = 0; i < nBuckets-1; ++i) {
                    BBox b0, b1;
                    int count0 = 0, count1 = 0;
                    for (int j = 0; j <= i; ++j) {
                        b0 = Union(b0, buckets[j].bounds);
                        count0 += buckets[j].count;
                    }
                    for (int j = i+1; j < nBuckets; ++j) {
                        b1 = Union(b1, buckets[j].bounds);
                        count1 += buckets[j].count;
                    }
                    cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
                              bbox.SurfaceArea();
                }

                // Find bucket to split at that minimizes SAH metric
                float minCost = cost[0];
                uint32_t minCostSplit = 0;
                for (int i = 1; i < nBuckets-1; ++i) {
                    if (cost[i] < minCost) {
                        minCost = cost[i];
                        minCostSplit = i;
                    }
                }

                // Either create leaf or split primitives at selected SAH bucket
                if (nPrimitives > maxPrimsInNode ||
                    minCost < nPrimitives) {
                    BVHPrimitiveInfo *pmid = std::partition(&buildData[start],
                        &buildData[end-1]+1,
                        CompareToBucket(minCostSplit, nBuckets, dim, centroidBounds));
                    mid = pmid - &buildData[0];
                }
                
                else {
                    // Create leaf _BVHBuildNode_
                    uint32_t firstPrimOffset = orderedPrims.size();
                    for (uint32_t i = start; i < end; ++i) {
                        uint32_t primNum = buildData[i].primitiveNumber;
                        orderedPrims.push_back(primitives[primNum]);
                    }
                    node->InitLeaf(firstPrimOffset, nPrimitives, bbox);
                    return node;
                }
            }
            break;
        }
        }
        node->InitInterior(dim,
                           recursiveBuild(buildArena, buildData, start, mid,
                                          totalNodes, orderedPrims),
                           recursiveBuild(buildArena, buildData, mid, end,
                                          totalNodes, orderedPrims));
    }
    return node;
}


uint32_t BVHContractAccel::flattenBVHTree(BVHBuildNode *node, uint32_t *offset, uint32_t parentOffset) {
    LinearBVHContractNode *linearNode = &nodes[*offset];
    linearNode->bounds = node->bounds;
	linearNode->childOffsetHead = linearNode->childOffsetTail = -1;
	linearNode->siblingOffsetNext = linearNode->siblingOffsetPrev = -1;
	linearNode->numChild = linearNode->visitCount = 0;
	linearNode->parentOffset = parentOffset;
    uint32_t myOffset = (*offset)++;
    if (node->nPrimitives > 0) {
        Assert(!node->children[0] && !node->children[1]);
        linearNode->primitivesOffset = node->firstPrimOffset;
        linearNode->nPrimitives = node->nPrimitives;
    }
    else {
        // Creater interior flattened BVH node
        linearNode->axis = node->splitAxis;
        linearNode->nPrimitives = 0;
		linearNode->numChild = 2;
        linearNode->childOffsetHead = flattenBVHTree(node->children[0], offset, myOffset);
		linearNode->childOffsetTail = flattenBVHTree(node->children[1], offset, myOffset);
		LinearBVHContractNode *prevSiblingNode = &nodes[linearNode->childOffsetHead];
		LinearBVHContractNode *nextSiblingNode = &nodes[linearNode->childOffsetTail];
		prevSiblingNode->siblingOffsetNext = linearNode->childOffsetTail;
		nextSiblingNode->siblingOffsetPrev = linearNode->childOffsetHead;
    }
    return myOffset;
}


BVHContractAccel::~BVHContractAccel() {
    FreeAligned(nodes);
}

void BVHContractAccel::tick() {
	if (!nodes)	return ;

	uint32_t test = nodes[0].visitCount;
	if ((test>>24) > intersectTest && intersectTest == 0) {
		intersectTest = test>>24;
		fprintf(stderr, "BVHContractAccel tick()\n");
		recursiveContractSA(0);
		fprintf(stderr, "Ray-Distribution Guided Contraction Remain Node %d\n", realNodes);
		fflush(stderr);
	}
	return ;
}
bool BVHContractAccel::recIntersect(uint32_t offset, const Ray &ray, Intersection *isect, 
									const Vector &invDir, const uint32_t dirIsNeg[]) const {
	LinearBVHContractNode *node = &nodes[offset];
	(node->visitCount)++;
	bool hit = false;
	if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
		if (node->nPrimitives > 0) {
			// Intersect ray with primitives in leaf BVH node
			PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHContractNode *>(node));
			for (uint32_t i = 0; i < node->nPrimitives; ++i)
			{
				PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
				if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
				{
					PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
					hit = true;
				}
				else {
					PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
				}
			}
		} else {
			if (dirIsNeg[node->axis]) {
				uint32_t ch = node->childOffsetTail;
				LinearBVHContractNode *p;
				for (; ch != -1; ch = p->siblingOffsetPrev) {
					p = &nodes[ch];
					hit |= recIntersect(ch, ray, isect, invDir, dirIsNeg);
				}
			} else {
				uint32_t ch = node->childOffsetHead;
				LinearBVHContractNode *p;
				for (; ch != -1; ch = p->siblingOffsetNext) {
					p = &nodes[ch];
					hit |= recIntersect(ch, ray, isect, invDir, dirIsNeg);
				}
			}
		}
	}
	return hit;
}
bool BVHContractAccel::Intersect(const Ray &ray, Intersection *isect) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTION_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    bool hit = false;
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };
    // Follow ray through BVH nodes to find primitive intersections
#ifdef TRAVERSAL_REC
	hit = recIntersect(0, ray, isect, invDir, dirIsNeg);
#endif
#ifdef TRAVERSAL_LOOP
	bool process = true;
	uint32_t idxOffset = 0;
	while (idxOffset != -1) {
		LinearBVHContractNode *node = &nodes[idxOffset];
		(node->visitCount)++;
		if (process) {
			// Check ray against BVH node
			if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
				if (node->nPrimitives > 0) {
					// Intersect ray with primitives in leaf BVH node
					PBRT_BVH_INTERSECTION_TRAVERSED_LEAF_NODE(const_cast<LinearBVHContractNode *>(node));
					for (uint32_t i = 0; i < node->nPrimitives; ++i)
					{
						PBRT_BVH_INTERSECTION_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
						if (primitives[node->primitivesOffset+i]->Intersect(ray, isect))
						{
							PBRT_BVH_INTERSECTION_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
							hit = true;
						}
						else {
							PBRT_BVH_INTERSECTION_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
					   }
					}
				}
			} else {
				process = false;
			}
		}

		if (dirIsNeg[node->axis]) {
			if (node->childOffsetTail != -1 && process) {
				idxOffset = node->childOffsetTail;
				process = true;
			} else if (node->parentOffset != -1) {
				const LinearBVHContractNode *pNode = &nodes[node->parentOffset];
				if (dirIsNeg[pNode->axis] && node->siblingOffsetPrev != -1) {
					idxOffset = node->siblingOffsetPrev;
					process = true;
				} else if (!dirIsNeg[pNode->axis] && node->siblingOffsetNext != -1) {
					idxOffset = node->siblingOffsetNext;
					process = true;
				} else {
					idxOffset = node->parentOffset;
					process = false;
				}
			} else {
				idxOffset = node->parentOffset;
				process = false;
			}
		} else {
			if (node->childOffsetHead != -1 && process) {
				idxOffset = node->childOffsetHead;
				process = true;
			} else if (node->parentOffset != -1) {
				const LinearBVHContractNode *pNode = &nodes[node->parentOffset];
				if (dirIsNeg[pNode->axis] && node->siblingOffsetPrev != -1) {
					idxOffset = node->siblingOffsetPrev;
					process = true;
				} else if (!dirIsNeg[pNode->axis] && node->siblingOffsetNext != -1) {
					idxOffset = node->siblingOffsetNext;
					process = true;
				} else {
					idxOffset = node->parentOffset;
					process = false;
				}
			} else {
				idxOffset = node->parentOffset;
				process = false;
			}
		}
	}
#endif
    PBRT_BVH_INTERSECTION_FINISHED();
    return hit;
}

bool BVHContractAccel::recIntersectP(uint32_t offset, const Ray &ray, 
									 const Vector &invDir, const uint32_t dirIsNeg[]) const {
	LinearBVHContractNode *node = &nodes[offset];
	(node->visitCount)++;
	if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
		if (node->nPrimitives > 0) {
			PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHContractNode *>(node));
				for (uint32_t i = 0; i < node->nPrimitives; ++i) {
				PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
				if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
					PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
					return true;
				}
				else {
					PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
				}
			}
		} else {
			if (dirIsNeg[node->axis]) {
				uint32_t ch = node->childOffsetTail;
				LinearBVHContractNode *p;
				for (; ch != -1; ch = p->siblingOffsetPrev) {
					p = &nodes[ch];
					if (recIntersectP(ch, ray, invDir, dirIsNeg))
						return true;
				}
			} else {
				uint32_t ch = node->childOffsetHead;
				LinearBVHContractNode *p;
				for (; ch != -1; ch = p->siblingOffsetNext) {
					p = &nodes[ch];
					if (recIntersectP(ch, ray, invDir, dirIsNeg))
						return true;
				}
			}
		}
	}
	return false;
}
bool BVHContractAccel::IntersectP(const Ray &ray) const {
    if (!nodes) return false;
    PBRT_BVH_INTERSECTIONP_STARTED(const_cast<BVHAccel *>(this), const_cast<Ray *>(&ray));
    Vector invDir(1.f / ray.d.x, 1.f / ray.d.y, 1.f / ray.d.z);
    uint32_t dirIsNeg[3] = { invDir.x < 0, invDir.y < 0, invDir.z < 0 };

#ifdef TRAVERSAL_REC
	PBRT_BVH_INTERSECTIONP_FINISHED();
	return recIntersectP(0, ray, invDir, dirIsNeg);
#endif
#ifdef TRAVERSAL_LOOP
	bool process = true;
	uint32_t idxOffset = 0;
	while (idxOffset != -1) {
		LinearBVHContractNode *node = &nodes[idxOffset];
		(node->visitCount)++;
		if (process) {
			// Check ray against BVH node
			if (::IntersectP(node->bounds, ray, invDir, dirIsNeg)) {
				if (node->nPrimitives > 0) {
					PBRT_BVH_INTERSECTIONP_TRAVERSED_LEAF_NODE(const_cast<LinearBVHContractNode *>(node));
					  for (uint32_t i = 0; i < node->nPrimitives; ++i) {
						PBRT_BVH_INTERSECTIONP_PRIMITIVE_TEST(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
						if (primitives[node->primitivesOffset+i]->IntersectP(ray)) {
							PBRT_BVH_INTERSECTIONP_PRIMITIVE_HIT(const_cast<Primitive *>(primitives[node->primitivesOffset+i].GetPtr()));
							return true;
						}
						else {
							PBRT_BVH_INTERSECTIONP_PRIMITIVE_MISSED(const_cast<Primitive *>(primitives[node->primitivesOffset + i].GetPtr()));
						}
					}
				}
			} else {
				process = false;
			}
		}
		
		if (dirIsNeg[node->axis]) {
			if (node->childOffsetTail != -1 && process) {
				idxOffset = node->childOffsetTail;
				process = true;
			} else if (node->parentOffset != -1) {
				const LinearBVHContractNode *pNode = &nodes[node->parentOffset];
				if (dirIsNeg[pNode->axis] && node->siblingOffsetPrev != -1) {
					idxOffset = node->siblingOffsetPrev;
					process = true;
				} else if (!dirIsNeg[pNode->axis] && node->siblingOffsetNext != -1) {
					idxOffset = node->siblingOffsetNext;
					process = true;
				} else {
					idxOffset = node->parentOffset;
					process = false;
				}
			} else {
				idxOffset = node->parentOffset;
				process = false;
			}
		} else {
			if (node->childOffsetHead != -1 && process) {
				idxOffset = node->childOffsetHead;
				process = true;
			} else if (node->parentOffset != -1) {
				const LinearBVHContractNode *pNode = &nodes[node->parentOffset];
				if (dirIsNeg[pNode->axis] && node->siblingOffsetPrev != -1) {
					idxOffset = node->siblingOffsetPrev;
					process = true;
				} else if (!dirIsNeg[pNode->axis] && node->siblingOffsetNext != -1) {
					idxOffset = node->siblingOffsetNext;
					process = true;
				} else {
					idxOffset = node->parentOffset;
					process = false;
				}
			} else {
				idxOffset = node->parentOffset;
				process = false;
			}
		}
	}
    PBRT_BVH_INTERSECTIONP_FINISHED();
    return false;
#endif
}

struct CompareContractPoints {
    CompareContractPoints(int d) { dim = d; }
    int dim;
    bool operator()(const std::pair<uint32_t, Point> &a,
                    const std::pair<uint32_t, Point> &b) const {
        return a.second[dim] < b.second[dim];
    }
};

struct CompareContractToBucket {
    CompareContractToBucket(int split, int num, int d, const BBox &b)
        : centroidBounds(b)
    { splitBucket = split; nBuckets = num; dim = d; }
    bool operator()(const std::pair<uint32_t, Point> &p, 
					const std::pair<uint32_t, Point> &q) const {
		int b = nBuckets * ((p.second[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
		int c = nBuckets * ((q.second[dim] - centroidBounds.pMin[dim]) /
            (centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
		if (b == nBuckets) b = nBuckets-1;
		if (c == nBuckets) c = nBuckets-1;
		Assert(b >= 0 && b < nBuckets);
		Assert(c >= 0 && c < nBuckets);
		return b <= c;
	}

    int splitBucket, nBuckets, dim;
    const BBox &centroidBounds;
};


void BVHContractAccel::reorderChild(uint32_t uoffset, LinearBVHContractNode mem[]) {
	LinearBVHContractNode *node = &mem[uoffset];

	if (node->numChild == 0)	return ;

	BBox centroidBounds;
	for (uint32_t offset = node->childOffsetHead; offset != -1; ) {
		LinearBVHContractNode *currChild = &mem[offset];
		Point centroid = .5f * currChild->bounds.pMin + .5f * currChild->bounds.pMax;
		offset = currChild->siblingOffsetNext;
		centroidBounds = Union(centroidBounds, centroid);
	}
	int dim = centroidBounds.MaximumExtent();

	switch (splitMethod) {
	case SPLIT_MIDDLE: case SPLIT_EQUAL_COUNTS: {
		uint32_t offset = node->childOffsetHead;
		vector<std::pair<uint32_t, Point>> A;
		while (offset != -1) {
			LinearBVHContractNode *currChild = &mem[offset];
			A.push_back(std::make_pair(offset, .5f * currChild->bounds.pMin + .5f * currChild->bounds.pMax));
			offset = currChild->siblingOffsetNext;
		}

		std::sort(A.begin(), A.end(), CompareContractPoints(dim));
		node->childOffsetHead = A.front().first;
		node->childOffsetTail = A.back().first;
		mem[A.front().first].siblingOffsetPrev = -1;
		mem[A.back().first].siblingOffsetNext = -1;
		for (uint32_t it = 1; it < A.size(); it++) {
			LinearBVHContractNode *a = &mem[A[it-1].first];
			LinearBVHContractNode *b = &mem[A[it].first];
			a->siblingOffsetNext = A[it].first;
			b->siblingOffsetPrev = A[it-1].first;
		}
	}
	case SPLIT_SAH: default: {
		// Allocate _BucketInfo_ for SAH partition buckets
		const int nBuckets = 12;
		struct BucketInfo {
			BucketInfo() { count = 0; }
			int count;
			BBox bounds;
		};
		BucketInfo buckets[nBuckets];

		uint32_t offset = node->childOffsetHead;
		vector<std::pair<uint32_t, Point>> A;
		while (offset != -1) {
			LinearBVHContractNode *currChild = &mem[offset];
			Point centroid = .5f * currChild->bounds.pMin + .5f * currChild->bounds.pMax;
			int b = nBuckets * 
				((centroid[dim] - centroidBounds.pMin[dim]) /
				(centroidBounds.pMax[dim] - centroidBounds.pMin[dim]));
			if (b == nBuckets) b = nBuckets-1;
			Assert(b >= 0 && b < nBuckets);
			buckets[b].count++;
			buckets[b].bounds = Union(buckets[b].bounds, currChild->bounds);
			A.push_back(std::make_pair(offset, centroid));
			offset = currChild->siblingOffsetNext;
		}

		// Compute costs for splitting after each bucket
		float cost[nBuckets-1];
		for (int i = 0; i < nBuckets-1; ++i) {
			BBox b0, b1;
			int count0 = 0, count1 = 0;
			for (int j = 0; j <= i; ++j) {
				b0 = Union(b0, buckets[j].bounds);
				count0 += buckets[j].count;
			}
			for (int j = i+1; j < nBuckets; ++j) {
				b1 = Union(b1, buckets[j].bounds);
				count1 += buckets[j].count;
			}
			cost[i] = .125f + (count0*b0.SurfaceArea() + count1*b1.SurfaceArea()) /
								node->bounds.SurfaceArea();
		}

		// Find bucket to split at that minimizes SAH metric
		float minCost = cost[0];
		uint32_t minCostSplit = 0;
		for (int i = 1; i < nBuckets-1; ++i) {
			if (cost[i] < minCost) {
				minCost = cost[i];
				minCostSplit = i;
			}
		}

		std::sort(A.begin(), A.end(), CompareContractToBucket(minCostSplit, nBuckets, dim, centroidBounds));
		node->childOffsetHead = A.front().first;
		node->childOffsetTail = A.back().first;
		mem[A.front().first].siblingOffsetPrev = -1;
		mem[A.back().first].siblingOffsetNext = -1;
		for (uint32_t it = 1; it < A.size(); it++) {
			LinearBVHContractNode *a = &mem[A[it-1].first];
			LinearBVHContractNode *b = &mem[A[it].first];
			a->siblingOffsetNext = A[it].first;
			b->siblingOffsetPrev = A[it-1].first;
		}
	}
	}
}
bool BVHContractAccel::contractionCriterionSA(LinearBVHContractNode *node, LinearBVHContractNode *pnode) {
	if (pnode == NULL || node == NULL || node->numChild == 0)
		return false;
	if (pnode->numChild + node->numChild > 16)
		return false;
	float alpha = node->bounds.SurfaceArea() / pnode->bounds.SurfaceArea();
	return alpha > 1.f - 1.f / node->numChild;
}
#define BVH_CONTRACT_RD_THRESHOLD 512
bool BVHContractAccel::contractionCriterionRD(LinearBVHContractNode *node, LinearBVHContractNode *pnode) {
	if (pnode == NULL || node == NULL || node->numChild == 0)
		return false;
	if (pnode->numChild + node->numChild > 16 || node->visitCount < BVH_CONTRACT_RD_THRESHOLD)
		return false;
	float alpha = 1.f * node->visitCount / pnode->visitCount;
	return alpha > 1.f - 1.f / node->numChild;
}
uint32_t BVHContractAccel::flattenLinearBVHTree(uint32_t uoffset, LinearBVHContractNode mem[], uint32_t *offset, uint32_t parentOffset) {
	LinearBVHContractNode *node = &nodes[uoffset];
	LinearBVHContractNode *linearNode = &mem[*offset];
	uint32_t myOffset = (*offset)++;
	*linearNode = *node;
	linearNode->parentOffset = parentOffset;
	uint32_t it = node->childOffsetHead, prevChildOffset = -1, currChildOffset = -1;
	while (it != -1) {
		LinearBVHContractNode *currChild = &nodes[it];
		currChild->parentOffset = uoffset;
		currChildOffset = flattenLinearBVHTree(it, mem, offset, myOffset);
		if (prevChildOffset == -1) {
			linearNode->childOffsetHead = currChildOffset;
		} else {
			mem[prevChildOffset].siblingOffsetNext = currChildOffset;
			mem[currChildOffset].siblingOffsetPrev = prevChildOffset;
		}
		prevChildOffset = currChildOffset;
		it = currChild->siblingOffsetNext;
	}
	linearNode->childOffsetTail = prevChildOffset;
	return myOffset;
}
void BVHContractAccel::recursiveContractSA(uint32_t uoffset) {
	LinearBVHContractNode *node = &nodes[uoffset];
	if (node->nPrimitives > 0)	// is leaf
		return ;
	uint32_t offset = node->childOffsetHead;
	while (offset != -1) {
		LinearBVHContractNode *currChild = &nodes[offset];
		if (contractionCriterionSA(currChild, node)) {
			realNodes--;
			if (currChild->siblingOffsetNext != -1) {
				LinearBVHContractNode *ctail = &nodes[currChild->childOffsetTail];
				ctail->siblingOffsetNext = currChild->siblingOffsetNext;
				LinearBVHContractNode *cnext = &nodes[currChild->siblingOffsetNext];
				cnext->siblingOffsetPrev = currChild->childOffsetTail;
			} else {
				node->childOffsetTail = currChild->childOffsetTail;
			}

			if (currChild->siblingOffsetPrev != -1) {
				LinearBVHContractNode *chead = &nodes[currChild->childOffsetHead];
				chead->siblingOffsetPrev = currChild->siblingOffsetPrev;
				LinearBVHContractNode *cprev = &nodes[currChild->siblingOffsetPrev];
			 	cprev->siblingOffsetNext = currChild->childOffsetHead;
			} else {
				node->childOffsetHead = currChild->childOffsetHead;
			}
			node->numChild = node->numChild + currChild->numChild - 1;
			offset = currChild->childOffsetHead;
		} else {
			offset = currChild->siblingOffsetNext;
		}
	}

	// reorder 
	// reorderChild(uoffset, nodes);

	offset = node->childOffsetHead;
	while (offset != -1) {
		LinearBVHContractNode *currChild = &nodes[offset];
		currChild->parentOffset = uoffset;
		recursiveContractSA(offset);
		offset = currChild->siblingOffsetNext;
	}
}
void BVHContractAccel::recursiveContractRD(uint32_t uoffset) {
	LinearBVHContractNode *node = &nodes[uoffset];
	if (node->nPrimitives > 0)	// is leaf
		return ;
	if (node->visitCount < BVH_CONTRACT_RD_THRESHOLD)
		return ;
	uint32_t offset = node->childOffsetHead;
	while (offset != -1) {
		LinearBVHContractNode *currChild = &nodes[offset];
		if (contractionCriterionRD(currChild, node)) {
			realNodes--;
			if (currChild->siblingOffsetNext != -1) {
				LinearBVHContractNode *ctail = &nodes[currChild->childOffsetTail];
				ctail->siblingOffsetNext = currChild->siblingOffsetNext;
				LinearBVHContractNode *cnext = &nodes[currChild->siblingOffsetNext];
				cnext->siblingOffsetPrev = currChild->childOffsetTail;
			} else {
				node->childOffsetTail = currChild->childOffsetTail;
			}

			if (currChild->siblingOffsetPrev != -1) {
				LinearBVHContractNode *chead = &nodes[currChild->childOffsetHead];
				chead->siblingOffsetPrev = currChild->siblingOffsetPrev;
				LinearBVHContractNode *cprev = &nodes[currChild->siblingOffsetPrev];
			 	cprev->siblingOffsetNext = currChild->childOffsetHead;
			} else {
				node->childOffsetHead = currChild->childOffsetHead;
			}
			node->numChild = node->numChild + currChild->numChild - 1;
			offset = currChild->childOffsetHead;
		} else {
			offset = currChild->siblingOffsetNext;
		}
	}

	offset = node->childOffsetHead;
	while (offset != -1) {
		LinearBVHContractNode *currChild = &nodes[offset];
		currChild->parentOffset = uoffset;
		recursiveContractSA(offset);
		offset = currChild->siblingOffsetNext;
	}
}

BVHContractAccel *CreateBVHContractAccelerator(const vector<Reference<Primitive> > &prims,
        const ParamSet &ps) {
    string splitMethod = ps.FindOneString("splitmethod", "sah");
    uint32_t maxPrimsInNode = ps.FindOneInt("maxnodeprims", 4);
    return new BVHContractAccel(prims, maxPrimsInNode, splitMethod);
}


