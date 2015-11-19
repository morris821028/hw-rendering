
#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_CAMERAS_REALISTIC_H
#define PBRT_CAMERAS_REALISTIC_H

#include "camera.h"
#include "paramset.h"
#include "film.h"

struct Lens {
	float radius, axpos, n, aperture;
	bool end;
	float z;
	bool Refraction(Ray *ray, float n2) const;
	bool outOfRange(float x, float y) const {
		return x * x + y * y > aperture * aperture;
	}
};

// RealisticCamera Declarations
class RealisticCamera : public Camera {
public:
	// RealisticCamera Public Methods
	RealisticCamera(const AnimatedTransform &cam2world,
						float hither, float yon, float sopen,
						float sclose, float filmdistance, float aperture_diameter, string specfile,
						float filmdiag, Film *film);
	float GenerateRay(const CameraSample &sample, Ray *) const;
  
private:
	// RealisticCamera Private Member
	float hither, yon, aperture_diameter, filmdistance, filmdiag;
	float distToBack, backLensAperture;
	vector<Lens> lensgroup;
	Transform RasterToCamera;
	// RealisticCamera Private Methods
	float ReadSpecFile(string fileName);
};


RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film);


#endif	// PBRT_CAMERAS_REALISTIC_H
