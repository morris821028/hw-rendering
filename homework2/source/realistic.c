
#include "stdafx.h"
#include "cameras/realistic.h"
#include <iostream>
#include <fstream>

// Option: Whitted, Heckbert, Other Method for Refraction.
// #define REFRACTION_Whitted
#define REFRACTION_Heckbert
// #define REFRACTION_Other

RealisticCamera::RealisticCamera(const AnimatedTransform &cam2world,
				 float hither, float yon, 
				 float sopen, float sclose, 
				 float filmdistance, float aperture_diameter, string specfile, 
				 float filmdiag, Film *f)
	: Camera(cam2world, sopen, sclose, f) // pbrt-v2 doesnot specify hither and yon
{
	// YOUR CODE HERE -- build and store datastructures representing the given lens
	// and film placement.

	this->hither = hither, this->yon = yon, this->aperture_diameter = aperture_diameter;
	this->filmdistance = filmdistance, this->filmdiag = filmdiag;
	
	// read lens spec file
	this->distToBack = ReadSpecFile(specfile);

	// compute transform: raster to camera, film to camera
	float diag = hypotf(f->xResolution, f->yResolution);
	float scale = filmdiag / diag;
	float X = scale * f->xResolution * 0.5f;
	float Y = scale * f->yResolution * 0.5f;
	RasterToCamera = Translate(Vector(0.f, 0.f, -(this->filmdistance + this->distToBack))) * 
					 Translate(Vector(X, -Y, 0.f)) *
					 Scale(-scale, scale, 1);
}

float RealisticCamera::ReadSpecFile(string fileName) {
	const char comment = '#';

	std::ifstream fin(fileName);
	string line;
	
	lensgroup = vector<Lens>();

	float dist = 0;
	while (std::getline(fin, line)) {
		if (line[0] == comment)
			continue;
		float radius, axpos, N, aperture;
		if (sscanf(line.c_str(), "%f %f %f %f", &radius, &axpos, &N, &aperture) != 4) {
			std::cerr << "SpecFile " << fileName << std::endl;
			std::cerr << "Len format usage: radius, axpos, N, aperture" << std::endl;
			exit(1);
		}

		aperture /= 2.0f;	// diameter to radius
		Lens lens = {
			radius, axpos, (N == 0) ? 1 : N, aperture, 
			N == 0, -dist, Point(0.f, 0.f, -dist - radius)
		};

		lensgroup.push_back(lens);
		dist += axpos;
	}
	backLensAperture = lensgroup.back().aperture;
	printf("#Lens = %d\n", lensgroup.size());
	fin.close();
	return dist;
}

float RealisticCamera::GenerateRay(const CameraSample &sample, Ray *ray) const {
	// Generate raster and back lens samples
	Point Praster(sample.imageX, sample.imageY, 0.f);
	Point Pcamera;
	RasterToCamera(Praster, &Pcamera);

	const Lens &backLens = lensgroup.back();
	float lensU, lensV, lensZ;
	// Get Sample Pback (lensU, lensV, -distToBack)
	ConcentricSampleDisk(sample.lensU, sample.lensV, &lensU, &lensV);
	lensU *= backLensAperture;
	lensV *= backLensAperture;
	lensZ = -distToBack;

	// camera -> inner lens -> outer lens
	// sphere coordinate x^2 + y^2 + z^2 = r^2

	// hit point on the surface of the inner lens
	ray->o = Pcamera;
	ray->d = Normalize(Point(lensU, lensV, lensZ) - Pcamera);
	ray->time = Lerp(sample.time, shutterOpen, shutterClose);
	ray->mint = 0.f;
		
	// iterator over the lens components from inner to outer
	for (int i = lensgroup.size() - 1; i >= 0; i--) {
		if (lensgroup[i].end) {
			// ray run pass aperture (not lens)
			float deltaZ = lensgroup[i].z - ray->o.z;
			ray->o = ray->o + (deltaZ / ray->d.z) * ray->d;
			if (lensgroup[i].outOfRange(ray->o.x, ray->o.y))
				return 0.f;
			
		} else {
			float n = (i == 0) ? 1.f : lensgroup[i-1].n;
			if (!lensgroup[i].Refraction(ray, n))
				return 0.f;
		}
	}	
	
	// set exposure weight
	float cosTheta = ray->d.z;
	float Z = filmdistance + ray->o.z;
	float weight = (backLens.aperture * backLens.aperture * M_PI) / ( Z * Z);
	weight = weight * cosTheta * cosTheta * cosTheta * cosTheta;

	ray->maxt = (yon - hither) / ray->d.z;
	CameraToWorld(*ray, ray);
	ray->d = Normalize(ray->d);
	return weight;
}

bool Lens::Refraction(Ray *ray, float n2) const {
	// I: incident vector I (normalize), 
	// C: incident circle center
	// Assume after t unit time, the point will arrive the surface,
	// then get equation below:
	// \overrightarrow{OC} + \overrightarrow{I} \times t = \overrightarrow{OP}, and |\overrightarrow{OP}| = \text{radius}
	// |\overrightarrow{OC} + \overrightarrow{I} \times t| = \overrightarrow{OP}
	// |\overrightarrow{I}|^2 t^2 + 2 \overrightarrow{OC} \cdot \overrightarrow{I} t + |\overrightarrow{OC}|^2 - \text{radius}^2 = 0
	// quadratic equation about t.
	Vector I = ray->d;
	Vector oc = ray->o - o;
	float b = I.x * oc.x + I.y * oc.y + I.z * oc.z;
	float c = oc.x * oc.x + oc.y * oc.y + oc.z * oc.z - radius * radius;
	float discrim = b * b - c, t;
	if (discrim <= 0.f)	return false;
	float rootDiscrim = sqrtf(discrim);

	if (radius > 0.f)
		t = -b + rootDiscrim;
	else
		t = -b - rootDiscrim;

	// move hit point to the surface of refraction indices
	ray->o = ray->o + t * I;
	if (outOfRange(ray->o.x, ray->o.y))
		return false;

	Vector N = Normalize(o - ray->o);
	if (radius < 0.f)
		N = -N;
	// Whitted's Method
#ifdef REFRACTION_Whitted
	float eta = n2 / n;
	Vector I2 = I / (- Dot(I, N));
	Vector J = I2 + N;
	float alpha = eta * eta * Dot(I2, I2) - Dot(J, J);
	if (alpha < 0.f)
		return false;
	alpha = 1.f / sqrt(alpha);
	Vector T2 = alpha * J - N;
	ray->d = T2 / T2.Length();
	return true;
#endif
	// Heckbert's Method
#ifdef REFRACTION_Heckbert
	float eta = n / n2;
	float c1 = - Dot(I, N);
	float c2 = 1.f - eta * eta * (1.f - c1 * c1);
	if (c2 < 0.f)
		return false;
	c2 = sqrtf(c2);
	ray->d = eta * I + (eta * c1 - c2) * N;
	return true;
#endif
	// Other Method
#ifdef REFRACTION_Other
	float eta = n2 / n;
	float c1 = - Dot(I, N);
	float beta = eta * eta - 1 + c1 * c1;
	if (beta < 0.f)
		return false;
	beta = c1 - sqrtf(beta);
	ray->d = (I + beta * N) / eta;
	return true;
#endif
}

RealisticCamera *CreateRealisticCamera(const ParamSet &params,
        const AnimatedTransform &cam2world, Film *film) {
	// Extract common camera parameters from \use{ParamSet}
	float hither = params.FindOneFloat("hither", -1);
	float yon = params.FindOneFloat("yon", -1);
	float shutteropen = params.FindOneFloat("shutteropen", -1);
	float shutterclose = params.FindOneFloat("shutterclose", -1);

	// Realistic camera-specific parameters
	string specfile = params.FindOneString("specfile", "");
	float filmdistance = params.FindOneFloat("filmdistance", 70.0); // about 70 mm default to film
 	float fstop = params.FindOneFloat("aperture_diameter", 1.0);	
	float filmdiag = params.FindOneFloat("filmdiag", 35.0);

	Assert(hither != -1 && yon != -1 && shutteropen != -1 &&
		shutterclose != -1 && filmdistance!= -1);
	if (specfile == "") {
	    Severe( "No lens spec file supplied!\n" );
	}
	return new RealisticCamera(cam2world, hither, yon,
				   shutteropen, shutterclose, filmdistance, fstop, 
				   specfile, filmdiag, film);
}

