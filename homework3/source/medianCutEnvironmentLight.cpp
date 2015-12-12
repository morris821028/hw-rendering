
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


// lights/medianCutEnvironmentLight.cpp*
#include "stdafx.h"
#include "lights/medianCutEnvironmentLight.h"
#include "sh.h"
#include "montecarlo.h"
#include "paramset.h"
#include "imageio.h"

#define DEBUG
// MedianCutEnvironmentLight Utility Classes
struct MedianCutEnvironmentLightCube {
    // MedianCutEnvironmentLight Public Methods
    MedianCutEnvironmentLightCube(const MedianCutEnvironmentLight *l, const Scene *s,
                     float t, bool cv, float pe)
        : light(l), scene(s), time(t), pEpsilon(pe), computeVis(cv) { }
    Spectrum operator()(int, int, const Point &p, const Vector &w) {
        Ray ray(p, w, pEpsilon, INFINITY, time);
        if (!computeVis || !scene->IntersectP(ray))
            return light->Le(RayDifferential(ray));
        return 0.f;
    }
    const MedianCutEnvironmentLight *light;
    const Scene *scene;
    float time, pEpsilon;
    bool computeVis;
};



// InfiniteAreaLight Method Definitions
MedianCutEnvironmentLight::~MedianCutEnvironmentLight() {
    delete distribution;
    delete radianceMap;
}


MedianCutEnvironmentLight::MedianCutEnvironmentLight(const Transform &light2world,
        const Spectrum &L, int ns, const string &texmap)
    : Light(light2world, ns) {
    int width = 0, height = 0;
    RGBSpectrum *texels = NULL;
    // Read texel data from _texmap_ into _texels_
    if (texmap != "") {
        texels = ReadImage(texmap, &width, &height);
        if (texels)
            for (int i = 0; i < width * height; ++i)
                texels[i] *= L.ToRGBSpectrum();
    }
    if (!texels) {
        width = height = 1;
        texels = new RGBSpectrum[1];
        texels[0] = L.ToRGBSpectrum();
    }
    radianceMap = new MIPMap<RGBSpectrum>(width, height, texels);
	
    // Initialize sampling PDFs for environment light

	// Remember to scale the light intensity with the areas (solid angles)
	float solidAngleScale = ((2.f * M_PI) / (width)) * ((M_PI) / (height));
	for (int v = 0; v < height; v++) {
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
		for (int u = 0; u < width; u++)
			texels[u+v*width] = texels[u+v*width] * solidAngleScale * sinTheta;
	}

	// Compute scalar-valued image _img_ from environment map
    float filter = 1.f / max(width, height);
    float *img = new float[width*height];
    for (int v = 0; v < height; ++v) {
        float vp = (float)v / (float)height;
        float sinTheta = sinf(M_PI * float(v+.5f)/float(height));
        for (int u = 0; u < width; ++u) {
            float up = (float)u / (float)width;
            img[u+v*width] = radianceMap->Lookup(up, vp, filter).y();
            img[u+v*width] *= sinTheta;
        }
    }

#ifdef DEBUG
	printf("Morris working ! ------------------------------------\n");
	printf("width %d height %d\n", width, height);
	fflush(stdin);
#endif

	float *sumTable = new float[width*height];
	float *columnSum = new float[width];
	for (int u = 0; u < width; u++)
		columnSum[u] = 0;
#define VERT(u, v) ((u)+(v)*width)
	for (int v = 0; v < height; v++) {
		float leftsum = 0;
		for (int u = 0; u < width; u++) {
			leftsum += img[VERT(u, v)] * solidAngleScale;
			if (v > 0)
				sumTable[VERT(u, v)] = sumTable[VERT(u, v-1)] + leftsum;
			else
				sumTable[VERT(u, v)] = leftsum;
		}
	}
#ifdef DEBUG
	/*
	for (int v = 0; v < height; v++) {
		for (int u = 0; u < width; u++)
			printf("%f ", sumTable[VERT(u, v)]);
		puts("");
	}
	fflush(stdin);
	*/
#endif
	struct Region {
		int lu, lv, ru, rv;
		Region(int lu = 0, int lv = 0, int ru = 0, int rv = 0): 
			lu(lu), lv(lv), ru(ru), rv(rv) {}
		float getEnergy(float sumTable[], int width, int height) {
			float e = sumTable[VERT(ru, rv)];
			if (lu)	e -= sumTable[VERT(lu-1, rv)];
			if (lv)	e -= sumTable[VERT(ru, lv-1)];
			if (lu && lv)	e += sumTable[VERT(lu-1, lv-1)];
			return e;
		}
		int getAxis() {
			return ru - lu > rv - lv ? 0 : 1;
		}
	};
#undef VERT
	
	vector<Region> Rs;
	Rs.push_back(Region(0, 0, width-1, height-1));

	// A light probe image is subdivided into 64 equal energy regions
	int partitions_regions = 64;
	for (int it = 0; it < 8 && Rs.size() < partitions_regions; it ++) {
		vector<Region> nextRs;
		for (Region region : Rs) {
			float half_energy = region.getEnergy(sumTable, width, height) * 0.5f;
			int div_axis = region.getAxis();

			// printf("Region %d: %d %d %d %d %f\n", it, region.lu, region.lv, region.ru, region.rv, region.getEnergy(sumTable, width, height));

			int lu = region.lu, lv = region.lv, ru = region.ru, rv = region.rv;
			if (div_axis == 0) {
				int l = lu, r = ru, mid, ret = lu;
				while (l <= r) {
					mid = (l + r)/2;
					// printf("binary x %d %d %d %f\n", l, r, mid, Region(lu, lv, mid, rv).getEnergy(sumTable, width, height));
					if (Region(lu, lv, mid, rv).getEnergy(sumTable, width, height) < half_energy)
						l = mid + 1, ret = mid;
					else
						r = mid - 1;
				}
				nextRs.push_back(Region(lu, lv, ret, rv));
				if (ret+1 <= ru)
					nextRs.push_back(Region(ret+1, lv, ru, rv));
			} else {
				int l = lv, r = rv, mid, ret = lv;
				while (l <= r) {
					mid = (l + r)/2;
					// printf("binary x %d %d %d %f\n", l, r, mid, Region(lu, lv, ru, mid).getEnergy(sumTable, width, height));
					if (Region(lu, lv, ru, mid).getEnergy(sumTable, width, height) < half_energy)
						l = mid + 1, ret = mid;
					else
						r = mid - 1;
				}
				nextRs.push_back(Region(lu, lv, ru, ret));
				if (ret+1 <= rv)
					nextRs.push_back(Region(lu, ret+1, ru, rv));
			}
		}
		Rs = nextRs;
	}

	// computing : A point light is created for each region at its centroid
	float v_scale = 1.f / height, u_scale = 1.f / width;
	this->PDF = 1.f / Rs.size();
	this->VLRAs = vector<VLRA>(Rs.size());
	for (int it = 0; it < Rs.size(); it++) {
		RGBSpectrum spectrum = RGBSpectrum(0.f);
		Region region = Rs[it];
		float cv = 0.f, cu = 0.f, sumf = 0;

#define VERT(u, v) ((u)+(v)*width)
		for (int v = region.lv; v <= region.rv; v++) {
			for (int u = region.lu; u <= region.ru; u++) {
				spectrum += texels[VERT(u, v)];
				float f = img[VERT(u, v)];
				cv += v * f, cu += u * f, sumf += f;
			}
		}
#undef VERT
		this->VLRAs[it] = VLRA(cu / sumf * v_scale, cv / sumf * u_scale, spectrum);
	}

	for (int it = 0; it < Rs.size(); it++) {
		printf("Region %d: %d %d %d %d %f\n", it, Rs[it].lu, Rs[it].lv, Rs[it].ru, Rs[it].rv, Rs[it].getEnergy(sumTable, width, height));
	}
	fflush(stdin);
	
    delete[] texels;

    // Compute sampling distributions for rows and columns of image
    distribution = new Distribution2D(img, width, height);
    delete[] img;
}


Spectrum MedianCutEnvironmentLight::Power(const Scene *scene) const {
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    return M_PI * worldRadius * worldRadius *
        Spectrum(radianceMap->Lookup(.5f, .5f, .5f), SPECTRUM_ILLUMINANT);
}


Spectrum MedianCutEnvironmentLight::Le(const RayDifferential &r) const {
    Vector wh = Normalize(WorldToLight(r.d));
    float s = SphericalPhi(wh) * INV_TWOPI;
    float t = SphericalTheta(wh) * INV_PI;
    return Spectrum(radianceMap->Lookup(s, t), SPECTRUM_ILLUMINANT);
}


void MedianCutEnvironmentLight::SHProject(const Point &p, float pEpsilon,
        int lmax, const Scene *scene, bool computeLightVis,
        float time, RNG &rng, Spectrum *coeffs) const {
    // Project _InfiniteAreaLight_ to SH using Monte Carlo if visibility needed
    if (computeLightVis) {
        Light::SHProject(p, pEpsilon, lmax, scene, computeLightVis,
                         time, rng, coeffs);
        return;
    }
    for (int i = 0; i < SHTerms(lmax); ++i)
        coeffs[i] = 0.f;
    int ntheta = radianceMap->Height(), nphi = radianceMap->Width();
    if (min(ntheta, nphi) > 50) {
        // Project _InfiniteAreaLight_ to SH from lat-long representation

        // Precompute $\theta$ and $\phi$ values for lat-long map projection
        float *buf = new float[2*ntheta + 2*nphi];
        float *bufp = buf;
        float *sintheta = bufp;  bufp += ntheta;
        float *costheta = bufp;  bufp += ntheta;
        float *sinphi = bufp;    bufp += nphi;
        float *cosphi = bufp;
        for (int theta = 0; theta < ntheta; ++theta) {
            sintheta[theta] = sinf((theta + .5f)/ntheta * M_PI);
            costheta[theta] = cosf((theta + .5f)/ntheta * M_PI);
        }
        for (int phi = 0; phi < nphi; ++phi) {
            sinphi[phi] = sinf((phi + .5f)/nphi * 2.f * M_PI);
            cosphi[phi] = cosf((phi + .5f)/nphi * 2.f * M_PI);
        }
        float *Ylm = ALLOCA(float, SHTerms(lmax));
        for (int theta = 0; theta < ntheta; ++theta) {
            for (int phi = 0; phi < nphi; ++phi) {
                // Add _InfiniteAreaLight_ texel's contribution to SH coefficients
                Vector w = Vector(sintheta[theta] * cosphi[phi],
                                  sintheta[theta] * sinphi[phi],
                                  costheta[theta]);
                w = Normalize(LightToWorld(w));
                Spectrum Le = Spectrum(radianceMap->Texel(0, phi, theta),
                                       SPECTRUM_ILLUMINANT);
                SHEvaluate(w, lmax, Ylm);
                for (int i = 0; i < SHTerms(lmax); ++i)
                    coeffs[i] += Le * Ylm[i] * sintheta[theta] *
                        (M_PI / ntheta) * (2.f * M_PI / nphi);
            }
        }

        // Free memory used for lat-long theta and phi values
        delete[] buf;
    }
    else {
        // Project _InfiniteAreaLight_ to SH from cube map sampling
        SHProjectCube(MedianCutEnvironmentLightCube(this, scene, time, computeLightVis,
                                       pEpsilon),
                      p, 200, lmax, coeffs);
    }
}


MedianCutEnvironmentLight *CreateMedianCutEnvironmentLight(const Transform &light2world,
        const ParamSet &paramSet) {
    Spectrum L = paramSet.FindOneSpectrum("L", Spectrum(1.0));
    Spectrum sc = paramSet.FindOneSpectrum("scale", Spectrum(1.0));
    string texmap = paramSet.FindOneFilename("mapname", "");
    int nSamples = paramSet.FindOneInt("nsamples", 1);
    if (PbrtOptions.quickRender) nSamples = max(1, nSamples / 4);
    return new MedianCutEnvironmentLight(light2world, L * sc, nSamples, texmap);
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Point &p, float pEpsilon,
        const LightSample &ls, float time, Vector *wi, float *pdf,
        VisibilityTester *visibility) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();

	// When pbrt asks a sample from environment light, uniformly 
	// select one from all lights and return its direction, intensity and PDF
	VLRA vlra = VLRAs[Floor2Int(ls.uComponent * VLRAs.size())];
	
    // Convert infinite light sample point to direction
    float theta = vlra.u * M_PI, phi = vlra.v * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    *wi = LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                              costheta));

    // Compute PDF for sampled MedianCut environment light direction
    *pdf = this->PDF;

    // Return radiance value for MedianCut environment light direction
    visibility->SetRay(p, pEpsilon, *wi, time);
	Spectrum Ls = Spectrum(vlra.spectrum,
                           SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


float MedianCutEnvironmentLight::Pdf(const Point &, const Vector &w) const {
    PBRT_INFINITE_LIGHT_STARTED_PDF();
    Vector wi = WorldToLight(w);
    float theta = SphericalTheta(wi), phi = SphericalPhi(wi);
    float sintheta = sinf(theta);
    if (sintheta == 0.f) return 0.f;
    float p = distribution->Pdf(phi * INV_TWOPI, theta * INV_PI) /
           (2.f * M_PI * M_PI * sintheta);
    PBRT_INFINITE_LIGHT_FINISHED_PDF();
    return p;
}


Spectrum MedianCutEnvironmentLight::Sample_L(const Scene *scene,
        const LightSample &ls, float u1, float u2, float time,
        Ray *ray, Normal *Ns, float *pdf) const {
    PBRT_INFINITE_LIGHT_STARTED_SAMPLE();
    // Compute direction for infinite light sample ray

    // Find $(u,v)$ sample coordinates in infinite light texture
    float uv[2], mapPdf;
    distribution->SampleContinuous(ls.uPos[0], ls.uPos[1], uv, &mapPdf);
    if (mapPdf == 0.f) return Spectrum(0.f);

    float theta = uv[1] * M_PI, phi = uv[0] * 2.f * M_PI;
    float costheta = cosf(theta), sintheta = sinf(theta);
    float sinphi = sinf(phi), cosphi = cosf(phi);
    Vector d = -LightToWorld(Vector(sintheta * cosphi, sintheta * sinphi,
                                    costheta));
    *Ns = (Normal)d;

    // Compute origin for infinite light sample ray
    Point worldCenter;
    float worldRadius;
    scene->WorldBound().BoundingSphere(&worldCenter, &worldRadius);
    Vector v1, v2;
    CoordinateSystem(-d, &v1, &v2);
    float d1, d2;
    ConcentricSampleDisk(u1, u2, &d1, &d2);
    Point Pdisk = worldCenter + worldRadius * (d1 * v1 + d2 * v2);
    *ray = Ray(Pdisk + worldRadius * -d, d, 0., INFINITY, time);

    // Compute _InfiniteAreaLight_ ray PDF
    float directionPdf = mapPdf / (2.f * M_PI * M_PI * sintheta);
    float areaPdf = 1.f / (M_PI * worldRadius * worldRadius);
    *pdf = directionPdf * areaPdf;
    if (sintheta == 0.f) *pdf = 0.f;
    Spectrum Ls = (radianceMap->Lookup(uv[0], uv[1]), SPECTRUM_ILLUMINANT);
    PBRT_INFINITE_LIGHT_FINISHED_SAMPLE();
    return Ls;
}


