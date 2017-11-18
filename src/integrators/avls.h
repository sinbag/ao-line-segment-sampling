#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_AVLS_H
#define PBRT_INTEGRATORS_AVLS_H

// integrators/avls.h*
#include "pbrt.h"
#include "integrator.h"
#include "camera.h"
#include "primitive.h"
#include "shape.h"
#include "shapes/sphere.h"
#include "scene.h"

namespace pbrt {
    
    // DirectLightingIntegrator Declarations
    class AVLSIntegrator : public SamplerIntegrator {
    public:
        // DirectLightingIntegrator Public Methods
        AVLSIntegrator(  std::shared_ptr<const Camera> camera,
                       std::shared_ptr<Sampler> sampler,
                       const Bounds2i &pixelBounds,
                       int numPhiSamples,
                       int numThetaSamples,
                       bool cosineWeightedSampling,
                       bool SolidAngleUniformSampling,
                       bool SphericalCoordinatesUniformSampling,
                       bool visualizeIntegrand,
                       bool computeFourierSpectrum,
                       float lineSegmentLength)
        
        :   camera(camera),
        SamplerIntegrator(camera, sampler, pixelBounds),
        m_numPhiSamples(numPhiSamples),
        m_numThetaSamples(numThetaSamples),
        cosineWeightedSampling(cosineWeightedSampling),
        SolidAngleUniformSampling(SolidAngleUniformSampling),
        SphericalCoordinatesUniformSampling(SphericalCoordinatesUniformSampling),
        visualizeIntegrand(visualizeIntegrand),
        computeFourierSpectrum(computeFourierSpectrum),
        lineSegmentLength(lineSegmentLength)
        { m_totalSamples = m_numPhiSamples*m_numThetaSamples;};
        
        Spectrum Li(const RayDifferential &ray, const Scene &scene,
                    Sampler &sampler, MemoryArena &arena, int depth) const;
        
        /*Extract all the spheres from the scene and adds them to a vector*/
        void Preprocess(const Scene &scene, Sampler &sampler) {
            int num_spheres = 0;
            for(std::shared_ptr<Primitive> prim: scene.primitives) {
                Primitive *raw_prim = prim.get();
                if(typeid(*raw_prim) == typeid(GeometricPrimitive))
                {
                    GeometricPrimitive geo_prim = dynamic_cast<GeometricPrimitive&>(*raw_prim);
                    Shape *raw_shape = geo_prim.shape.get();
                    if(typeid(*raw_shape) == typeid(Sphere))
                    {
                        Sphere s = dynamic_cast<Sphere&>(*raw_shape);
                        spheres.push_back(s);
                        num_spheres += 1;
                    }
                }
            };
            
            sampler.Request2DArray(m_totalSamples);
        };
        
    private:
        std::shared_ptr<const Camera> camera;
        
        //AVLIntegrator Private Methods
        static Vector3f sphericalDirection(float theta, float phi);
        
        static Vector3f squareToCosineHemisphere(float u, float v);
        
        static Point2f squareToUniformDisk(float u, float v);
        
        static Vector3f squareToUniformHemisphere(float u, float v);
        
        static float distanceToRay(Point3f p, Ray r);
        
        static double getBlurredValue(double phi, double theta, double length,
                                      SurfaceInteraction *isect,
                                      std::vector<Sphere> spheres,
                                      bool cosineWeightedIntegrand=true);
        
        static bool compareIntervals(const std::pair<float,float> &i1, const std::pair<float,float> &is);
        
        static std::vector<std::pair<float,float>> combineIntervals(std::vector<std::pair<float,float>> intervals);
        
        static std::vector<std::pair<float,float>> getCappedRegions(
                                                                    float startTheta , float stopTheta,
                                                                    std::pair<float,float> integralInterval);
        
        static std::vector<std::pair<float,float>> getComplimentOfIntervals(
                                                                            std::vector<std::pair<float,float>> intervals);
        
        static double calculateVisibility(
                                          std::vector<std::pair<float,float>> occludedIntervals, double angleScale);
        
        // AVLIntegrator Private Data
        std::vector<Sphere> spheres;
        int m_numPhiSamples, m_numThetaSamples, m_totalSamples;
        bool cosineWeightedSampling, SolidAngleUniformSampling, SphericalCoordinatesUniformSampling;
        bool visualizeIntegrand, computeFourierSpectrum;
        float lineSegmentLength;
    };
    
    AVLSIntegrator *CreateAVLSIntegrator(
                                         const ParamSet &params, std::shared_ptr<Sampler> sampler,
                                         std::shared_ptr<const Camera> camera);
    
}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_AVLS_H
