// integrators/avl.cpp*
#include "api.h"
#include "core/reflection.h"
#include "integrators/avl.h"
#include "interaction.h"
#include "paramset.h"
#include "primitive.h"
#include "scene.h"
#include "shape.h"
#include "shapes/sphere.h"
#include "camera.h"
#include "film.h"
#include "stats.h"
#include <random>
#include <cmath>

namespace pbrt {
    
    // DirectLightingIntegrator Method Definitions
    Spectrum AVLIntegrator::Li(const RayDifferential &ray,
                               const Scene &scene, Sampler &sampler,
                               MemoryArena &arena, int depth) const {
        // Find closest ray intersection or return background radiance
        SurfaceInteraction isect;
        if (!scene.Intersect(ray, &isect)) {
            return Spectrum(1.0f);
        }
        
        if (typeid(*isect.shape) != typeid(Sphere))
        {
            std::cout<<"Scene should only contain spheres"<<std::endl;
            return Spectrum(0.0f);
        }
        
        ///Homogenization:
        //@sinbag=================================
        ///Code added by Gurprit Singh [2016-10-10]
        /// New random seed is created for each tile
        /// which changes with tile coordinates
        ///
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0, 1);
        float randomShift = dis(gen);
        
        //@sinbag=================================
        
        /****/
        ///Added by Gurprit Singh
        ///To sample Phi directly from the built-in sampler
        
        const Float *phiSampleArray = sampler.Get1DArray(m_numPhiSamples);
        
        //Get the hit sphere information to compare with spheres in scene
        const Sphere *hit_sphere = static_cast<const Sphere*> (isect.shape);
        Point3f hit_sphere_center = (*hit_sphere->ObjectToWorld)(Point3f(0,0,0));
        float hit_sphere_radius = hit_sphere->radius;
        
        float sum = 0.f;
        // constants to optionally account for cosine-weighted integration
        float angleScale = cosineWeighted ? 2.f : 1.f;
        float integralScale = cosineWeighted ? 0.25f : 1.0f;
        float maxIntegral = cosineWeighted ? 1.f : 2.f;
        
        for (int k = 0; k < m_numPhiSamples; k++) {
            ///Homogenizing the samples by adding a
            ///random shift vector to each sample.
            float u = phiSampleArray[k] + randomShift;
            
            if(u > 1)
                u = u - 1.f;
            else if(u < 0)
                u = 1.f + u;
            
            float phi = u*Pi;
            
            // For each line sample, calculate the intervals that are occluded
            std::vector<std::pair<float, float>> occludedRegions;
            std::vector<std::pair<float, float>> visibleRegions;
            
            // plane normal is the vector in the tangent plane perpendicular to phi
            BSDF *bsdf = new BSDF(isect);
            Vector3f planeNormal = sphericalDirection(PiOver2, phi + PiOver2);
            Vector3f bitangent = sphericalDirection(PiOver2, phi);
            planeNormal = bsdf->LocalToWorld(planeNormal);
            bitangent = bsdf->LocalToWorld(bitangent);
            
            for(Sphere s: spheres) {
                Point3f sphere_center = (*s.ObjectToWorld)(Point3f(0,0,0));
                float sphere_radius = s.radius;
                
                if ( sphere_center == hit_sphere_center
                    && sphere_radius == hit_sphere_radius) continue;
                
                // Check if an intersection occurs between the plane and the sphere
                float plane_to_sphere_center_distance = Dot((sphere_center - isect.p),planeNormal);
                if(fabsf(plane_to_sphere_center_distance) > sphere_radius-MachineEpsilon) continue;
                
                // calculate the circle formed by the intersection
                float circle_radius = sqrt(sphere_radius*sphere_radius - plane_to_sphere_center_distance*plane_to_sphere_center_distance);
                Point3f circle_center = sphere_center - plane_to_sphere_center_distance*planeNormal;
                Vector3f hp_to_circle = circle_center - isect.p;
                float hp_to_circle_center_distance = hp_to_circle.Length();
                hp_to_circle /= hp_to_circle_center_distance;
                
                // shouldn't fail unless we are inside another sphere
                assert(circle_radius <= hp_to_circle_center_distance);
                
                float thetaRange = std::asin(std::min(1.0f,circle_radius/hp_to_circle_center_distance));
                
                // Calculate the angle between the shading normal and hp to circle
                float centerTheta = acosf(std::max(-1.0f, std::min(1.0f,Dot(hp_to_circle, isect.n))));
                
                float startTheta = centerTheta - thetaRange;
                float stopTheta = std::min(PiOver2, centerTheta + thetaRange);
                
                if (thetaRange != thetaRange) std::cout << "thetaRange: " << circle_radius << "," << hp_to_circle_center_distance <<std::endl;
                if (centerTheta != centerTheta) std::cout << "centerTheta: " << Dot(hp_to_circle, isect.n) << std::endl;
                if (startTheta != startTheta) std::cout << "startTheta" <<std::endl;
                if (stopTheta != stopTheta) std::cout << "stopTheta: " <<std::endl;
                
                if (startTheta > PiOver2) {
                    // occluded is completely below horizon
                    continue;
                }
                
                if (Dot(hp_to_circle, bitangent) < 0.0f) {
                    float temp = startTheta;
                    startTheta = -stopTheta;
                    stopTheta = -temp;
                }
                
                // This is necessary for sorting
                if ( stopTheta < startTheta) std::cout << "startTheta should be less than stopTheta" << std::endl;
                
                std::pair<float, float> region(startTheta, stopTheta);
                occludedRegions.push_back(region);
            }
            
            occludedRegions = AVLIntegrator::combineIntervals(occludedRegions);
            visibleRegions = AVLIntegrator::getComplimentOfIntervals(occludedRegions);
            double perSampleRadiance = AVLIntegrator::calculateVisibility(visibleRegions,
                                                                          integralScale,
                                                                          angleScale);
            
            //perSampleColorData[k] = perSampleRadiance/maxIntegral;
            sum += perSampleRadiance;
        }
        //    for(int i = 0; i < m_numPhiSamples; i++)
        //        std::cout << i << " "<< perSampleColorData[i] << std::endl;
        return Spectrum((sum/m_numPhiSamples/maxIntegral));
    }
    
    // Utility Functions
    Vector3f AVLIntegrator::sphericalDirection(float theta, float phi) {
        float sinTheta, cosTheta, sinPhi, cosPhi;
        
        //__sincosf(theta,&sinTheta, &cosTheta);
        //__sincosf(phi, &sinPhi, &cosPhi);
        sinTheta = std::sin(theta);
        cosTheta = std::cos(theta);
        sinPhi = std::sin(phi);
        cosPhi = std::cos(phi);
        
        return Vector3f(
                        sinTheta * cosPhi,
                        sinTheta * sinPhi,
                        cosTheta );}
    
    Vector3f AVLIntegrator::squareToCosineHemisphere(float u, float v) {
        Point2f p = AVLIntegrator::squareToUniformDisk(u,v);
        float z = std::sqrt(std::max((float) 0,
                                     1.0f - p.x*p.x- p.y*p.y));
        
        return Vector3f(p.x, p.y, z);
    }
    
    Point2f AVLIntegrator::squareToUniformDisk(float u, float v) {
        float r = std::sqrt(v);
        float sinPhi, cosPhi;
        //__sincosf(2.0f * Pi * u, &sinPhi, &cosPhi);
        sinPhi = std::sin(2.0f * Pi * u);
        cosPhi = std::cos(2.0f * Pi * u);
        
        return Point2f(
                       cosPhi * r,
                       sinPhi * r
                       );
    }
    
    Vector3f AVLIntegrator::squareToUniformHemisphere(float u, float v) {
        float cosTheta = v;
        float sinTheta = std::sqrt(std::max((float) 0, 1-cosTheta*cosTheta));
        
        float sinPhi, cosPhi;
        //__sincosf(2.0f * Pi * u, &sinPhi, &cosPhi);
        sinPhi = std::sin(2.0f * Pi * u);
        cosPhi = std::cos(2.0f * Pi * u);
        
        
        return Vector3f(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
    }
    
    bool AVLIntegrator::compareIntervals(const std::pair<float,float> &i1, const std::pair<float,float> &i2) {
        return i1.first < i2.first;
    }
    
    std::vector<std::pair<float,float>> AVLIntegrator::combineIntervals(std::vector<std::pair<float,float>> intervals) {
                                                                        
        
        if (intervals.size() <= 1) return intervals;
        
        // sort step
        std::sort(intervals.begin(), intervals.end(), compareIntervals);
        
        // merge step
        std::vector<std::pair<float,float>> disjointIntervals;
        disjointIntervals.push_back(intervals[0]);
        for (int i=1; i<intervals.size(); i++) {
            std::pair<float,float> nextInterval = intervals[i];
            std::pair<float,float> currInterval = disjointIntervals.back();
            
            if (currInterval.second < nextInterval.first) {
                // intervals don't overlap, add next interval
                disjointIntervals.push_back(nextInterval);
            } else if (currInterval.second < nextInterval.second) {
                // intervals overlap and the current interval doesn't contain the next.
                // Remove current and add a combined interval
                currInterval.second = nextInterval.second;
                disjointIntervals.pop_back();
                disjointIntervals.push_back(currInterval);
            }
        }
        return disjointIntervals;
    }
    
    /* getComplimentOfIntervals -- returns the compliment of intervals on the range -pi/2 to pi/2. For example
     * given the vector of intervals [(-pi/4,0), (pi/4,pi/2)] would return [(-pi/2,-pi/4),(0,pi4)]
     *
     */
    std::vector<std::pair<float,float>> AVLIntegrator::getComplimentOfIntervals(std::vector<std::pair<float,float>> intervals) {
        std::vector<std::pair<float,float>> complimentIntervals;
        //make sure there is at least one interval
        if (intervals.size() <= 0) {
            std::pair<float,float> new_interval(-PiOver2,PiOver2);
            complimentIntervals.push_back(new_interval);
            return complimentIntervals;
        }
        
        //Add interval from -pi/2 to start of first interval
        if (intervals.front().first > -PiOver2) {
            std::pair<float,float> new_interval(-PiOver2,intervals.front().first);
            complimentIntervals.push_back(new_interval);
        }
        
        //Add interval from end of last interval to pi/2
        if (intervals.back().second < PiOver2) {
            std::pair<float,float> new_interval(intervals.back().second,PiOver2);
            complimentIntervals.push_back(new_interval);
        }
        
        for(int i=1; i<intervals.size(); i++) {
            std::pair<float,float> new_interval(intervals.at(i-1).second, intervals.at(i).first);
            complimentIntervals.push_back(new_interval);
        }
        return complimentIntervals;
    }
    
    float AVLIntegrator::calculateVisibility(std::vector<std::pair<float,float>> occludedIntervals,
                                             float integralScale,
                                             float angleScale) {
        float amountOccluded = 0;
        for (std::pair<float,float> interval : occludedIntervals) {
            float start = interval.first;
            float stop = interval.second;
            
            if (start > 0.0f) {
                // occlusion occurs on one side of the normtal
                // unweighted:
                // Integrate[Sin[t], {t, startTheta, stopTheta}]
                //
                // Cosine-Weighted:
                // Integrate[sin(t)cos(t),{t, stopTheta, startTheta]
                amountOccluded += cos(angleScale*start) - cos(angleScale*stop);
            } else if (stop < 0.0f) {
                amountOccluded += cos(angleScale*stop) - cos(angleScale*start);
            } else {
                // occlusion straddles normal
                // unweighted:
                // Integrate[Sin[t], {t, startTheta, stopTheta}]
                //
                // Cosine-Weighted:
                // Integrate[sin(t)cos(t),  {t, 0, -start}] +
                // Integrate[sint(t)cos(t), {t, 0, stop}]
                amountOccluded += 2 - cos(angleScale*start) - cos(angleScale*stop);
            }
        }
        return integralScale*amountOccluded;
    }
    
    AVLIntegrator *CreateAVLIntegrator(const ParamSet &params, std::shared_ptr<Sampler> sampler,
                                       std::shared_ptr<const Camera> camera) {
        float numPhiSamples;
        numPhiSamples = params.FindOneInt("numPhiSamples",1);
        
        bool cosineWeighted;
        cosineWeighted = params.FindOneBool("cosineWeighted",false);
        
        int np;
        const int *pb = params.FindInt("pixelbounds", &np);
        Bounds2i pixelBounds = camera->film->GetSampleBounds();
        if (pb) {
            if (np != 4)
                Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                      np);
            else {
                pixelBounds = Intersect(pixelBounds,
                                        Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
                if (pixelBounds.Area() == 0)
                    Error("Degenerate \"pixelbounds\" specified.");
            }
        }
        return new AVLIntegrator(camera, sampler,
                                 pixelBounds,
                                 numPhiSamples,
                                 cosineWeighted);
    }
    
}  // namespace pbrt
