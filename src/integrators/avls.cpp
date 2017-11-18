
// integrators/avls.cpp*
#include "api.h"
#include "core/reflection.h"
#include "integrators/avls.h"
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
#include "imageio.h"
#include <cmath>
#include <cstdlib>

#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace pbrt {
    
    // DirectLightingIntegrator Method Definitions
    Spectrum AVLSIntegrator::Li(const RayDifferential &ray,
                                const Scene &scene, Sampler &sampler,
                                MemoryArena &arena, int depth) const {
        // Find closest ray intersection or return background radiance
        SurfaceInteraction isect;
        if (!scene.Intersect(ray, &isect)) {
            return Spectrum(1.0f);
        }
        
        ///HOmogenization:
        //@sinbag=================================
        ///Code added by Gurprit Singh [2016-10-10]
        /// New random seed is created for each tile
        /// which changes with tile coordinates
        ///
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(0, 1);
        double randomShift = dis(gen);
        
        //@sinbag=================================
        
        /****/
        ///Added by Gurprit Singh
        ///To sample Phi directly from the built-in sampler
        
        const Point2f *uvSampleArray = sampler.Get2DArray(m_totalSamples);
        
        if (computeFourierSpectrum) {
            /****/
            ///Added by Gurprit Singh
            ///To visualize the integrand within a pixel
            /// This is not complete at the moment!!!
            
            std::string prefix = "ao-segments";
            
            int width = 512, height = 512;
            Bounds2i outputBounds = {Point2i(0,0), Point2i(width,height)};
            Point2i totalResolution(width, height);
            std::vector<Point2f> sampleLocations(m_totalSamples, Point2f(0.,0.));
            std::vector<double> functionValue(m_totalSamples, 0.);
            /****/
            
            
            // constants to optionally account for cosine-weighted integration
            double length = std::max(0.0f, std::min(Pi, lineSegmentLength));
            double pdf = 1.0;
            double weight = InvPi;
            for (int k = 0; k < m_totalSamples; k++) {
                ///Homogenizing the samples by adding a
                ///random shift vector to each sample.
                double u = uvSampleArray[k].x + randomShift;
                double v = uvSampleArray[k].y + randomShift;
                
                if(u > 1)
                    u = u - 1.f;
                else if(u < 0)
                    u = 1.f + u;
                
                if(v > 1)
                    v = v - 1.f;
                else if(v < 0)
                    v = 1.f + v;
                
                if((u > 1 || u < 0) || (v > 1 || v < 0)){
                    std::cerr << "Abort! AVLSIntegrator random samples not in unit square domain: " << u <<" " << v << std::endl;
                    exit(-1);
                }
                
                Point2f unitSample(u,v);
                sampleLocations[k] = unitSample;
                
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2.0f*v-1.0f);
                    pdf = InvPi*InvPi;
                }
                else if(SolidAngleUniformSampling){
                    if ( v <= 0.5) {
                        theta = -acos(2.0f*v);
                    } else {
                        theta = acos(2.0f*(1.0f-v));
                    }
                    pdf = std::abs(std::sin(theta)) * Inv2Pi;
                }
                else if(cosineWeightedSampling){
                    if ( v <= 0.5 ) {
                        theta = -0.5*acos(4.0f*v-1.0f);
                    } else {
                        theta = 0.5*acos(3.0f-4.0f*v);
                    }
                    pdf = std::abs(std::sin(theta)) * std::cos(theta) * InvPi;
                }
                double pointValue = AVLSIntegrator::getBlurredValue(phi, theta, length,
                                                                    &isect, spheres);
                pointValue = weight*pointValue/pdf;
                functionValue[k] = pointValue;
            }
            std::stringstream ss;
            
            //        int half_xRes = width * 0.5;
            //        int half_yRes = height * 0.5;
            int npoints = sampleLocations.size();
            double fracPts = 1.0/double(npoints);
            double _frequencyStep = 1.0;
            float *powerSpectrum = new float[3*width*height]();
            float *realFourierSpectrum = new float[3*width*height]();
            float *imagFourierSpectrum = new float[3*width*height]();
            
            std::cerr << "Fourier computation..." << std::endl;
            
            ///
            /// Compute the Power spectrum
            ///
            //tbb::tick_count t0 = tbb::tick_count::now();
            tbb::task_scheduler_init init(8);
            tbb::parallel_for(
                              tbb::blocked_range2d<int>(0,width, 16, 0, height, 16),
                              [=](const tbb::blocked_range2d<int>& imgblock ) {
                                  for( int row = imgblock.rows().begin(); row != imgblock.rows().end(); ++row ){
                                      for( int col = imgblock.cols().begin(); col != imgblock.cols().end(); ++col ) {
                                          double fx = 0.f, fy = 0.f;
                                          //double w1 = (col-half_xRes)*_frequencyStep;
                                          //double w2 = (row-half_yRes)*_frequencyStep;
                                          double w1 = (col)*_frequencyStep;
                                          double w2 = (row)*_frequencyStep;
                                          double exp = 0.0;
                                          for (int i = 0; i < npoints; ++i) {
                                              exp = -2 * M_PI * (w1 * sampleLocations[i].x + w2 * sampleLocations[i].y);
                                              fx += functionValue[i] * cos(exp);
                                              fy += functionValue[i] * sin(exp);
                                          }
                                          int index = row * width + col;
                                          
                                          realFourierSpectrum[3*index+0] = fx * fracPts;
                                          realFourierSpectrum[3*index+1] = fx * fracPts;
                                          realFourierSpectrum[3*index+2] = fx * fracPts;
                                          
                                          imagFourierSpectrum[3*index+0] = fy * fracPts;
                                          imagFourierSpectrum[3*index+1] = fy * fracPts;
                                          imagFourierSpectrum[3*index+2] = fy * fracPts;
                                          
                                          double powerVal = (fx*fx + fy*fy) * fracPts;
                                          powerSpectrum[3*index+0] = powerVal;
                                          powerSpectrum[3*index+1] = powerVal;//(gfx*gfx + gfy*gfy) * fracPts;
                                          powerSpectrum[3*index+2] = powerVal;//(bfx*bfx + bfy*bfy) * fracPts;
                                          
                                          //std::cerr << row << " " << col << " " << fx << " "<< fy  <<" " << powerVal << std::endl;
                                      }
                                  }
                              }
                              );
            ///
            /// Write the powerspectrum to an image
            ///
            ss.str(std::string());
            ss << prefix << "-ao-fourier-real.exr";
            WriteImage(ss.str(), realFourierSpectrum, outputBounds, totalResolution);
            
            ss.str(std::string());
            ss << prefix << "-ao-fourier-imag.exr";
            WriteImage(ss.str(), imagFourierSpectrum, outputBounds, totalResolution);
            
            ss.str(std::string());
            ss << prefix << "-ao-fourier-power.exr";
            WriteImage(ss.str(), powerSpectrum, outputBounds, totalResolution);
            
        }
        
        
        double sum = 0.f;
        if (visualizeIntegrand) {
            /****/
            ///Added by Gurprit Singh
            ///To visualize the integrand within a pixel
            
            int width = 512, height = 512;
            Bounds2i outputBounds = {Point2i(0,0), Point2i(width,height)};
            Point2i totalResolution(width, height);
            std::unique_ptr<Float[]> sampleColorData(new Float[3 * width * height]);
            /****/
            
            
            // constants to optionally account for cosine-weighted integration
            double length = std::max(0.0f, std::min(Pi, lineSegmentLength));
            double pdf = 1.0;
            double weight = InvPi;
            for (int k = 0; k < m_totalSamples; k++) {
                
                double u = uvSampleArray[k].x;
                double v = uvSampleArray[k].y;
                
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2.0f*v-1.0f);
                    pdf = InvPi*InvPi;
                }
                else if(SolidAngleUniformSampling){
                    if ( v <= 0.5) {
                        theta = -acos(2.0f*v);
                    } else {
                        theta = acos(2.0f*(1.0f-v));
                    }
                    pdf = std::abs(std::sin(theta)) * Inv2Pi;
                }
                else if(cosineWeightedSampling){
                    if ( v <= 0.5 ) {
                        theta = -0.5*acos(4.0f*v-1.0f);
                    } else {
                        theta = 0.5*acos(3.0f-4.0f*v);
                    }
                    pdf = std::abs(std::sin(theta)) * std::cos(theta) * InvPi;
                }
                double pointValue = AVLSIntegrator::getBlurredValue(phi, theta, length,
                                                                    &isect, spheres);
                pointValue = weight*pointValue/pdf;
                int col = u * width;
                int row = v * height;
                int index = row*width+col;
                sampleColorData[3*index+0] = pointValue;
                sampleColorData[3*index+1] = pointValue;
                sampleColorData[3*index+2] = pointValue;
            }
            std::stringstream ss;
            
            if (cosineWeightedSampling) {
                ss << "ao-segments-cws-l" << length << ".exr";
            }
            else if (SphericalCoordinatesUniformSampling){
                ss << "ao-segments-usscoords-l" << length << ".exr";
            }
            else if (SolidAngleUniformSampling) {
                ss << "ao-segments-ussa-l" << length << ".exr";
            }
            
            WriteImage(ss.str(), &sampleColorData[0], outputBounds, totalResolution);
        }
        else {
            // constants to optionally account for cosine-weighted integration
            double length = std::max(0.0f, std::min(Pi, lineSegmentLength));
            double pdf = 1.0;
            double weight = InvPi;
            for (int k = 0; k < m_totalSamples; k++) {
                
                double u = uvSampleArray[k].x + randomShift;
                double v = uvSampleArray[k].y + randomShift;
                
                if(u > 1)
                    u = u - 1.f;
                else if(u < 0)
                    u = 1.f + u;
                
                if(v > 1)
                    v = v - 1.f;
                else if(v < 0)
                    v = 1.f + v;
                
                if((u > 1 || u < 0) || (v > 1 || v < 0)){
                    std::cerr << "Abort! AVLSIntegrator random samples not in unit square domain: " << u <<" " << v << std::endl;
                    exit(-1);
                }
                
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2.0f*v-1.0f);
                    pdf = InvPi*InvPi;
                }
                else if(SolidAngleUniformSampling){
                    if ( v <= 0.5) {
                        theta = -acos(2.0f*v);
                    } else {
                        theta = acos(2.0f*(1.0f-v));
                    }
                    pdf = std::abs(std::sin(theta)) * Inv2Pi;
                }
                else if(cosineWeightedSampling){
                    if ( v <= 0.5 ) {
                        theta = -0.5*acos(4.0f*v-1.0f);
                    } else {
                        theta = 0.5*acos(3.0f-4.0f*v);
                    }
                    pdf = std::abs(std::sin(theta)) * std::cos(theta) * InvPi;
                }
                double pointValue = weight * AVLSIntegrator::getBlurredValue(phi, theta, length,
                                                                             &isect, spheres);
                sum += pointValue / pdf;
            }
        }
        double visibility = sum/double(m_totalSamples);
        return Spectrum(visibility);
    }
    
    double AVLSIntegrator::getBlurredValue(double phi, double theta, double length,
                                           SurfaceInteraction *isect_ptr,
                                           std::vector<Sphere> spheres,
                                           bool cosineWeightedIntegrand) {
        SurfaceInteraction isect = *isect_ptr;
        if (typeid(*isect.shape) != typeid(Sphere)) {
            std::cout<<"Scene should only contain spheres"<<std::endl;
            std::exit(-1);
        }
        //Get the hit sphere information to compare with spheres in scene
        const Sphere *hit_sphere = static_cast<const Sphere*> (isect.shape);
        Point3f hit_sphere_center = (*hit_sphere->ObjectToWorld)(Point3f(0,0,0));
        float hit_sphere_radius = hit_sphere->radius;
        BSDF *bsdf = new BSDF(isect);
        
        double angleScale = 2;// cosineWeightedIntegrand ? 2 : 1;
        double integralScale = 0.25;// cosineWeightedIntegrand ? 0.25 : 1;
//        double maxIntegral = cosineWeightedIntegrand ? 1 : 0.5;
        
        float start = theta - length/2.f;
        float stop = theta + length/2.f;
        
        if (-PiOver2 - start > MachineEpsilon) start += Pi;
        if (stop - PiOver2 > MachineEpsilon) stop -= Pi;
        std::pair<float, float> integralInterval(start, stop);
        
        // Given theta, determine a random range to integrate over
        //plane normal is the vector in the tangent plane perpendicular to phi
        Vector3f planeNormal = sphericalDirection(PiOver2, phi + PiOver2);
        Vector3f bitangent = sphericalDirection(PiOver2, phi);
        planeNormal = bsdf->LocalToWorld(planeNormal);
        bitangent = bsdf->LocalToWorld(bitangent);
        
        // Determine the occlusion contributed by each sphere
        std::vector<std::pair<float, float>> occludedRegions;
        for(Sphere s: spheres) {
            Point3f sphere_center = (*s.ObjectToWorld)(Point3f(0,0,0));
            float sphere_radius = s.radius;
            
            if ( sphere_center == hit_sphere_center
                && sphere_radius == hit_sphere_radius) {
                continue;
            }
            
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
            std::pair<float,float> occludedRegion(startTheta,stopTheta);
            occludedRegions.push_back(occludedRegion);
        }
        
        //Calculate the amount of visbility given the occluded regions found above
        occludedRegions = AVLSIntegrator::combineIntervals(occludedRegions);
        
        std::vector<std::pair<float, float>> visibleRegions;
        std::vector<std::pair<float, float>> cappedVisibleRegions;
        visibleRegions = AVLSIntegrator::getComplimentOfIntervals(occludedRegions);
        for(auto region : visibleRegions) {
            for(auto cappedRegion: getCappedRegions(region.first, region.second, integralInterval)) {
                cappedVisibleRegions.push_back(cappedRegion);
            }
        }
        
        double integralSum = AVLSIntegrator::calculateVisibility(cappedVisibleRegions, angleScale);
        
        //std::cout << integralSum + occludedIntegralSum << std::endl;
        //std::cout << "Phi: " << phi << std::endl;
        //std::cout << "Integral Sum: " << integralSum << std::endl;
        return integralScale*integralSum/length;
    }
    
    // Utility Functions
    Vector3f AVLSIntegrator::sphericalDirection(float theta, float phi) {
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
                        cosTheta );
    }
    
    Vector3f AVLSIntegrator::squareToCosineHemisphere(float u, float v) {
        Point2f p = AVLSIntegrator::squareToUniformDisk(u,v);
        float z = std::sqrt(std::max((float) 0,
                                     1.0f - p.x*p.x- p.y*p.y));
        
        return Vector3f(p.x, p.y, z);
    }
    
    Point2f AVLSIntegrator::squareToUniformDisk(float u, float v) {
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
    
    Vector3f AVLSIntegrator::squareToUniformHemisphere(float u, float v) {
        float cosTheta = v;
        float sinTheta = std::sqrt(std::max((float) 0, 1-cosTheta*cosTheta));
        
        float sinPhi, cosPhi;
        //__sincosf(2.0f * Pi * u, &sinPhi, &cosPhi);
        sinPhi = std::sin(2.0f * Pi * u);
        cosPhi = std::cos(2.0f * Pi * u);
        
        return Vector3f(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
    }
    
    bool AVLSIntegrator::compareIntervals(const std::pair<float,float> &i1, const std::pair<float,float> &i2) {
        return i1.first < i2.first;
    }
    
    std::vector<std::pair<float,float>> AVLSIntegrator::combineIntervals(std::vector<std::pair<float,float>> intervals) {
        
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
    
    /*  caps the occluded regions so that only the portion intersection the integral
     *  interval is preserved. If hemisphere wrapping has occured, may have two intervals
     *  generated from startTheta and stopTheta
     *
     *  @param: startTheta - the starting theta value for the occluded region
     *	@param: stopTheta - the stopping theta value for the occluded region
     *  @param: integralInterval (start, stop) - pair of values representing the range of
     *		integration for the line segment
     *  @returns: cappedRegions - the capped occluded regions
     */
    std::vector<std::pair<float,float>> AVLSIntegrator::getCappedRegions(
                                                                         float startTheta, float stopTheta,
                                                                         std::pair<float,float> integralInterval) {
        float start = integralInterval.first;
        float stop = integralInterval.second;
        std::vector<std::pair<float,float>> cappedRegions;
        
        if (std::fabs(stop - start) < ShadowEpsilon) {
            std::pair<float,float> interval(startTheta, stopTheta);
            cappedRegions.push_back(interval);
            return cappedRegions;
        }
        
        if (start < stop) {
            // the integral interval doesn't wrap around the hemisphere
            float cappedStart = std::max(startTheta, start);
            float cappedStop = std::min(stopTheta, stop);
            
            if (cappedStart < cappedStop) {
                // integral interval and the occluded region actually intersect
                std::pair<float,float> interval(cappedStart, cappedStop);
                cappedRegions.push_back(interval);
            }
        } else if (startTheta > start) {
            // integral interval wraps around horizon and occluded region
            // falls entirely within (start, PiOver2)
            float cappedStop = std::min(PiOver2, stopTheta);
            std::pair<float,float> interval(startTheta, cappedStop);
            cappedRegions.push_back(interval);
        } else if (stopTheta < stop) {
            // integral interval wraps around horizon and occluded region
            // falls entirely within (-PiOver2, stop)
            float cappedStart = std::max(-PiOver2, startTheta);
            std::pair<float,float> interval(cappedStart, stopTheta);
            cappedRegions.push_back(interval);
        } else {
            // integral interval wraps around horizon
            if (startTheta < stop) {
                // part of occluded region falls within (-PiOver2, stop)
                std::pair<float,float> interval(startTheta, stop);
                cappedRegions.push_back(interval);
            }
            if (stopTheta > start) {
                // part of occluded region falls within (start, PiOver2)
                std::pair<float,float> interval(start, stopTheta);
                cappedRegions.push_back(interval);
            }
        }
        return cappedRegions;
    }
    
    /* getComplimentOfIntervals -- returns the compliment of intervals on the range -pi/2 to pi/2. For example
     * given the vector of intervals [(-pi/4,0), (pi/4,pi/2)] would return [(-pi/2,-pi/4),(0,pi4)]
     *
     */
    std::vector<std::pair<float,float>> AVLSIntegrator::getComplimentOfIntervals(
                                                                                 std::vector<std::pair<float,float>> intervals) {
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
    
    double AVLSIntegrator::calculateVisibility(
                                               std::vector<std::pair<float,float>> visibleIntervals, double angleScale) {
        // calculate the integral of the occluded region
        double visibility = 0.f;
        for (std::pair<float,float> interval : visibleIntervals) {
            double start = interval.first;
            double stop = interval.second;
            if (start > 0.0f) {
                // occlusion occurs on one side of the normal
                // unweighted:
                // Integrate[Sin[t], {t, startTheta, stopTheta}]
                //
                // Cosine-Weighted:
                // Integrate[sin(t)cos(t),{t, startTheta, stopTheta]
                visibility += cosf(angleScale*start) - cosf(angleScale*stop);
            } else if (stop < 0.0f) {
                visibility += cosf(angleScale*stop) - cosf(angleScale*start);
            } else {
                // occlusion straddles normal
                // unweighted:
                // Integrate[Sin[t], {t, startTheta, stopTheta}]
                //
                // Cosine-Weighted:
                // Integrate[sin(t)cos(t),  {t, 0, -start}] +
                // Integrate[sint(t)cos(t), {t, 0, stop}]
                visibility += 2.f - cosf(angleScale*start) - cosf(angleScale*stop);
            }
        }
        return visibility;
    }
    
    AVLSIntegrator *CreateAVLSIntegrator(
                                         const ParamSet &params, std::shared_ptr<Sampler> sampler,
                                         std::shared_ptr<const Camera> camera) {
        
        int numPhiSamples, numThetaSamples;
        numPhiSamples = params.FindOneInt("numPhiSamples",1);
        numThetaSamples = params.FindOneInt("numThetaSamples", 1);
        
        bool cosineWeightedSampling;
        cosineWeightedSampling = params.FindOneBool("cosineWeightedSampling",false);
        
        bool SolidAngleUniformSampling;
        SolidAngleUniformSampling = params.FindOneBool("SolidAngleUniformSampling",false);
        
        bool SphericalCoordinatesUniformSampling;
        SphericalCoordinatesUniformSampling = params.FindOneBool("SphericalCoordinatesUniformSampling",false);
        
        bool visualizeIntegrand;
        visualizeIntegrand = params.FindOneBool("visualizeIntegrand", false);
        
        bool computeFourierSpectrum;
        computeFourierSpectrum = params.FindOneBool("computeFourierSpectrum", false);
        
        float lineSegmentLength;
        lineSegmentLength = params.FindOneFloat("lineSegmentLength", PiOver2);
        
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
        return new AVLSIntegrator(camera,
                                  sampler,
                                  pixelBounds,
                                  numPhiSamples,
                                  numThetaSamples,
                                  cosineWeightedSampling,
                                  SolidAngleUniformSampling,
                                  SphericalCoordinatesUniformSampling,
                                  visualizeIntegrand,
                                  computeFourierSpectrum,
                                  lineSegmentLength);
    }
    
}  // namespace pbrt
