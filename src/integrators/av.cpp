// integrators/av.cpp*
#include "api.h"
#include "core/reflection.h"
#include "integrators/av.h"
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

#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

namespace pbrt {
    
    Spectrum AVIntegrator::Li(const RayDifferential &ray,
                              const Scene &scene, Sampler &sampler,
                              MemoryArena &arena, int depth) const {
        
        // Find closest ray intersection or return background radiance
        SurfaceInteraction isect;
        
        if (!scene.Intersect(ray, &isect)) {
            return Spectrum(1.0f);
        }
        
        Spectrum L(0.0f);
        float sum = 0.0;
        /****/
        ///To sample Phi directly from the built-in sampler
        const Point2f *uvSampleArray = sampler.Get2DArray(m_numSecondarySamples);
        
        
        ///HOmogenization:
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
        
        if (computeFourierSpectrum) {
            /****/
            ///Added by Gurprit Singh
            ///To visualize the integrand within a pixel
            /// This is not complete at the moment!!!
            
            std::string prefix = "ao-segments";
            
            int width = 512, height = 512;
            Bounds2i outputBounds = {Point2i(0,0), Point2i(width,height)};
            Point2i totalResolution(width, height);
            std::vector<Point2f> sampleLocations(m_numSecondarySamples, Point2f(0.,0.));
            std::vector<double> functionValue(m_numSecondarySamples, 0.);
            /****/
            
            
            // constants to optionally account for cosine-weighted integration
            for (int k = 0; k < m_numSecondarySamples; k++) {
                
                //double u = std::fmod((uvSampleArray[k].x + m_jitter*sampler.Get1D()), 1.0f);
                //double v = std::fmod((uvSampleArray[k].y + m_jitter*sampler.Get1D()), 1.0f);
                
                double u = uvSampleArray[k].x + randomShift;
                double v = uvSampleArray[k].y + randomShift;
                
                //std::cerr << u << " " << v << " "<< randomShift << std::endl;
                
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
                
                double pdf = 1.f;
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2*v-1);
                    pdf = InvPi * InvPi;
                }
                if(SolidAngleUniformSampling){
                    if(v <= 0.5)
                        theta = -acos(2*v);
                    else
                        theta = acos(2*(1-v));
                    pdf = fabs(sin(theta)) * Inv2Pi;
                }
                else if(cosineWeightedSampling){
                    if(v <= 0.5)
                        theta = - 0.5 * acos(4*v-1);
                    else
                        theta = 0.5 * acos(3 - 4*v);
                    pdf = fabs(sin(theta))*cos(theta) * InvPi;
                }
                
                Vector3f localDir = sphericalDirection(theta, phi);
                
                BSDF *bsdf = new BSDF(isect);
                Vector3f worldDir = bsdf->LocalToWorld(localDir);
                Ray r(isect.p,worldDir,lineLength,0);
                SurfaceInteraction isect2;
                if(!scene.Intersect(r,&isect2)){
                    
                    float radiance = (InvPi*cos(theta)*fabs(sin(theta))) / pdf;
                    functionValue[k] = radiance;
                }
            }
            std::stringstream ss;
            
//            int half_xRes = width * 0.5;
//            int half_yRes = height * 0.5;
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
            
            delete [] powerSpectrum;
            delete [] realFourierSpectrum;
            delete [] imagFourierSpectrum;
        }
        
        if(visualizeIntegrand){
            /****/
            ///Added by Gurprit Singh
            ///To visualize the integrand within a pixel
            
            int width = 512, height = 512;
            Bounds2i outputBounds = {Point2i(0,0), Point2i(width,height)};
            Point2i totalResolution(width, height);
            std::unique_ptr<Float[]> sampleColorData(new Float[3 * width * height]);
            /****/
            
            for (int k = 0; k < m_numSecondarySamples; k++) {
                float u = uvSampleArray[k].x+randomShift;
                float v = uvSampleArray[k].y+randomShift;
                
                //std::cerr << u << " "<< v << " " << randomShift << std::endl;
                
                if(u > 1)
                    u = u - 1.f;
                else if(u < 0)
                    u = 1.f + u;
                
                if(v > 1)
                    v = v - 1.f;
                else if(v < 0)
                    v = 1.f + v;
                
                if((u > 1 || u < 0) || (v > 1 || v < 0)){
                    std::cerr << "Abort! AVIntegrator random samples not in unit square domain: " << u <<" " << v << std::endl;
                    exit(-1);
                }
                
                double pdf = 1.f;
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2*v-1);
                    pdf = InvPi * InvPi;
                }
                if(SolidAngleUniformSampling){
                    if(v <= 0.5)
                        theta = -acos(2*v);
                    else
                        theta = acos(2*(1-v));
                    pdf = fabs(sin(theta)) * Inv2Pi;
                }
                else if(cosineWeightedSampling){
                    if(v <= 0.5)
                        theta = - 0.5 * acos(4*v-1);
                    else
                        theta = 0.5 * acos(3 - 4*v);
                    pdf = fabs(sin(theta))*cos(theta) * InvPi;
                }
                
                
                Vector3f localDir = sphericalDirection(theta, phi);
                
                BSDF *bsdf = new BSDF(isect);
                Vector3f worldDir = bsdf->LocalToWorld(localDir);
                Ray r(isect.p,worldDir,lineLength,0);
                SurfaceInteraction isect2;
                if(!scene.Intersect(r,&isect2)){
                    
                    float radiance = (InvPi*cos(theta)*fabs(sin(theta))) / pdf;
                    
                    int col = u * width;
                    int row = v * height;
                    int index = row*width+col;
                    sampleColorData[3*index+0] = radiance;
                    sampleColorData[3*index+1] = radiance;
                    sampleColorData[3*index+2] = radiance;
                }
            }
            std::stringstream ss;
            ss << "ao-points.exr";
            WriteImage(ss.str(), &sampleColorData[0], outputBounds, totalResolution);
        }
        else{
            for (int k = 0; k < m_numSecondarySamples; k++) {
                float u = uvSampleArray[k].x+randomShift;
                float v = uvSampleArray[k].y+randomShift;
                
                //std::cerr << u << " "<< v << " " << randomShift << std::endl;
                
                if(u > 1)
                    u = u - 1.f;
                else if(u < 0)
                    u = 1.f + u;
                
                if(v > 1)
                    v = v - 1.f;
                else if(v < 0)
                    v = 1.f + v;
                
                if((u > 1 || u < 0) || (v > 1 || v < 0)){
                    std::cerr << "Abort! AVIntegrator random samples not in unit square domain: " << u <<" " << v << std::endl;
                    exit(-1);
                }
                
                double pdf = 1.f;
                double phi = Pi*u;
                double theta = 0;
                if(SphericalCoordinatesUniformSampling){
                    theta = PiOver2 * (2*v-1);
                    pdf = InvPi * InvPi;
                }
                if(SolidAngleUniformSampling){
                    if(v <= 0.5)
                        theta = -acos(2*v);
                    else
                        theta = acos(2*(1-v));
                    pdf = fabs(sin(theta)) * Inv2Pi;
                }
                else { //if(cosineWeightedSampling)
                    if(v <= 0.5)
                        theta = - 0.5 * acos(4*v-1);
                    else
                        theta = 0.5 * acos(3 - 4*v);
                    pdf = fabs(sin(theta))*cos(theta) * InvPi;
                }
                
                Vector3f localDir = sphericalDirection(theta, phi);
                
                //			Vector3f localDir = cosineWeighted
                //				? AVIntegrator::squareToCosineHemisphere(u,v)
                //				: AVIntegrator::squareToUniformHemisphere(u,v);
                BSDF *bsdf = new BSDF(isect);
                Vector3f worldDir = bsdf->LocalToWorld(localDir);
                Ray r(isect.p,worldDir,lineLength,0);
                SurfaceInteraction isect2;
                if(!scene.Intersect(r,&isect2)){
                    
                    double perSampleRadiance = (InvPi*cos(theta)*fabs(sin(theta))) / pdf;
                    sum += (perSampleRadiance);
                }
            }
        }
        return Spectrum((sum/m_numSecondarySamples));
    }
    
    // Utility Functions
    Vector3f AVIntegrator::sphericalDirection(float theta, float phi) {
        float sinTheta, cosTheta, sinPhi, cosPhi;
        //__sincosf(theta,&sinTheta, &cosTheta);
        //__sincosf(phi, &sinPhi, &cosPhi);
        sinTheta = std::sin(theta);
        cosTheta = std::cos(theta);
        sinPhi = std::sin(phi);
        cosPhi = std::cos(phi);
        
        return Vector3f(
                        (sinTheta) * cosPhi,
                        (sinTheta) * sinPhi,
                        cosTheta );
    }
    
    Vector3f AVIntegrator::squareToCosineHemisphere(float u, float v) {
        Point2f p = AVIntegrator::squareToUniformDisk(u,v);
        float z = std::sqrt(std::max((float) 0,
                                     1.0f - p.x*p.x- p.y*p.y));
        
        return Vector3f(p.x, p.y, z);
    }
    
    Point2f AVIntegrator::squareToUniformDisk(float u, float v) {
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
    
    Vector3f AVIntegrator::squareToUniformHemisphere(float u, float v) {
        float cosTheta = v;
        float sinTheta = std::sqrt(std::max((float) 0, 1-cosTheta*cosTheta));
        
        float sinPhi, cosPhi;
        //__sincosf(2.0f * Pi * u, &sinPhi, &cosPhi);
        sinPhi = std::sin(2.0f * Pi * u);
        cosPhi = std::cos(2.0f * Pi * u);
        
        return Vector3f(cosPhi * sinTheta, sinPhi * sinTheta, cosTheta);
    }
    
    
    AVIntegrator *CreateAVIntegrator(
                                     const ParamSet &params, std::shared_ptr<Sampler> sampler,
                                     std::shared_ptr<const Camera> camera) {
        
        int numSecondarySamples;
        numSecondarySamples = params.FindOneInt("numSecondarySamples",1);
        
        bool cosineWeightedSampling;
        cosineWeightedSampling = params.FindOneBool("cosineWeightedSampling",true);
        
        bool SphericalCoordinatesUniformSampling;
        SphericalCoordinatesUniformSampling = params.FindOneBool("SphericalCoordinatesUniformSampling",false);
        
        bool SolidAngleUniformSampling;
        SolidAngleUniformSampling = params.FindOneBool("SolidAngleUniformSampling",false);
        
        bool visualizeIntegrand;
        visualizeIntegrand = params.FindOneBool("visualizeIntegrand",false);
        
        bool computeFourierSpectrum;
        computeFourierSpectrum = params.FindOneBool("computeFourierSpectrum",false);
        
        float lineLength;
        lineLength = params.FindOneFloat("lineLength",500.0);
        
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
        return new AVIntegrator(camera, sampler,
                                pixelBounds,
                                numSecondarySamples,
                                cosineWeightedSampling,
                                SphericalCoordinatesUniformSampling,
                                SolidAngleUniformSampling,
                                visualizeIntegrand,
                                computeFourierSpectrum,
                                lineLength);
    }
    
}  // namespace pbrt
