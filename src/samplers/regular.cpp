//
//  regular.cpp
//  PBRT-V3
//
//  Created by Gurprit Singh on 11/3/16.
//
//


/*
 pbrt source code is Copyright(c) 1998-2015
 Matt Pharr, Greg Humphreys, and Wenzel Jakob.
 
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



// samplers/Regular.cpp*
#include "samplers/regular.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt {

// RegularSampler Method Definitions
void RegularSampler::StartPixel(const Point2i &p) {
    // Generate single Regular samples for the pixel
    for (size_t i = 0; i < samples1D.size(); ++i) {
        RegularSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
                           jitterSamples);
    }
    for (size_t i = 0; i < samples2D.size(); ++i) {
        RegularSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng,
                           jitterSamples);
    }
    
    // Generate arrays of Regular samples for the pixel
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples1DArraySizes[i];
            RegularSample1D(&sampleArray1D[i][j * count], count, rng,
                               jitterSamples);
        }
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            
            int xSamples = sqrtf(samples2DArraySizes[i]);
            int ySamples = xSamples;
            int count = xSamples * ySamples;
            RegularSample2D(&sampleArray2D[i][j * count], xSamples, ySamples,
                               rng, jitterSamples);
            //Shuffle(&sampleArray2D[i][j * count], xSamples * ySamples, 1, rng);
            //LatinHypercube(&sampleArray2D[i][j * count].x, count, 2, rng);
        }
    PixelSampler::StartPixel(p);
}

std::unique_ptr<Sampler> RegularSampler::Clone(int seed) {
    RegularSampler *ss = new RegularSampler(*this);
    //ss->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(ss);
}

RegularSampler *CreateRegularSampler(const ParamSet &params) {
                                     
    //bool jitter = params.FindOneBool("jitter", false);
    int xsamp = params.FindOneInt("xsamples", 1);
    int ysamp = params.FindOneInt("ysamples", 1);
    int sd = params.FindOneInt("dimensions", 2);
    bool jitter = false;
    int offset = params.FindOneInt("offset", 0);
    
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new RegularSampler(xsamp, ysamp, offset, jitter, sd);
}

}  // namespace pbrt

