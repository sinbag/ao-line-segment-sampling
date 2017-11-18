//
//  multijitter.h
//  PBRT-V3
//
//  Created by Gurpreet Singh Bagga on 23/12/16.
//
//

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SAMPLERS_MultiJitter_H
#define PBRT_SAMPLERS_MultiJitter_H


// samplers/MultiJitter.h*
#include "sampler.h"
#include "rng.h"


namespace pbrt {


// MultiJitterSampler Declarations
class MultiJitterSampler : public PixelSampler {
public:
    // MultiJitterSampler Public Methods
    MultiJitterSampler(int xPixelSamples, int yPixelSamples, bool jitterSamples,
                      int nSampledDimensions)
    : PixelSampler(xPixelSamples * yPixelSamples, nSampledDimensions),
    xPixelSamples(xPixelSamples),
    yPixelSamples(yPixelSamples),
    jitterSamples(jitterSamples) {}
    void StartPixel(const Point2i &);
    std::unique_ptr<Sampler> Clone(int seed);
    
private:
    // MultiJitterSampler Private Data
    const int xPixelSamples, yPixelSamples;
    const bool jitterSamples;
};

MultiJitterSampler *CreateMultiJitterSampler(const ParamSet &params);
    
} //namespace pbrt

#endif  // PBRT_SAMPLERS_MultiJitter_H
