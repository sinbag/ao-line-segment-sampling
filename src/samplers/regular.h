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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_SAMPLERS_Regular_H
#define PBRT_SAMPLERS_Regular_H

// samplers/Regular.h*
#include "sampler.h"
#include "rng.h"

namespace pbrt {

// RegularSampler Declarations
class RegularSampler : public PixelSampler {
public:
    // RegularSampler Public Methods
    RegularSampler(int xPixelSamples, int yPixelSamples, uint64_t rngOffset,
                   bool jitterSamples, int nSampledDimensions)
    : PixelSampler(xPixelSamples*yPixelSamples, nSampledDimensions),
    xPixelSamples(xPixelSamples),
    yPixelSamples(yPixelSamples),
    jitterSamples(jitterSamples) {}
    void StartPixel(const Point2i &);
    std::unique_ptr<Sampler> Clone(int seed);
    
private:
    // RegularSampler Private Data
    const int xPixelSamples, yPixelSamples;
    const bool jitterSamples;
};

RegularSampler *CreateRegularSampler(const ParamSet &params);

}  // namespace pbrt
#endif  // PBRT_SAMPLERS_Regular_H
