//
//  multijitter.cpp
//  PBRT-V3
//
//  Created by Gurpreet Singh Bagga on 23/12/16.
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

// samplers/MultiJitter.cpp*
#include "samplers/multijitter.h"
#include "paramset.h"
#include "sampling.h"
#include "stats.h"

namespace pbrt{

// MultiJitterSampler Method Definitions
void MultiJitterSampler::StartPixel(const Point2i &p) {
    // Generate single MultiJitter samples for the pixel
    for (size_t i = 0; i < samples1D.size(); ++i) {
        StratifiedSample1D(&samples1D[i][0], xPixelSamples * yPixelSamples, rng,
                           jitterSamples);
        Shuffle(&samples1D[i][0], xPixelSamples * yPixelSamples, 1, rng);
    }
        
    for (size_t i = 0; i < samples2D.size(); ++i) {
        MultiJitterSample2D(&samples2D[i][0], xPixelSamples, yPixelSamples, rng,
                           jitterSamples);
        
        for (int col = 0; col < xPixelSamples; col++){
            for (int row = yPixelSamples-1; row >= 1; row--) {
                int other = rng.UniformUInt32(row);
                int index1 = row * xPixelSamples + col;
                int index2 = other * yPixelSamples + col;
                std::swap(samples2D[i][index1].x, samples2D[i][index2].x);
            }
        }

        for (int row = 0; row < yPixelSamples; row++){
            for (int col = xPixelSamples-1; col >= 1; col--) {
                int other = rng.UniformUInt32(col);
                int index1 = row * xPixelSamples + col;
                int index2 = row * yPixelSamples + other;
                std::swap(samples2D[i][index1].y, samples2D[i][index2].y);
            }
        }
        Shuffle(&samples2D[i][0], xPixelSamples * yPixelSamples, 1, rng);
        
//        // Shuffle the columns
//        for (int col = 0; col < xPixelSamples; ++col){
//            for (int row = 0; row < yPixelSamples; ++row) {
//                int other = col + rng.UniformUInt32(xPixelSamples-col);;
//                int index1 = col * xPixelSamples + row;
//                int index2 = other * yPixelSamples + row;
//                std::swap(samples2D[i][index1].x, samples2D[i][index2].x);
//            }
//        }
//        
//        // Shuffle the rows
//        for (int col = 0; col < xPixelSamples; ++col){
//            for (int row = 0; row < yPixelSamples; ++row) {
//                int other = row + rng.UniformUInt32(xPixelSamples-row);
//                int index1 = col * xPixelSamples + row;
//                int index2 = col * yPixelSamples + other;
//                std::swap(samples2D[i][index1].y, samples2D[i][index2].y);
//            }
//        }
    }
    
    // Generate arrays of MultiJitter samples for the pixel
    for (size_t i = 0; i < samples1DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int count = samples1DArraySizes[i];
            StratifiedSample1D(&sampleArray1D[i][j * count], count, rng,
                               jitterSamples);
            Shuffle(&sampleArray1D[i][j * count], count, 1, rng);
        }
    for (size_t i = 0; i < samples2DArraySizes.size(); ++i)
        for (int64_t j = 0; j < samplesPerPixel; ++j) {
            int xSamples = sqrtf(samples2DArraySizes[i]);
            int ySamples = xSamples;
            int count = xSamples * ySamples;
            MultiJitterSample2D(&sampleArray2D[i][j * count], xSamples, ySamples,
                               rng, jitterSamples);
            //int count = samples2DArraySizes[i];
            //LatinHypercube(&sampleArray2D[i][j * count].x, count, 2, rng);
            
            for (int col = 0; col < xSamples; col++){
                for (int row = ySamples-1; row >= 1; row--) {
                    int other = rng.UniformUInt32(row);
                    int index1 = row * xSamples + col;
                    int index2 = other * ySamples + col;
                    std::swap(sampleArray2D[i][index1].x, sampleArray2D[i][index2].x);
                }
            }
            
            for (int row = 0; row < ySamples; row++){
                for (int col = xSamples-1; col >= 1; col--) {
                    int other = rng.UniformUInt32(col);
                    int index1 = row * xSamples + col;
                    int index2 = row * ySamples + other;
                    std::swap(sampleArray2D[i][index1].y, sampleArray2D[i][index2].y);
                }
            }
            Shuffle(&sampleArray2D[i][0], xSamples * ySamples, 1, rng);

        }
    
    PixelSampler::StartPixel(p);
}

std::unique_ptr<Sampler> MultiJitterSampler::Clone(int seed) {
    MultiJitterSampler *ss = new MultiJitterSampler(*this);
    ss->rng.SetSequence(seed);
    return std::unique_ptr<Sampler>(ss);
}

MultiJitterSampler *CreateMultiJitterSampler(const ParamSet &params) {
    bool jitter = params.FindOneBool("jitter", true);
    int xsamp = params.FindOneInt("xsamples", 4);
    int ysamp = params.FindOneInt("ysamples", 4);
    int sd = params.FindOneInt("dimensions", 4);
    if (PbrtOptions.quickRender) xsamp = ysamp = 1;
    return new MultiJitterSampler(xsamp, ysamp, jitter, sd);
}

} //namespace pbrt
