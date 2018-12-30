// Copyright 2015 Lior Kogan (koganlior1@gmail.com)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is 
// distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and limitations under the License.

// ==========================================================================
#include "StreamSamplerTest.h"

#include "StreamSampler.h"

#include <chrono>
#include <iostream>    // cout
#include <vector>
#include <numeric>
#include <functional>
#include <math.h>

using namespace std;
using namespace StreamSampler;

// ==========================================================================
// Example
// ==========================================================================

// a simple stream of incremental integers: 0..nStreamSize-1
class SimpleStream
{
public:
    SimpleStream(size_t nStreamSize)
        : m_nStreamSize(nStreamSize), m_nNextElement(0)
    {}

    bool GetNextElement(size_t& Element)
    {
        if (m_nNextElement == m_nStreamSize)
            return false; // end of stream

        Element = m_nNextElement++;
        return true;
    }

private:
    size_t m_nStreamSize ;
    size_t m_nNextElement;
};

// example: using a StreamSampler
bool StreamSamplerExample()
{
    SimpleStream                Stream(1000)         ; // 0..999
    CStreamSamplerWOR_Z<size_t> StreamSampler(10, 10); // 10 sample sets of 10 samples each

    uint64_t nSkip = 0;                                // number of future stream elements to skip
    size_t   Element  ;                                // stream element

    while (Stream.GetNextElement(Element))             // while not end of stream
        nSkip ? --nSkip : nSkip = StreamSampler.AddElement(Element);

    auto SampleSets = StreamSampler.GetSampleSets();   // get sample sets

    for (const auto& SampleSet : SampleSets)           // for each sample set
    {
        for (const auto& Sample : SampleSet)           // for each sample
            cout << Sample << " ";                     // print it

        cout << "\n";
    }

    return true;
}

// static auto a = StreamSamplerExample();

// ==========================================================================
// Chi-squared test (for uniformity test)
// code based on "df.c" and "funct.c" from the Marsaglia's Diehard Battery of Tests of Randomness
// (http://www.stat.fsu.edu/pub/diehard/)
// ==========================================================================

// gamma(z) when 2*fZ is integer
double Gamma(double fZ)
{
    uint64_t nTmp = (uint64_t)(2. * fZ);
    if (nTmp != 2. * fZ || fZ == 0)
        throw invalid_argument("Gamma: invalid argument"); 

         if (nTmp == 1) return 1.7724538509055160272981674833411; // sqrt(PI)
    else if (nTmp == 2) return 1.                               ;
    else                return (fZ - 1.) * Gamma(fZ - 1.)       ;
}

// CDF of Standard Normal
double Phi(double fX)
{
    double fTmp = fX / sqrt(2.);
    fTmp        = 1. + erf(fTmp); 
    return fTmp / 2.;
}

// PDF of Chi-square
double PDFchisq(double fV, uint64_t nDF)
{
    return pow(fV / 2., (nDF - 2) / 2.) * exp(-fV / 2.) / (2. * Gamma(nDF / 2.));
}

// CDF of Chi-square
// equivalent to Excel's CHISQ.DIST(fV,df,TRUE)
double CDFchisq(double fV, uint64_t nDF)
{
    switch (nDF)
    {
        case 1: return 2. * Phi(sqrt(fV)) - 1.;
        case 2: return 1. - exp(-fV / 2.);
        default: break;
    }

    return CDFchisq(fV, nDF - 2) - 2*PDFchisq(fV, nDF);
}

// Return value: p(observed fV or higher) can happen by chance
// The value from the chi-squared distribution for the statistic and the appropriate degrees of freedom. 
double ChiSqTest(double fV, uint64_t nDF)
{
    return 1. - CDFchisq(fV, nDF);
}

// ==========================================================================
// Uniformity test for StreamSamplers
// ==========================================================================

template <typename Alg> 
bool StreamSamplerTestUnif()
{
    const size_t   nTrailSets  =  1000; // number of trail sets
    const uint64_t nTrails     =  1000; // number of trails per trail set
    const uint64_t nVals       = 10000; // length of stream
    const uint64_t nSampleSize =    10; // size of sample set
    const uint64_t nBins       =   100; // vnCounts supposed to be splitted uniformally among bins

    vector<double> vfX1(nTrailSets);       // for each trail set: p(observed fV or higher) can happen by chance (with the assumption of independence) (for large nTrails*nSampleSize/nBin)
    vector<double> vfX2(nTrailSets);

    for (size_t nTrailSet = 0; nTrailSet < nTrailSets; ++nTrailSet) // for each trail set
    {
        vector<uint64_t> vnCounts1(nBins); // count sum(instances per bin) for all trails. (sum(vnCounts) = nTrails * nSampleSize)
        vector<uint64_t> vnCounts2(nBins);

        Alg Sampler1(1, nSampleSize);      // one sample set of nSampleSize samples
        Alg Sampler2(1, nSampleSize);

        for (uint64_t nTrail = 0; nTrail < nTrails; ++nTrail)       // for each trail
        {
            // ----- type 1 test
            uint64_t i1 = 0;

            while (i1 < nVals)
                i1 += 1 + Sampler1.AddElement(i1 % nBins        );

            auto SS1 = Sampler1.GetSampleSets(); // also reset the samples
            for (auto& Set : SS1)                // for each sample set
                for (auto& n : Set)              // for each sample
                    ++vnCounts1[(size_t)n];      // inc count of this value

            // ----- type 2 test
            uint64_t i2 = 0;

            while (i2 < nVals)
                i2 += 1 + Sampler2.AddElement(i2 / (nVals / nBins));

            auto SS2 = Sampler2.GetSampleSets(); // also reset the samples
            for (auto& Set : SS2)                // for each sample set
                for (auto& n : Set)              // for each sample
                    ++vnCounts2[(size_t)n];      // inc count of this value
        }

        // check for uniformity (same count per bin) using chi-square test
        // see "The Art of Computer Programming" [Knuth] Vol.2, 3.3.1.
        // fV is calculated according to page 43 formula (8)
        double   fV1 = inner_product(begin(vnCounts1), end(vnCounts1), begin(vnCounts1), 0.) / (nTrails * nSampleSize / nBins) - (nTrails * nSampleSize);
        double   fV2 = inner_product(begin(vnCounts2), end(vnCounts2), begin(vnCounts2), 0.) / (nTrails * nSampleSize / nBins) - (nTrails * nSampleSize);
        uint64_t nDF = nBins - 1;            // degrees of freedom
        
        // p(observed fV or higher) can happen by chance (with assumption of independence) (for large nTrails*nSampleSize/nBin)
        vfX1[nTrailSet] = ChiSqTest(fV1, nDF);
        vfX2[nTrailSet] = ChiSqTest(fV2, nDF);
    }

    // avrg vfX for all trail sets. should be around 0.5
    double fAvrgP1 = accumulate(vfX1.begin(), vfX1.end(), 0.) / vfX1.size();
    cout << "Test 1: Avrg p(observed fV or higher): " << fAvrgP1 << "     Should be around 0.5 ---> ";
    cout << ((fabs(fAvrgP1 - 0.5) < 0.05) ? "PASS\n" : "FAIL\n");

    double fAvrgP2 = accumulate(vfX2.begin(), vfX2.end(), 0.) / vfX2.size();
    cout << "Test 2: Avrg p(observed fV or higher): " << fAvrgP2 << "     Should be around 0.5 ---> ";
    cout << ((fabs(fAvrgP2 - 0.5) < 0.05) ? "PASS\n" : "FAIL\n");
    return true;
}

// --------------------------------------------------------------------------
bool StreamSamplerTestUniformity()
{
    cout << "Uniformity Test\n=============\n";
    cout << "Algorithm M: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_M <uint64_t>>();
    cout << "Algorithm L: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_L <uint64_t>>();
    cout << "Algorithm K: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_K <uint64_t>>();
    cout << "Algorithm Z: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_Z <uint64_t>>();
    cout << "Algorithm Y: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_Y <uint64_t>>();
    cout << "Algorithm X: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_X <uint64_t>>();
    cout << "Algorithm R: \n"; StreamSamplerTestUnif<CStreamSamplerWOR_R <uint64_t>>();
    cout << "Algorithm R0:\n"; StreamSamplerTestUnif<CStreamSamplerWOR_R0<uint64_t>>();
    cout << "\n";
    return true;
}

// ==========================================================================
// Performance benchmark for StreamSamplers
// ==========================================================================

template <typename Alg> 
bool StreamSamplerTestPer(size_t nVals, size_t nSampleSize) // length of stream, size of sample set
{
    uint64_t nTrail   = 0                                   ;
    auto     base     = chrono::high_resolution_clock::now();
    auto     TimeDiff = 0ns                                 ;
    
    while (TimeDiff < 10s) // 10 seconds
    {
        ++nTrail;
        Alg Sampler(1, nSampleSize, nTrail); // one sample set of nSampleSize samples, different seed for each sample set
        uint64_t i = 0;
        while (i < nVals)
            i += 1 + Sampler.AddElement(i);

        TimeDiff = chrono::high_resolution_clock::now() - base;
    }

    cout << (double)nTrail / chrono::duration_cast<chrono::nanoseconds>(TimeDiff).count() * 1e9 << " trails/sec\n";
    return true;
}

// --------------------------------------------------------------------------
bool StreamSamplerPerformanceBenchmark()
{
    const size_t nVals       = 100000000; // length of stream
    const size_t nSampleSize =      1000; // size of sample set

    cout << "Performance Test (" << nSampleSize << " / " << nVals << ")\n===========================\n";
    cout << "Algorithm M:  "; StreamSamplerTestPer<CStreamSamplerWOR_M <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm L:  "; StreamSamplerTestPer<CStreamSamplerWOR_L <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm K:  "; StreamSamplerTestPer<CStreamSamplerWOR_K <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm Z:  "; StreamSamplerTestPer<CStreamSamplerWOR_Z <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm Y:  "; StreamSamplerTestPer<CStreamSamplerWOR_Y <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm X:  "; StreamSamplerTestPer<CStreamSamplerWOR_X <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm R:  "; StreamSamplerTestPer<CStreamSamplerWOR_R <uint64_t>>(nVals, nSampleSize);
    cout << "Algorithm R0: "; StreamSamplerTestPer<CStreamSamplerWOR_R0<uint64_t>>(nVals, nSampleSize);
    cout << "\n";
    return true;
}

// ==========================================================================
// static auto b = StreamSamplerTestUniformity() && StreamSamplerPerformanceBenchmark();
