// Copyright 2015 Lior Kogan (koganlior1@gmail.com)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software distributed under the License is 
// distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
// See the License for the specific language governing permissions and limitations under the License.

#pragma once

#include <random>
#include <vector>
#include <algorithm>  // min_element

using namespace std;

// ==========================================================================
// CStreamSampler
// ==========================================================================

// Abstract base class for stream samplers
// A stream sampler extract a given number of independent sample sets, each with a given number of elements, from a stream.
// The number of elements in the stream is usually very large and unknown. Only one pass is possible.
template <typename ElementType, typename RNE = mt19937_64> // Random Number Engine
class CStreamSampler
{
public:
    CStreamSampler(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                   size_t                    nSetSize            ,   // [i] size of each sample set
                   typename RNE::result_type nSeed = RandomSeed() ); // [i] Random Number Engine seed

    // return the number of future stream elements the caller should skip before calling AddElement again
    virtual uint64_t AddElement(const ElementType&  Sample) = 0;

    // return the number of future stream elements the caller should skip before calling AddElement again
    // if an element is sampled into more then one set, it will be moved first and then copied
    // caller should call AddElement(move(s)) or simply AddElement(s) if s is a temporary object
    virtual uint64_t AddElement(      ElementType&& Sample) = 0; 
    
    // get the SampleSets and reset it (implemented using move semantics)
    vector<vector<ElementType>> GetSampleSets();

protected:
    bool                        m_bValid     ; // invalidated when calling GetSampleSets(). validated on costruction, and on AddElement() if invalid
    vector<RNE>                 m_vRndGen    ; // random number engine. see note in constructor's impl.

    const size_t                m_nSampleSets; // number of independent sample sets
    const size_t                m_nSetSize   ; // size of each sample set
    uint64_t                    m_nElements  ; // number of stream elements seen so far
    vector<vector<ElementType>> m_vSampleSets; // for each sample set: vector of samples (reservoir)

    virtual void Reset();

    static typename RNE::result_type RandomSeed()
    {
        random_device rd;

        if (is_same<RNE::result_type, uint32_t>::value)      // if RNE::result_type is uint32_t: return 32-bit seed
            return rd();
        else if (is_same<RNE::result_type, uint64_t>::value) // if RNE::result_type is uint64_t: return 64-bit seed
            return ((uint64_t)rd() << 32) | rd();
        else
            throw invalid_argument("please modify CStreamSampler::RandomSeed() according to the RNE used");
    }
};

// ==========================================================================
// CStreamSamplerWOR_R0
// ==========================================================================

// sample without replacement (WOR)
// each tuple has the same probability for being selected
// Based on "The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
// In this implementation AddElement always return 0 (no skip)
template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_R0 : public CStreamSampler<ElementType, RNE>
{
public:
    CStreamSamplerWOR_R0(size_t                    nSampleSets         ,  // [i] number of independent sample sets
                         size_t                    nSetSize            ,  // [i] size of each sample set
                         typename RNE::result_type nSeed = RandomSeed() ) // [i] RNE seed
        : CStreamSampler(nSampleSets, nSetSize, nSeed) {}

    // return 0 (no skip)
    virtual uint64_t AddElement(const ElementType&  Sample) override;
    virtual uint64_t AddElement(      ElementType&& Sample) override;
};

// ==========================================================================
// CStreamSamplerWOR
// ==========================================================================

// sample without replacement (WOR)
// each tuple has the same probability for being selected
template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR : public CStreamSampler<ElementType, RNE>
{
public:
    CStreamSamplerWOR(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                      size_t                    nSetSize            ,   // [i] size of each sample set
                      typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

    // return the number of future stream elements the caller should skip before calling AddElement again
    virtual uint64_t AddElement(const ElementType&  Sample) override;
    virtual uint64_t AddElement(      ElementType&& Sample) override;

protected:
    vector<uint64_t> m_vnSkip   ; // for each sample set: number of next elements to skip
    vector<size_t  > m_vnNextIdx; // for each sample set: next index to fill in reservoir
    uint64_t         m_nNextSkip; // next stream skip size = min(m_vnSkip)

    virtual void Reset() override;

    // draw number of next elements to skip and next index to replace
    virtual void DrawNext(size_t nSampleSetIdx) = 0; // [i] idx of sample set
};

// ==========================================================================
// Algorithms:
//
// R    : Based on "The Art of Computer Programming" [Knuth] 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman
//        Modified according to Ex.10, and such that AddElement calculates the number of future stream elements to skip
// X,Y,Z: Based on "Random Sampling with a Reservoir" [Jeferey Scott Vitter, 1985]
// K,L,M: Based on "Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))" [Kim-Hung Li ,1994]

// ==========================================================================
// Algorithm R
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_R : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_R(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm X
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_X : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_X(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm Y
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_Y : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_Y(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm Z
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_Z : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_Z(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm K
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_K : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_K(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    const double   m_fHs;
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm L
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_L : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_L(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
// Algorithm M
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_M : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_M(size_t                    nSampleSets         ,   // [i] number of independent sample sets
                        size_t                    nSetSize            ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed() ); // [i] RNE seed

private:
    uint64_t         nR                 ;
    vector<bool    > m_vbStep2          ;
    vector<double  > m_vfU, m_vfW, m_vfQ;
    vector<uint64_t> m_vnT, m_vnCount   ;

    void DrawNext(size_t nSampleSetIdx) override;    // [i] idx of sample set
};

// ==========================================================================
#include "StreamSampler.tpp"
