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
#include <stdexcept>  // invalid_argument

// ==========================================================================
namespace StreamSampler {

// A stream is a sequence of data elements made available over time.
// The number of elements in the stream is usually large and unknown a priori.

// A stream sampler extracts a sample set with a given size from a stream.
// Each possible sample set (of the given size) has an equal probability of being extracted. 
// A stream sampler is an online algorithm: The size of the input is unknown, and only one pass over the stream is possible.

using namespace std;

// ==========================================================================
// RandomSeed
// ==========================================================================

template <typename RNE> // Random Number Engine
static typename RNE::result_type RandomSeed()
{
    random_device rd;

    if (is_same<typename RNE::result_type, uint32_t>::value)      // if RNE::result_type is uint32_t: return 32-bit seed
        return rd();
    else if (is_same<typename RNE::result_type, uint64_t>::value) // if RNE::result_type is uint64_t: return 64-bit seed
        return ((uint64_t)rd() << 32) | rd();
    else
        throw invalid_argument("please modify StreamSampler::RandomSeed() according to the RNE used");
}

// ==========================================================================
// CStreamSampler
// ==========================================================================

// Abstract base class for stream samplers
template <typename ElementType, typename RNE = mt19937_64> // Random Number Engine
class CStreamSampler
{
public:
    CStreamSampler(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                   size_t                    nSetSize                 ,   // [i] size of each sample set
                   typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] Random Number Engine seed

    // return the number of future stream elements the caller should skip before calling AddElement again
    virtual uint64_t AddElement(const ElementType&  Element) = 0;

    // return the number of future stream elements the caller should skip before calling AddElement again
    // if an element is sampled into more then one set, it will be moved first and then copied
    // caller should call AddElement(move(s)) or simply AddElement(s) if Sample is a temporary object
    virtual uint64_t AddElement(      ElementType&& Element) = 0;
    
    // get the SampleSets and reset it (implemented using move semantics)
    vector<vector<ElementType>> GetSampleSets();

    virtual void Reset();

protected:
    bool                        m_bValid     ; // invalidated when calling GetSampleSets(). validated on construction, and on AddElement() if invalid
    const size_t                m_nSampleSets; // number of independent sample sets
    const size_t                m_nSetSize   ; // size of each sample set
    vector<vector<ElementType>> m_vSampleSets; // for each sample set: vector of samples (reservoir)
    uint64_t                    m_nElements  ; // number of stream elements seen so far
    vector<RNE>                 m_vRndGen    ; // for each sample set: random number engine. see note in constructor's impl.
};

// ==========================================================================
// CStreamSamplerWOR_R0
// ==========================================================================

// Basic stream sampler without replacement (WOR)
// Based on "The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_R0 : public CStreamSampler<ElementType, RNE>
{
public:
    CStreamSamplerWOR_R0(size_t                    nSampleSets              ,  // [i] number of independent sample sets
                         size_t                    nSetSize                 ,  // [i] size of each sample set
                         typename RNE::result_type nSeed = RandomSeed<RNE>() ) // [i] RNE seed
        : CStreamSampler<ElementType, RNE>(nSampleSets, nSetSize, nSeed) {}

    // always return 0 (no skip)
    virtual uint64_t AddElement(const ElementType&  Element) override;
    virtual uint64_t AddElement(      ElementType&& Element) override;
};

// ==========================================================================
// CStreamSamplerWOR
// ==========================================================================

// Abstract base class for more advanced stream samplers without replacement (WOR)
template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR : public CStreamSampler<ElementType, RNE>
{
public:
    CStreamSamplerWOR(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                      size_t                    nSetSize                 ,   // [i] size of each sample set
                      typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

    // return the number of future stream elements the caller should skip before calling AddElement again
    virtual uint64_t AddElement(const ElementType&  Element) override;
    virtual uint64_t AddElement(      ElementType&& Element) override;

    virtual void Reset() override;

protected:
    vector<uint64_t> m_vnSkip   ; // for each sample set: number of next elements to skip
    vector<size_t  > m_vnNextIdx; // for each sample set: next index to fill in reservoir
    uint64_t         m_nNextSkip; // next stream skip size = min(m_vnSkip)

    // draw number of next elements to skip and next index to replace
    virtual void DrawNext(size_t nSampleSetIdx) = 0; // [i] idx of sample set
};

// ==========================================================================
// Algorithms:
//
// R    : Based on "The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
//        Implemented such that AddElement calculates the number of future stream elements to skip
// X,Y,Z: Based on "Random Sampling with a Reservoir" [Jeferey Scott Vitter, 1985]
// K,L,M: Based on "Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))" [Kim-Hung Li, 1994]

// ==========================================================================
// Algorithm R
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_R : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_R(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm X
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_X : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_X(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm Y
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_Y : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_Y(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

private:
    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm Z
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_Z : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_Z(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

    virtual void Reset() override;

private:
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm K
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_K : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_K(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

    virtual void Reset() override;

private:
    const double   m_fHs;
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm L
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_L : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_L(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

    virtual void Reset() override;

private:
    vector<double> m_vfW;

    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};

// ==========================================================================
// Algorithm M
// ==========================================================================

template <typename ElementType, typename RNE = mt19937_64>
class CStreamSamplerWOR_M : public CStreamSamplerWOR<ElementType, RNE>
{
public:
    CStreamSamplerWOR_M(size_t                    nSampleSets              ,   // [i] number of independent sample sets
                        size_t                    nSetSize                 ,   // [i] size of each sample set
                        typename RNE::result_type nSeed = RandomSeed<RNE>() ); // [i] RNE seed

    virtual void Reset() override;

private:
    const uint64_t   nR                 ;
    vector<bool    > m_vbStep2          ;
    vector<double  > m_vfU, m_vfW, m_vfQ;
    vector<uint64_t> m_vnT, m_vnCount   ;

    void DrawNext(size_t nSampleSetIdx) override; // [i] idx of sample set
};



// **************************************************************************
// implementation
// **************************************************************************



// ==========================================================================
// CStreamSampler
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSampler<ElementType, RNE>::CStreamSampler(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                 size_t                    nSetSize   ,  // [i] size of each sample set
                                                 typename RNE::result_type nSeed       ) // [i] type of seed: as type of 1st class template parameter of mt19937_64
    : m_bValid(true), 
      m_nSampleSets(nSampleSets), m_nSetSize(nSetSize),
      m_vSampleSets(nSampleSets, vector<ElementType>(nSetSize)),
      m_nElements(0)
{
    if (0 == nSampleSets) throw invalid_argument("Stream Sampler: 0 sample sets");
    if (0 == nSetSize   ) throw invalid_argument("Stream Sampler: sample size 0");

    // Note: why using many RNEs?
    // When using many sample sets, it is important to use an arbitrary k-dimensional equidistribution RNE
    // (where every possible k-tuple of RNs will occur, and they will all occur the same number of times)
    // Mersenne Twister 19337 is 623-dimensionally 32-bit / 312-dimensionally 64-bit equidistributed
    // (see "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator" [Matsumoto, Nishimura 1998])
    // If there are more than 312 sample sets, using the same RNE for all of them may result in biased samples. 
    // Therefore we prefer to use a different RNE (different seed) for each sample set.
    m_vRndGen.resize(nSampleSets);
    for (auto& RndGen : m_vRndGen)
        RndGen.seed(nSeed++);
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE> 
void CStreamSampler<ElementType, RNE>::Reset()
{
    if (m_vSampleSets.size() == 0)
        m_vSampleSets.assign(m_nSampleSets, vector<ElementType>(m_nSetSize));

    m_nElements = 0;
    m_bValid    = true;
}

// --------------------------------------------------------------------------
// get the SampleSets and reset it (implemented using move semantics)
template <typename ElementType, typename RNE> 
vector<vector<ElementType>> CStreamSampler<ElementType, RNE>::GetSampleSets()
{
    m_bValid = false;
    return move(m_vSampleSets);
}

// ==========================================================================
// CStreamSamplerWOR_R0
// ==========================================================================

template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR_R0<ElementType, RNE>::AddElement(const ElementType& Element)
{
    if (!this->m_bValid)
        this->Reset();

    if (this->m_nElements < this->m_nSetSize)                              // first m_nSetSize elements
        for (auto& vSampleSet : this->m_vSampleSets)
            vSampleSet[this->m_nElements] = Element;                       // copy element into each SampleSet
    else
        for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)             // for each sample set
        {
            auto r = uniform_int_distribution<uint64_t>(0, this->m_nElements)(this->m_vRndGen[nSS]); // inclusive range
            if (r < this->m_nSetSize)
                this->m_vSampleSets[nSS][r] = Element;                     // copy element
        }

    ++this->m_nElements;
    return 0;
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR_R0<ElementType, RNE>::AddElement(ElementType&& Element)
{
    if (!this->m_bValid)
        this->Reset();

    if (this->m_nElements < this->m_nSetSize)                              // first m_nSetSize elements
    {   
        for (size_t nSS = 1; nSS < this->m_nSampleSets; ++nSS)             // copy element into each SampleSet except the 1st
            this->m_vSampleSets[nSS][(size_t)this->m_nElements] = Element;

        this->m_vSampleSets[0][(size_t)this->m_nElements] = move(Element); // move element into 1st SampleSet
    }
    else
    {
        ElementType* pS = nullptr;
        for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)             // for each sample set
        {
            auto r = uniform_int_distribution<uint64_t>(0, this->m_nElements)(this->m_vRndGen[nSS]); // inclusive range
            if (r < this->m_nSetSize)
            {
                if (pS)
                    this->m_vSampleSets[nSS][r] = *pS;                     // copy element
                else
                {
                    pS  = &this->m_vSampleSets[nSS][r];
                    *pS = move(Element);                                   // move element
                }
            }
        }
    }
         
    ++this->m_nElements;
    return 0;
}

// ==========================================================================
// CStreamSamplerWOR
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR<ElementType, RNE>::CStreamSamplerWOR(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                       size_t                    nSetSize   ,  // [i] size of each sample set
                                                       typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSampler<ElementType, RNE>(nSampleSets, nSetSize, nSeed),
      m_vnSkip(nSampleSets), m_vnNextIdx(nSampleSets, -1), m_nNextSkip(0)
    {}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR<ElementType, RNE>::Reset()
{
    CStreamSampler<ElementType, RNE>::Reset();

    m_nNextSkip = 0;
    m_vnSkip     .assign(this->m_nSampleSets,  0);
    m_vnNextIdx  .assign(this->m_nSampleSets, -1);
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR<ElementType, RNE>::AddElement(const ElementType& Element)
{
    if (!this->m_bValid)
        Reset();

    if (this->m_nElements < this->m_nSetSize)                             // first m_nSetSize elements: fill reservoir
        for (auto& vSampleSet : this->m_vSampleSets)
            vSampleSet[(size_t)this->m_nElements] = Element;              // copy element into each SampleSet

    if (++this->m_nElements >= this->m_nSetSize)                          // after initial fill of reservoir
    {
        this->m_nElements += m_nNextSkip;                                 // update according to number of elements skipped

        for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)            // for each sample set
            if (m_vnSkip[nSS] -= m_nNextSkip)                             // decrease by number of elements skipped. still >0 ?
                --m_vnSkip[nSS];                                          // skip current element
            else
            {
                if (m_vnNextIdx[nSS] != (size_t)-1)                       // not on the (m_nSetSize-1)'th element
                    this->m_vSampleSets[nSS][m_vnNextIdx[nSS]] = Element; // copy element

                DrawNext(nSS);                                            // draw number of next elements to skip and next index to replace
            }
    }

    m_nNextSkip = *min_element(begin(m_vnSkip), end(m_vnSkip));
    return m_nNextSkip;
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR<ElementType, RNE>::AddElement(ElementType&& Element)
{
    if (!this->m_bValid)
        Reset();

    if (this->m_nElements < this->m_nSetSize)                              // first m_nSetSize elements: fill reservoir
    {
        for (size_t nSS = 1; nSS < this->m_nSampleSets; ++nSS)             // copy element into each SampleSet except the 1st
            this->m_vSampleSets[nSS][(size_t)this->m_nElements] = Element;

        this->m_vSampleSets[0][(size_t)this->m_nElements] = move(Element); // move element into 1st SampleSet
    }

    if (++this->m_nElements >= this->m_nSetSize)                           // after initial fill of reservoir
    {
        this->m_nElements += m_nNextSkip;                                  // update according to number of elements skipped

        ElementType* pS = nullptr;
        for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)             // for each sample set
            if (m_vnSkip[nSS] -= m_nNextSkip)                              // decrease by number of elements skipped. still >0 ?
                --m_vnSkip[nSS];                                           // skip current element
            else
            {
                if (m_vnNextIdx[nSS] != (size_t)-1)                        // not on the (m_nSetSize-1)'th element 
                {
                    if (pS)
                        this->m_vSampleSets[nSS][m_vnNextIdx[nSS]] = *pS;  // copy element
                    else
                    {
                        pS  = &this->m_vSampleSets[nSS][m_vnNextIdx[nSS]];
                        *pS = move(Element);                               // move element
                    }
                }

                DrawNext(nSS);                                             // draw number of next elements to skip and next index to replace
            }
    }

    m_nNextSkip = *min_element(begin(m_vnSkip), end(m_vnSkip));
    return m_nNextSkip;
}

// ==========================================================================
// Algorithm R
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_R<ElementType, RNE>::CStreamSamplerWOR_R(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_R<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    while (1)
    {
        uint64_t nNextIdx = uniform_int_distribution<uint64_t>(0, this->m_nElements + this->m_vnSkip[nSS])(this->m_vRndGen[nSS]); // inclusive range
        if (nNextIdx < this->m_nSetSize)
        {
            this->m_vnNextIdx[nSS] = (size_t)nNextIdx;
            break;
        }

        ++this->m_vnSkip[nSS];
    }
}

// ==========================================================================
// Algorithm X
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_X<ElementType, RNE>::CStreamSamplerWOR_X(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_X<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    double fHs = (double)(this->m_nElements + 1 - this->m_nSetSize) / (this->m_nElements + 1);
    double fV  = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);

    while (fHs > fV) // increase skip till fHs <= fV
    {
        ++this->m_vnSkip[nSS];
        fHs *= (double)(this->m_nElements + 1 - this->m_nSetSize + this->m_vnSkip[nSS]) / (this->m_nElements + 1 + this->m_vnSkip[nSS]);
    }

    this->m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]); // draw next index to replace
}

// ==========================================================================
// Algorithm Y
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_Y<ElementType, RNE>::CStreamSamplerWOR_Y(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_Y<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    double fHs = (double)(this->m_nElements + 1 - this->m_nSetSize) / (this->m_nElements + 1);
    double fV  = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);

    while (fHs > fV) // increase skip till fHs <= fV
    {
        // naive
        // double fHsPlus1 = fHs * (m_nElements + 1 - m_nSetSize + m_vnSkip[n] + 1) / (m_nElements + 1 + m_vnSkip[n] + 1); // H(s + 1)
        // double fDeltaS  = - (fHs - fV) / (fHsPlus1 - fHs); // Newton interpolation with discrete derivative
 
        // optimized
        const uint64_t nZ     = this->m_nElements + 1 + this->m_vnSkip[nSS];
        const double  fDeltaS = (fHs - fV) * (nZ + 1) / (fHs * this->m_nSetSize);
        const uint64_t nS     = (uint64_t)ceil(fDeltaS);

        // naive (no advantage over Algorithm X)
        // for (uint64_t i = 1; i <= nS; ++i)
        //     fHs *= (double)(nZ - m_nSetSize + i) / (nZ + i);
        
        // optimizated: when nS >= m_nSetSize we can cancel equal terms in the numerator and in the denominator
        // so the loop size is always <= m_nSetSize
        if (nS < this->m_nSetSize)
            for (uint64_t i = 1; i <= nS; ++i)
                fHs *= (double)(nZ - this->m_nSetSize + i) / (nZ + i);
        else
            for (uint64_t i = 0; i < this->m_nSetSize; ++i)
                fHs *= (double)(nZ - i) / (nZ + nS - i);

        this->m_vnSkip[nSS]+= nS;
    }

    this->m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]); // draw next index to replace
}

// ==========================================================================
// Algorithm Z
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_Z<ElementType, RNE>::CStreamSamplerWOR_Z(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed),
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS < nSampleSets; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / nSetSize); // initial W = exp(–log(random()) / nSetSize)
    }

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR_Z<ElementType, RNE>::Reset()
{
    CStreamSamplerWOR<ElementType, RNE>::Reset();

    for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)
        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / this->m_nSetSize); // initial W = exp(–log(random()) / nSetSize)
}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_Z<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    uint64_t& nS    = this->m_vnSkip[nSS];
    uint64_t  nTerm = this->m_nElements - this->m_nSetSize + 1;

    if (this->m_nElements <= 22 * this->m_nSetSize) // execute algorithm X
    {
        double fHs = (double)nTerm / (this->m_nElements + 1);
        double fV  = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);

        while (fHs > fV) // increase skip till fHs <= fV
        {
            ++nS;
            fHs *= (double)(nTerm+nS) / (this->m_nElements + 1 + nS);
        }
    }
    else while (1)
    {
        double fU = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);
        double fX = (m_vfW[nSS] - 1.) * this->m_nElements;
        nS = (uint64_t)floor(fX);

        // test if u<=h(s)/cg(x)
        double fRHS  = (this->m_nElements + fX) / (nTerm+nS) * nTerm;
        double fLHS  = pow(fU * (this->m_nElements + 1.) / fRHS * (this->m_nElements + 1.) / nTerm, 1. / this->m_nSetSize);
               fRHS /= this->m_nElements;

        if (fLHS <= fRHS)
        {
            m_vfW[nSS] = fRHS / fLHS;
            break;
        }

        // test if u<=f(s)/cg(x)
        double   fY     = fU * (this->m_nElements + 1) / nTerm * (this->m_nElements + nS + 1) / (this->m_nElements + fX);
        uint64_t nDenom = nTerm;
        uint64_t nNumer = (this->m_nSetSize < nS) ? nTerm+nS : this->m_nElements + 1;
        while (nNumer <= this->m_nElements + nS)
            fY *= (double)nNumer++ / nDenom++;

        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / this->m_nSetSize); // generate W in advance
        if (pow(fY, 1. / this->m_nSetSize) <= (this->m_nElements + fX) / this->m_nElements)
            break;
    }

    this->m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]);    // draw next index to replace
}

// ==========================================================================
// Algorithm K
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_K<ElementType, RNE>::CStreamSamplerWOR_K(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed), 
      m_fHs((double)nSetSize/2.), // a,b in Kim-Hung Li's paper
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS < nSampleSets; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / nSetSize); // initial W = exp(–log(random()) / nSetSize)
    }

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR_K<ElementType, RNE>::Reset()
{
    CStreamSamplerWOR<ElementType, RNE>::Reset();

    for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)
        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / this->m_nSetSize); // initial W = exp(–log(random()) / nSetSize)
}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_K<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    uint64_t& nS    = this->m_vnSkip[nSS];
    uint64_t  nTerm = this->m_nElements - this->m_nSetSize + 1;

    if (this->m_nElements <= 14 * this->m_nSetSize) // execute algorithm X
    {
        double fHs = (double)nTerm / (this->m_nElements + 1);
        double fV  = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);

        while (fHs > fV) // increase skip till fHs <= fV
        {
            ++nS;
            fHs *= (double)(nTerm+nS) / (this->m_nElements + 1 + nS);
        }
    }
    else  while (1)
    {
        double fU = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);
        double fX = (m_vfW[nSS] - 1.) * (this->m_nElements - m_fHs);
        nS = (uint64_t)floor(fX);

        // test if u<=h(s)/cg(x)
        double fRHS = (this->m_nElements + fX - m_fHs) / (nTerm + nS) * nTerm / (this->m_nElements - m_fHs);
        double fLHS = pow(fU * (this->m_nElements + 1.) / (nTerm - 1) / fRHS, 1. / this->m_nSetSize);

        if (fLHS <= fRHS)
        {
            m_vfW[nSS] = fRHS / fLHS;
            break;
        }

        // test if u<=f(s)/cg(x)
        double   fY     = fU * (this->m_nElements - m_fHs) / (nTerm- 1 ) * (this->m_nElements + nS + 1) / (this->m_nElements - m_fHs + fX);
        uint64_t nDenom = nTerm;
        uint64_t nNumer = (this->m_nSetSize < nS) ? nTerm + nS : this->m_nElements + 1;
        while (nNumer <= this->m_nElements + nS)
            fY *= (double)nNumer++ / nDenom++;

        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), -1. / this->m_nSetSize); // generate W in advance
        if (pow(fY, 1. / this->m_nSetSize) <= (this->m_nElements - m_fHs + fX) / (this->m_nElements - m_fHs))
            break;
    }

    this->m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]);    // draw next index to replace
}

// ==========================================================================
// Algorithm L
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_L<ElementType, RNE>::CStreamSamplerWOR_L(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed),
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS < nSampleSets; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), 1. / nSetSize); // initial W = exp(–log(random()) / nSetSize)
    }

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR_L<ElementType, RNE>::Reset()
{
    CStreamSamplerWOR<ElementType, RNE>::Reset();

    for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)
        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), 1. / this->m_nSetSize); // initial W = exp(–log(random()) / nSetSize)
}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_L<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    this->m_vnSkip   [nSS]  = (uint64_t)floor(log(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS])) / log(1. - m_vfW[nSS]));
    this->m_vnNextIdx[nSS]  = uniform_int_distribution<size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]);           // draw next index to replace
          m_vfW      [nSS] *= pow(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]), 1. / this->m_nSetSize); // generate W in advance
}

// ==========================================================================
// Algorithm M
// ==========================================================================

// produce a random floating-point value according to a beta distribution
// Based on "Generating Beta Variates with Nonintegral Shape Parameters" [R. C. H. Cheng, 1978] algorithm BB
// handles only cases where min(a0, b0) > 1
template <typename RNE>
inline double beta_distribution(RNE& engine, double a0, double b0)
{
    double a     = min(a0, b0);
    double b     = max(a0, b0);
    double alpha = a + b;
    double beta  = sqrt((alpha - 2) / (2. * a * b - alpha));
    double gamma = a + 1 / beta;
    double w;

    while (1)
    {
        double u1 = uniform_real_distribution<double>(0, 1)(engine);
        double u2 = uniform_real_distribution<double>(0, 1)(engine);
        double v  = beta * log(u1 / (1. - u1));
               w  = a * exp(v);
        double z  = u1 * u1 * u2;
        double r  = gamma * v - 1.3862943611198906188344642429164; // const = log(4)
        double s  = a + r - w;

        if (s + 2.6094379124341003746007593332262 >= 5. * z) // const = log(5) + 1
            break;

        double t  = log(z);
        if (s >= t)
            break;

        if (r + alpha * log(alpha / (b + w)) >= t)
            break;
    }

    return (a == a0)  ?  w / (b + w) : b / (b + w);
}

// --------------------------------------------------------------------------
constexpr double fStreamSamplerWOR_M_Theta = 10.5; // see table 2 in Li's paper
constexpr double fStreamSamplerWOR_M_Tau   = 2.07; 

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
CStreamSamplerWOR_M<ElementType, RNE>::CStreamSamplerWOR_M(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR<ElementType, RNE>(nSampleSets, nSetSize, nSeed),
      nR((uint64_t)floor(fStreamSamplerWOR_M_Tau * sqrt(nSetSize))),
      m_vbStep2(nSampleSets), m_vfU(nSampleSets), m_vfW(nSampleSets), m_vfQ(nSampleSets), m_vnT(nSampleSets), m_vnCount(nSampleSets)
    {
        double fH = 0; // sum(1 / u); u = nSetSize ... nSetSize+nR
        for (uint64_t i = nSetSize; i <= nSetSize + nR; ++i)
            fH += 1. / i;
    
        uint64_t nC = (uint64_t)floor(fStreamSamplerWOR_M_Theta * (fStreamSamplerWOR_M_Tau * fStreamSamplerWOR_M_Tau / 2. + 1. + nR) / fH - nSetSize);

        for (size_t nSS = 0; nSS < nSampleSets; ++nSS)
        {           
            m_vfW    [nSS] = beta_distribution(this->m_vRndGen[nSS], nSetSize + 0., nC + 1.);
            m_vfU    [nSS] = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);
            m_vnCount[nSS] = nC;
        }
    }

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR_M<ElementType, RNE>::Reset()
{
    CStreamSamplerWOR<ElementType, RNE>::Reset();

    m_vbStep2.assign(this->m_nSampleSets, false);
    m_vfQ    .assign(this->m_nSampleSets, 0.   );
    m_vnT    .assign(this->m_nSampleSets, 0    );

    double fH = 0; // sum(1 / u); u = nSetSize ... nSetSize+nR
    for (uint64_t i = this->m_nSetSize; i <= this->m_nSetSize + nR; ++i)
        fH += 1. / i;
    
    uint64_t nC = (uint64_t)floor(fStreamSamplerWOR_M_Theta * (fStreamSamplerWOR_M_Tau * fStreamSamplerWOR_M_Tau / 2. + 1. + nR) / fH - this->m_nSetSize);

    for (size_t nSS = 0; nSS < this->m_nSampleSets; ++nSS)
    {           
        m_vfW    [nSS] = beta_distribution(this->m_vRndGen[nSS], this->m_nSetSize + 0., nC + 1.);
        m_vfU    [nSS] = uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS]);
        m_vnCount[nSS] = nC;
    }
}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_M<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    while (m_vnCount[nSS] > 0)
    {
        if (m_vbStep2[nSS])
            this->m_vnSkip[nSS] += (uint64_t)floor(log(uniform_real_distribution<double>(0, 1)(this->m_vRndGen[nSS])) / m_vfQ[nSS]);

        --m_vnCount [nSS];
        m_vfU       [nSS] *= (1. + (double)this->m_nSetSize / ++m_vnT[nSS]);

        if (1 <= m_vfU[nSS])
        {
                  m_vfU      [nSS] = uniform_real_distribution<double>(0, 1                   )(this->m_vRndGen[nSS]);
            this->m_vnNextIdx[nSS] = uniform_int_distribution<size_t >(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]); // draw next index to replace
            return;
        }
            
        ++this->m_vnSkip[nSS];
    }

          m_vbStep2  [nSS]  = true;
          m_vnCount  [nSS]  = nR;
          m_vfQ      [nSS]  = log(1. - m_vfW[nSS]);
    this->m_vnSkip   [nSS] += (uint64_t)floor(log(m_vfU[nSS]) / m_vfQ[nSS]);
          m_vnT      [nSS]  = 0;
          m_vfW      [nSS] *= beta_distribution(this->m_vRndGen[nSS], this->m_nSetSize + 0., nR + 1.);
          m_vfU      [nSS]  = uniform_real_distribution<double>(0, 1                   )(this->m_vRndGen[nSS]);    
    this->m_vnNextIdx[nSS]  = uniform_int_distribution <size_t>(0, this->m_nSetSize - 1)(this->m_vRndGen[nSS]); // draw next index to replace
}

// ==========================================================================
} // namespace StreamSampler
