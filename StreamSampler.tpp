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
// CStreamSampler
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSampler<ElementType, RNE>::CStreamSampler(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                 size_t                    nSetSize   ,  // [i] size of each sample set
                                                 typename RNE::result_type nSeed       ) // [i] type of seed: as type of 1st class template parameter of mt19937_64
    : m_bValid(true),
        m_nSampleSets(nSampleSets), m_nSetSize(nSetSize), m_nElements(0), 
        m_vSampleSets(nSampleSets, vector<ElementType>(nSetSize))
{
    if (0 == nSampleSets) throw invalid_argument("Stream Sampler: 0 sample sets"); 
    if (0 == nSetSize   ) throw invalid_argument("Stream Sampler: sample size 0"); 

    // Note: why using many RNEs?
    // When many sample sets may be used, it is important to use an arbitrary k-dimensional equidistribution RNE
    // (where every possible k-tuple will occur, and they will all occur the same number of times)
    // Mersenne Twister 19337 is 623-dimensionally 32-bit / 312-dimensionally 64-bit equidistributed
    // (see "Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator" [Matsumoto, Nishimura 1998])
    // If there are >312 sample sets, using the same RNE for all of them may result in biased samples. 
    // Therefore we use a different RNE for each sample set.
    m_vRndGen.resize(nSampleSets);
    for (auto& RndGen : m_vRndGen)
        RndGen.seed(nSeed++);
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE> 
void CStreamSampler<ElementType, RNE>::Reset()
{
    m_nElements = 0;
    m_vSampleSets.assign(m_nSampleSets, vector<ElementType>(m_nSetSize));
    m_bValid    = true;
}

// --------------------------------------------------------------------------
// get the SampleSets and reset it (implemented using move semantics)
template <typename ElementType, typename RNE> 
vector<vector<ElementType>> CStreamSampler<ElementType, RNE>::GetSampleSets()
{
    auto v = move(m_vSampleSets);
    m_bValid = false;
    return v;
}

// ==========================================================================
// CStreamSamplerWOR_R0
// ==========================================================================

template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR_R0<ElementType, RNE>::AddElement(const ElementType& Sample)
{
    if (!m_bValid)
        Reset();

    if (m_nElements < m_nSetSize)                            // first m_nSetSize elements
        for (auto& vSampleSet : m_vSampleSets)
            vSampleSet[m_nElements] = Sample;                // copy element into each SampleSet
    else
        for (size_t nSS = 0; nSS < m_nSampleSets; ++nSS)     // for each sample set
        {
            auto r = uniform_int_distribution<uint64_t>(0, m_nElements)(m_vRndGen[nSS]); // inclusive range
            if (r < m_nSetSize)
                    m_vSampleSets[nSS][r] = Sample;          // copy element
        }

    ++m_nElements;
    return 0;
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR_R0<ElementType, RNE>::AddElement(ElementType&& Sample)
{
    if (!m_bValid)
        Reset();

    if (m_nElements < m_nSetSize)                            // first m_nSetSize elements
    {   
        for (size_t nSS = 1; nSS<m_nSampleSets; ++nSS)       // copy element into each SampleSet except the 1st
            m_vSampleSets[ss][m_nElements] = Sample;

        m_vSampleSets[0][m_nElements] = move(Sample);        // move element into 1st SampleSet
    }
    else
    {
        ElementType* pS = nullptr;
        for (size_t nSS = 0; nSS < m_nSampleSets; ++nSS)     // for each sample set
        {
            auto r = uniform_int_distribution<uint64_t>(0, m_nElements)(m_vRndGen[nSS]); // inclusive range
            if (r < m_nSetSize)
                if (pS)
                    m_vSampleSets[nSS][r] = *pS;             // copy element
                else
                {
                    pS  = &m_vSampleSets[nSS][r];
                    *pS = move(Sample);                      // move element
                }
        }
    }
         
    ++m_nElements;
    return 0;
}

// ==========================================================================
// CStreamSamplerWOR
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR<ElementType, RNE>::CStreamSamplerWOR(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                       size_t                    nSetSize   ,  // [i] size of each sample set
                                                       typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSampler(nSampleSets, nSetSize, nSeed),
      m_vnSkip(nSampleSets), m_vnNextIdx(nSampleSets, -1), m_nNextSkip(0)
    {}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR<ElementType, RNE>::AddElement(const ElementType& Sample)
{
    if (!m_bValid)
        Reset();

    if (m_nElements < m_nSetSize)                                  // first m_nSetSize elements: fill reservoir
        for (auto& vSampleSet : m_vSampleSets)
            vSampleSet[(size_t)m_nElements] = Sample;              // copy element into each SampleSet

    if (++m_nElements >= m_nSetSize)                               // after initial fill of reservoir
    {
        m_nElements += m_nNextSkip;                                // update according to number of elements skipped

        for (size_t nSS = 0; nSS < m_nSampleSets; ++nSS)           // for each sample set
            if (m_vnSkip[nSS] -= m_nNextSkip)                      // decrease by number of elements skipped. still >0 ?
                --m_vnSkip[nSS];                                   // skip current element
            else
            {
                if (m_vnNextIdx[nSS] != -1)                        // not on the (m_nSetSize-1)'th element
                    m_vSampleSets[nSS][m_vnNextIdx[nSS]] = Sample; // copy element

                DrawNext(nSS);                                     // draw number of next elements to skip and next index to replace
            }
    }

    m_nNextSkip = *min_element(begin(m_vnSkip), end(m_vnSkip));
    return m_nNextSkip;
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
uint64_t CStreamSamplerWOR<ElementType, RNE>::AddElement(ElementType&& Sample)
{
    if (!m_bValid)
        Reset();

    if (m_nElements < m_nSetSize)                                    // first m_nSetSize elements: fill reservoir
    {
        for (size_t nSS = 1; nSS<m_nSampleSets; ++nSS)               // copy element into each SampleSet except the 1st
            m_vSampleSets[nSS][(size_t)m_nElements] = Sample;

        m_vSampleSets[0][(size_t)m_nElements] = move(Sample);        // move element into 1st SampleSet
    }

    if (++m_nElements >= m_nSetSize)                                 // after initial fill of reservoir
    {
        m_nElements += m_nNextSkip;                                  // update according to number of elements skipped

        ElementType* pS = nullptr;
        for (size_t nSS = 0; nSS < m_nSampleSets; ++nSS)             // for each sample set
            if (m_vnSkip[nSS] -= m_nNextSkip)                        // decrease by number of elements skipped. still >0 ?
                --m_vnSkip[nSS];                                     // skip current element
            else
            {
                if (m_vnNextIdx[nSS] != -1)                          // not on the (m_nSetSize-1)'th element 
                    if (pS)
                        m_vSampleSets[nSS][m_vnNextIdx[nSS]] = *pS;  // copy element
                    else
                    {
                        pS  = &m_vSampleSets[nSS][m_vnNextIdx[nSS]];
                        *pS = move(Sample);                          // move element
                    }

                DrawNext(nSS);                                       // draw number of next elements to skip and next index to replace
            }
    }

    m_nNextSkip = *min_element(begin(m_vnSkip), end(m_vnSkip));
    return m_nNextSkip;
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
void CStreamSamplerWOR<ElementType, RNE>::Reset()
{
    CStreamSampler::Reset();

    m_nNextSkip = 0;
    m_vnSkip     .assign(m_nSampleSets,  0);
    m_vnNextIdx  .assign(m_nSampleSets, -1);

    // todo: reset for alg Z,K
}

// ==========================================================================
// Algorithm R
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_R<ElementType, RNE>::CStreamSamplerWOR_R(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_R<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    while (1)
    {
        uint64_t nNextIdx = uniform_int_distribution<uint64_t>(0, m_nElements+m_vnSkip[nSS])(m_vRndGen[nSS]); // inclusive range
        if (nNextIdx < m_nSetSize)
        {
            m_vnNextIdx[nSS] = (size_t)nNextIdx;
            break;
        }

        ++m_vnSkip[nSS];
    }
}

// ==========================================================================
// Algorithm X
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_X<ElementType, RNE>::CStreamSamplerWOR_X(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_X<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    double fHs = (double)(m_nElements+1-m_nSetSize) / (m_nElements+1);
    double fV  = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);

    while (fHs > fV) // increase skip till fHs <= fV
    {
        ++m_vnSkip[nSS];
        fHs *= (double)(m_nElements+1-m_nSetSize+m_vnSkip[nSS]) / (m_nElements+1+m_vnSkip[nSS]);
    }

    m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]); // draw next index to replace
}

// ==========================================================================
// Algorithm Y
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_Y<ElementType, RNE>::CStreamSamplerWOR_Y(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed)
    {}

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_Y<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    double fHs = (double)(m_nElements+1-m_nSetSize) / (m_nElements+1);
    double fV  = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);

    while (fHs > fV) // increase skip till fHs <= fV
    {
        // naive
        // double fHsPlus1 = fHs*(m_nElements+1-m_nSetSize+m_vnSkip[n]+1) / (m_nElements+1+m_vnSkip[n]+1); // H(s+1)
        // double fDeltaS  = - (fHs-fV) / (fHsPlus1-fHs); // Newton interpolation with discrete derivative
 
        // optimized
        const uint64_t nZ     = m_nElements+1+m_vnSkip[nSS];
        const double  fDeltaS = (fHs-fV) * (nZ+1) / (fHs*m_nSetSize);
        const uint64_t nS     = (uint64_t)ceil(fDeltaS);

        // naive (no advantage over Algorithm X)
        // for (uint64_t i = 1; i<=nS; ++i)
        //     fHs *= (double)(nZ-m_nSetSize+i) / (nZ+i);
        
        // optimizated: when nS>=m_nSetSize we can cancel equal terms in the numerator and in the denominator
        // so the loop size is always <= m_nSetSize
        if (nS < m_nSetSize)
            for (uint64_t i = 1; i<=nS; ++i)
                fHs *= (double)(nZ-m_nSetSize+i) / (nZ+i);
        else
            for (uint64_t i = 0; i<m_nSetSize; ++i)
                fHs *= (double)(nZ-i) / (nZ+nS-i);

        m_vnSkip[nSS]+= nS;
    }

    m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]); // draw next index to replace
}

// ==========================================================================
// Algorithm Z
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_Z<ElementType, RNE>::CStreamSamplerWOR_Z(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed),
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS<nSampleSets ; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]), -1./nSetSize); // initial W
    }

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_Z<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    uint64_t& nS    = m_vnSkip[nSS];
    uint64_t  nTerm = m_nElements-m_nSetSize+1;

    if (m_nElements <= 22*m_nSetSize) // execute algorithm X
    {
        double fHs = (double)nTerm / (m_nElements+1);
        double fV  = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);

        while (fHs > fV) // increase skip till fHs <= fV
        {
            ++nS;
            fHs *= (double)(nTerm+nS) / (m_nElements+1+nS);
        }
    }
    else while (1)
    {
        double fU = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);
        double fX = (m_vfW[nSS]-1.)*m_nElements;
        nS = (uint64_t)floor(fX);

        // test if u<=h(s)/cg(x)
        double fRHS = (m_nElements+fX) / (nTerm+nS) * nTerm;
        double fLHS = pow(fU * (m_nElements+1.) / fRHS * (m_nElements+1.) / nTerm, 1./m_nSetSize);
               fRHS/= m_nElements;

        if (fLHS <= fRHS)
        {
            m_vfW[nSS] = fRHS/fLHS;
            break;
        }

        // test if u<=f(s)/cg(x)
        double   fY     = fU * (m_nElements+1) / nTerm * (m_nElements+nS+1) / (m_nElements+fX);
        uint64_t nDenom = nTerm;
        uint64_t nNumer = (m_nSetSize < nS) ? nTerm+nS : m_nElements+1;
        while (nNumer <= m_nElements+nS)
            fY *= (double)nNumer++/nDenom++; 

        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]), -1./m_nSetSize); // generate W in advance
        if (pow(fY, 1./m_nSetSize) <= (m_nElements+fX) / m_nElements)
            break;
    }

    m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]);          // draw next index to replace
}

// ==========================================================================
// Algorithm K
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_K<ElementType, RNE>::CStreamSamplerWOR_K(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed), 
      m_fHs((double)nSetSize/2.), // a,b in Kim-Hung Li's paper
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS<nSampleSets ; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]), -1./nSetSize); // initial W
    }

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_K<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    uint64_t& nS    = m_vnSkip[nSS];
    uint64_t  nTerm = m_nElements-m_nSetSize+1;

    if (m_nElements <= 14*m_nSetSize) // execute algorithm X
    {
        double fHs = (double)nTerm / (m_nElements+1);
        double fV  = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);

        while (fHs > fV) // increase skip till fHs <= fV
        {
            ++nS;
            fHs *= (double)(nTerm+nS) / (m_nElements+1+nS);
        }
    }
    else  while (1)
    {
        double fU = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);
        double fX = (m_vfW[nSS]-1.)*(m_nElements-m_fHs);
        nS = (uint64_t)floor(fX);

        // test if u<=h(s)/cg(x)
        double fRHS = (m_nElements+fX-m_fHs) / (nTerm+nS) * nTerm / (m_nElements-m_fHs);
        double fLHS = pow(fU * (m_nElements+1.) / (nTerm-1) / fRHS, 1./m_nSetSize);

        if (fLHS <= fRHS)
        {
            m_vfW[nSS] = fRHS/fLHS;
            break;
        }

        // test if u<=f(s)/cg(x)
        double   fY     = fU * (m_nElements-m_fHs) / (nTerm-1) * (m_nElements+nS+1) / (m_nElements-m_fHs+fX);
        uint64_t nDenom = nTerm;
        uint64_t nNumer = (m_nSetSize < nS) ? nTerm+nS : m_nElements+1;
        while (nNumer <= m_nElements+nS)
            fY *= (double)nNumer++/nDenom++; 

        m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]), -1./m_nSetSize); // generate W in advance
        if (pow(fY, 1./m_nSetSize) <= (m_nElements-m_fHs+fX) / (m_nElements-m_fHs))
            break;
    }

    m_vnNextIdx[nSS] = uniform_int_distribution<size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]);          // draw next index to replace
}

// ==========================================================================
// Algorithm L
// ==========================================================================

template <typename ElementType, typename RNE>
CStreamSamplerWOR_L<ElementType, RNE>::CStreamSamplerWOR_L(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed),
      m_vfW(nSampleSets)
    {
        for (size_t nSS = 0; nSS<nSampleSets ; ++nSS)
            m_vfW[nSS] = pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]),  1./nSetSize); // initial W
    }

// --------------------------------------------------------------------------
// draw number of next elements to skip and next index to replace
template <typename ElementType, typename RNE>
inline void CStreamSamplerWOR_L<ElementType, RNE>::DrawNext(size_t nSS) // [i] idx of sample set
{
    m_vnSkip   [nSS]  = (uint64_t)floor(log(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS])) / log(1.-m_vfW[nSS]));
    m_vnNextIdx[nSS]  = uniform_int_distribution<size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]);           // draw next index to replace
    m_vfW      [nSS] *= pow(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]), 1./m_nSetSize); // generate W in advance
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
    double a     = __min(a0, b0);
    double b     = __max(a0, b0);
    double alpha = a+b;
    double beta  = sqrt((alpha-2)/(2.*a*b-alpha));
    double gamma = a + 1/beta;
    double w;

    while (1)
    {
        double u1 = uniform_real_distribution<double>(0, 1)(engine);
        double u2 = uniform_real_distribution<double>(0, 1)(engine);
        double v  = beta * log(u1/(1.-u1));
               w  = a*exp(v);
        double z  = u1*u1*u2;
        double r  = gamma*v - 1.3862943611198906188344642429164; // const = log(4)
        double s  = a+r-w;

        if (s + 2.6094379124341003746007593332262 >= 5.*z) // const = log(5)+1
            break;

        double t  = log(z);
        if (s >= t)
            break;

        if (r + alpha * log(alpha/(b+w)) >= t)
            break;
    }

    return (a == a0)  ?  w / (b+w) : b / (b+w);
}

// --------------------------------------------------------------------------
template <typename ElementType, typename RNE>
CStreamSamplerWOR_M<ElementType, RNE>::CStreamSamplerWOR_M(size_t                    nSampleSets,  // [i] number of independent sample sets
                                                           size_t                    nSetSize   ,  // [i] size of each sample set
                                                           typename RNE::result_type nSeed       ) // [i] random seed
    : CStreamSamplerWOR(nSampleSets, nSetSize, nSeed),
      m_vfW(nSampleSets), m_vfU(nSampleSets), m_vfQ(nSampleSets), m_vnT(nSampleSets), m_vnCount(nSampleSets), m_vbStep2(nSampleSets)
    {
        double fTheta = 10.5; // see table 2 in Li's paper
        double fTau   = 2.07; 

        nR = (uint64_t)floor(fTau * sqrt(nSetSize));
        double fH = 0; // sum(1/u); u= nSetSize ... nSetSize+nR
        for (uint64_t i = nSetSize; i <= nSetSize+nR; ++i)
            fH += 1./i;
    
        uint64_t nC = (uint64_t)floor(fTheta*(fTau*fTau/2.+1.+nR)/fH - nSetSize);

        for (size_t nSS = 0; nSS<nSampleSets ; ++nSS)
        {           
            m_vfW    [nSS] = beta_distribution(m_vRndGen[nSS], nSetSize+0., nC+1.);
            m_vfU    [nSS] = uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS]);
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
            m_vnSkip[nSS] += (uint64_t)floor(log(uniform_real_distribution<double>(0, 1)(m_vRndGen[nSS])) / m_vfQ[nSS]);

        --m_vnCount [nSS];
        m_vfU       [nSS] *= (1. + (double)m_nSetSize / ++m_vnT[nSS]);

        if (1 <= m_vfU[nSS])
        {
            m_vfU      [nSS] = uniform_real_distribution<double>(0, 1           )(m_vRndGen[nSS]);
            m_vnNextIdx[nSS] = uniform_int_distribution<size_t >(0, m_nSetSize-1)(m_vRndGen[nSS]); // draw next index to replace
            return;
        }
            
        ++m_vnSkip[nSS];
    }

    m_vbStep2  [nSS]  = true;
    m_vnCount  [nSS]  = nR;
    m_vfQ      [nSS]  = log(1. - m_vfW[nSS]);
    m_vnSkip   [nSS] += (uint64_t)floor(log(m_vfU[nSS]) / m_vfQ[nSS]);
    m_vnT      [nSS]  = 0;
    m_vfW      [nSS] *= beta_distribution(m_vRndGen[nSS], m_nSetSize+0., nR+1.);
    m_vfU      [nSS]  = uniform_real_distribution<double>(0, 1           )(m_vRndGen[nSS]);    
    m_vnNextIdx[nSS]  = uniform_int_distribution <size_t>(0, m_nSetSize-1)(m_vRndGen[nSS]); // draw next index to replace
}
