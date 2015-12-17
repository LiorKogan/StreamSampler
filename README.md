# StreamSampler
A stream is a sequence of data elements made available over time. The number of elements in the stream is usually very large and unknown a priori. 

A stream sampler extracts a given number of independent sample sets, each with a given number of elements, from the stream. Only one pass over the stream is possible. Each possible sample set (of a given size) has an equal probability of being extracted.

The following sampling-without-replacement algorithms are implemented:

 - R    : Based on "The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
 - X,Y,Z: Based on "Random Sampling with a Reservoir" [Jeferey Scott Vitter, 1985]
 - K,L,M: Based on "Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))" [Kim-Hung Li, 1994]
