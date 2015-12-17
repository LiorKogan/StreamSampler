# StreamSampler
A [stream](https://en.wikipedia.org/wiki/Stream_(computing)) is a sequence of data elements made available over time. The number of elements in the stream is usually large and unknown a priori. 

A stream sampler extracts a sample set with a given number of elements from a stream. Each possible sample set (of a given size) has an equal probability of being extracted. A stream sampler is an [online algorithm](https://en.wikipedia.org/wiki/Online_algorithm): The size of the input is unknown, and only [one pass](https://en.wikipedia.org/wiki/One-pass_algorithm) over the stream is possible. 

The following [sampling without replacement](https://en.wikipedia.org/wiki/Simple_random_sample) [reservoir algorithms](https://en.wikipedia.org/wiki/Reservoir_sampling) are implemented:

 - R    : Based on "The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
 - X,Y,Z: Based on ["Random Sampling with a Reservoir"](http://www.cs.umd.edu/~samir/498/vitter.pdf) [Jeferey Scott Vitter, 1985]
 - K,L,M: Based on "Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))" [Kim-Hung Li, 1994]

The implementation supports simultaneous extraction of an arbitrary number of independent sample sets.

