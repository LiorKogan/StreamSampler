# StreamSampler

Header-only C++11 library

Copyright 2015 Lior Kogan (koganlior1@gmail.com)

Released under the Apache License, Version 2.0

--

A [stream](https://en.wikipedia.org/wiki/Stream_(computing)) is a sequence of data elements made available over time. The number of elements in the stream is usually unknown a priori and can be large or even infinite.

A **stream sampler** maintains one or more up-to-date sample sets, each with a fixed number of elements.

We concentrate here on [simple random samples](https://www.scribbr.com/methodology/simple-random-sampling/): for each sample set, each stream element (from the start of the sampling till the latest element available) has an equal chance of being a member.

Up-to-dateness means that each sample set remains a simple random sample as stream element becomes available.

Stream samplers are implemented using [online algorithms](https://en.wikipedia.org/wiki/Online_algorithm): The size of the stream is unknown, and only [one pass](https://en.wikipedia.org/wiki/One-pass_algorithm) over the stream is possible. The memory size required by each stream sampler is constant and very small.

The following seven unweighted [sampling without replacement](https://en.wikipedia.org/wiki/Simple_random_sample) [reservoir](https://en.wikipedia.org/wiki/Reservoir_sampling) [randomized](https://en.wikipedia.org/wiki/Randomized_algorithm) algorithms are implemented:

 - R    : Presented in ["The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R](https://books.google.co.il/books?id=Zu-HAwAAQBAJ&printsec=frontcover&hl=iw&source=gbs_ge_summary_r&cad=0#v=onepage&q&f=false) (Reservoir Sampling) [attributed to Waterman](https://markkm.com/blog/reservoir-sampling/), modified according to Ex.10
 - X,Y,Z: Presented in ["Random Sampling with a Reservoir"](http://www.cs.umd.edu/~samir/498/vitter.pdf) [Jeferey Scott Vitter, 1985]
 - K,L,M: Presented in ["Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))"](http://dl.acm.org/citation.cfm?id=198435) [Kim-Hung Li, 1994]

Algorithm R is the standard 'textbook algorithm'. Algorithms X, Y, Z, K, L, and M offer huge performance improvement by drawing the number of stream elements to skip at each stage, so much less random numbers are generated, especially for large streams. Z, K, L, and M are typically 100's of times faster than R, while M is usually the most performant.

In all these papers, the algorithms were formulated such that the algorithm controls elements fetching from the stream (An external function, *GetNextElement()*, is called from within the algorithms). Such flow control is generally less suitable for real-world scenarios. In this implementation, the algorithms were reformulated such that a process can fetch elements from the stream, and a member function of the stream sampler class (*AddElement*) should be called. *AddElement* returns the number of future stream elements the caller should skip before calling *AddElement* again (hence the sublinear complexity).

Two versions of *AddElement* are implemented: one using copy semantics (*AddElement(const ElementType& Element)*) and one using move semantics (*AddElement(ElementType&& Element)*).

This implementation also extends the algorithms by supporting simultaneous construction of multiple independent sample sets.

*StreamSamplerTest* contains a usage example: *StreamSamplerExample()*, a comparative performance benchmark function *StreamSamplerPerformanceBenchmark()* and a uniformity test function *StreamSamplerTestUniformity()*.

