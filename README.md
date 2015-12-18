# StreamSampler
A [stream](https://en.wikipedia.org/wiki/Stream_(computing)) is a sequence of data elements made available over time. The number of elements in the stream is usually large and unknown a priori. 

A stream sampler extracts a sample set with a given size from a stream. Each possible sample set (of the given size) has an equal probability of being extracted. A stream sampler is an [online algorithm](https://en.wikipedia.org/wiki/Online_algorithm): The size of the input is unknown, and only [one pass](https://en.wikipedia.org/wiki/One-pass_algorithm) over the stream is possible. 

The following 7 unweighted [sampling without replacement](https://en.wikipedia.org/wiki/Simple_random_sample) [reservoir algorithms](https://en.wikipedia.org/wiki/Reservoir_sampling) are implemented:

 - R    : Presented in ["The Art of Computer Programming" [Knuth] Vol.2, 3.4.2 Algorithm R](https://books.google.co.il/books?id=Zu-HAwAAQBAJ&printsec=frontcover&hl=iw&source=gbs_ge_summary_r&cad=0#v=onepage&q&f=false) (Reservoir Sampling) attributed to Waterman, modified according to Ex.10
 - X,Y,Z: Presented in ["Random Sampling with a Reservoir"](http://www.cs.umd.edu/~samir/498/vitter.pdf) [Jeferey Scott Vitter, 1985]
 - K,L,M: Presented in ["Reservoir-Sampling Algorithms of Time Complexity O(n(1+log(N)-log(n)))"](http://dl.acm.org/citation.cfm?id=198435) [Kim-Hung Li, 1994]

In all these papers, the algorithms were drafted such that the fetching of elements from the stream is controlled by the algorithm (A  *GetNextElement()* function is called from within the algorithms). Such flow-control is generally less suitable for real-world scenarios. In this implementation, the algorithms were rearranged such that a process can fetch elements from the stream, and a member-function of the stream sampler class (*AddElement()*) should be called.

This implementation also extend the algorithms by supporting simultaneous extraction of any given number of independent sample sets.
