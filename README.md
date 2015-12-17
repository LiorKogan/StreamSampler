# StreamSampler
A stream is a sequence of data elements made available over time. 
The number of elements in the stream is usually very large and unknown a priori. 
A stream sampler extracts nSampleSets independent sample sets, each with nSetSize elements, from the stream.
Only one pass over the stream is possible.
Each possible selection of nSetSize elements items has an equal probability of occurring.
