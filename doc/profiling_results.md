# Profiling: GPU

# SIRT

Reconstructing a 1024x1024 image using SIRT with 20 iterations.

> 6.372s

To compare: CPU takes 9.346s, ASTRA takes 0.44s

Made the mistake of writing to global memory in inner loop of FP. (no effect, gets optimized)

IDEAS:
- First, no longer fill but do it in a kernel, no effect.
- Allocate using pitch. Logical width of 2d array is called 'pitch'
- Trying to precompute the number of steps, did not help, not sure if quality deteriorates..

# SART

Reconstructing a 1024x1024 image using SART with 10 iterations.

> 14.117 sec

First low hanging fruit. Updating the volume, and comparing with zero instead of epsilon

> [11.036, 10.6, 10.58, 10.58] sec

Remove branching from projetor (`if index < allowed`), also 'nudged' the origin towards delta

> [10.2, 10.2]

- Doing the FP in parallel, instead of per block (better parallelization)
- Precomputing the number of steps, putting it in the line

# Profiling: ASTRA vs Galactica

Test: Running 10 iterations of ART on a 128x128 SL-phantom, measuring time in seconds. 64 angles and 64 detectors.

## Galactica

Make sure to compile with `-O3` and without `-pg -g`. Before any optimizations:

> 1.893 sec

When running a profiler, I noticed that the if statement / filling of vector in `linear_projector::operator++` was (clearly) in a critical path, which may lead to poor performance.
Attempted fix: construct the entire queue in constructor.

> 0.874 sec

Made sqrt(2.0) constexpr constant, and merged `beta` and `C` in SIRT.

> 0.809 sec

Use volume centered around the origin instead of left bottom corner at the origin. I assume this will save time primarily in the line iterator. (saves one vector addition)

> \pm 0.809 (no effect)

Remove vector allocation in `interpolate` function, instead use a C array

> 0.382 sec

Next we move 'projector' from the inside loop to a reset function to construct

> 0.212 sec

Also cleared image instead of reallocate for SIRT.

Other ideas: multithreading and explicit sse
