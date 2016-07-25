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

## ASTRA
