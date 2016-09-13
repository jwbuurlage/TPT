# Galactica

Galactica is a library for developing generic tomographic reconstruction methods. It also comes with a number of reference implementations of common algorithms and common experimental setups. As an introduction to the library, and to fix terminology, we will briefly discuss tomographic reconstruction -- the target application of this library.

## Tomographic reconstruction

Tomography is a field that deals with image reconstruction based on various physical phenomena which are generally a type of _wave_ that penetrate an object that is to be reconstructed. A typical picture of this is something like this:

                                    _____________
                                    |            |
                                    |   o   o    |
      + ----------------------------|-o -- - - - |- - - - - - - - - - - - - - - |
                                    |     o  o   |
                                    |____________|

    source           ray/line          object(s)                              detector
                                    inside a volume

There can be multiple source, and multiple detectors. Furthermore, there can be different _views_ or _projection angles_ from which the object is being ''illuminated'' because the source and detector move around the object (or, equivalently, the object itself moves and/or rotates).

The complete collection of _measurements_ made by the detector is called a **sinogram**. The (discretized) reconstruction of the object is called the **image**.

In Galactica, this general experimental setup is modelled. For this, a number of constructs are used in the library.

## Components and building blocks of Galactica

### Describing the setup

Because of the focus on genericity, galactica supports arbitrary dimensional problem. Since we are measuring an experimental setup, in which we want to reconstruct a chosen volume as an image, we first have to describe this setup.

In Galactica, the *volume* is modelled as a box of size $X_1 \times X_2 \times \ldots \times X_D$. This box consists of voxels of unit volume, and is constructed along the positive axes. For details, see [`tomo::volume`](tomo::volume).

Next, the *geometry*, the collection of rays that are sent through the object, has to be specified. A number of common options are available:

 * Parallel geometry, see [`tomo::parallel_geometry`](tomo::parallel_geometry).
 * List of lines, see [`tomo::list_geometry`](tomo::list_geometry).

The ''list geometry'' is the most generic one, and allows image reconstruction for arbitrary experiments. However, when using a specific geometry some algorithms may be able to make use of the symmetry and structure to make optimizations.

An important ingredient for iterative reconstruction methods is the *projector* (possibly a better name would be 'ray integrator'). The projector decides how the ray is discretized, there is a trade off between more accurate discretization and computation efficiency. Available projectors are:

 * Closest, see [`tomo::closest_projector`](tomo::closest_projector).
 * Interpolation, see [`tomo::linear_projector`](tomo::linear_projector).
 * Joseph, see [`tomo::joseph_projector`](tomo::joseph_projector).

If there is no experimental data available, but instead Galactica is being as a testing bed on top of which algorithms are developed, then we need to simulate an experiment in the code. This is done by means of a *phantom*, which represents the 'actual' image. Available phantoms are

 * Shepp-logan, see [`tomo::shepp_logan_phantom`](tomo::shepp_logan_phantom).

To obtain 'experimental results', run `tomo::forward_projection` on the phantom to obtain a *sinogram*. Alternatively, the *sinogram* can be constructed from actual experimental data.

### Performing the reconstruction

We assume we have a sinogram which matches the geometry of the setup. We then have a number of reconstruction methods available:

 * ART, see [`tomo::art`](tomo::art).
 * SART, see [`tomo::sart`](tomo::sart).
 * SIRT, see [`tomo::sirt`](tomo::sirt).

### Visualizing the result

There is as of yet little support for this. The python module has a plot function for 2D images, and there is also a 2D ascii-based plotter as a utility.

### Choosing implementations

CUDA / MPI
