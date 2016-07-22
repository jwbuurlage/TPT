from py_galactica import *

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

from math import sqrt


def plot_helper_(f, ax):
    data = f.data()
    dim = f.dimensions()
    ax.imshow(np.split(np.array(data),
                       len(data) // dim[0]),
              interpolation="none",
              cmap='bone',
              extent=[0, dim[0], 0, dim[1]],
              aspect='auto')


def plot(f):
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    plot_helper_(f, ax)


# returns the l1 norm of vector like objects
def l1(reconstruction, phantom):
    result = 0
    r = reconstruction.data()
    p = phantom.data()
    assert(len(r) == len(p))
    for i in range(0, len(r)):
        result += abs(r[i] - p[i])
    return result


# returns the l2 norm of vector-like objects
def l2(reconstruction, phantom):
    result = 0
    r = reconstruction.data()
    p = phantom.data()
    assert(len(r) == len(p))
    for i in range(0, len(r)):
        result += (r[i] - p[i]) * (r[i] - p[i])
    return sqrt(result)


def test_algorithm(alg, volume, geometry, sino, betas, iterations, phantom):
    # make a betas x iterations figure for reconstructions
    fig, axes = plt.subplots(len(betas), len(
        iterations), figsize=(4 * len(iterations), 4 * len(betas)))

    errors = []
    # for every combination, run the reconstruction algorithm (use initial
    # guess to minimize work)
    for i in range(0, len(betas)):
        errors.append(([], []))
        f = image(volume)
        k = 0
        j = 0
        while len(iterations) > 0 and k <= iterations[-1]:
            if k == iterations[j]:
                plot_helper_(f, axes[i, j])
                axes[i, j].set_title(
                    'beta = {:.2f}, iters = {}'.format(betas[i], k))
                j += 1
            f = alg(volume, geometry, sino, betas[i], f)
            k += 1
            errors[i][0].append(k)
            errors[i][1].append(l2(f, phantom))

    fig.tight_layout()

    # make a figure for norms? plot all together in a single (large) graph
    fig2 = plt.figure(figsize=(5, 5))
    ax = fig2.add_axes([0.0, 0.0, 1.0, 1.0])
    for i in range(0, len(betas)):
        ax.plot(errors[i][0], errors[i][1],
                label="beta = {:.2f}".format(betas[i]))
    ax.set_title("errors")
    ax.set_yscale("log")
    ax.autoscale()
    ax.legend()


def timer(fn, *args):
    import time

    class Timer:

        def __enter__(self):
            self.start = time.clock()
            return self

        def __exit__(self, *args):
            self.end = time.clock()
            self.interval = self.end - self.start

    with Timer() as t:
        r = fn(*args)

    print('Function took: %.03f sec.' % t.interval)

    return r
