from . import py_tomos
from tomos.py_tomos import *

import matplotlib
import matplotlib.pyplot as plt

import astra

import numpy as np

from math import sqrt


def plot_helper_(f, ax):
    data = f.data()
    dim = f.shape()
    ax.imshow(np.split(np.array(data),
                       len(data) // dim.x),
              interpolation="none",
              cmap='bone',
              extent=[0, dim.x, 0, dim.y])
    plt.show()


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


def partition_bisection_astra(astra_proj_geom, astra_vol_geom, p, eps):
    def convert_cone(cone_vec):
        src, d, u, v = np.split(cone_vec, 4)
        return projection(
            vec3f(src[0], src[1], src[2]),
            vec3f(d[0], d[1], d[2]),
            vec2f(np.linalg.norm(v) * shape[0], np.linalg.norm(u) * shape[1]),
            [vec3f(v[0], v[1], v[2]), vec3f(u[0], u[1], u[2])], 
            vec2i(shape[0], shape[1]),
            False
            )

    def convert_par(par_vec):
        ray, d, u, v = np.split(par_vec, 4)
        # FIXME this hack should depend on the volume size which is unknown here..
        # I think this is a safe size, but should have an optional warning here
        ray *= 10000.0
        return projection(
            vec3f(d[0] - ray[0], d[1] - ray[1], d[2] - ray[2]),
            vec3f(d[0], d[1], d[2]),
            vec2f(np.linalg.norm(v) * shape[0], np.linalg.norm(u) * shape[1]),
            [vec3f(v[0], v[1], v[2]), vec3f(u[0], u[1], u[2])], 
            vec2i(shape[0], shape[1]),
            True
            )

    # convert astra geometry
    geom_type = astra_proj_geom['type']
    if geom_type == 'cone':
        vec_geom = astra.functions.geom_2vec(astra_proj_geom)
        shape = np.array([vec_geom['DetectorColCount'], vec_geom['DetectorRowCount']])
        projection_list = list(map(convert_cone, vec_geom['Vectors']))
        custom_geom = custom_geometry(projection_list)
    elif geom_type =='parallel3d':
        vec_geom = astra.functions.geom_2vec(astra_proj_geom)
        shape = np.array([vec_geom['DetectorColCount'], vec_geom['DetectorRowCount']])
        projection_list = list(map(convert_par, vec_geom['Vectors']))
        custom_geom = custom_geometry(projection_list)
    elif geom_type =='cone_vec':
        vec_geom = astra_proj_geom
        shape = np.array([vec_geom['DetectorColCount'], vec_geom['DetectorRowCount']])
        projection_list = list(map(convert_cone, vec_geom['Vectors']))
        custom_geom = custom_geometry(projection_list)
    elif geom_type =='parallel3d_vec':
        vec_geom = astra_proj_geom
        shape = np.array([vec_geom['DetectorColCount'], vec_geom['DetectorRowCount']])
        projection_list = list(map(convert_par, vec_geom['Vectors']))
        custom_geom = custom_geometry(projection_list)
    else:
        print('ERROR running GRCB Partitioner: unknown ASTRA geometry type.')
        return

    # convert volume
    opts = astra_vol_geom['option']
    vol_shape = [astra_vol_geom['GridRowCount'], astra_vol_geom['GridColCount'], astra_vol_geom['GridSliceCount']]
    vol_origin = [opts['WindowMinX'], opts['WindowMinY'], opts['WindowMinZ']]
    vol_size = [opts['WindowMaxX'] - opts['WindowMinX'],
                opts['WindowMaxY'] - opts['WindowMinY'],
                opts['WindowMaxZ'] - opts['WindowMinZ']]
    tomos_vol = volume_3d(
            vec3i(vol_shape[0], vol_shape[1], vol_shape[2]),
            vec3f(vol_origin[0], vol_origin[1], vol_origin[2]),
            vec3f(vol_size[0], vol_size[1], vol_size[2]))

    # partition
    return partition_bisection(custom_geom, tomos_vol, p, eps)


def print_partitioning(tree):
    node = tree.root()
    q = [(node, 0)]
    while q:
        c, d = q.pop()

        s = "{} ({}, {})".format('>' * (d + 1), c.value.d, c.value.a)
        print(s)

        if c.left():
            q.append((c.left(), d + 1))
        if c.right():
            q.append((c.right(), d + 1))
