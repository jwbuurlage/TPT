from py_galactica import *

import matplotlib
import matplotlib.pyplot as plt

import numpy as np

def plot(image):
    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    data = image.data()
    dim = image.dimensions()
    ax.imshow(np.split(np.array(data), len(
        data) // dim[0]), interpolation="none", cmap='Greys', extent=[0, dim[0], 0, dim[1]], aspect='auto')
