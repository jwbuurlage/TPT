import tomo
from vpython import *

import ipywidgets as wd
from ipywidgets import interact, interactive, fixed

from threading import *


# Globals
stop = False
steps = 40
k = 4
detector_width = 4
v = tomo.volume_3d(k, k, k)

# Lock for the visualization
vis_lock = RLock()

geometries = {
    'Cone beam': tomo.cone_beam_geometry(v, steps, k / (detector_width - 1), tomo.vec2i(detector_width, detector_width), 1.0, 1.0)
    #    'Cone beam': tomo.cone_beam_geometry(v, steps, 3.33, tomo.vec2i(4, 4), 1.0, 1.0),
    #    'Cone beam': tomo.cone_beam_geometry(v, steps, 3.33, tomo.vec2i(4, 4), 1.0, 1.0),
    #    'Cone beam': tomo.cone_beam_geometry(v, steps, 3.33, tomo.vec2i(4, 4), 1.0, 1.0),
}

# FIXME
# WANT TO USE ARBITRARY SIZES

#@interact(k_set=16)
#def g(k_set):
#    global k
#    k = k_set
# END FIXME

# FIRST:


def clear():
    for obj in scene.objects:
        obj.visible = False
        del obj


def prepare_scene(volume):
    vol_center = vector(volume.x() / 2, volume.y() / 2, volume.z() / 2)
    scene.center = vol_center
    box(pos=vol_center, length=volume.x(), height=volume.z(), width=volume.y(), opacity=0.5)

    w = wd.Dropdown(
        options=geometries,
        description='Geometry:',
    )
    display(w)

    def handle_start(s):
        show_geometry(w.value, v, steps, detector_width ** 2)

    b_start = wd.Button(description='Start')
    b_start.on_click(handle_start)
    display(b_start)

    def handle_stop(s):
        global stop
        stop = True

    b = wd.Button(description='Stop')
    b.on_click(handle_stop)
    display(b)

prepare_scene(v)

def show_geometry(geom, volume, steps, detector_count):
    global stop
    objects = []
    try:
        vis_lock.acquire()
        ball = sphere(pos=vector(geom.get_line(0).origin.x, geom.get_line(
            0).origin.y, geom.get_line(0).origin.z), radius=0.3, color=color.red)
        objects.append(ball)

        ls = []
        while not stop:
            trail = curve(color=color.magenta, radius=0.1)
            for j in range(0, steps):
                rate(steps // 4)

                for ar in ls:
                    ar.visible = False
                    del ar
                ls.clear()

                if stop:
                    break

                o = geom.source_location(j)
                for i in range(detector_count * j, detector_count * (j + 1)):
                    if stop:
                        break

                    l = geom.get_line(i)

                    delta = vector(l.length * l.delta.x,
                                   l.length * l.delta.y,
                                   l.length * l.delta.z)

                    origin = vector(l.origin.x, l.origin.y, l.origin.z)

                    ls.append(arrow(pos=origin, axis=delta,
                                    shaftwidth=0.1, color=color.yellow))

                ball.pos = vector(o.x, o.y, o.z)
                trail.append(ball.pos, retain=(steps // 2))
            trail.visible = False
            del trail
    finally:
        for obj in objects:
            obj.visible = False
            del obj
        stop = False
        vis_lock.release()


def show_partitioning(partitioning):
    colors = [color.red, color.blue, color.green, color.yellow, color.orange, color.white, color.blue, color.cyan, color.magenta]
    boxes = []
    for x in range(0, v.x()):
        for y in range(0, v.y()):
            for z in range(0, v.z()):
                col = colors[partitioning.owner(v.index(x, y, z))]
                boxes.append(box(pos=vector(0.5 + x, 0.5 + y, 0.5 + z), length=1, height=1, width=1, opacity=0.5, color=col))
