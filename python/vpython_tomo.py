import tomo
from vpython import *

import ipywidgets as wd

stop = False

def clear():
    for obj in scene.objects:
        obj.visible=False
        del obj

def show_geometry(geom, volume, steps, detector_count):
    global stop

    vol_center = vector(volume.x() / 2, volume.y() / 2, volume.z() / 2)
    scene.center = vol_center
    box(pos=vol_center, length=volume.x(), height=volume.z(), width=volume.y(), opacity=0.5)
    ball = sphere(pos=vector(geom.get_line(0).origin.x, geom.get_line(0).origin.y, geom.get_line(0).origin.z), radius=0.3, color=color.red)
    ls = []
    while not stop:
        trail = curve(color=color.magenta, radius=0.1)
        for j in range(0, steps):
            rate(steps // 4)
            for ar in ls:
                    ar.visible = False
                    del ar
            ls.clear()
            o = geom.source_location(j);
            for i in range(detector_count * j, detector_count * (j + 1)):
                l = geom.get_line(i)

                delta =  vector(l.length * l.delta.x,
                                   l.length * l.delta.y,
                                   l.length * l.delta.z)

                origin = vector(l.origin.x, l.origin.y, l.origin.z)

                ls.append(arrow(pos = origin, axis = delta, shaftwidth=0.1, color = color.yellow))

            ball.pos = vector(o.x, o.y, o.z);
            trail.append(ball.pos, retain=(steps // 2))
        trail.visible = False
        del trail

    clear()
    stop = False

def show_parititoning():
    # also show geometry at the same time, would be perfect
    pass

def handle(s):
    global stop
    stop = True

b = wd.Button(description='Stop')
b.on_click(handle)
display(b)
