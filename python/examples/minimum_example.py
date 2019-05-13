import tpt

n = 256

v = tpt.volume(n)
f = tpt.modified_sl_phantom(v)
tpt.plot(f)

k = tpt.closest(v);
g = tpt.parallel(v, n)
p = tpt.forward_project(f, g, k)
tpt.plot(p)

r = tpt.sirt(v, g, k, p, beta = 1.0, iterations = 10)
tpt.plot(r)

