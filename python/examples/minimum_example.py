import tomos

n = 256

v = tomos.volume(n)
f = tomos.modified_sl_phantom(v)
tomos.plot(f)

k = tomos.closest(v);
g = tomos.parallel(v, n)
p = tomos.forward_project(f, g, k)
tomos.plot(p)

r = tomos.sirt(v, g, k, p, beta = 1.0, iterations = 10)
tomos.plot(r)

