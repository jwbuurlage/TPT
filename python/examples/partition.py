import tomos

n = 64

big_vol = tomos.volume_3d(tomos.vec3i(n, n, n))
cone = tomos.cone_beam(big_vol, 100, tomos.vec2f(1, 1), tomos.vec2i(n, n), 1.0, 1.0)

tomos.partition_bisection(cone, big_vol, 4, 0.05)
