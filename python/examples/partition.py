import tpt

n = 64

big_vol = tpt.volume_3d(tpt.vec3i(n, n, n))
cone = tpt.cone_beam(big_vol, 100, tpt.vec2f(1, 1), tpt.vec2i(n, n), 1.0, 1.0)

tree = tpt.partition_bisection(cone, big_vol, 8, 0.05)
tpt.print_partitioning(tree)
