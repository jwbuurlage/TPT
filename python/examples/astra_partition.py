import tpt
import numpy as np
import astra

angles = np.linspace(0, 2*np.pi, 48, False)
vol_geom = astra.create_vol_geom(64, 64, 64)

proj_geom = astra.create_proj_geom('cone', 1.0, 1.0, 32, 64, angles, 1000, 0)
cone_tree = tpt.partition_bisection_astra(proj_geom, vol_geom, 4, 0.05)
tpt.print_partitioning(cone_tree)

par_proj_geom = astra.create_proj_geom('parallel3d', 1.0, 1.0, 32, 64, angles)
parallel_tree = tpt.partition_bisection_astra(par_proj_geom, vol_geom, 4, 0.05)
tpt.print_partitioning(parallel_tree)

cone_vec_geom = astra.functions.geom_2vec(proj_geom)
cone_vec_tree = tpt.partition_bisection_astra(cone_vec_geom, vol_geom, 4, 0.05)
tpt.print_partitioning(cone_vec_tree)

par_vec_geom = astra.functions.geom_2vec(par_proj_geom)
par_vec_tree = tpt.partition_bisection_astra(par_vec_geom, vol_geom, 4, 0.05)
tpt.print_partitioning(par_vec_tree)
