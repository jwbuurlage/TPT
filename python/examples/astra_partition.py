import tomos
import numpy as np
import astra

angles = np.linspace(0, 2*np.pi, 48, False)
proj_geom = astra.create_proj_geom('cone', 1.0, 1.0, 32, 64, angles, 1000, 0)
vol_geom = astra.create_vol_geom(64, 64, 64)

tomos.partition_bisection_astra(proj_geom, vol_geom, 4, 0.05)

par_proj_geom = astra.create_proj_geom('parallel3d', 1.0, 1.0, 32, 64, angles)
tomos.partition_bisection_astra(par_proj_geom, vol_geom, 4, 0.05)

cone_vec_geom = astra.functions.geom_2vec(proj_geom)
tomos.partition_bisection_astra(cone_vec_geom, vol_geom, 4, 0.05)

par_vec_geom = astra.functions.geom_2vec(par_proj_geom)
tomos.partition_bisection_astra(par_vec_geom, vol_geom, 4, 0.05)
