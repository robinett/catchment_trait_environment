from scipy.stats import qmc
import numpy as np

upper_bounds = [2,2,5,5]
lower_bounds = [0,0,-5,-5]
num_particles = 10
num_params = 4

sampler = qmc.LatinHypercube(d=num_params)
sample = sampler.random(n=num_particles)

positions = qmc.scale(sample,lower_bounds,upper_bounds)
positions = np.array(positions)
print(positions)
