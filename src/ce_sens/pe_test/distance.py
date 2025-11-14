import numpy as np
import argparse
from pycbc.inference.io import loadfile  

def check_distance():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")
    parser.add_argument("--config_override")
    parser.add_argument("--distance", type=float)

    opts = parser.parse_args()
    path = opts.results
    override_path = opts.config_override
    dist = opts.distance

    data = loadfile(path, 'r')
    d = data.read_samples(parameters=['distance'])['distance']

    frac = 0.01

    min_d = min(d)
    max_d = max(d)

    w = max_d - min_d
    eps = frac * w

    low_dist = 3
    up_dist = 3

    if np.mean(d >= max_d - eps) > 0.001:
        print(path, 'max distance bound is not converged')
        triggered = True
        up_dist = 4
    if np.mean(d <= min_d + eps) > 0.001:
        print(path, 'min distance bound is not converged')
        triggered = True
        low_dist = 4
    if triggered:
        f = open(override_path, 'w')
        f.write(f"""
[prior-distance]
name = uniform_radius
min-distance = {dist / low_dist}
max-distance = {dist * up_dist}
    """)
        f.close()
    else:
        f = None
    return f

check_distance()
