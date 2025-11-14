import numpy as np
import argparse
from pycbc.inference.io import loadfile  

def check_mchirp():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")
    parser.add_argument("--config_override")
    parser.add_argument("--mchirp", type=float)

    opts = parser.parse_args()
    path = opts.results
    override_path = opts.config_override
    mchirp = opts.mchirp

    data = loadfile(path, 'r')
    m = data.read_samples(parameters=['mchirp'])['mchirp']

    frac = 0.01

    min_m = min(m)
    max_m = max(m)

    w = max_m - min_m
    eps = frac * w

    up_mchirp = 150
    low_mchirp = 150

    print(np.mean(m >= max_m - eps))
    print(np.mean(m <= min_m + eps))
    if np.mean(m >= max_m - eps) > 0.001:
        print(path, 'max mchirp bound is not converged')
        triggered = True
        up_mchirp = up_mchirp + 50
    if np.mean(m <= min_m + eps) > 0.001:
        print(path, 'min mchirp bound is not converged')
        triggered = True
        low_mchirp = low_mchirp + 50
    if triggered:
        f = open(override_path, 'w')
        f.write(f"""
[prior-mchirp]
name = uniform
min-mchirp = {max(mchirp - low_mchirp, 0.1)}
max-mchirp = {mchirp + up_mchirp}
    """)
        f.close()
    else:
        f = None
    return f

check_mchirp()
