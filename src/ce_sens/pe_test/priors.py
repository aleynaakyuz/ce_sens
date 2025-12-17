import numpy as np
import argparse
from pycbc.inference.io import loadfile  

def posterior_quantile_bounds(samples, q_low=0.005, q_high=0.995):
    low = np.quantile(samples, q_low)
    high = np.quantile(samples, q_high)
    return low, high

def additive_margin(low, high, pct=0.20):
    width = high - low
    low_m = low - pct * width
    high_m = high + pct * width
    return low_m, high_m

def write_priors():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")
    parser.add_argument("--priors")

    opts = parser.parse_args()
    path = opts.results
    priors = opts.priors

    data = loadfile(path, 'r')
    d = data.read_samples(parameters=['distance'])['distance']
    q = data.read_samples(parameters=['q'])['q']
    mchirp = data.read_samples(parameters=['mchirp'])['mchirp']

    low_d, high_d = posterior_quantile_bounds(d)
    low_d, high_d = additive_margin(low_d, high_d, pct=0.10)

    low_q, high_q = posterior_quantile_bounds(q)
    low_q, high_q = additive_margin(low_q, high_q, pct=0.10)

    low_mchirp, high_mchirp = posterior_quantile_bounds(mchirp)
    low_mchirp, high_mchirp = additive_margin(low_mchirp, high_mchirp, pct=0.10)

    f = open(priors, 'w')
    f.write(f"""
[prior-distance]
name = uniform_radius
min-distance = {max(10, low_d)}
max-distance = {high_d}

[prior-q]
name = uniform
min-q = {max(1, low_q)}
max-q = {high_q}

[prior-mchirp]
name = uniform
min-mchirp = {max(0.1, low_mchirp)}
max-mchirp = {high_mchirp}
""")
    f.close()