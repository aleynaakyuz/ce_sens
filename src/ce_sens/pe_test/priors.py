import numpy as np
import argparse
from pycbc.inference.io import loadfile  

def posterior_quantile_bounds(samples, q_low=0.005, q_high=0.995):
    low = np.quantile(samples, q_low)
    high = np.quantile(samples, q_high)
    return low, high

def additive_margin(low, high, inj, pct=0.20):
    add_inj = False
    if inj < low:
        low = inj
        add_inj = True
    if inj > high:
        high = inj
        add_inj = True
    width = high - low
    if not add_inj:
        low_m = low - pct * width
        high_m = high + pct * width
    else:
        pct = pct * 2
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

    inj_distance = data['injections'].attrs['distance']
    inj_q = data['injections']['q'][:][0]
    inj_mchirp = data['injections']['mchirp'][:][0]

    low_d, high_d = posterior_quantile_bounds(d)
    low_d, high_d = additive_margin(low_d, high_d, inj_distance, pct=0.35)

    low_q, high_q = posterior_quantile_bounds(q)
    low_q, high_q = additive_margin(low_q, high_q, inj_q, pct=0.35)

    low_mchirp, high_mchirp = posterior_quantile_bounds(mchirp)
    low_mchirp, high_mchirp = additive_margin(low_mchirp, high_mchirp, inj_mchirp, pct=0.35)

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
min-mchirp = {max(10, low_mchirp)}
max-mchirp = {high_mchirp}
""")
    f.close()