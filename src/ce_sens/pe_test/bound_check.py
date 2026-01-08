import numpy as np
import argparse
import csv
from pycbc.inference.io import loadfile  

def bound_check():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")

    parser.add_argument(
        "--prior-bound",
        nargs=3,
        action="append",
        metavar=("PARAM", "LOW", "HIGH"),
        help="Prior bounds as: PARAM LOW HIGH",
    )

    parser.add_argument("--csv-path")    

    opts = parser.parse_args()
    path = opts.results
    csv_path = opts.csv_path

    prior_bounds = {
        param: [float(low), float(high)]
        for param, low, high in (opts.prior_bound or [])
    }

    data = loadfile(path, 'r')
    d = data.read_samples(parameters=['distance'])['distance']
    q = data.read_samples(parameters=['q'])['q']
    mchirp = data.read_samples(parameters=['mchirp'])['mchirp']


    if np.mean(d > 0.99 * prior_bounds['distance'][1]) > 0.01: # Checks if more than 1% of samples are within 1% of the max prior bound
        max_dist = "Max distance bound is not converged"
    else:
        max_dist = "Converged"
    if np.mean(d < 1.01 * prior_bounds['distance'][0]) > 0.01: # Checks if more than 1% of samples are within 1% of the min prior bound
        min_dist = "Min distance bound is not converged"
    else:
        min_dist = "Converged"

    if np.mean(q > 0.99 * prior_bounds['q'][1]) > 0.01: # Checks if more than 1% of samples are within 1% of the max prior bound
        max_q = "Max q bound is not converged"
    else:
        max_q = "Converged"
    if np.mean(q < 1.01 * prior_bounds['q'][0]) > 0.01: # Checks if more than 1% of samples are within 1% of the min prior bound
        min_q = "Min q bound is not converged"
    else: 
        min_q = "Converged"

    if np.mean(mchirp > 0.99 * prior_bounds['mchirp'][1]) > 0.01: # Checks if more than 1% of samples are within 1% of the max prior bound
        max_mchirp = "Max mchirp bound is not converged"
    else:
        max_mchirp = "Converged"
    if np.mean(mchirp < 1.01 * prior_bounds['mchirp'][0]) > 0.01: # Checks if more than 1% of samples are within 1% of the min prior bound
        min_mchirp = "Min mchirp bound is not converged"
    else: 
        min_mchirp = "Converged"


    row = [
    path,
    max_dist,
    min_dist,
    max_q,
    min_q,
    max_mchirp,
    min_mchirp
    ]

    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)