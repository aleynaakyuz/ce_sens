import numpy as np
import argparse
import csv
from pycbc.inference.io import loadfile  
from pycbc.workflow.configuration import WorkflowConfigParser


def load_prior_bounds(config_file):
	config = WorkflowConfigParser()
	config.read(config_file)

	prior_bounds = {}

	for section in config.sections():
		if section.startswith("prior-"):
			param = section.replace("prior-", "")
			min_key = f"min-{param}"
			max_key = f"max-{param}"

			prior_bounds[param] = [
				config.getfloat(section, min_key),
				config.getfloat(section, max_key),
			]
	return prior_bounds

def bound_check():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")

    parser.add_argument("--distance-bound", nargs=2, type=float)
    parser.add_argument("--q-bound", nargs=2, type=float)
    parser.add_argument("--mchirp-bound", nargs=2, type=float)

    parser.add_argument("--csv-path")
    parser.add_argument("--priors")
    

    opts = parser.parse_args()
    path = opts.results
    csv_path = opts.csv_path
    priors = opts.priors

    if priors:
        prior_bounds = load_prior_bounds(priors)

    else:
         prior_bounds = {}
         prior_bounds["distance"] = opts.distance_bound
         prior_bounds["q"] = opts.q_bound
         prior_bounds["mchirp"] = opts.mchirp_bound
         

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
        ## TODO: make sure this is not triggered when q_min = 1.
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

    cfg_str = data['config_file']['0'][:].tobytes().decode('latin-1')

    for line in cfg_str.splitlines():
        if line.strip().startswith("injection-file"):
            inj_file = line
    

    for line in cfg_str.splitlines():
        if line.strip().startswith("nlive"):
            nlive = line

    row = [
    inj_file,
    nlive,
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