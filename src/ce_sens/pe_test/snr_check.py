import argparse
import csv
from pycbc.inference.io import loadfile

def snr_check():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results")
    parser.add_argument("--csv-path")
    parser.add_argument("--injection-snr", type=float)



    opts = parser.parse_args()
    path = opts.results
    csv_path = opts.csv_path
    injection_snr = opts.injection_snr    

    data = loadfile(path, 'r')

    estimated_snr = (2 * max(data.read_samples(parameters='loglr')['loglr']))**0.5

    difference_snr = abs(estimated_snr - injection_snr) / injection_snr * 100
    difference_flag = difference_snr > 5.0

    cfg_str = data['config_file']['0'][:].tobytes().decode('latin-1')

    for line in cfg_str.splitlines():
        if line.strip().startswith("injection-file"):
            inj_file = line
    

    for line in cfg_str.splitlines():
        if line.strip().startswith("nlive"):
            nlive = line

    row = [
    path,
    inj_file,
    nlive,
	injection_snr,
	estimated_snr,
	difference_snr,
	difference_flag
    ]


    with open(csv_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(row)    
    
