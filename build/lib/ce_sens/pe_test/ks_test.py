from pycbc.inference.burn_in import ks_test
import numpy as np
import argparse
from pycbc.inference.io import loadfile  
import sys

def ks_test():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data")
    parser.add_argument("--data-retry")

    opts = parser.parse_args()
    data_path = opts.data
    data_retry_path = opts.data_retry

    data_retry = loadfile(data_retry_path, 'r')
    data = loadfile(data_path, 'r')

    d = data.read_samples(parameters=['distance'])['distance']
    d_retry = data_retry.read_samples(parameters=['distance'])['distance']

    d_dic = {'distance':d}
    d_retry_dic = {'distance':d_retry}

    test = ks_test(d_dic, d_retry_dic)

    if test:
        sys.exit(1)
    else:
        sys.exit(0)
