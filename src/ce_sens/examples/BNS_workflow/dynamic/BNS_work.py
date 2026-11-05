#!/usr/bin/env

import h5py
import re
import pycbc.workflow as wf
import pycbc
import argparse, os
import configparser
import numpy as np
from pycbc.workflow.configuration import WorkflowConfigParser
import pycbc.workflow.pegasus_workflow as wdax
from ce_sens.utils import calculate_snr, get_parameter_list
from pycbc.psd.read import from_txt
from ce_sens.early_warning.early_warning import stitching_psds
from scipy.signal.windows import tukey
from pycbc.types.frequencyseries import FrequencySeries

parser = argparse.ArgumentParser()
parser.add_argument("--injections", nargs='*')
parser.add_argument("--inference-config")
parser.add_argument("--detectors", nargs='*')
parser.add_argument("--psd-paths", nargs='*')

pycbc.init_logging(True)
wf.add_workflow_command_line_group(parser)
wf.add_workflow_settings_cli(parser)
opts = parser.parse_args()

workflow = wf.Workflow(opts, "gw")
wf.makedir(opts.output_dir)

def get_net_snr(params):
	df = 0.001
	net_snr_sq = 0
	for k, det in enumerate(opts.detectors):
		psd = from_txt(opts.psd_paths[k], low_freq_cutoff=5.1, length=int(4000/df), delta_f=df)
		snr = calculate_snr(det, psd, params, 5.1)
		net_snr_sq += snr**2
	return np.sqrt(net_snr_sq)

infhand = wf.resolve_url_to_file(opts.inference_config)

inference_exe = wf.Executable(workflow.cp, "inference", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

psd_lst = [from_txt(path, length=int(5000/0.01), low_freq_cutoff=5.1, delta_f=0.01, is_asd_file=False) for path in opts.psd_paths]

lag = 0
switch_duration = 0.0001

df = 0.001
for i, inj in enumerate(opts.injections):
	params2 = h5py.File(inj)
	params = {key:params2.attrs[key] for key in list(params2.attrs.keys())}
	params.update({'delta_f':df})
	net_snr = get_net_snr(params)
	match = re.search(r'CE40-INJECTION_(.+)\-0-10.hdf', inj)
	i = match.group(1)

	params.update({'phase_order' : -1})
	params.update({'delta_f':0.01})
	new_psd, start, end = stitching_psds(psd_lst[0], psd_lst[-1], lag, switch_duration, params)
	new_psd.save(opts.output_dir + f"/psd{i}.txt")

	path = opts.output_dir + '/inj_%s.ini' % i
	asd_path = f"{opts.detectors[0]}:" + opts.output_dir + f"/psd{i}.txt" + " " + f"{opts.detectors[1]}:" + f"{opts.psd_paths[1]}"
	print('ASD PATH:', asd_path)

	tc = params['tc']
	ra = params['ra']
	dec = params['dec']
	dist = params['distance']
	mchirp = params['mchirp']
	eta = params['eta']
	lambda_tilde = params['lambda_tilde']
	delta_lambda_tilde = params['delta_lambda_tilde']
	lambda1 = params['lambda1']
	lambda2 = params['lambda2']
	mass1 = params['mass1']
	mass2 = params['mass2']

	f = open(path, 'w')
	f.write("""
[model]
tc_ref = %s

mass1_ref = %s
mass2_ref = %s

marginalize_distance_snr_range = %s, %s

lambda1_ref = %s
lambda2_ref = %s

[data]
injection-file = %s
trigger-time = %s
asd-file = %s
fake-strain-from-file = %s

[static_params]
ra = %s
dec = %s

[prior-tc]
name = uniform
min-tc = %s
max-tc = %s

[prior-mchirp]
name = uniform
min-mchirp = %s
max-mchirp = %s

[prior-eta]
name = uniform
min-eta = %s
max-eta = %s

[prior-distance]
name = uniform_radius
min-distance = %s
max-distance = %s

[prior-delta_lambda_tilde]
name = uniform
min-delta_lambda_tilde = %s
max-delta_lambda_tilde = %s

[prior-lambda_tilde]
name = uniform
min-lambda_tilde = %s
max-lambda_tilde = %s
""" % (tc, mass1, mass2, min(5, net_snr - 50), net_snr + 50, lambda1, lambda2, inj, tc, asd_path, asd_path, ra, dec, tc - 0.05, tc + 0.05, mchirp - 0.001, mchirp + 0.001, eta - 0.001, min(0.25, eta + 0.001), max(10, dist - 150), dist + 150, delta_lambda_tilde - 25, delta_lambda_tilde + 25, lambda_tilde - 100, lambda_tilde + 100))
	f.close()
	fhand = wf.resolve_url_to_file(path)
	node = inference_exe.create_node()

	node.add_input_list_opt("--config-file", [infhand, fhand])

	inference_file = node.new_output_file_opt(workflow.analysis_time, '.hdf', '--output-file', tags=[str(inj)])

	workflow += node

workflow.save()
