#!/usr/bin/env

import h5py
import numpy as np
import pycbc.workflow as wf
import pycbc
import argparse
import csv
from pycbc.workflow.configuration import WorkflowConfigParser
import pycbc.workflow.pegasus_workflow as wdax
from pycbc.conversions import mchirp_from_mass1_mass2, q_from_mass1_mass2
from pycbc.psd.read import from_txt
from ce_sens.utils import calculate_snr, get_parameter_list
from pycbc.waveform import get_waveform_filter_length_in_time

def get_net_snr(params):
	df = 0.001
	net_snr_sq = 0
	for k, det in enumerate(opts.detectors):
		psd = from_txt(opts.psd_paths[k], low_freq_cutoff=5.1, length=int(4000/df), delta_f=df)
		snr = calculate_snr(det, psd, params, 5.1)
		net_snr_sq += snr**2
	return np.sqrt(net_snr_sq)

def create_constant_config(params, path):
	net_snr = get_net_snr(params)
	time = get_waveform_filter_length_in_time(**params)

	tc = params['tc']

	f = open(path, 'w')
	f.write(f"""
[model]
marginalize_distance_snr_range = {max(net_snr - 20, 5)}, {net_snr + 20}

[data]
trigger-time = {tc}
injection-file = {inj}
analysis-start-time = {-(int(time) + 5)}
analysis-end-time = 5

[prior-tc]
name = uniform
min-tc = {tc - 0.05}
max-tc = {tc + 0.05}
""")
	f.close()
	fhand = wf.resolve_url_to_file(path)
	return fhand

def create_initial_priors(params, path):
	dist = params['distance']
	mass1 = params['mass1']
	mass2 = params['mass2']
	mchirp = mchirp_from_mass1_mass2(mass1, mass2)
	q = q_from_mass1_mass2(mass1, mass2)

	params.update({'mchirp': mchirp, 'q': q})

	f = open(path, 'w')
	f.write(f"""
[prior-mchirp]
name = uniform
min-mchirp = {max(mchirp - 400, 10)}
max-mchirp = {mchirp + 400}

[prior-q]
name = uniform
min-q = {max(q - 6, 1)}
max-q = {(q + 6)}

[prior-distance]
name = uniform_radius
min-distance = {dist/7}
max-distance = {dist*8}

""")
	f.close()
	fhand = wf.resolve_url_to_file(path)
	return fhand

def create_nlive_file(nlive, path):
	f = open(path, 'w')
	f.write(f"""
[sampler]
nlive = {nlive}
""")
	f.close()
	fhand = wf.resolve_url_to_file(path)
	return fhand

parser = argparse.ArgumentParser()
parser.add_argument("--injections", nargs='*')
parser.add_argument("--inference-config")
parser.add_argument("--detectors", nargs='*')
parser.add_argument("--psd-paths", nargs='*')

pycbc.init_logging(True)
wf.add_workflow_command_line_group(parser)
wf.add_workflow_settings_cli(parser)
opts = parser.parse_args()

wf.makedir(opts.output_dir)

bound_csv_file = opts.output_dir + "/bound_checks.csv"

header_bound = [
	"path",
    "injection_file",
	"nlive",
    "max_distance_bound_converged",
    "min_distance_bound_converged",
    "max_q_bound_converged",
    "min_q_bound_converged",
	"max_mchirp_bound_converged",
	"min_mchirp_bound_converged"
]

with open(bound_csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header_bound)

bound_csv_file = wf.resolve_url_to_file(bound_csv_file)

snr_csv_file = "snr_checks.csv"

header_snr = [
	"path",
    "injection_file",
	"nlive",
	"injection_snr",
	"estimated_snr",
	"difference_snr",
	"snr_difference_more_than_5_percent"
]

with open(snr_csv_file, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(header_snr)	

snr_csv_file = wf.resolve_url_to_file(snr_csv_file)

workflow = wf.Workflow(opts, "gw")
infhand = wf.resolve_url_to_file(opts.inference_config)


inference_exe_first = wf.Executable(workflow.cp, "inference_first", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

inference_exe_second = wf.Executable(workflow.cp, "inference_second", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

demarg_exe = wf.Executable(workflow.cp, "demarg", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

prior_exe = wf.Executable(workflow.cp, "priors", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

bound_check_exe = wf.Executable(workflow.cp, "bound_check", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

snr_check_exe = wf.Executable(workflow.cp, "snr_check", ifos=workflow.ifos,
                              out_dir=opts.output_dir)

df = 0.001
for inj_idx, inj in enumerate(opts.injections):
	params2 = h5py.File(inj)
	params = {key:params2[key][:][0] for key in list(params2.keys())}
	att = {key:params2.attrs[key] for key in list(params2.attrs.keys())}
	params.update(att)
	params.update({'delta_f':df})

	initial_prior_path = opts.output_dir + f'/initial_priors_{inj_idx}.ini'
	nlive_path = opts.output_dir + f'/initial_nlive_{inj_idx}.ini'
	const_path = opts.output_dir + f'/constant_config_{inj_idx}.ini'

	constant_file = create_constant_config(params, const_path)
	fhand = create_initial_priors(params, initial_prior_path)
	nlive_file = create_nlive_file(100, nlive_path)

	node = inference_exe_first.create_node()
	node.add_input_list_opt("--config-file", [infhand, fhand, nlive_file, constant_file])
	inference_file = node.new_output_file_opt(workflow.analysis_time, '.hdf', 
										   '--output-file', 
										   tags=[str(inj)])
	workflow += node

	node = demarg_exe.create_node()
	node.add_input_opt("--input-file", inference_file)
	demarg_file_initial = node.new_output_file_opt(workflow.analysis_time, ".hdf",
                                              "--output-file",
                                              tags=[str(inj) + 'demarg'])
	workflow += node

	node = bound_check_exe.create_node()
	node.add_input_opt("--results", demarg_file_initial)
	node.add_input_opt("--csv-path", bound_csv_file)
	node.add_list_opt("--distance-bound", [f"{params['distance']/7}", f"{params['distance']*7}"])
	node.add_list_opt("--q-bound", [f"{params['q'] - 6}", f"{params['q'] + 6}"])
	node.add_list_opt("--mchirp-bound", [f"{max(params['mchirp'] - 400, 10)}", f"{params['mchirp'] + 400}"])

	workflow += node

	node = snr_check_exe.create_node()
	node.add_input_opt("--results", demarg_file_initial)
	node.add_input_opt("--csv-path", snr_csv_file)
	node.add_opt("--injection-snr", f"{get_net_snr(params)}")

	workflow += node

	node = prior_exe.create_node()
	node.add_input_opt("--results", demarg_file_initial)
	prior_file = node.new_output_file_opt(workflow.analysis_time, ".ini",
                                              "--priors",
                                              tags=[f'priors_{inj_idx}'])
	workflow += node

	nlive2 = 500

	tag = f'inj{inj_idx}_nlive{nlive2}'
	nlive_path = opts.output_dir + f'/nlive_{nlive2}_inj_{inj_idx}.ini'

	nlive_file = create_nlive_file(nlive2, nlive_path)

	node = inference_exe_second.create_node()
	node.add_input_list_opt("--config-file", [infhand, prior_file, nlive_file, constant_file])
	inference_file = node.new_output_file_opt(workflow.analysis_time, '.hdf', 
										'--output-file', 
										tags=[tag])
	workflow += node

	node = demarg_exe.create_node()
	node.add_input_opt("--input-file", inference_file)
	demarg_file = node.new_output_file_opt(workflow.analysis_time, ".hdf",
											  "--output-file",
											  tags=[f'{tag}_demarg'])
	workflow += node

	node = bound_check_exe.create_node()
	node.add_input_opt("--results", demarg_file)
	node.add_input_opt("--csv-path", bound_csv_file)
	node.add_input_opt("--priors", prior_file)

	workflow += node

	node = snr_check_exe.create_node()
	node.add_input_opt("--results", demarg_file)
	node.add_input_opt("--csv-path", snr_csv_file)
	node.add_opt("--injection-snr", f"{get_net_snr(params)}")

	workflow += node

	nlive3 = 1000

	tag = f'inj{inj_idx}_nlive{nlive3}'
	nlive_path = opts.output_dir + f'/nlive_{nlive3}_inj_{inj_idx}.ini'

	nlive_file = create_nlive_file(nlive3, nlive_path)

	node = inference_exe_second.create_node()
	node.add_input_list_opt("--config-file", [infhand, prior_file, nlive_file, constant_file])
	inference_file = node.new_output_file_opt(workflow.analysis_time, '.hdf', 
										'--output-file', 
										tags=[tag])
	workflow += node

	node = demarg_exe.create_node()
	node.add_input_opt("--input-file", inference_file)
	demarg_file = node.new_output_file_opt(workflow.analysis_time, ".hdf",
											  "--output-file",
											  tags=[f'{tag}_demarg'])
	workflow += node

	node = snr_check_exe.create_node()
	node.add_input_opt("--results", demarg_file)
	node.add_input_opt("--csv-path", snr_csv_file)
	node.add_opt("--injection-snr", f"{get_net_snr(params)}")

	workflow += node

	nlive4 = 3000

	tag = f'inj{inj_idx}_nlive{nlive4}'
	nlive_path = opts.output_dir + f'/nlive_{nlive4}_inj_{inj_idx}.ini'
		
	nlive_file = create_nlive_file(nlive4, nlive_path)

	node = inference_exe_second.create_node()
	node.add_input_list_opt("--config-file", [infhand, prior_file, nlive_file, constant_file])
	inference_file = node.new_output_file_opt(workflow.analysis_time, '.hdf', 
										'--output-file', 
										tags=[tag])
	workflow += node

	node = demarg_exe.create_node()
	node.add_input_opt("--input-file", inference_file)
	demarg_file = node.new_output_file_opt(workflow.analysis_time, ".hdf",
											  "--output-file",
											  tags=[f'{tag}_demarg'])
	workflow += node

	node = snr_check_exe.create_node()
	node.add_input_opt("--results", demarg_file)
	node.add_input_opt("--csv-path", snr_csv_file)
	node.add_opt("--injection-snr", f"{get_net_snr(params)}")

	workflow += node	
workflow.save()				
