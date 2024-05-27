from ce_sens.early_warning.early_warning import early_warning
from ce_sens.utils import get_parameter_list
from tqdm import tqdm
import numpy as np
import argparse
import h5py

def early_warning_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("parameter_path", type=str, help="File that has parameters")
    parser.add_argument("approximant", help="approximant to used in the simulation")
    parser.add_argument("output_path", help="output path")
    parser.add_argument("psd_path_1", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("psd_path_2", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("time", nargs=2, type=int, help="Early warning time in seconds. If there are more than one time requested, give as a list")
    parser.add_argument("snr_path_1", type=str, help="First detector that has snr calculation")
    parser.add_argument("snr_path_2", nargs='?', help="First detector that has snr calculation")
    parser.add_argument("dynamic_psd", nargs='?', help="First detector that has snr calculation")
    parser.add_argument("lag", type=float, nargs='?', help="If the psd is dynamic, switch time before the intersection")
    parser.add_argument("switch_duration", type=float, nargs='?', help="If psd is dynamic, duration of the switch")
    parser.add_argument("snr_path_3", nargs='?', help="First detector that has snr calculation")
    parser.add_argument("psd_path_3", nargs='?', help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")

    args = parser.parse_args()
    param_path = args.parameter_path
    apx = args.approximant
    out_path = args.output_path
    time = args.time

    snr_path_1 = args.snr_path_1
    snr_path_2 = args.snr_path_2
    snr_path_3 = args.snr_path_3

    lag = args.lag
    switch_duration = args.switch_duration

    psd_path_1 = args.psd_path_1
    psd_path_2 = args.psd_path_2
    psd_path_3 = args.psd_path_3
    dynamic_psd = args.dynamic_psd

    low_freq_dic = {'J1': 5.2, 'J2':5.2, 'E1':1, 'I1': 3}
    df = 0.1 ## REVISIT THIS

    psd_path_lst = [psd_path_1, psd_path_2, psd_path_3]
    psd_paths = [p for p in psd_path_lst if p is not None]

    snr_path_lst = [snr_path_1, snr_path_2, snr_path_3]
    snr_paths = [x for x in snr_path_lst if x is not None]

    snr_vals_dic = {}
    for path in snr_paths:
        snr = h5py.File(path, 'r')
        network_key = list(snr.keys())
        snr_vals_dic[network_key[0]] = snr[network_key[0]][:]

    det_list = list(snr_vals_dic.keys())

    net_snr_shape =  len(snr_vals_dic[det_list[0]])
    net_snr = np.zeros(net_snr_shape)
    for det_snr in snr_vals_dic.values():
        sum_snr = net_snr + det_snr**2
    det_snr = np.sqrt(sum_snr)

    psd_dic = {key: psd for key in det_list for psd in psd_paths}

    detections = det_snr > 10

    hf = h5py.File(out_path, 'w')
    for det in det_list:
        snr_l = []
        param_list = get_parameter_list(param_path, apx, df, low_freq_dic[det])
        detected_params = np.array(param_list)[detections]
        for t in time:
            for param in tqdm(detected_params):
                ew_snrs = early_warning(t, param, det, psd_dic[det], dynamic_psd, lag, switch_duration)
            hf.create_dataset(str(t), data=ew_snrs)
        hf.create_dataset(str(det), data=snr_l)
    hf.close()
