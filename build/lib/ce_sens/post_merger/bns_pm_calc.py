import numpy as np
import argparse
import h5py
from pycbc.psd.read import from_txt
from ce_sens.utils import get_dic
from ce_sens.post_merger.bns_postmerger import get_temp, create_td_data, match_data, align_normalize, to_freq, pm_snr

def bns_pm_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("output_path", help="output path")
    parser.add_argument("snr_path_1", type=str, help="First detector that has snr calculation")
    parser.add_argument("snr_path_2", type=str, help="First detector that has snr calculation")
    parser.add_argument("psd_path_1", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("psd_path_2", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("psd_path_3", nargs='?', help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("dynamic_psd", nargs='?', help="First detector that has snr calculation")
    parser.add_argument("snr_path_3", nargs='?', help="First detector that has snr calculation")

    args = parser.parse_args()

    input_path = args.path
    out_path = args.output_path
    
    snr_path_1 = args.snr_path_1
    snr_path_2 = args.snr_path_2
    snr_path_3 = args.snr_path_3

    psd_path_1 = args.psd_path_1
    psd_path_2 = args.psd_path_2
    psd_path_3 = args.psd_path_3
    dynamic_psd = args.dynamic_psd

    end = args.end
    start = args.start

    slen = 10
    srate = 10000
    tlen = slen * srate
    flen = tlen // 2 + 1
    df = 1.0 / slen
    dt = 1.0 / srate

    hp, hc = get_temp()
    flow = {'J1': 5.2, 'J2':5.2, 'E1':1, 'I1': 3}

    snr_path_lst = [snr_path_1, snr_path_2, snr_path_3]
    snr_paths = [x for x in snr_path_lst if x is not None]

    psd_path_lst = [psd_path_1, psd_path_2, psd_path_3]
    psd_paths = [p for p in psd_path_lst if p is not None]

    snr_vals_dic = {}
    for path in snr_paths:
        snr = h5py.File(path, 'r')
        network_key = list(snr.keys())
        snr_vals_dic[network_key[0]] = snr[network_key[0]][:]

    det_list = list(snr_vals_dic.keys())

    psd_dic = {key: psd for key in det_list for psd in psd_paths}
    psds = {key:from_txt(psd_dic[key], length=flen, delta_f=df, low_freq_cutoff=flow[key]) for key in psd_dic}
    
    net_snr_shape =  len(snr_vals_dic[det_list[0]])
    net_snr = np.zeros(net_snr_shape)
    for det_snr in snr_vals_dic.values():
        sum_snr = net_snr + det_snr**2
    det_snr = np.sqrt(sum_snr)

    detections = det_snr[start:end] > 10
    
    lenn = sum(detections)

    data_dic = get_dic(input_path)

    parameters2 = {
    "approximant": "TaylorF2",
    "delta_t": hp.delta_t,
    "f_lower": 100,
    "phase_order": -1
    }

    temp_data = {key: value[start:end][detections] for key, value in data_dic.items()}

    hr_l, mlen_l = create_td_data(temp_data, parameters2, hp, lenn)

    s_lst, i_lst, hp_lst, hc_lst, hr_lst, snr_l = match_data(hp, hc, hr_l, mlen_l, lenn)
    hp_lst, hc_lst = align_normalize(hp_lst, hc_lst, hr_lst, s_lst, i_lst, snr_l, lenn)

    z = data_dic['redshift'][:][start:end]

    php_l, phc_l = to_freq(hp_lst, hc_lst, lenn, z)

    snrs = pm_snr(psds, temp_data, php_l, phc_l, lenn)

    hf = h5py.File(out_path, 'w')
    for k in snrs:
        snrs[k] = np.array(snrs[k])
        hf.create_dataset(str(k), data=snrs[k])
    hf.close()

