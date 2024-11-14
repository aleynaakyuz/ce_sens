import numpy as np
import argparse
import h5py
from tqdm import tqdm
from pycbc.detector import Detector
from ce_sens.utils import get_dic, calculate_snr
from pycbc.psd.read import from_txt
#from ce_sens.snr_calc.optimal_df import opt_df_static
from ce_sens.post_merger.bns_postmerger import get_temp, create_td_data, match_data, align_normalize, to_freq, post_merger_snr

def bns_pm_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("output_path", help="output path")
    parser.add_argument("det", type=str, help="detector")
    #parser.add_argument("snr_path_1", type=str, help="First detector that has snr calculation")
    #parser.add_argument("snr_path_2", type=str, help="First detector that has snr calculation")
    parser.add_argument("psd_path", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    #parser.add_argument("psd_path_2", help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    #parser.add_argument("psd_path_3", nargs='?', help="Path to psd(s). Gave them as lists. If there is two psds to be sticted gave them in the same list")
    #parser.add_argument("snr_path_3", nargs='?', help="First detector that has snr calculation")
    #parser.add_argument("dynamic_psd", nargs='?', help="First detector that has snr calculation")


    args = parser.parse_args()

    input_path = args.path
    out_path = args.output_path

    #snr_path_1 = args.snr_path_1
    #snr_path_2 = args.snr_path_2
    #snr_path_3 = args.snr_path_3

    psd_path = args.psd_path
    #psd_path_2 = args.psd_path_2
    #psd_path_3 = args.psd_path_3
    #dynamic_psd = args.dynamic_psd

    end = args.end
    start = args.start

    slen = 10
    srate = 10000
    tlen = slen * srate
    flen = tlen // 2 + 1
    df = 1.0 / slen
    dt = 1.0 / srate

    hp, hc = get_temp()

    #snr_path_lst = [snr_path_1, snr_path_2, snr_path_3]
    #snr_paths = [x for x in snr_path_lst if x is not None]

    #psd_path_lst = [psd_path_1, psd_path_2, psd_path_3]
    #psd_paths = [p for p in psd_path_lst if p is not None]

    #snr_vals_dic = {}
    #for path in snr_paths:
    #    snr = h5py.File(path, 'r')
    #    network_key = list(snr.keys())
    #    snr_vals_dic[network_key[0]] = snr[network_key[0]][:]

    #det_list = list(snr_vals_dic.keys())

    #psd_dic = {key: psd for key in det_list for psd in psd_paths}
    #psds = {key:from_txt(psd_dic[key], length=flen, delta_f=df, low_freq_cutoff=flow[key]) for key in psd_dic}

    #net_snr_shape =  len(snr_vals_dic[det_list[0]])
    #net_snr = np.zeros(net_snr_shape)
    #for det_snr in snr_vals_dic.values():
    #    sum_snr = net_snr + det_snr**2
    #det_snr = np.sqrt(sum_snr)

    #detections = det_snr[start:end] > 10

    #lenn = sum(detections)

    data_dic = get_dic(input_path)
    det = Detector(args.det)

    parameters2 = {
    "approximant": "TaylorF2",
    "delta_t": hp.delta_t,
    "f_lower": 100,
    "phase_order": -1
    }

    parameters3 = {
    "approximant": "TaylorF2",
    "f_lower": 100,
    "phase_order": -1,
    "delta_f" : df
    }

    psd = from_txt(psd_path, length=flen, delta_f= df, low_freq_cutoff=5.2)

    pm_snr_lst = []
    hf = h5py.File(out_path, 'w')
    for i in tqdm(range(start, end)):
        #try:
        temp_data = {key: value[i] for key, value in data_dic.items()}
        param = {**temp_data, **parameters3}
        snr = calculate_snr(args.det, psd, param, 5.2)
        #except Exception as e:
        #    print(f"An error occurred: {e}")
        #    snr = 0

        if snr > 10:
            pm_param = {**temp_data, **parameters2}

            hr_l, mlen_l = create_td_data(pm_param, hp)

            s_lst, i_lst, hp_lst, hc_lst, hr_lst, snr_l = match_data(hp, hc, hr_l, mlen_l)
            hp_lst, hc_lst = align_normalize(hp_lst, hc_lst, hr_lst, s_lst, i_lst, snr_l)

            z = data_dic['redshift'][:][i]

            php_l, phc_l = to_freq(hp_lst, hc_lst, z)

            pm_snr = post_merger_snr(args.psd_path, det, temp_data, php_l, phc_l, z)

        else:
            pm_snr = 0

        pm_snr_lst.append(pm_snr)
    hf.create_dataset('pm_snrs', data=pm_snr_lst)
    hf.close()

