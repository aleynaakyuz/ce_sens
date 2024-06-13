from pycbc.filter import sigma
from pycbc.detector import Detector
from pycbc.waveform import get_td_waveform
from tqdm import trange
from pycbc.psd.read import from_txt
import h5py
import numpy as np
import argparse
from ce_sens.utils import get_dic

def bbh_pm_calc():
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

    time = 1697205750
    slen = 10
    srate = 10000
    tlen = slen * srate
    flen = tlen // 2 + 1
    df = 1.0 / slen
    dt = 1.0 / srate
    flow = 5.2

    parameters2 = {
    "approximant": "IMRPhenomXPHM",
    "delta_t": dt,
    "f_lower": 100,
    "phase_order": -1
    }

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
    psds = {key:from_txt(psd_dic[key], length=flen, delta_f=df, low_freq_cutoff=flow) for key in psd_dic}

    net_snr_shape =  len(snr_vals_dic[det_list[0]])
    net_snr = np.zeros(net_snr_shape)
    for det_snr in snr_vals_dic.values():
        sum_snr = net_snr + det_snr**2
    det_snr = np.sqrt(sum_snr)

    detections = det_snr[start:end] > 10

    lenn_det = sum(detections)

    data_dic = get_dic(input_path)
    temp_data = {key: value[start:end][detections] for key, value in data_dic.items()}

    ra = data_dic['ra']
    dec = data_dic['dec']
    pol = data_dic['polarization']

    hf = h5py.File(out_path, 'w')
    for det_str in det_list:
        det = Detector(det_str)
        snr_l = []
        for i in trange(lenn_det):
            temp_param = {key: temp_data[key][i] for key in temp_data.keys()}
            param = {**temp_param, **parameters2}
            hp, hc = get_td_waveform(**param)
            _, loc_hp = hp.abs_max_loc()
            _, loc_hc = hc.abs_max_loc()
            hp_pm = hp[loc_hp:]
            hc_pm = hc[loc_hc:]
            lenn = max(len(hp_pm), len(hc_pm))
            hp_pm.resize(lenn)
            hc_pm.resize(lenn)
            hp_pm = hp_pm.to_frequencyseries(delta_f=df)
            hc_pm = hc_pm.to_frequencyseries(delta_f=df)
            fp, fc = det.antenna_pattern(ra[i], dec[i], pol[i], time)
            proj_strain = hp_pm * fp + hc_pm * fc
            snr = sigma(proj_strain, psd=psds[det_str], low_frequency_cutoff=1000,
                            high_frequency_cutoff=4800)
            snr_l.append(snr)
        hf.create_dataset(str(det_str), data=snr_l)
    hf.create_dataset('detections', data=detections)
    hf.close()