from pycbc.filter import sigma
from pycbc.detector import Detector
from pycbc.waveform import get_td_waveform
from tqdm import trange
from pycbc.psd.read import from_txt
import h5py
import numpy as np
import argparse
from tqdm import tqdm
from ce_sens.utils import get_dic
from ce_sens.snr_calc.optimal_df import opt_df_static
from ce_sens.snr_calc.snr import read_psds

def bbh_pm_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("det", type=str, help="detector")
    parser.add_argument("psd_path", help="Path to psd. Gave them as lists. If there is two psds to be sticted gave them in the same list")
    parser.add_argument("output_path", help="output path")

    args = parser.parse_args()

    input_path = args.path
    out_path = args.output_path

    psd_path = args.psd_path

    end = args.end
    start = args.start

    time = 1697205750
    slen = 10
    srate = 10000
    tlen = slen * srate
    flen = tlen // 2 + 1
    df = 1.0 / slen
    dt = 1.0 / srate
    flow = 5.1

    parameters2 = {
    "approximant": "IMRPhenomXPHM",
    "delta_t": dt,
    "f_lower": 50,
    "phase_order": -1
    }

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
    #psds = {key:from_txt(psd_dic[key], length=flen, delta_f=df, low_freq_cutoff=flow) for key in psd_dic}

    #net_snr_shape =  len(snr_vals_dic[det_list[0]])
    #net_snr = np.zeros(net_snr_shape)
    #for det_snr in snr_vals_dic.values():
    #    sum_snr = net_snr + det_snr**2
    #det_snr = np.sqrt(sum_snr)

    #detections = det_snr[start:end] > 10

    #lenn_det = sum(detections)
    det = Detector(args.det)

    data_dic = get_dic(input_path)
    #temp_data = {key: value[start:end][detections] for key, value in data_dic.items()}

    ra = data_dic['ra']
    dec = data_dic['dec']
    pol = data_dic['polarization']

    snr_l = []
    hf = h5py.File(out_path, 'w')
    for i in tqdm(range(start, end)):
        try:
            temp_data = {key: value[i] for key, value in data_dic.items()}
            param = {**temp_data, **parameters2}
            snr = opt_df_static(param, det, psd)
        except Exception as e:
            print(f"An error occurred: {e}")
            snr = 0

        if snr > 10:
            hp, hc = get_td_waveform(**param)
            _, loc_hp = hp.abs_max_loc()
            _, loc_hc = hc.abs_max_loc()
            hp_pm = hp[loc_hp:]
            hc_pm = hc[loc_hc:]
            lenn = max(len(hp_pm), len(hc_pm))
            hp_pm.resize(lenn)
            hc_pm.resize(lenn)
            hp_pm = hp_pm.to_frequencyseries(delta_f= df * (1 + param['redshift'][:]))
            hc_pm = hc_pm.to_frequencyseries(delta_f= df * (1 + param['redshift'][:]))
            fp, fc = det.antenna_pattern(ra[i], dec[i], pol[i], time)
            proj_strain = hp_pm * fp + hc_pm * fc
            psd = from_txt(psd_path, length=flen, delta_f= df * (1 + param['redshift'][:]), low_freq_cutoff=flow)
            snr = sigma(proj_strain, psd=psd, low_frequency_cutoff=flow)
            snr_l.append(snr)
    hf.create_dataset('snrs', data=snr_l)
    hf.close()