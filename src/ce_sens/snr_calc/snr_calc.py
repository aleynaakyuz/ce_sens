import h5py
import argparse
from tqdm import tqdm
from ce_sens.early_warning.early_warning import merger_time, early_warning
from ce_sens.utils import get_dic
from ce_sens.snr_calc.snr import read_psds, opt_df_static, opt_df_dynamic, calculate_snr

def snr_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("det", type=str, help="detector names")
    parser.add_argument("low_freq_cutoff", type=float, help="low_freq_cutoff")
    parser.add_argument("df_max", type=float, help="maximum delta_f")
    parser.add_argument("approximant", type=str, help="approximant to used in the simulation")
    parser.add_argument("psd_path", type=str, help="paths of the psd")
    parser.add_argument("out_path", type=str, help="Output path")
    parser.add_argument("dynamic_psd_path", type=str, nargs='?', help="paths of the second psd")
    parser.add_argument("lag",  type=float, nargs='?', help="time that switch ends before intersection")
    parser.add_argument("switch_duration", type=float, nargs='?', help="Time required for switch to occur")


    args = parser.parse_args()
    input_path = args.path
    start = args.start 
    end = args.end
    det = args.det
    low_freq = args.low_freq_cutoff
    df_max = args.df_max
    apx = args.approximant
    psd_path = args.psd_path
    dynamic_psd = args.dynamic_psd_path
    out_path = args.out_path
    lag = args.lag
    switch_duration = args.switch_duration

    hf = h5py.File(out_path, 'w')

    data_dic = get_dic(input_path)

    parameters2 = {
    "approximant": apx,
    "delta_f": df_max,
    "f_lower": low_freq,
    "phase_order": -1
    }

    if dynamic_psd:
        snr_list = []
        sf_list = []
        ef_list = []
        psd_dic = read_psds(psd_path, df_max, low_freq)
        dynamic_psd_dic = read_psds(dynamic_psd, df_max, low_freq)
        for i in tqdm(range(start, end)):
            temp_data = {key: value[i] for key, value in data_dic.items()}
            param = {**temp_data, **parameters2}
            mer_t = merger_time(param)
            time = mer_t + switch_duration + lag
            ew_snr, sf, ef = early_warning(time, param, det, psd_dic, dynamic_psd_dic, lag, switch_duration)
            if (mer_t > 0) and (ew_snr > 10):
                try:
                    snr, sf, ef = opt_df_dynamic(param, det, psd_dic, dynamic_psd_dic, lag, switch_duration)
                except:
                    snr = 0
                    sf = 0
                    ef = 0
            else:
                try:
                    snr = opt_df_static(param, det, psd_dic)
                    sf = 0
                    ef = 0
                except:
                    snr = 0
                    sf = 0
                    ef = 0
            snr_list.append(snr)
            sf_list.append(sf)
            ef_list.append(ef)

        hf.create_dataset('start_freq', data=sf_list)
        hf.create_dataset('end_freq', data=ef_list)
    else:
        snr_list = []
        psd = read_psds(psd_path, df_max, low_freq)
        for i in tqdm(range(start, end)):
            try:
                temp_data = {key: value[i] for key, value in data_dic.items()}
                param = {**temp_data, **parameters2}
                snr = opt_df_static(param, det, psd)
            except Exception as e:
                print(f"An error occurred: {e}")
                snr = 0
            snr_list.append(snr)

    hf.create_dataset(str(det), data=snr_list)
    hf.close()
