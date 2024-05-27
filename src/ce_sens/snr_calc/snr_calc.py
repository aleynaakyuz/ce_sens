import h5py
import argparse
from tqdm import tqdm
from ce_sens.early_warning.early_warning import merger_time
from ce_sens.utils import get_dic
from ce_sens.snr_calc.snr import read_psds, opt_df_static, opt_df_dynamic

def snr_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("det", help="detector names")
    parser.add_argument("low_freq_cutoff", type=float, help="low_freq_cutoff")
    parser.add_argument("df_max", type=float, help="maximum delta_f")
    parser.add_argument("approximant", help="approximant to used in the simulation")
    parser.add_argument("psd_path", help="paths of the psd")
    parser.add_argument("out_path", help="Output path")
    parser.add_argument("dynamic_psd_path", nargs='?', help="paths of the second psd")
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
        psd_dic = read_psds(psd_path, df_max, low_freq)
        dynamic_psd_dic = read_psds(dynamic_psd, df_max, low_freq)
        for i in tqdm(range(start, end)):
            temp_data = {key: value[i] for key, value in data_dic.items()}
            param = {**temp_data, **parameters2}
            mer_t = merger_time(param)
            if mer_t > 0:
                snr = opt_df_dynamic(param, det, psd_dic, dynamic_psd_dic, lag, switch_duration)
            else:
                snr = opt_df_static(param, det, psd_dic)
            snr_list.append(snr)


    else:
        snr_list = []
        psd = read_psds(psd_path, df_max, low_freq)
        for i in tqdm(range(start, end)):
            temp_data = {key: value[i] for key, value in data_dic.items()}
            param = {**temp_data, **parameters2}
            snr = opt_df_static(param, det, psd)
            snr_list.append(snr)

    hf.create_dataset(str(det), data=snr_list)
    hf.close()
