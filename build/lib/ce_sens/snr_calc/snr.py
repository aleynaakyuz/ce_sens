import numpy as np
from pycbc.psd.read import from_txt
from ce_sens.utils import calculate_snr

df_min = 0.0001
f_max = 4000

def read_psds(psd_path, df_max, low_freq): 
    dfs = [df_max]
    
    while dfs[-1] > df_min:
        dfs.append(dfs[-1]/2)
    df_arr = np.array(dfs)
    
    psds = {}
    for i in df_arr:
        lenn = int(f_max / i) 
        psd = from_txt(psd_path, low_freq_cutoff=low_freq, 
                        length=lenn, delta_f=i)
        psds.update({i: psd})
    return psds
    
def calc_diff(param, det, psd1, psd2, sf, ef=None):
    snr_f = calculate_snr(det, psd1, param, sf, high_freq=ef)
    snr_i = calculate_snr(det, psd2, param, sf, high_freq=ef)
    return snr_f - snr_i

def opt_df_static(final_data, det, psd):
    df = final_data['delta_f']
    low_freq = final_data['f_lower']
    snr_l = calculate_snr(det, psd[df], final_data, low_freq) 
    while df > df_min:
        df = df/2
        final_data.update({"delta_f": df}) 
        snr_s = calculate_snr(det, psd[df], final_data, low_freq)
        if abs(snr_l - snr_s) / snr_l < 0.01:
            break
        else:
            snr_l = snr_s
            continue
    return snr_l

