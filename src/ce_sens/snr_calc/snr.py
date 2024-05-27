import numpy as np
from pycbc.psd.read import from_txt
from pycbc.psd.read import from_txt
from pycbc.filter import sigma
from ce_sens.utils import get_proj_strain
from ce_sens.early_warning.early_warning import stitching_psds

df_min = 0.125
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

def calculate_snr(det, psd, param, low_freq):
    proj_strain = get_proj_strain(det, param)
    amp = sigma(proj_strain, psd=psd, low_frequency_cutoff=low_freq)
    return amp  

def opt_df_static(final_data, det, psd):
    df = final_data['delta_f']
    low_freq = final_data['f_lower']
    snr_l = calculate_snr(det, psd[df], final_data, low_freq) 
    while df > df_min:
        df = df/2
        final_data.update({"delta_f": df}) 
        snr_s = calculate_snr(det, psd[df], final_data, low_freq)
        if abs(snr_l - snr_s) < 0.1:
            break
        else:
            snr_l = snr_s
            continue
    return snr_l

def opt_df_dynamic(param, det, psd, dynamic_psd, lag, switch_duration):
    df = param['delta_f']
    low_freq = param['f_lower']
    new_psd = stitching_psds(psd[df], dynamic_psd[df], lag, switch_duration, param)
    snr_l = calculate_snr(det, new_psd, param, low_freq) 
    while df > df_min:
        df = df/2
        param.update({"delta_f": df})
        new_psd = stitching_psds(psd[df], dynamic_psd[df], lag, switch_duration, param)
        snr_s = calculate_snr(det, new_psd, param, low_freq)
        if abs(snr_l - snr_s) < 0.1:
            break
        else:
            snr_l = snr_s
            continue
    return snr_l
