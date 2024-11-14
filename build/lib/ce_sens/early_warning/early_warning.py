import numpy as np
from pycbc.waveform.spa_tmplt import spa_length_in_time
from pycbc.types.frequencyseries import FrequencySeries
from scipy.optimize import minimize
from ce_sens.utils import calculate_snr

df_min = 0.0001
f_max = 4000

def find_freq(f_lower, data, target):
    temp_data = data.copy()
    temp_data.update({'f_lower':f_lower})
    guess = spa_length_in_time(**temp_data)
    return abs(target - guess)

def calculate_freqs(time, param):
    soln_1 = minimize(find_freq, x0=70, args=(param,time), method='Nelder-Mead')
    freq = soln_1.x[0]
    return freq

def merger_time(param):
    temp_param = param.copy()
    temp_param.update({'f_lower':297.5})
    mer_t = spa_length_in_time(**temp_param)
    return mer_t

def early_warning(time, param, det, psd_1, dynamic_psd, lag=None, switch_duration=None):
    soln = minimize(find_freq, x0=70, args=(param,time), method='Nelder-Mead')
    x = soln.x
    snr, sf, ef = opt_df_dynamic(param, det, psd_1, dynamic_psd, lag, switch_duration, high_freq=max(x, 5.5))
    return snr, sf, ef

def opt_df_dynamic(param, det, psd, dynamic_psd, lag, switch_duration, high_freq=None, save_psd=False):
    df = param['delta_f']
    low_freq = param['f_lower']
    new_psd, sf, ef = stitching_psds(psd[df], dynamic_psd[df], lag, switch_duration, param)
    snr_l = calculate_snr(det, new_psd, param, low_freq, high_freq)
    while df > df_min:
        df = df/2
        param.update({"delta_f": df})
        new_psd, sf, ef = stitching_psds(psd[df], dynamic_psd[df], lag, switch_duration, param)
        snr_s = calculate_snr(det, new_psd, param, low_freq, high_freq)
        if abs(snr_l - snr_s) / snr_l < 0.01:
            break
        else:
            snr_l = snr_s
            continue
    if save_psd:
        return snr_l, sf, ef, new_psd
    else:
        return snr_l, sf, ef

def stitching_psds(psd_1, psd_2, lag, switch_duration, param):
    df = param['delta_f']
    merger_time_after_intersection = merger_time(param)
    end_time = merger_time_after_intersection + lag
    start_time = merger_time_after_intersection + switch_duration + lag
    sf = calculate_freqs(start_time, param)
    ef = calculate_freqs(end_time, param)
    end_inx = round(ef * (1/df))
    start_inx = round(sf * (1/df))
    first_arr = psd_1.data[:start_inx]
    second_arr = psd_2.data[end_inx:]
    margin = abs(end_inx - start_inx)
    gap = np.ones(margin) #* 10**(-47)
    new_psd = np.concatenate((first_arr, gap, second_arr))
    new_psd_fs = FrequencySeries(new_psd, delta_f=df)
    return new_psd_fs, sf, ef
