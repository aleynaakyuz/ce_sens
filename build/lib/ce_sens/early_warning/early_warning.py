import numpy as np
from pycbc.waveform.spa_tmplt import spa_length_in_time
from pycbc.types.frequencyseries import FrequencySeries
from pycbc.filter import sigma
from scipy.optimize import minimize
from pycbc.psd.read import from_txt
from ce_sens.utils import get_proj_strain

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
    temp_param.update({'f_lower':312})
    mer_t = spa_length_in_time(**temp_param)
    return mer_t

def stitching_psds(psd_1, psd_2, lag, switch_duration, param):
    df = param['delta_f']
    merger_time_after_intersection = merger_time(param)
    end_time = merger_time_after_intersection + lag
    start_time = merger_time_after_intersection + switch_duration + lag 
    sf = calculate_freqs(start_time, param)
    ef = calculate_freqs(end_time, param)
    end_inx = round(ef) * 10
    start_inx = round(sf) * 10
    first_arr = psd_1.data[:start_inx]
    second_arr = psd_2.data[end_inx:]
    margin = abs(end_inx - start_inx)
    gap = np.ones(margin)  #* 10**(-47)
    new_psd = np.concatenate((first_arr, gap, second_arr))
    new_psd_fs = FrequencySeries(new_psd, delta_f=df)
    return new_psd_fs, sf, ef

def early_warning(time, param, det, psd_path_1, dynamic_psd, lag=None, switch_duration=None):
    low_f = param['f_lower']
    df = param['delta_f']
    soln = minimize(find_freq, x0=70, args=(param,time), method='Nelder-Mead')
    x = soln.x
    proj_strain = get_proj_strain(det, param)
    if dynamic_psd:
        print('here')
        psd_1 = from_txt(psd_path_1, low_freq_cutoff=low_f, length=int(4000/df), delta_f=df)
        psd_2 = from_txt(dynamic_psd, low_freq_cutoff=low_f, length=int(4000/df), delta_f=df)
        psd, sf, ef = stitching_psds(psd_1, psd_2, lag, switch_duration, param)
        snr = sigma(proj_strain, psd=psd, low_frequency_cutoff=low_f, high_frequency_cutoff=max(x, 5.5))
        return snr, sf, ef
    else:
        #print('here')
        psd = from_txt(psd_path_1, low_freq_cutoff=low_f, length=int(4000/df), delta_f=df)
        snr = sigma(proj_strain, psd=psd, low_frequency_cutoff=low_f, high_frequency_cutoff=max(x, 5.5))
        return snr
