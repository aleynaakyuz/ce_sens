from ce_sens.early_warning.early_warning import stitching_psds
from ce_sens.snr_calc.snr import calculate_snr

df_min = 0.0001
f_max = 4000

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