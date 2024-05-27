from pycbc.detector import Detector
from pycbc.waveform import get_fd_waveform
import numpy as np
import h5py

def check_length(hp, hc, df, f_max=4000):
    if len(hp) * df > f_max:
        inx = int(f_max//df)
        hp = hp[:inx]
        hc = hc[:inx]
    return hp, hc

def get_proj_strain(det, param, freq=False):
    ra = param['ra']
    dec = param['dec']
    pol = param['polarization']
    df = param['delta_f']
    
    t_gps=1697205750
    detec = Detector(det)
    hp, hc = get_fd_waveform(**param)
    hp, hc = check_length(hp, hc, df)
    if freq:
        freq = hp.sample_frequencies.data
        cross_l = []
        plus_l = []
        for f in freq:
            f_plus_samp, f_cross_samp = detec.antenna_pattern(ra, dec, pol, t_gps, frequency=f)
            cross_l.append(f_cross_samp)
            plus_l.append(f_plus_samp)
        f_cross = np.array(cross_l)
        f_plus = np.array(plus_l)
    else:
        f_plus, f_cross = detec.antenna_pattern(ra, dec, pol, t_gps)
    
    proj_strain = f_plus * hp + f_cross * hc
    return proj_strain


def get_parameter_list(data_path, apx, df, low_freq):
    data = h5py.File(data_path, 'r')
    lenn = len(data['mass1'][:])
    temp_params = list(data.keys())
    parameters = {key: key for key in temp_params}

    parameters2 = {
    "approximant": apx,
    "delta_f": df,
    "f_lower": low_freq,
    "phase_order": -1
    }

    data_dic = {key: data[dataset_key][:] for key, dataset_key in parameters.items()}
    temp_data = [{key: value[i] for key, value in data_dic.items()} for i in range(lenn)]
    final_data = [{**temp_data[i], **parameters2} for i in range(lenn)]
    return final_data
