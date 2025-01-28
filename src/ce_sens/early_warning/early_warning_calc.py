from ce_sens.early_warning.early_warning import early_warning
from ce_sens.utils import get_parameter_list, get_proj
from pycbc.waveform.spa_tmplt import spa_length_in_time
from pycbc.psd.read import from_txt
from pycbc.waveform import get_fd_waveform
from pycbc.filter import sigma
from pycbc.detector import Detector
from scipy.optimize import minimize
from tqdm import trange
import numpy as np
import argparse
import h5py

def early_warning_calc():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="path to parameters")
    parser.add_argument("start", type=int, help="start of the parameter index")
    parser.add_argument("end", type=int, help="end of the parameter index")
    parser.add_argument("det1", type=str, help="detector names")
    parser.add_argument("det2", type=str, help="detector names")
    parser.add_argument("psd_path1", type=str, help="paths of the psd")
    parser.add_argument("psd_path2", type=str, help="paths of the psd")
    parser.add_argument("out_path", type=str, help="Output path")

    args = parser.parse_args()
    input_path = args.path
    start = args.start
    end = args.end
    det1 = args.det1
    det2 = args.det2
    psd_path1 = args.psd_path1
    psd_path2 = args.psd_path2
    out_path = args.out_path

    data = get_parameter_list(input_path, 'TaylorF2', 0.1, 5.1)


    psd1 = from_txt(psd_path1, length=int(4000/0.1), low_freq_cutoff=5.1, delta_f=0.1)
    psd2 = from_txt(psd_path2, length=int(4000/0.1), low_freq_cutoff=5.1, delta_f=0.1)


    def find_snr(f_higher, data, psd1, psd2, det1, det2):  
        temp = data.copy()
        if f_higher <= temp['f_lower']:
            return float('inf')
        temp.update({'f_higher':f_higher})
        
        hp, hc = get_fd_waveform(**temp)
        
        f_plus1, f_cross1 = get_proj(det1, data)
        f_plus2, f_cross2 = get_proj(det2, data)

        proj_strain1 = hp * f_plus1 + hc * f_cross1
        proj_strain2 = hp * f_plus2 + hc * f_cross2
        
        amp1 = sigma(proj_strain1, psd=psd1, low_frequency_cutoff=temp['f_lower'], high_frequency_cutoff=temp['f_higher'])
        amp2 = sigma(proj_strain2, psd=psd2, low_frequency_cutoff=temp['f_lower'], high_frequency_cutoff=temp['f_higher'])
        
        amp = np.sqrt(amp1**2 + amp2**2)
        return abs(amp - 10)


    hf = h5py.File(out_path, 'w')
    full_snr_lst = []
    ew_snr_lst = []
    time_lst = []
    f_higher_lst = []
    for i in trange(start, end):
        hp, hc = get_fd_waveform(**data[i])
        strain = hp + hc

        f_plus1, f_cross1 = get_proj(det1, data[i])
        f_plus2, f_cross2 = get_proj(det2, data[i])

        proj_strain1 = hp * f_plus1 + hc * f_cross1
        proj_strain2 = hp * f_plus2 + hc * f_cross2
        
        snr1 = sigma(proj_strain1, psd=psd1, low_frequency_cutoff=5.1)
        snr2 = sigma(proj_strain2, psd=psd2, low_frequency_cutoff=5.1)
        full_snr = np.sqrt(snr1**2 + snr2**2)
        if full_snr > 10:
            soln = minimize(find_snr, x0=20, args=(data[i],psd1, psd2, det1, det2), method='Nelder-Mead')
            
            amp1 = sigma(proj_strain1, psd=psd1, low_frequency_cutoff=5.1, high_frequency_cutoff=soln.x[0])
            amp2 = sigma(proj_strain2, psd=psd2, low_frequency_cutoff=5.1, high_frequency_cutoff=soln.x[0])

            amp = np.sqrt(amp1**2 + amp2**2)
            
            data[i].update({'f_lower':soln.x[0]})
            time = spa_length_in_time(**data[i])

            full_snr_lst.append(full_snr)
            ew_snr_lst.append(amp)
            time_lst.append(time)
            f_higher_lst.append(soln.x[0])
        else:
            amp = 0
            time = 0
            f_higher = 0
            full_snr_lst.append(full_snr)
            ew_snr_lst.append(amp)
            time_lst.append(time)
            f_higher_lst.append(f_higher)

    hf.create_dataset('full_snr', data=full_snr_lst)
    hf.create_dataset('early_warning', data=ew_snr_lst)
    hf.create_dataset('time', data=time_lst)
    hf.create_dataset('f_higher', data=f_higher_lst)
    hf.close()
