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
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    parser.add_argument("--dets", nargs="+", required=True,
                        help="Detector names")
    parser.add_argument("--psds", nargs="+", required=True,
                        help="PSD paths")
    parser.add_argument("out_path", type=str)

    args = parser.parse_args()

    if len(args.dets) != len(args.psds):
        raise ValueError("Number of detectors must match number of PSDs")

    data = get_parameter_list(args.path, 'TaylorF2', 0.1, 5.1)

    psds = [
        from_txt(psd_path,
                 length=int(4000/0.1),
                 low_freq_cutoff=5.1,
                 delta_f=0.1)
        for psd_path in args.psds
    ]

    def find_snr(f_higher, data_i):
        temp = data_i.copy()

        if f_higher <= temp['f_lower']:
            return float('inf')

        temp.update({'f_higher': f_higher})
        hp, hc = get_fd_waveform(**temp)

        network_snr_sq = 0

        for det, psd in zip(args.dets, psds):
            f_plus, f_cross = get_proj(det, data_i)
            proj_strain = hp * f_plus + hc * f_cross

            amp = sigma(proj_strain,
                        psd=psd,
                        low_frequency_cutoff=temp['f_lower'],
                        high_frequency_cutoff=temp['f_higher'])

            network_snr_sq += amp**2

        network_snr = np.sqrt(network_snr_sq)
        return abs(network_snr - 10)

    hf = h5py.File(args.out_path, 'w')

    full_snr_lst = []
    ew_snr_lst = []
    time_lst = []
    f_higher_lst = []

    for i in trange(args.start, args.end):

        hp, hc = get_fd_waveform(**data[i])

        network_snr_sq = 0

        for det, psd in zip(args.dets, psds):
            f_plus, f_cross = get_proj(det, data[i])
            proj_strain = hp * f_plus + hc * f_cross

            snr = sigma(proj_strain,
                        psd=psd,
                        low_frequency_cutoff=5.1)

            network_snr_sq += snr**2

        full_snr = np.sqrt(network_snr_sq)

        if full_snr > 10:

            soln = minimize(find_snr,
                            x0=20,
                            args=(data[i],),
                            method='Nelder-Mead')

            f_cut = soln.x[0]

            network_snr_sq = 0

            for det, psd in zip(args.dets, psds):
                f_plus, f_cross = get_proj(det, data[i])
                proj_strain = hp * f_plus + hc * f_cross

                amp = sigma(proj_strain,
                            psd=psd,
                            low_frequency_cutoff=5.1,
                            high_frequency_cutoff=f_cut)

                network_snr_sq += amp**2

            ew_snr = np.sqrt(network_snr_sq)

            data[i].update({'f_lower': f_cut})
            time = spa_length_in_time(**data[i])

        else:
            ew_snr = 0
            time = 0
            f_cut = 0

        full_snr_lst.append(full_snr)
        ew_snr_lst.append(ew_snr)
        time_lst.append(time)
        f_higher_lst.append(f_cut)

    hf.create_dataset('full_snr', data=full_snr_lst)
    hf.create_dataset('early_warning', data=ew_snr_lst)
    hf.create_dataset('time', data=time_lst)
    hf.create_dataset('f_higher', data=f_higher_lst)

    hf.close()
