from pycbc.filter import sigma, matched_filter
from pycbc.detector import Detector
from pycbc.types import TimeSeries
from pycbc.waveform import get_td_waveform
from tqdm import trange
from pycbc.psd.read import from_txt
import matplotlib.pyplot as plt
import h5py, lal
import numpy as np

slen = 10
srate = 10000
tlen = slen * srate
flen = tlen // 2 + 1
df = 1.0 / slen
dt = 1.0 / srate
flow = {'CE40': 5.2, 'CE40_LF':5.2, 'CE20': 5.2, 'CE20_PM':5.2, 'E1':1, 'I1': 3}

def get_temp():
    f = h5py.File('/home/aakyuz/runs5/pm/data.h5', 'r')
    d = f['rh_22/Rh_l2_m2_r01000.txt'][:]
    t = d[:,0] * lal.MTSUN_SI * 2.7
    dx = t[1] - t[0]
    hp = TimeSeries(d[:, 1], delta_t=dx, epoch=t[0])
    hc = TimeSeries(d[:, 2], delta_t=dx, epoch=t[0])
    hp = hp.time_slice(-1, hp.end_time)
    hc = hc.time_slice(-1, hc.end_time)
    return hp, hc

def create_td_data(param, hp):
    hr, _ = get_td_waveform(**param)
    mlen = max(len(hp), len(hr))
    return hr, mlen

def match_data(hp, hc, hr_l, mlen_l):
    hp_cp = hp.copy()
    hc_cp = hc.copy()

    hr_l.resize(mlen_l)
    hp_cp.resize(mlen_l)
    hc_cp.resize(mlen_l)

    snr = matched_filter(hp_cp, hr_l, low_frequency_cutoff=300, high_frequency_cutoff=700)
    s, i = snr.abs_max_loc()
    return s, i, hp, hc, hr_l, snr

def align_normalize(hp_lst, hc_lst, hr_lst, s_lst, i_lst, snr_l):
    n = sigma(hp_lst, low_frequency_cutoff=300, high_frequency_cutoff=700)
    sdif = hp_lst.start_time - hr_lst.start_time
    hp_lst.start_time += i_lst * snr_l.delta_t  - sdif
    hc_lst.start_time += i_lst * snr_l.delta_t  - sdif
    hp_lst *= s_lst / n
    hc_lst *= s_lst / n
    return hp_lst, hc_lst

def to_freq(hp_lst, hc_lst, z):
    hp_cp = hp_lst.copy()
    hc_cp = hc_lst.copy()

    php = hp_cp.time_slice(-.0001, 0.01)
    phc = hc_cp.time_slice(-.0001, 0.01)

    php = php.to_frequencyseries(delta_f=df * (1 + z))
    phc = phc.to_frequencyseries(delta_f=df * (1 + z))
    return php, phc

def post_merger_snr(psd, det, temp_data, php_l, phc_l, z):
    time = temp_data['tc']

    fp, fc = det.antenna_pattern(temp_data['ra'], temp_data['dec'],
                            temp_data['polarization'], time)

    pm = fp * php_l + fc * phc_l
    psd = from_txt(psd, length=flen, delta_f=(df* (1 + z)), low_freq_cutoff=5.2)
    snr = sigma(pm, psd=psd, low_frequency_cutoff=1000,
                high_frequency_cutoff=4800)
    return snr