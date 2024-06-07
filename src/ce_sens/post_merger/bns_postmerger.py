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

def get_temp():
    f = h5py.File('/home/aakyuz/runs3/pm/data.h5', 'r')
    d = f['rh_22/Rh_l2_m2_r01000.txt'][:]
    t = d[:,0] * lal.MTSUN_SI * 2.7
    dx = t[1] - t[0]
    hp = TimeSeries(d[:, 1], delta_t=dx, epoch=t[0])
    hc = TimeSeries(d[:, 2], delta_t=dx, epoch=t[0])
    hp = hp.time_slice(-1, hp.end_time)
    hc = hc.time_slice(-1, hc.end_time)
    return hp, hc

def create_td_data(temp_data, param2, hp, lenn):
    hr_l = []
    mlen_l = []
    for i in trange(lenn):
        temp_param = {key: temp_data[key][i] for key in temp_data.keys()}
        param = {**temp_param, **param2}
        hr, _ = get_td_waveform(**param)
        mlen = max(len(hp), len(hr))
        hr_l.append(hr)
        mlen_l.append(mlen)
    return hr_l, mlen_l

def match_data(hp, hc, hr_l, mlen_l, lenn):
    s_lst = []
    i_lst = []
    hp_lst = []
    hc_lst = []
    hr_lst = []
    snr_l = []
    for inx in trange(lenn):
        hp_cp = hp.copy()
        hc_cp = hc.copy()
        
        hr_l[inx].resize(mlen_l[inx])
        hp_cp.resize(mlen_l[inx])
        hc_cp.resize(mlen_l[inx])
        
        hr_lst.append(hr_l[inx])
        hc_lst.append(hc_cp)
        hp_lst.append(hp_cp)
        
        snr = matched_filter(hp_cp, hr_l[inx], low_frequency_cutoff=300, high_frequency_cutoff=700)
        s, i = snr.abs_max_loc()
        s_lst.append(s)
        i_lst.append(i)
        snr_l.append(snr)
    return s_lst, i_lst, hp_lst, hc_lst, hr_lst, snr_l

def align_normalize(hp_lst, hc_lst, hr_lst, s_lst, i_lst, snr_l, lenn):
    for inx in trange(lenn):
        n = sigma(hp_lst[inx], low_frequency_cutoff=300, high_frequency_cutoff=700)
        sdif = hp_lst[inx].start_time - hr_lst[inx].start_time
        hp_lst[inx].start_time += i_lst[inx] * snr_l[inx].delta_t  - sdif
        hc_lst[inx].start_time += i_lst[inx] * snr_l[inx].delta_t  - sdif
        hp_lst[inx] *= s_lst[inx] / n
        hc_lst[inx] *= s_lst[inx] / n
    return hp_lst, hc_lst

def to_freq(hp_lst, hc_lst, lenn):
    php_l = []
    phc_l = []
    for inx in trange(lenn):
        hp_cp = hp_lst[inx].copy()
        hc_cp = hc_lst[inx].copy()

        php = hp_cp.time_slice(-.0001, 0.01)
        phc = hc_cp.time_slice(-.0001, 0.01)

        php = php.to_frequencyseries(delta_f=df)
        phc = phc.to_frequencyseries(delta_f=df)
        php_l.append(php)
        phc_l.append(phc)
    return php_l, phc_l

def pm_snr(psds, temp_data, php_l, phc_l, lenn):
    time = 1697205750
    det = {k: Detector(k) for k in psds.keys()}
    snrs = {k: [] for k in psds.keys()}
    for i in trange(lenn):   
        ant = {}
        for ifo in psds.keys():
            fp, fc = det[ifo].antenna_pattern(temp_data['ra'][i], temp_data['dec'][i],
                                    temp_data['pol'][i], time)
            ant[ifo] = fp, fc
            
        for ifo in psds.keys():
            fp, fc = ant[ifo]
            pm = fp * php_l[i] + fc * phc_l[i]
            snr = sigma(pm, psd=psds[ifo], low_frequency_cutoff=1000,
                        high_frequency_cutoff=4800)
            snrs[ifo].append(snr)
    return snrs