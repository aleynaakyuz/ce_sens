import numpy as np
import h5py
from pycbc.cosmology import redshift
import pycbc.workflow.configuration as wfc
from pycbc.distributions.utils import prior_from_config
from pycbc.population.population_models import distance_from_rate, merger_rate_density, coalescence_rate, sfr_madau_dickinson_2014, total_rate_upto_redshift
from tqdm import tqdm

nps = 100000

def normalization_const(rho, time):
    d = 0.0001
    z = redshift(d)
    V = (4/3)*np.pi*d**3
    loc_events = rho*time*V
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z, npoints=nps)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z, npoints=nps)
    tot_rate = total_rate_upto_redshift(z, coa_rate)
    c = loc_events / tot_rate 
    return c

def number_of_samples(rho, c, z_max):
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z_max, npoints=nps)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z_max, npoints=nps)
    tot_rate = total_rate_upto_redshift(z_max, coa_rate)
    number = tot_rate * c
    return int(number)

def distance_calc(samples, coa_rate, z_max):
    params = {}
    z_lst = []
    ld_l = []

    for tot_rate in tqdm(samples['total_rate']):
        dist = distance_from_rate(tot_rate, coa_rate, maxz=z_max, npoints=nps)
        z = redshift(dist)
        ld_l.append(dist)
        z_lst.append(z)
    params.update({'distance':np.array(ld_l)})
    params.update({'redshift':np.array(z_lst)})
    return params

def other_params(samples, params, type):

    z = params['redshift']
    if type=='BBH':
        print('IN BBH')
        srcmass2 = samples['srcmass1'][:] * samples['q'][:]
        mass2 = srcmass2 * (1+z)
        params.update({'srcmass2':srcmass2})
        params.update({'mass2':mass2})    
    else:
        print('IN BNS')
        mass2 = samples['srcmass2'][:] * (1+z)
        params.update({'mass2':mass2})

    mass1 = samples['srcmass1'][:] * (1+z)
    params.update({'mass1':mass1})
    return params

def make_hdf5(samples, output_path, params):
    with h5py.File(output_path, 'w') as hf:
        for i in samples.fieldnames:
            hf.create_dataset(str(i), data=samples[str(i)])
        for key in params.keys():
            hf.create_dataset(key, data=params[key])
    hf.close()

def create_data(inp_path, out_path, rho, time, z_max, type):
    c = normalization_const(rho, time)
    num = number_of_samples(rho, c, z_max)
    print('number of parameters:', num)
    cp = wfc.WorkflowConfigParser(inp_path)
    joint_dist = prior_from_config(cp, prior_section='prior')
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z_max, npoints=nps)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z_max, npoints=nps)

    samples = joint_dist.rvs(num)

    params = distance_calc(samples, coa_rate, z_max)
    params = other_params(samples, params, type)

    data = make_hdf5(samples, out_path, params)

    return data
