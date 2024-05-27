import numpy as np
import h5py
from pycbc.cosmology import redshift
import pycbc.workflow.configuration as wfc
from pycbc.distributions.utils import prior_from_config
from pycbc.population.population_models import distance_from_rate, merger_rate_density, coalescence_rate, sfr_madau_dickinson_2014, total_rate_upto_redshift


def normalization_const(rho, time):
    d = 0.0001
    z = redshift(d)
    V = (4/3)*np.pi*d**3
    loc_events = rho*time*V
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z, npoints=1000)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z, npoints=1000)
    tot_rate = total_rate_upto_redshift(z, coa_rate)
    c = loc_events / tot_rate 
    return c

def number_of_samples(rho, c, z_max):
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z_max, npoints=1000)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z_max, npoints=1000)
    tot_rate = total_rate_upto_redshift(z_max, coa_rate)
    number = tot_rate * c
    return int(number)

def other_params(samples, tot_rate, coa_rate, z_max):
    dist = distance_from_rate(tot_rate, coa_rate, maxz=z_max, npoints=1000)
    z = redshift(dist)
    mass1 = samples['srcmass1'] * (1+z)
    mass2 = samples['srcmass2'] * (1+z)
    return dist, z, mass1, mass2

def make_hdf5(samples, output_path, dist, z, mass1, mass2):
    with h5py.File(output_path, 'w') as hf:
        for i in samples.fieldnames:
            hf.create_dataset(str(i), data=samples[str(i)])
        hf.create_dataset('distance', data=dist)
        hf.create_dataset('redshift', data=z)
        hf.create_dataset('mass1', data=mass1)
        hf.create_dataset('mass2', data=mass2)
    hf.close()

def create_data(inp_path, out_path, rho, time, z_max):
    c = normalization_const(rho, time)
    num = number_of_samples(rho, c, z_max)
    
    cp = wfc.WorkflowConfigParser(inp_path)
    joint_dist = prior_from_config(cp, prior_section='prior')
    merger_rate_dens = merger_rate_density(sfr_madau_dickinson_2014, 'inverse', 
                                           rho_local=rho, maxz=z_max, npoints=1000)
    coa_rate = coalescence_rate(merger_rate_dens, maxz=z_max, npoints=1000)

    
    samples = joint_dist.rvs(num)
    tot_rate = samples['total_rate']

    dist, z, mass1, mass2 = other_params(samples, tot_rate, coa_rate, z_max)

    data = make_hdf5(samples, out_path, dist, z, mass1, mass2)

    return data
