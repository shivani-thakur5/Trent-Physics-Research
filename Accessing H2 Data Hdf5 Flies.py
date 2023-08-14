# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 10:59:46 2022

@author: shivani thakur 
""" 

import h5py
H2_filename = "hih2_galaxy_067.hdf5"

with h5py.File(H2_filename, "r") as f:
    # Print all root level object names (aka keys) 
    # these can be group or dataset names 
    print("Keys: %s" % f.keys())
    # get first object name/key; may or may NOT be a group
    a_group_key = list(f.keys())[0]

    # get the object type for a_group_key: usually group or dataset
    print(type(f[a_group_key])) 

    # If a_group_key is a group name, 
    # this gets the object names in the group and returns as a list
    data = list(f[a_group_key])

    # If a_group_key is a dataset name, 
    # this gets the dataset values and returns as a list
    data = list(f[a_group_key])
    # preferred methods to get dataset values:
    ds_obj = f[a_group_key]      # returns as a h5py dataset object
    ds_arr = f[a_group_key][()]  # returns as a numpy array
    
    id_group = f['id_group'][()]
    id_subhalo = f['id_subhalo'][()]
    m_h2_GD14_map = f['m_h2_S14_map'][()]
    m_h2_GD14_map = f['m_h2_S14_map'][()]
    gas_mass_profile = f['profile_gas_rho_map'][()]
    stellar_mass_profile = f['profile_star_rho_map'][()]
    sfr_profile = f['profile_sfr_map'][()]
    

#Loading to a txt file - open() function with "w" for creating a new file or writing over an existing one
fname5 = "D:\Shivani\Summer 2022 Internship\TNG Illustrus API\TNG 100-1\Output\orbit txt files\Orbit_"+ str(orbit) +"_" + subtxt +".txt"
with open(fname5, "w") as g:
    np.savetxt(g, np.c_[np.asarray(snaps), np.asarray(subhaloIDs),np.asarray(MassFracss), np.asarray(lst_vx_weightedAvg), np.asarray(lst_vy_weightedAvg), np.asarray(lst_vz_weightedAvg), np.asarray(lst_totvel_weightedavg)], fmt = '%.5f',  delimiter=',')
