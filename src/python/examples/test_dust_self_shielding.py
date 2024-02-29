########################################################################
#
# yt fields for grackle functions
#
#
# Copyright (c) Enzo/Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import os
import yt
from yt.utilities.physical_constants import kboltz, mass_hydrogen_cgs, gravitational_constant_cgs
import unyt
yt.set_log_level("error")

import caesar as cs
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5

from tqdm import tqdm

from pygrackle import \
    add_grackle_fields_simba
    
     
constants = {"dust_optical_depth_LW_small": unyt.unyt_quantity(8.2e-21, ""),
             "dust_optical_depth_LW_large": unyt.unyt_quantity(2.3e-21, ""),
             "Gamma":5/3}
    
#* >>> Functions to check if string represents a float or an int <<<
def check_type(input: str) -> str:
    input = input.strip()
    input.replace(".", "a", 1)
    input.replace("e", "a", 1)
    input.replace("E", "a")
    
    if input.isdigit():
        input_type = "int"
    else:
        input_type = "float"
        
    return input_type
    

#* >>> Function to load grackle parameters into dictionary from Gizmo text output <<<
def load_gizmo_parameters(grackle_info_filepath:str) -> dict:
    """
    INPUTS:
        + grackle_info_filepath -> The path to the GRACKLE_INFO output from Gizmo
    """
    parameters = {}
    known_strings = ["grackle_data_file"]
    known_arrays = {"SolarAbundances":[None for i in range(10)]}
    with open(grackle_info_filepath, "r") as fptr:
        skip = True
        for line in fptr.readlines():
            if not skip:
                par_name = line.split(" = ")[0]
                par_val = line.split(" = ")[1]
                
                # Handle known string values
                if par_name.strip() in known_strings:
                    par_val = line.split(" = ")[1]
                    par_val = str(os.path.join(os.path.dirname("/home/ejones/codes/kiara-stack/grackle/input/"), par_val.strip()))
                # Handle known array values
                elif par_name.split("[")[0].strip() in known_arrays:
                    index = int(par_name.split("[")[1].split("]")[0].strip())
                    known_arrays[par_name.split("[")[0].strip()][index] = float(par_val.strip())
                    continue
                else:
                    # Convert numeric values to floats and ints accordingly
                    if check_type(par_val) == "float":
                        par_val = float(par_val)
                    elif check_type(par_val) == "int":
                        par_val = int(par_val)
                    
                parameters[par_name.strip()] = par_val
            
            if "grackle run-time parameters:" in line.lower():
                skip = False
                
    # Add arrays to the parameters dictionary
    #! Arrays cannot be added to the parameter dictonary at the moment so ignore them
    # for par in known_arrays.keys():
    #     parameters[par] = known_arrays[par]
    
    return parameters

def save_dust_shielding_data(snapdir, snapnamebase, snapnumber):
    fptr = h5.File(snapdir + f"dust_shielding/Tdust_{snapnamebase}_{str(snapnumber).zfill(3)}.hdf5", "a")

    for dust_shielding in ["on", "off"]:
        
        groupname = "dust_shielding_" + dust_shielding
        if groupname not in fptr:
            fptr_group = fptr.create_group(groupname)
        else:
            fptr_group = fptr[groupname]
            
        ds = yt.load(snapdir + f"snap_{snapnamebase}_{str(snapnumber).zfill(3)}.hdf5")
        grp = cs.load(f"{snapdir}caesar_files/caesar_{snapnamebase}_{str(snapnumber).zfill(3)}.hdf5")
        particleIDs = ds.all_data()["PartType0", "ParticleIDs"]
        
        grackle_pars = load_gizmo_parameters(snapdir + "GRACKLE_INFO")
        grackle_pars["dust_self_shielding"] = 1 if dust_shielding == "on" else 0
        if grackle_pars["use_isrf_field"] == 0:
            print("Error! use_isrf_field has been turned off! Exiting...")
            exit(1)
        add_grackle_fields_simba(ds, parameters=grackle_pars)
        
        all_gals = [gal for gal in grp.galaxies]
            
        for gal_index in tqdm(range(len(all_gals))):
            
            gal = all_gals[gal_index]

            if f"galID_{gal.GroupID}" in fptr_group:
                continue
                
            gal_pos = gal.pos.to("code_length")
            gal_radii = gal.radii["gas_r80"].to("kpc").d * 2
            gal_partIDs = particleIDs[gal.glist]
            
            sp = ds.sphere(gal_pos, gal_radii)
            
            dust_temperature = sp["gas", "grackle_dust_temperature"].d
            galaxy_filter = np.array([partID in gal_partIDs for partID in sp["PartType0", "ParticleIDs"]], dtype=bool)
            dust_temperature = dust_temperature[galaxy_filter]
                
            fptr_dataset = fptr_group.create_dataset(f"galID_{gal.GroupID}", data=dust_temperature)
            
def calculate_hydrogen_column_density(gas_temperature, gas_density, mean_molecular_weight, HI_density, HII_density,
                                        HM_density, redshift, parameters):
    
    a = 1 / (1 + redshift)
    cl_jeans = np.sqrt((parameters["Gamma"] * np.pi * kboltz) / \
                       (gravitational_constant_cgs * mass_hydrogen_cgs))
    
    ldustShield = cl_jeans * np.sqrt(gas_temperature / (gas_density * mean_molecular_weight))
    
    N_H =  ldustShield / mass_hydrogen_cgs * (HI_density + HII_density + HM_density).to("g/cm**3")
    
    return N_H
            
def calculate_dust_shielding_factor(dust_to_gas, hydrogen_column_density):
    dust_optical_depth_small = constants["dust_optical_depth_LW_small"] * \
                            (dust_to_gas / 0.01) * hydrogen_column_density
    dust_optical_depth_large = constants["dust_optical_depth_LW_large"] * \
                                (dust_to_gas / 0.01) * hydrogen_column_density
    dust_shielding_factor = np.exp(-1 * (dust_optical_depth_small + dust_optical_depth_large))
    dust_shielding_factor[dust_shielding_factor > 1.] = 1.
    
    return dust_shielding_factor
    
    
snapnamebase = "m50n512"
snapnumber = 97
snapdir = f"/Backup00/ejones/{snapnamebase}/"

#save_dust_shielding_data(snapdir, snapnamebase, snapnumber)

ds = yt.load(snapdir + f"snap_{snapnamebase}_{str(snapnumber).zfill(3)}.hdf5")

grackle_pars = load_gizmo_parameters(snapdir + "GRACKLE_INFO")
grackle_pars["dust_self_shielding"] = 1
add_grackle_fields_simba(ds, parameters=grackle_pars)

sp = ds.sphere(ds.domain_center, (1000, "kpc"))
filter = np.logical_and(sp["gas", "dust_density"] > 0., sp["gas", "isrf_habing"] > 0.)

mean_molecular_weight = sp["gas", "grackle_mean_molecular_weight"][filter]
gas_temperature = sp["gas", "temperature"][filter]
HI_density = sp["gas", "H_p0_density"][filter]
HII_density = sp["gas", "H_p1_density"][filter]
HM_density = sp["gas", "H_m1_density"][filter]
dust_density = sp["gas", "dust_density"][filter]
gas_density = sp["gas", "density"][filter]
dust_to_gas = dust_density / gas_density
isrf = sp["gas", "isrf_habing"][filter]

hydrogen_column_density = calculate_hydrogen_column_density(gas_temperature, gas_density,
                            mean_molecular_weight, HI_density, HII_density, HM_density,
                            ds.current_redshift, grackle_pars)

dust_shielding_factor = calculate_dust_shielding_factor(dust_to_gas, hydrogen_column_density)



isrf_shielded = isrf * dust_shielding_factor

fig, ax = plt.subplots(ncols=2, sharey=True)
ax[0].hist(np.log10(isrf), bins=50, density=True, log=True, histtype="step", label="No shielding", color="red")
ax[0].hist(np.log10(isrf_shielded), bins=100, density=True, log=True, histtype="step", label="Shielding", color="blue")
ax[0].legend()
_, bins, _ = ax[1].hist(np.log10(isrf), bins=50, density=True, log=True, histtype="step", label="No shielding", color="red")
ax[1].hist(np.log10(isrf_shielded), bins=bins, density=True, log=True, histtype="step", label="Shielding", color="blue")
for a in ax: a.grid(color="grey", alpha=0.2)
for a in ax: a.set_xlabel(r"$\log_{10}G_0$ [Habing]")
ax[0].set_ylabel("Probability Density")
# ax.set_xlim(-20, 20)
fig.tight_layout()
fig.savefig("isrf.png")