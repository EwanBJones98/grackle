########################################################################
#
# Cooling cell example script
#
#  This will initialize a single cell at a given temperature,
#  iterate the cooling solver for a fixed time, and output the
#  temperature vs. time.
#
#
# Copyright (c) 2015-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import yt

from pygrackle import \
    FluidContainer, \
    chemistry_data, \
    evolve_constant_density

from pygrackle.one_zone import \
    ConstantPressureModel

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from read_debug import read_debug

tiny_number = 1e-20

# Constants for unit conversions
KPC_IN_CM = 3.086e21
SOLARMASS_IN_GRAMS = 1.989e33
BOLTZMANN = 1.38064852e-23

if __name__ == "__main__":
    

    # ---- #
    # Ordering in file: ID, a, time, density, energy, metallicity, fH2, G0, note
    # Set initial values
    a                   = 0.0862079
    density             = 4.52669e-20 # g / cm^3
    initial_energy      = 1.34674e+13 # cm^2 / s^2 ``
    metallicity         = 0.00148776
    heating_rate        = -0.160125
    #initial_temperature = 1.e4 # K
    final_time          = 100.0 # Myr
    current_redshift    = 1/a - 1
    # ---- #

    # Read external fields from debug file
    externalFieldNames = ("specific_heating_rate", "time")
    externalFieldDict = read_debug("ISM_DEBUG.txt", externalFieldNames)
    
    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 1
    my_chemistry.primordial_chemistry = 2
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.self_shielding_method = 3
    my_chemistry.H2_self_shielding = 3
    my_chemistry.h2_on_dust = 1
    my_chemistry.interstellar_radiation_field = 1.0
    my_chemistry.use_specific_heating_rate = 1
    my_dir = os.path.dirname(os.path.abspath(__file__))
    grackle_data_file = os.path.join(
        my_dir, "..", "..", "..", "input", "CloudyData_UVB=HM2012.h5")
    my_chemistry.grackle_data_file = grackle_data_file

    #Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.a_value = 1. / (1. + current_redshift) / \
        my_chemistry.a_units
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Myr          # 1 Myr in s
    my_chemistry.velocity_units = my_chemistry.a_units * \
        (my_chemistry.length_units / my_chemistry.a_value) / \
        my_chemistry.time_units

    # my_chemistry.comoving_coordinates = 0 # proper units
    # my_chemistry.a_units = 1.0
    # # my_chemistry.a_value = 1. / (1. + current_redshift) / \
    # #     my_chemistry.a_units
    # my_chemistry.a_value = initial_a
    # my_chemistry.density_units = 1.0
    # my_chemistry.length_units = 1.0
    # my_chemistry.time_units = sec_per_Myr          # 1 Myr in s
    # my_chemistry.velocity_units = 1.0

    rval = my_chemistry.initialize()

    fc = FluidContainer(my_chemistry, 1)
    fc["density"][:] = density / my_chemistry.density_units
    fc["energy"][:] = initial_energy / (my_chemistry.velocity_units**2)

    if my_chemistry.primordial_chemistry > 0:
        fc["HI"][:] = 0.76 * fc["density"]
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HeI"][:] = (1.0 - 0.76) * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 1:
        fc["H2I"][:] = tiny_number * fc["density"]
        fc["H2II"][:] = tiny_number * fc["density"]
        fc["HM"][:] = tiny_number * fc["density"]
        fc["de"][:] = tiny_number * fc["density"]
    if my_chemistry.primordial_chemistry > 2:
        fc["DI"][:] = 2.0 * 3.4e-5 * fc["density"]
        fc["DII"][:] = tiny_number * fc["density"]
        fc["HDI"][:] = tiny_number * fc["density"]
    if my_chemistry.metal_cooling == 1:
        fc["metal"][:] = metallicity * fc["density"] * \
          my_chemistry.SolarMetalFractionByMass

    fc["x-velocity"][:] = 0.0
    fc["y-velocity"][:] = 0.0
    fc["z-velocity"][:] = 0.0

    #fc["specific_heating_rate"][:] = heating_rate

    # fc["energy"][:] = initial_temperature / \
    #     fc.chemistry_data.temperature_units
    # fc.calculate_temperature()
    # fc["energy"][:] *= initial_temperature / fc["temperature"]


    model = ConstantPressureModel(
        fc, safety_factor=0.01,
        final_time=final_time, external_data=externalFieldDict)
    model.evolve()
    data = model.finalize_data()

    p1, = pyplot.loglog(data["time"].to("Myr"), data["temperature"],
                        color="black", label="T")
    p2, = pyplot.loglog(data["time"].to("Myr"), data["dust_temperature"],
                        color="blue", label="T$_{dust}$")
    pyplot.xlabel("Time [Myr]")
    pyplot.ylabel("T [K]")

    pyplot.twinx()
    p3, = pyplot.loglog(data["time"].to("Myr"), data["H2I"]/data["density"],
                        color="red", label="$f_{H2}$")
    pyplot.ylabel("$f_{H2}$")
    pyplot.legend()
    pyplot.legend([p1,p2, p3], ["T$_{gas}$", "T$_{dust}$", "$f_{H2}$"],
                  fancybox=True, loc="center left")
    pyplot.ylim(1e-8, 1)
    pyplot.tight_layout()
    name = model.name
    pyplot.savefig(f"{name}.png")

    # save data arrays as a yt dataset
    yt.save_as_dataset({}, f"{name}.h5", data)
