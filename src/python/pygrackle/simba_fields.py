"""
Written by: Ewan Jones
Date of Creation: 15/02/2024

This file is analgous to 'yt_fields.py' except it works
on gizmo datasets rather than enzo datasets.
"""
import numpy as np

from yt.fields.field_detector import \
    FieldDetector
from yt.frontends.gizmo.data_structures import \
    GizmoDataset
from yt.utilities.physical_constants import \
    me, mp
from pygrackle import \
    FluidContainer, \
    chemistry_data

_parameter_map = {}
_parameter_map[GizmoDataset] = {
    "use_grackle": "use_grackle",
    "Gamma": "Gamma",
    "primordial_chemistry": "MultiSpecies",
    "metal_cooling": "MetalCooling",
    "h2_on_dust": "H2FormationOnDust",
    "cmb_temperature_floor": "CMBTemperatureFloor",
    "three_body_rate": "ThreeBodyRate",
    "cie_cooling": "CIECooling",
    "h2_optical_depth_approximation": "H2OpticalDepthApproximation",
    "photoelectric_heating": "PhotoelectricHeating",
    "photoelectric_heating_rate": "PhotoelectricHeatingRate",
    "NumberOfTemperatureBins": "NumberOfTemperatureBins",
    "CaseBRecombination": "CaseBRecombination",
    "TemperatureStart": "TemperatureStart",
    "TemperatureEnd": "TemperatureEnd",
    "NumberOfDustTemperatureBins": "NumberOfDustTemperatureBins",
    "DustTemperatureStart": "DustTemperatureStart",
    "DustTemperatureEnd": "DustTemperatureEnd",
    "HydrogenFractionByMass": "HydrogenFractionByMass",
    "DeuteriumToHydrogenRatio": "DeuteriumToHydrogenRatio",
    "SolarMetalFractionByMass": "SolarMetalFractionByMass",
    "UVbackground_redshift_on": "RadiationRedshiftOn",
    "UVbackground_redshift_off": "RadiationRedshiftOff",
    "UVbackground_redshift_fullon": "RadiationRedshiftFullOn",
    "UVbackground_redshift_drop": "RadiationRedshiftDropOff",
    "use_radiative_transfer": "RadiativeTransfer",
    "radiative_transfer_coupled_rate_solver": "RadiativeTransferCoupledRateSolver",
    "radiative_transfer_hydrogen_only": "RadiativeTransferHydrogenOnly",
    "with_radiative_cooling": "with_radiative_cooling",
    "use_volumetric_heating_rate": "use_volumetric_heating_rate",
    "use_specific_heating_rate": "use_specific_heating_rate",
    "self_shielding_method": "self_shielding_method",
    "dust_self_shielding": "dust_self_shielding",
    "H2_self_shielding": "H2_self_shielding",
    "grackle_data_file": "grackle_data_file",
    "UVbackground": "UVbackground",
    "Compton_xray_heating": "Compton_xray_heating",
    "LWbackground_intensity": "LWbackground_intensity",
    "LWbackground_sawtooth_suppression": "LWbackground_sawtooth_suppression",
    "use_isrf_field": "use_isrf_field"
}

_field_map = {
    'density': (('gas', 'density'), 'code_mass / code_length**3'),
    'HI': (('gas', 'H_p0_density'), 'code_mass / code_length**3'),
    'HII': (('gas', 'H_p1_density'), 'code_mass / code_length**3'),
    'HM': (('gas', 'H_m1_density'), 'code_mass / code_length**3'),
    'HeI': (('gas', 'He_p0_density'), 'code_mass / code_length**3'),
    'HeII': (('gas', 'He_p1_density'), 'code_mass / code_length**3'),
    'HeIII': (('gas', 'He_p2_density'), 'code_mass / code_length**3'),
    'H2I': (('gas', 'H2_p0_density'), 'code_mass / code_length**3'),
    'H2II': (('gas', 'H2_p1_density'), 'code_mass / code_length**3'),
    'DI': (('gas', 'D_p0_density'), 'code_mass / code_length**3'),
    'DII': (('gas', 'D_p1_density'), 'code_mass / code_length**3'),
    'HDI': (('gas', 'HD_p0_density'), 'code_mass / code_length**3'),
    'de': (('gas', 'El_density'), 'code_mass / code_length**3'),
    'metal': (('gas', 'total_metal_density'), 'code_mass / code_length**3'),
    'dust': (('gas', 'dust_density'), 'code_mass / code_length**3'),
    'x-velocity': (('gas', 'velocity_x'), 'code_velocity'),
    'y-velocity': (('gas', 'velocity_y'), 'code_velocity'),
    'z-velocity': (('gas', 'velocity_z'), 'code_velocity'),
    'energy': (('gas', 'specific_thermal_energy'), 'code_velocity**2'),
    'RT_heating_rate': (('gas', 'photo_gamma'), 'erg/s'),
    'isrf_habing': (("gas", "isrf_habing"), "dimensionless")
}

def _get_needed_fields(my_chemistry):
    fields = ['density', 'energy'] + \
      ['%s-velocity' % ax for ax in 'xyz']
    if my_chemistry.primordial_chemistry > 0:
        fields += ['HI', 'HII', 'HeI', 'HeII', 'HeIII', 'de']
    if my_chemistry.primordial_chemistry > 1:
        fields += ['HM', 'H2I', 'H2II']
    if my_chemistry.primordial_chemistry > 2:
        fields += ['DI', 'DII', 'HDI']
    if my_chemistry.metal_cooling == 1:
        fields += ['metal']
    if my_chemistry.use_dust_density_field == 1:
        fields += ['dust']
    if my_chemistry.use_radiative_transfer == 1:
        fields += ['RT_heating_rate']
    if my_chemistry.use_isrf_field == 1:
        fields += ['isrf_habing']
    return fields

def _data_to_fc(data, size=None, fc=None):
    if size is None:
        size = data['gas', 'density'].size
    if fc is None:
        fc = FluidContainer(data.ds.grackle_data, size)

    flatten = len(data['gas', 'density'].shape) > 1

    fields = _get_needed_fields(fc.chemistry_data)
    
    for gfield in fields:
        yfield, units = _field_map[gfield]
        fdata = data[yfield].to(units)

        if flatten:
            fdata = fdata.flatten()
        
        fc[gfield][:] = fdata
        
    if 'de' in fc:
        fc['de'] *= (mp/me)
        
    return fc

def prepare_grackle_data(ds, parameters=None):
    sim_type = type(ds)
    par_map = _parameter_map.get(sim_type)
    if par_map is None and parameters is None:
        raise RuntimeError(
            "Simulation type not supported: %s." % sim_type)
    elif par_map is None and parameters is not None:
        print(f"Unknown simulation type {sim_type}." +
              "\nAttempting calculation with supplied parameters...\n")
        
        all_parameters = dict([(gpar, ds.parameters[dpar])
                               for gpar, dpar in parameters.items()
                               if dpar in ds.parameters])
        all_parameters.update(parameters)
    else:

        all_parameters = \
        dict([(gpar, ds.parameters[dpar])
                for gpar, dpar in par_map.items()
                if dpar in ds.parameters])
        all_parameters['use_grackle'] = 1

        if parameters is None:
            parameters = {}
        all_parameters.update(parameters)

    my_chemistry = chemistry_data()
    for gpar, val in all_parameters.items():
        if val is None:
            continue
        if isinstance(val, str):
            sval = bytes(val, 'utf-8')
            setattr(my_chemistry, gpar, sval)
        else:
            setattr(my_chemistry, gpar, val)

    my_chemistry.comoving_coordinates = ds.cosmological_simulation
    my_chemistry.density_units = (ds.mass_unit / ds.length_unit**3).in_cgs().d
    my_chemistry.length_units = ds.length_unit.in_cgs().d
    my_chemistry.time_units = ds.time_unit.in_cgs().d
    my_chemistry.a_units = 1 / (1 + ds.parameters.get('CosmologyInitialRedshift', 0))
    my_chemistry.a_value = 1 / (1 + ds.current_redshift) / my_chemistry.a_units
    my_chemistry.velocity_units = ds.velocity_unit.in_cgs().d
    my_chemistry.initialize()
    
    ds.grackle_data = my_chemistry

_grackle_fields = {
    'cooling_time': 'code_time',
    'dust_temperature': 'K',
    'gamma': '',
    'mean_molecular_weight': '',
    'pressure': 'code_mass * code_velocity**2 / code_length**3',
    'temperature': 'K',
    }

def _grackle_field(field, data):
    gfield = field.name[1][len("grackle_"):]
    units = _grackle_fields[gfield]

    if not hasattr(data.ds, "grackle_data"):
        raise RuntimeError("Grackle has not been initialized.")

    fc = _data_to_fc(data)
    if not isinstance(data, FieldDetector):
        func = "calculate_%s" % gfield
        # print("calling ", func, "with isrf_habing = ", fc["isrf_habing"][:10])
        getattr(fc, func)()

    fdata = fc[gfield]
    
    return fdata * data.ds.quan(1, units).in_cgs()

def _specific_thermal_energy(field, data):
    return data.ds.arr(data["PartType0", "specific_thermal_energy"].d, "code_velocity**2")

def _total_metal_density(field, data):
    metallicity = data["PartType0", "metallicity"]
    density = data["PartType0", "Density"]
    return (metallicity * density).to("code_mass/code_length**3")

def _dust_density(field, data):
    gas_density = data.ds.arr(data["PartType0", "Density"].d, "code_mass/code_length**3")
    gas_mass = data.ds.arr(data["PartType0", "Masses"].d, "code_mass")
    dust_mass = data.ds.arr(data["PartType0", "Dust_Masses"].d, "code_mass")
    dust_density =  dust_mass * gas_density / gas_mass
    return dust_density.to("code_mass/code_length**3")

def _primordial_species_density(field, data):
    """
    This function converts Simba's ['PartType0','Grackle<species>'] primordial fields into
     the ['gas', '<species>_<ionisation>_density'] fields expected for the yt_fields
     functions.
    """
    
    # This dictionary converts between the two naming conventions
    primordial_species = {"H_p0_density":'HI', "H_p1_density":'HII',
                        "He_p0_density":'HeI', "He_p1_density":'HeII',
                        "He_p2_density":'HeIII',
                        "H_m1_density":'HM', "H2_p0_density":'H2I',
                        "H2_p1_density":'H2II', "D_p0_density":"DI",
                        "D_p1_density":"DII", "HD_p0_density":"HDI"}
    
    fieldName = f"Grackle{primordial_species[field.name[1]]}"
    gasDensity = data["PartType0", "Density"]
    grackleSpeciesFraction = data[("PartType0", fieldName)]
    
    return (grackleSpeciesFraction * gasDensity).to("code_mass/code_length**3")

def _electron_density_nine_species_network(field, data):
    electron_density = data["gas", "H_p1_density"] + \
                       0.25 * data["gas", "He_p1_density"] + \
                       0.5 * data["gas", "He_p2_density"] + \
                       0.5 * data["gas", "H2_p1_density"] - \
                       data["gas", "H_m1_density"]
    return electron_density.to("code_mass/code_length**3")

def _isrf_habing(field, data):
    isrf_habing = data["PartType0", "LocalG0"].d
    return data.ds.arr(isrf_habing, "dimensionless")

def add_grackle_fields_simba(ds, parameters=None):
    
    prepare_grackle_data(ds, parameters=parameters)
    
    # Add specific thermal energy
    ds.add_field(("gas", "specific_thermal_energy"),
                 function=_specific_thermal_energy,
                 units="code_velocity**2",
                 sampling_type="local",
                 force_override=True)
    
    # Add metal density
    ds.add_field(("gas", "total_metal_density"),
                 function=_total_metal_density,
                 units="code_mass/code_length**3",
                 sampling_type="local",
                 force_override=True)
    
    # Add dust density
    ds.add_field(("gas", "dust_density"),
                 function=_dust_density,
                 units="code_mass/code_length**3",
                 sampling_type="local",
                 force_override=True)
    
    # Add isrf field
    ds.add_field(("gas", "isrf_habing"),
                 function=_isrf_habing,
                 units="dimensionless",
                 sampling_type="local",
                 force_override=True)
    
    # Add densities of all primordial chemical species
    primordial_species = []
    if parameters["primordial_chemistry"] > 0:
        primordial_species += ['H_p0', 'H_p1', 'He_p0', 'He_p1', 'He_p2']
    if parameters["primordial_chemistry"] > 1:
        primordial_species += ['H_m1', 'H2_p0', 'H2_p1']
    if parameters["primordial_chemistry"] > 2:
        primordial_species += ['D_p0', 'D_p1', 'HD_p0']
        
    for species in primordial_species:
        ds.add_field(("gas", species+"_density"),
                     function=_primordial_species_density,
                     units="code_mass/code_length**3",
                     sampling_type="local",
                     force_override=True)
        
    # Add electron density
    if parameters["primordial_chemistry"] == 2:
        ds.add_field(("gas", "El_density"),
                    function=_electron_density_nine_species_network,
                    units="code_mass/code_length**3",
                    sampling_type="local",
                    force_override=True)
        
    for field, units in _grackle_fields.items():
        fname = f"grackle_{field}"
        funits = str(ds.quan(1, units).in_cgs().units)
        ds.add_field(('gas', fname), function=_grackle_field,
                     sampling_type="local", units=funits,
                     force_override=True)    