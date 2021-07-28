import numpy as np

from yt.fields.field_detector import \
    FieldDetector
from yt.frontends.enzo.data_structures import \
    EnzoDataset
from yt.utilities.physical_constants import \
    me, mp
from pygrackle import \
    FluidContainer, \
    chemistry_data

_parameter_map = {}
_parameter_map[EnzoDataset] = {
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
    "H2_self_shielding": "H2_self_shielding",
    "grackle_data_file": "grackle_data_file",
    "UVbackground": "UVbackground",
    "Compton_xray_heating": "Compton_xray_heating",
    "LWbackground_intensity": "LWbackground_intensity",
    "LWbackground_sawtooth_suppression": "LWbackground_sawtooth_suppression"
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
    'energy': (('gas', 'thermal_energy'), 'code_velocity**2'),
    'RT_heating_rate': (('gas', 'photo_gamma'), 'erg/s'),
    'RT_HI_ionization_rate': (('gas', 'H_p0_ionization_rate'), '1 / code_time'),
    'RT_HeI_ionization_rate': (('gas', 'He_p0_ionization_rate'), '1 / code_time'),
    'RT_HeII_ionization_rate': (('gas', 'He_p1_ionization_rate'), '1 / code_time'),
    'RT_H2_dissociation_rate': (('gas', 'H2_p0_dissociation_rate'), '1 / code_time'),
    'dark_matter': (('data', 'dark_matter_density'), 'code_mass / code_length**3'),
    'external_pressure': (('data', 'external_pressure'), 'code_mass * code_velocity**2 / code_length**3'),
    'metallicity': (('data', 'metallicity3'), ''),
}

def _get_needed_fields(my_chemistry):
    fields = \
      ['density', 'energy'] + \
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
        fields += ['RT_heating_rate',
                   'RT_HI_ionization_rate',
                   'RT_HeI_ionization_rate',
                   'RT_HeII_ionization_rate',
                   'RT_H2_dissociation_rate']
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

def prepare_grackle_data(ds, parameters=None, sim_type=None, initialize=True):
    if sim_type is None:
        sim_type = type(ds)
    par_map = _parameter_map.get(sim_type)
    if par_map is None:
        raise RuntimeError(
            "Simulation type not supported: %s." % sim_type)

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

    if initialize:
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
        getattr(fc, func)()

    fdata = fc[gfield]
    if hasattr(data, 'ActiveDimensions'):
        fdata = fdata.reshape(data.ActiveDimensions)

    return fdata * data.ds.quan(1, units).in_cgs()

def _total_metal_density(field, data):
    field_data = data.ds.arr(np.zeros(data['index', 'ones'].shape),
                             'code_mass / code_length**3')
    fields = [
        ("enzo", "Metal_Density"),
        ("enzo", "SN_Colour")]

    for field in fields:
        if field not in data.ds.field_list:
            continue
        field_data += data[field]
    return field_data

def add_grackle_fields(ds, parameters=None):
    ds.add_field(('gas', 'total_metal_density'),
                 function=_total_metal_density,
                 units='g/cm**3',
                 sampling_type='cell')

    prepare_grackle_data(ds, parameters=parameters)
    for field, units in _grackle_fields.items():
        fname = "grackle_%s" % field
        funits = str(ds.quan(1, units).in_cgs().units)
        ds.add_field(fname, function=_grackle_field,
                     sampling_type="cell", units=funits)
