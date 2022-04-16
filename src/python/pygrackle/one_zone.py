"""
OneZoneModel class


"""

import abc
from collections import defaultdict
import functools
import numpy as np
from unyt import unyt_array, unyt_quantity
from unyt.unit_registry import UnitRegistry

from pygrackle.utilities.physical_constants import \
    gravitational_constant_cgs, \
    sec_per_year

class OneZoneModel(abc.ABC):
    """
    Class for running one-zone models with a pygrackle fluid container.
    """

    name = None
    verbose = True
    stopping_criteria = ()

    def __init__(self, fc, data=None, external_data=None,
                 unit_registry=None, event_trigger_fields=None):

        self.fc = fc
        self.data = data

        if unit_registry is None:
            unit_registry = UnitRegistry()
        self.unit_registry = unit_registry

        if external_data is None:
            external_data = []
        self.external_data = external_data

        for field in self.external_data:
            if field not in self.fc.fields:
                self.fc._setup_fluid(field)

        self.check_stopping_criteria()
        self.prepare_event_times(event_trigger_fields)

    def prepare_event_times(self, fields):
        """
        Create a list of times to make sure we hit.

        These fields need to be in the external data.
        """

        self.events = []
        if fields is None:
            return

        for field in fields:
            data = self.external_data[field]
            onoff = np.where(~((data[:-1] > 0) ^ (data[1:] <= 0)))[0] + 1
            for t in self.external_data["time"][onoff]:
                if t not in self.events:
                    self.events.append(t)

        self.events = np.array(self.events)

    _arr = None
    @property
    def arr(self):
        """
        Create a unyt_array using the Arbor's unit registry.
        """
        if self._arr is not None:
            return self._arr
        self._arr = functools.partial(unyt_array,
                                      registry=self.unit_registry)
        return self._arr

    _quan = None
    @property
    def quan(self):
        """
        Create a unyt_quantity using the Arbor's unit registry.
        """
        if self._quan is not None:
            return self._quan
        self._quan = functools.partial(unyt_quantity,
                                       registry=self.unit_registry)
        return self._quan

    def check_stopping_criteria(self):
        canstop = False
        for cri in self.stopping_criteria:
            if getattr(self, cri, None) is not None:
                canstop = True
                break
        if not canstop:
            raise RuntimeError(
                f"Must specify at least one stopping criteria: {self.stopping_criteria}.")

    def calculate_timestep(self):
        self.fc.calculate_cooling_time()
        dt = self.safety_factor * np.abs(self.fc["cooling_time"][0])
        dt = min(dt, self.remaining_time)
        return dt

    @property
    def external_fields(self):
        return [field for field in self.external_data if field != "time"]

    @property
    def remaining_time(self):
        """
        Return time to nearest event or final time.
        """
        if len(self.events):
            upcoming = self.events > self.current_time
            if not upcoming.any():
                t_next = np.inf
            else:
                t_next = (self.events[upcoming] - self.current_time).min()
        else:
            t_next = np.inf

        if self.final_time is None:
            final_time = np.inf
        else:
            final_time = self.final_time

        rt = final_time - self.current_time

        return min(t_next, rt)

    @abc.abstractproperty
    def finished(self):
        pass

    @abc.abstractmethod
    def update_quantities(self):
        pass

    def update_external_fields(self):
        if not self.external_data:
            return

        edata = self.external_data
        time = edata["time"]
        itime = np.digitize(self.current_time, time) - 1
        if itime < 0 or itime >= time.size - 1:
            return

        new_fields = {}
        for field in self.external_fields:
            fdata = edata[field]
            if fdata[itime] < 0 or fdata[itime+1] < 0:
                new_fields[field] = fdata[itime]
            else:
                fdata = np.log(np.clip(fdata, a_min=1e-50, a_max=np.inf))
                slope = (fdata[itime+1] - fdata[itime]) / (time[itime+1] - time[itime])
                new_fields[field] = \
                  np.exp(slope * (self.current_time - time[itime]) + fdata[itime])

        for field in new_fields:
            self.fc[field][:] = new_fields[field]

    def add_to_data(self):
        """
        Add current fluid container values to the data structure.
        """

        fc = self.fc
        if self.data is None:
            self.data = defaultdict(list)
        data = self.data

        data["time"].append(self.current_time)

        for field in fc.density_fields:
            data[field].append(fc[field][0])
        data["energy"].append(fc["energy"][0])

        fc.calculate_gamma()
        data["gamma"].append(fc["gamma"][0])

        fc.calculate_temperature()
        data["temperature"].append(fc["temperature"][0])

        fc.calculate_pressure()
        data["pressure"].append(fc["pressure"][0])

        fc.calculate_mean_molecular_weight()
        data["mean_molecular_weight"].append(fc["mean_molecular_weight"][0])

        if fc.chemistry_data.h2_on_dust:
            fc.calculate_dust_temperature()
            data["dust_temperature"].append(fc["dust_temperature"][0])

        for field in self.external_fields:
            data[field].append(fc[field][0])

    def finalize_data(self):
        """
        Turn lists of values into array with proper cgs units.
        """

        fc = self.fc
        my_chemistry = fc.chemistry_data
        data = self.data

        for field in data:
            if field in fc.density_fields:
                data[field] = my_chemistry.density_units * \
                  self.arr(data[field], "g/cm**3")
            elif field == "energy":
                data[field] = my_chemistry.energy_units * \
                  self.arr(data[field], "erg/g")
            elif "time" in field:
                data[field] = my_chemistry.time_units * \
                  self.arr(data[field], "s")
            elif "temperature" in field:
                data[field] = self.arr(data[field], "K")
            elif "pressure" in field:
                data[field] = my_chemistry.pressure_units * \
                  self.arr(data[field], "dyne/cm**2")
            elif "radius" in field:
                data[field] = my_chemistry.length_units * \
                  self.arr(data[field], "cm")
            elif "mass" in field:
                data[field] = my_chemistry.density_units * my_chemistry.length_units**3 * \
                  self.arr(data[field], "g")
            elif "velocity" in field:
                data[field] = my_chemistry.velocity_units * \
                  self.arr(data[field], "cm/s")
            else:
                data[field] = np.array(data[field])
        return data

    def print_status(self):
        if not self.verbose:
            return

        fc = self.fc
        my_chemistry = fc.chemistry_data
        ctime = self.current_time * my_chemistry.time_units / sec_per_year
        cdensity = fc["density"][0] * my_chemistry.density_units
        ctemperature = fc["temperature"][0]

        status = f"Evolve {self.name} - t: {ctime:e} yr, " + \
          f"rho: {cdensity:e} g/cm^3, T: {ctemperature:e} K."
        print (status, flush=True)

    def scale_density_fields(self, factor, exclude=None):
        fc = self.fc

        if exclude is None:
            exclude = ["dark_matter"]

        for field in fc.density_fields:
            if field in exclude:
                continue
            fc[field][:] *= factor

    def before_solve_chemistry(self):
        pass

    def evolve(self):
        if self.data is None:
            self.current_time = 0
            self.update_external_fields()
            self.add_to_data()
        else:
            self.current_time = self.data["time"][-1]

        self.print_status()

        while not self.finished:
            self.update_external_fields()
            self.before_solve_chemistry()

            dt = self.calculate_timestep()

            self.fc.solve_chemistry(dt)

            self.current_time += dt

            self.update_quantities()
            self.print_status()
            self.add_to_data()

class CoolingModel(OneZoneModel):
    stopping_criteria = ("final_time", "final_temperature")

    def __init__(self, fc, data=None, safety_factor=0.01,
                 final_time=None, final_temperature=None,
                 external_data=None):

        self.safety_factor = safety_factor
        self.final_time = final_time
        self.final_temperature = final_temperature
        super().__init__(fc, data=data, external_data=external_data)

    @property
    def finished(self):
        if self.final_temperature is not None and \
          self.fc["temperature"][0] <= self.final_temperature:
            return True

        if self.final_time is not None and \
          self.current_time >= self.final_time:
            return True

        return False

class ConstantDensityModel(CoolingModel):
    name = "constant density"

    def update_quantities(self):
        pass

class ConstantPressureModel(CoolingModel):
    name = "constant pressure"

    def update_quantities(self):
        fc = self.fc
        my_chemistry = fc.chemistry_data

        T1 = self.data["temperature"][-1]
        mu1 = self.data["mean_molecular_weight"][-1]
        fc.calculate_temperature()
        T2 = fc["temperature"][0]
        fc.calculate_mean_molecular_weight()
        mu2 = fc["mean_molecular_weight"][0]
        factor = (T1 * mu2) / (T2 * mu1)

        self.scale_density_fields(factor)

class ConstantEntropyModel(ConstantPressureModel):
    name = "constant entropy"

    def update_quantities(self):
        fc = self.fc
        my_chemistry = fc.chemistry_data

        T1 = self.data["temperature"][-1]
        mu1 = self.data["mean_molecular_weight"][-1]
        rho1 = self.data["density"][-1]
        g1 = self.data["gamma"][-1]

        fc.calculate_temperature()
        T2 = fc["temperature"][0]
        fc.calculate_mean_molecular_weight()
        mu2 = fc["mean_molecular_weight"][0]
        fc.calculate_gamma()
        g2 = fc['gamma'][0]
        factor = (mu2 / mu1) * (rho1 / mu1)**((g1-g2)/(g2-1)) * \
          (T2 / T1)**(1 / (g2 - 1))

        self.scale_density_fields(factor)

class FreeFallModel(OneZoneModel):
    name = "free-fall"
    stopping_criteria = ("final_time", "final_density")

    def __init__(self, fc, data=None,
                 external_data=None, unit_registry=None,
                 safety_factor=0.01, include_pressure=True,
                 final_time=None, final_density=None,
                 event_trigger_fields=None):

        self.include_pressure = include_pressure
        self.safety_factor = safety_factor
        self.final_time = final_time
        self.final_density = final_density
        super().__init__(fc, data=data, external_data=external_data,
                         unit_registry=unit_registry,
                         event_trigger_fields=event_trigger_fields)

    @property
    def gravitational_constant(self):
        """
        Gravitational constant in internal units.
        """
        my_chemistry = self.fc.chemistry_data
        val = gravitational_constant_cgs * \
          my_chemistry.density_units * my_chemistry.time_units**2
        return val

    @property
    def freefall_constant(self):
        """
        Constant used in analytical free-fall solution.
        """
        return np.sqrt((32 * self.gravitational_constant) / (3. * np.pi))

    @property
    def finished(self):
        if self.final_density is not None and \
          self.fc["density"][0] >= self.final_density:
            return True

        if self.final_time is not None and \
          self.current_time >= self.final_time:
            return True

        return False

    def calculate_timestep(self):
        fc = self.fc
        dt_ff = self.safety_factor / self.freefall_constant / \
          np.sqrt(fc["density"][0])

        fc.calculate_cooling_time()
        dt_cool = self.safety_factor * \
          np.abs(fc["cooling_time"][0])

        dt = min(dt_ff, dt_cool)
        dt = min(dt, self.remaining_time)
        self.dt = dt
        return dt

    def update_quantities(self):
        fc = self.fc
        my_chemistry = fc.chemistry_data

        self.calculate_collapse_factor()
        force_factor = fc["force_factor"][0]

        self.fc.calculate_pressure()
        pressure = self.fc["pressure"][0]
        density = self.fc["density"][0]
        factor = np.sqrt(1 - force_factor) * \
          self.freefall_constant * np.sqrt(density) * self.dt + 1

        # update energy assuming un-altered free-fall collapse
        e_factor = self.freefall_constant * np.sqrt(density) * self.dt + 1

        self.scale_density_fields(factor)

        de = - pressure * (1 - e_factor) / (e_factor * density)
        fc["energy"][0] += de

    def calculate_collapse_factor(self):
        """
        Compute the new density using the modified
        free-fall collapse as per Omukai et al. (2005).
        """

        data = self.data

        if not self.include_pressure:
            self.fc["force_factor"] = np.array([0])
            return

        # Calculate the effective adiabatic index, dlog(p)/dlog(rho).
        density = data["density"]
        pressure = data["pressure"]

        if len(pressure) < 3:
            self.fc["force_factor"] = np.array([0])
            return

        # compute dlog(p) / dlog(rho) using last two timesteps
        gamma_eff = np.log10(pressure[-1] / pressure[-2]) / \
            np.log10(density[-1] / density[-2])

        # compute a higher order derivative if more than two points available
        if len(pressure) > 2:
            gamma_eff += 0.5 * ((np.log10(pressure[-2] / pressure[-3]) /
                                 np.log10(density[-2] / density[-3])) - gamma_eff)

        gamma_eff = min(gamma_eff, 4./3.)

        # Equation 9 of Omukai et al. (2005)
        if gamma_eff < 0.83:
            force_factor = 0.0
        elif gamma_eff < 1.0:
            X = gamma_eff - 1
            force_factor = 0.6 + 2.5 * X - 6.0 * np.power(X, 2.)
        else:
            X = gamma_eff - 4/3
            force_factor = 1.0 + 0.2 * X - 2.9 * np.power(X, 2.)
        force_factor = max(force_factor, 0.0)
        force_factor = min(force_factor, 0.95)

        self.fc["force_factor"] = np.array([force_factor])

class MinihaloModel(FreeFallModel):
    name = "minihalo"
    stopping_criteria = ("final_time", "final_density", "gas_mass")
    use_dark_matter = False

    def __init__(self, fc, data=None,
                 external_data=None, unit_registry=None,
                 safety_factor=0.01, include_pressure=True,
                 final_time=None, final_density=None,
                 initial_radius=None, gas_mass=None,
                 include_turbulence=True,
                 event_trigger_fields=None, cosmology=None):

        self.initial_radius = initial_radius
        self.gas_mass = gas_mass
        self.include_turbulence = include_turbulence

        super().__init__(fc, data=data,
                 external_data=external_data,
                 unit_registry=unit_registry,
                 safety_factor=safety_factor,
                 include_pressure=include_pressure,
                 final_time=final_time,
                 final_density=final_density,
                 event_trigger_fields=event_trigger_fields)

        self.cosmology = cosmology
        self.initialize_cosmology()

    def initialize_cosmology(self):
        if self.cosmology is None:
            self.initial_redshift = None
            self.initial_cosmo_time = None
            return

        my_chemistry = self.fc.chemistry_data
        a_initial = my_chemistry.a_value * my_chemistry.a_units
        self.initial_redshift = 1 / a_initial - 1
        self.initial_cosmo_time = \
          self.cosmology.t_from_z(self.initial_redshift).to('s').d / \
          my_chemistry.time_units

    @property
    def current_redshift(self):
        if self.cosmology is None:
            return None
        cosmo_time = self.fc.chemistry_data.time_units * \
          (self.current_time + self.initial_cosmo_time)
        return self.cosmology.z_from_t(cosmo_time)

    @property
    def finished(self):
        if super().finished:
            return True

        if self.gas_mass is not None:
            m_BE = self.calculate_bonnor_ebert_mass()
            if self.gas_mass >= m_BE:
                return True

        return False

    def before_solve_chemistry(self):
        if "metallicity" in self.external_fields:
            self.fc["metal"][0] = self.fc["density"][0] * self.fc["metallicity"][0]
        self.fc.chemistry_data.override_redshift = self.current_redshift

    @property
    def current_radius(self):
        return self.initial_radius * (self.data["density"][0] /
                                      self.fc["density"][0])**(1/3)

    def update_quantities(self):
        fc = self.fc
        my_chemistry = fc.chemistry_data

        tcs = self.calculate_sound_crossing_time()
        tff = self.calculate_free_fall_time()

        fc.calculate_pressure()
        pressure = fc["pressure"][0]

        # free-fall
        if (tff < tcs):
            self.calculate_collapse_factor()
            force_factor = fc["force_factor"][0]

            total_density = fc["density"][0]
            if self.use_dark_matter:
                total_density += fc["dark_matter"][0]

            factor = np.sqrt(1 - force_factor) * \
              self.freefall_constant * np.sqrt(total_density) * self.dt + 1

            # update energy assuming un-altered free-fall collapse
            e_factor = self.freefall_constant * np.sqrt(total_density) * self.dt + 1

        # pressure-dominated
        else:
            external_pressure = self.data["external_pressure"][-1]

            P1 = self.data["pressure"][-1]
            T1 = self.data["temperature"][-1]
            mu1 = self.data["mean_molecular_weight"][-1]
            P2 = max(pressure, external_pressure)
            fc.calculate_temperature()
            T2 = fc["temperature"][0]
            fc.calculate_mean_molecular_weight()
            mu2 = fc["mean_molecular_weight"][0]
            factor = (P2 * T1 * mu2) / (T2 * mu1 * P1)

            e_factor = factor

        self.scale_density_fields(factor)

        de = - pressure * (1 - e_factor) / (e_factor * self.fc["density"][0])
        fc["energy"][0] += de

    def calculate_sound_speed(self):
        fc = self.fc
        my_chemistry = fc.chemistry_data

        density = fc["density"][0]
        fc.calculate_pressure()
        pressure = fc["pressure"][0]
        fc.calculate_gamma()
        gamma = fc["gamma"][0]

        cs = np.sqrt(gamma * pressure / density)
        if self.include_turbulence:
            v_turb = self.data["turbulent_velocity"][-1]
            cs = np.sqrt(cs**2 + v_turb**2)

        return cs

    def calculate_sound_crossing_time(self):
        cs = self.calculate_sound_speed()
        return 2 * self.current_radius * \
          self.fc.chemistry_data.a_value / cs

    def calculate_free_fall_time(self):
        density = self.fc["density"][0]
        if self.use_dark_matter:
            density += self.data["dark_matter"][-1]
        return 1 / (self.freefall_constant * np.sqrt(density))

    def calculate_bonnor_ebert_mass(self):
        ### Bonnor-Ebert Mass constant
        a = 1.67
        b = (225 / (32 * np.sqrt(5 * np.pi))) * a**-1.5

        fc = self.fc
        my_chemistry = fc.chemistry_data

        fc.calculate_pressure()
        pressure = fc["pressure"][0]

        # Convert to CGS because I am tired of
        # messing up cosmological units.
        # Ignore turbulent velocity for cs in BE mass.
        include_turbulence = self.include_turbulence
        self.include_turbulence = False
        cs_cgs = self.calculate_sound_speed() * \
          my_chemistry.velocity_units
        self.include_turbulence = include_turbulence

        G = gravitational_constant_cgs
        pressure_cgs = pressure * my_chemistry.pressure_units
        m_BE = (b * (cs_cgs**4 / G**1.5) * pressure_cgs**-0.5)

        mass_units = my_chemistry.density_units * \
          my_chemistry.length_units**3
        return m_BE / mass_units

    def add_to_data(self):
        super().add_to_data()

        fc = self.fc

        t_ff = self.calculate_free_fall_time()
        self.data["free_fall_time"].append(t_ff)

        t_cs = self.calculate_sound_crossing_time()
        self.data["sound_crossing_time"].append(t_cs)

        # Turn off CMB floor to calculate cooling time
        cmb = fc.chemistry_data.cmb_temperature_floor
        fc.chemistry_data.cmb_temperature_floor = 0
        fc.calculate_cooling_time()
        self.data["cooling_time"].append(fc["cooling_time"][0])
        fc.chemistry_data.cmb_temperature_floor = cmb

        m_BE = self.calculate_bonnor_ebert_mass()
        self.data["mass_BE"].append(m_BE)
