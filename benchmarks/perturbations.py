import functools

import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle, solar_system_ephemeris

from poliastro.twobody.propagation import cowell
from poliastro.core.elements import rv2coe
from poliastro.ephem import build_ephem_interpolant

from poliastro.core.perturbations import (
    J2_perturbation,
    J3_perturbation,
    third_body,
    radiation_pressure,
    atmospheric_drag_exponential,
)
from poliastro.bodies import Earth, Moon, Sun
from poliastro.constants import rho0_earth, H0_earth


class J2_Propagation:
    def setup(self):
        from poliastro.twobody import Orbit

        self.r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
        self.v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s
        self.k = Earth.k.to(u.km ** 3 / u.s ** 2).value

        self.orbit = Orbit.from_vectors(Earth, self.r0 * u.km, self.v0 * u.km / u.s)

        self.tof = [1.0] * u.min

    def time_J2_propagation(self):
        cowell(
            self.orbit.attractor.k,
            self.orbit.r,
            self.orbit.v,
            self.tof,
            ad=J2_perturbation,
            J2=Earth.J2.value,
            R=Earth.R.to(u.km).value,
        )


class J3_Propagation:
    def setup(self):
        from poliastro.twobody import Orbit

        self.a_ini = 8970.667 * u.km
        self.ecc_ini = 0.25 * u.one
        self.raan_ini = 1.047 * u.rad
        self.nu_ini = 0.0 * u.rad
        self.argp_ini = 1.0 * u.rad
        self.inc_ini = 0.2618 * u.rad

        self.k = Earth.k.to(u.km ** 3 / u.s ** 2).value

        self.orbit = Orbit.from_classical(
            Earth,
            self.a_ini,
            self.ecc_ini,
            self.inc_ini,
            self.raan_ini,
            self.argp_ini,
            self.nu_ini,
        )
        self.tof = [1.0] * u.min
        self.a_J2J3 = lambda t0, u_, k_: J2_perturbation(
            t0, u_, k_, J2=Earth.J2.value, R=Earth.R.to(u.km).value
        ) + J3_perturbation(t0, u_, k_, J3=Earth.J3.value, R=Earth.R.to(u.km).value)

    def time_J3_propagation(self):
        cowell(
            self.orbit.attractor.k,
            self.orbit.r,
            self.orbit.v,
            self.tof,
            ad=self.a_J2J3,
            rtol=1e-8,
        )


class AtmosphericDrag:
    def setup(self):
        from poliastro.twobody import Orbit

        self.R = Earth.R.to(u.km).value
        self.k = Earth.k.to(u.km ** 3 / u.s ** 2).value

        self.orbit = Orbit.circular(Earth, 250 * u.km)

        # parameters of a body
        self.C_D = 2.2  # dimentionless (any value would do)
        self.A_over_m = ((np.pi / 4.0) * (u.m ** 2) / (100 * u.kg)).to_value(
            u.km ** 2 / u.kg
        )  # km^2/kg
        self.B = self.C_D * self.A_over_m

        # parameters of the atmosphere
        self.rho0 = rho0_earth.to(u.kg / u.km ** 3)  # kg/km^3
        self.H0 = H0_earth.to(u.km).value
        self.tof = [60.0] * u.s

    def time_atmospheric_drag_exponential(self):
        cowell(
            self.orbit.attractor.k,
            self.orbit.r,
            self.orbit.v,
            self.tof,
            ad=atmospheric_drag_exponential,
            R=self.R,
            C_D=self.C_D,
            A_over_m=self.A_over_m,
            H0=self.H0,
            rho0=self.rho0,
        )
