import functools

import numpy as np

from astropy import units as u
from astropy.time import Time
from astropy.coordinates import Angle, solar_system_ephemeris

from poliastro.twobody.propagation import cowell
from poliastro.core.elements import rv2coe
from poliastro.ephem import build_ephem_interpolant

from poliastro.core.util import norm
from poliastro.core.perturbations import (
    J2_perturbation, J3_perturbation, atmospheric_drag, third_body, radiation_pressure
)
from poliastro.bodies import Earth, Moon, Sun
from poliastro.twobody import Orbit


class J2_propagation():
    def setup(self):
        self.r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
        self.v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s
        self.k = Earth.k.to(u.km**3 / u.s**2).value

        self.orbit = Orbit.from_vectors(Earth, self.r0 * u.km, self.v0 * u.km / u.s)

        self.tof = (1.0 * u.min).to(u.s).value
    def time_J2_propagation(self):
        cowell(self.orbit, self.tof, ad=J2_perturbation, J2=Earth.J2.value, R=Earth.R.to(u.km).value)


class J3_propagation():
    def setup(self):
        self.a_ini = 8970.667 * u.km
        self.ecc_ini = 0.25 * u.one
        self.raan_ini = 1.047 * u.rad
        self.nu_ini = 0.0 * u.rad
        self.argp_ini = 1.0 * u.rad
        self.inc_ini = 0.2618 * u.rad

        self.k = Earth.k.to(u.km**3 / u.s**2).value

        self.orbit = Orbit.from_classical(Earth, self.a_ini, self.ecc_ini, self.inc_ini, self.raan_ini, self.argp_ini, self.nu_ini)
        self.tof = (1.0 * u.min).to(u.s).value
        self.a_J2J3 = lambda t0, u_, k_: J2_perturbation(t0, u_, k_, J2=Earth.J2.value, R=Earth.R.to(u.km).value) + \
                                         J3_perturbation(t0, u_, k_, J3=Earth.J3.value, R=Earth.R.to(u.km).value)
    def time_J3_propagation(self):
        cowell(self.orbit, self.tof, ad=self.a_J2J3, rtol=1e-8)


class drag():
    def setup(self):
        self.R = Earth.R.to(u.km).value
        self.k = Earth.k.to(u.km**3 / u.s**2).value

        self.orbit = Orbit.circular(Earth, 250 * u.km)

        # parameters of a body
        self.C_D = 2.2  # dimentionless (any value would do)
        self.A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2).value  # km^2
        self.m = 100  # kg
        self.B = self.C_D * self.A / self.m

        # parameters of the atmosphere
        self.rho0 = Earth.rho0.to(u.kg / u.km**3).value  # kg/km^3
        self.H0 = Earth.H0.to(u.km).value
        self.tof = 60  # s
    def time_atmospheric_drag(self):
        cowell(self.orbit, self.tof, ad=atmospheric_drag, R=self.R, C_D=self.C_D, A=self.A, m=self.m, H0=self.H0, rho0=self.rho0)
