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


def time_J2_propagation_Earth():
    r0 = np.array([-2384.46, 5729.01, 3050.46])  # km
    v0 = np.array([-7.36138, -2.98997, 1.64354])  # km/s
    k = Earth.k.to(u.km**3 / u.s**2).value

    orbit = Orbit.from_vectors(Earth, r0 * u.km, v0 * u.km / u.s)

    tof = (1.0 * u.min).to(u.s).value
    r, v = cowell(orbit, tof, ad=J2_perturbation, J2=Earth.J2.value, R=Earth.R.to(u.km).value)


def time_J3_propagation_Earth():
    a_ini = 8970.667 * u.km
    ecc_ini = 0.25 * u.one
    raan_ini = 1.047 * u.rad
    nu_ini = 0.0 * u.rad
    argp_ini = 1.0 * u.rad
    inc_ini = 0.2618 * u.rad

    k = Earth.k.to(u.km**3 / u.s**2).value

    orbit = Orbit.from_classical(Earth, a_ini, ecc_ini, inc_ini, raan_ini, argp_ini, nu_ini)

    tof = (1.0 * u.min).to(u.s).value
    r_J2, v_J2 = cowell(orbit, tof, ad=J2_perturbation,
                        J2=Earth.J2.value, R=Earth.R.to(u.km).value, rtol=1e-8)


def time_atmospheric_drag():
    R = Earth.R.to(u.km).value
    k = Earth.k.to(u.km**3 / u.s**2).value

    # parameters of a circular orbit with h = 250 km (any value would do, but not too small)
    orbit = Orbit.circular(Earth, 250 * u.km)

    # parameters of a body
    C_D = 2.2  # dimentionless (any value would do)
    A = ((np.pi / 4.0) * (u.m**2)).to(u.km**2).value  # km^2
    m = 100  # kg
    B = C_D * A / m

    # parameters of the atmosphere
    rho0 = Earth.rho0.to(u.kg / u.km**3).value  # kg/km^3
    H0 = Earth.H0.to(u.km).value
    tof = 60  # s
    r, v = cowell(orbit, tof, ad=atmospheric_drag, R=R, C_D=C_D, A=A, m=m, H0=H0, rho0=rho0)
