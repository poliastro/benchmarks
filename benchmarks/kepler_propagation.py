from poliastro.bodies import Sun, Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import cowell, kepler, mean_motion

from astropy import units as u

def time_propagation(method, ecc):
    a = 7000 * u.km
    _a = 0.0 * u.rad
    if ecc < 1.0:
        orbit = Orbit.from_classical(Earth, a, ecc * u.one, _a, _a, _a, _a)
    else:
        orbit = Orbit.from_classical(Earth, -a, ecc * u.one, _a, _a, _a, _a)
    orbit.propagate(1.0 * u.min, method=method)

time_propagation.params = ([cowell, kepler, mean_motion], [0.0, 0.5, 0.995, 1.5, 10.0, 100.0])
time_propagation.param_names = ['method', 'ecc']
