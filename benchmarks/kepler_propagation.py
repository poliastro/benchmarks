from astropy import units as u

from poliastro.bodies import Earth
from poliastro.twobody.propagation import cowell, vallado, farnocchia


class KeplerPropagation:

    params = (["cowell", "vallado", "farnocchia"], [0.0, 0.5, 0.995, 1.5, 10.0, 100.0])
    param_names = ["method", "ecc"]

    def setup(self, *args):
        from poliastro.twobody import Orbit

        self.Orbit = Orbit
        self.method = {
            "cowell": cowell,
            "vallado": vallado,
            "farnocchia": farnocchia,
        }

    def time_propagation(self, method, ecc):
        a = 7000 * u.km
        _a = 0.0 * u.rad
        if ecc < 1.0:
            orbit = self.Orbit.from_classical(Earth, a, ecc * u.one, _a, _a, _a, _a)
        else:
            orbit = self.Orbit.from_classical(Earth, -a, ecc * u.one, _a, _a, _a, _a)
        orbit.propagate(1.0 * u.min, method=self.method[method])
