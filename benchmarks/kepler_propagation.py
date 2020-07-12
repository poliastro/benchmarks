from astropy import units as u


class PropagationSuite:

    params = (['cowell', 'vallado', 'farnocchia'], [0.0, 0.5, 0.995, 1.5, 10.0, 100.0])
    param_names = ['method', 'ecc']

    def setup(self, *args):
        from poliastro.bodies import Earth
        from poliastro.twobody import Orbit
        from poliastro.twobody.propagation import cowell, vallado, farnocchia
        self.Earth = Earth
        self.Orbit = Orbit
        self.method = {
            'cowell': cowell,
            'vallado': vallado,
            'farnocchia': farnocchia,
        }
        
    def time_propagation(self, method, ecc): 
        a = 7000 * u.km
        _a = 0.0 * u.rad
        if ecc < 1.0:
            orbit = self.Orbit.from_classical(self.Earth, a, ecc * u.one, _a, _a, _a, _a)
        else:
            orbit = self.Orbit.from_classical(self.Earth, -a, ecc * u.one, _a, _a, _a, _a)
        orbit.propagate(1.0 * u.min, method=self.method[method])


if __name__=="__main__":
    from itertools import product
    suite = Suite()
    for param in suite.params:
        suite.setup()
        suite.time_propagation(method=method, ecc=ecc)
    