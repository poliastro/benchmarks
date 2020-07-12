from astropy import units as u

# def time_propagation(method, ecc):

class Suite:

    def setup(self):
        from poliastro.bodies import Earth
        from poliastro.twobody import Orbit
        from poliastro.twobody.propagation import cowell
        self.Earth = Earth
        self.Orbit = Orbit
        self.cowell = cowell

    def time_propagation(self):
        method = self.cowell
        ecc = 0.5
        a = 7000 * u.km
        _a = 0.0 * u.rad
        if ecc < 1.0:
            orbit = self.Orbit.from_classical(self.Earth, a, ecc * u.one, _a, _a, _a, _a)
        else:
            orbit = self.Orbit.from_classical(self.Earth, -a, ecc * u.one, _a, _a, _a, _a)
        orbit.propagate(1.0 * u.min, method=method)

# time_propagation.params = ([cowell, kepler, mean_motion], [0.0, 0.5, 0.995, 1.5, 10.0, 100.0])
# time_propagation.param_names = ['method', 'ecc']
# time_propagation.params = ([cowell], [0.0, 0.5, 0.995, 1.5, 10.0, 100.0])
# time_propagation.param_names = ['method', 'ecc']


if __name__=="__main__":
    from poliastro.bodies import Earth
    from poliastro.twobody import Orbit
    from poliastro.twobody.propagation import cowell
    suite = Suite()
    suite.time_propagation()
    