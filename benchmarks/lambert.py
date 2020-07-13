from astropy import units as u

from poliastro.bodies import Earth
from poliastro.iod import izzo, vallado


class Lambert:

    params = ["vallado.lambert", "izzo.lambert"]
    param_names = ["lambert"]

    def setup(self, *args):
        self.lambert = {
            "vallado.lambert": vallado.lambert,
            "izzo.lambert": izzo.lambert,
        }

    def time_lambert(self, lambert):
        k = Earth.k
        r0 = [15945.34, 0.0, 0.0] * u.km
        r = [12214.83399, 10249.46731, 0.0] * u.km
        tof = 76.0 * u.min
        next(self.lambert[lambert](k, r0, r, tof))
