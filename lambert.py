from astropy import units as u

from poliastro.bodies import Earth

from poliastro.iod import izzo, vallado


def time_lambert(lambert):
    k = Earth.k
    r0 = [15945.34, 0.0, 0.0] * u.km
    r = [12214.83399, 10249.46731, 0.0] * u.km
    tof = 76.0 * u.min

    next(lambert(k, r0, r, tof))

time_lambert.params = ([vallado.lambert, izzo.lambert])
time_lambert.param_names = ['lambert']
