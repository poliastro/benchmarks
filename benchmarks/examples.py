from poliastro.examples import iss


def time_propagate_iss_one_period():
    iss.propagate(iss.period)

time_propagate_iss_one_period.repeat = 7
