from poliastro.examples import iss


def time_propagate_iss_one_period():
    iss.propagate(iss.period)
