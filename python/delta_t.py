c_special_units = 29.9792458


def vertex_time(sc_time, sc_pathlength, relativistic_beta):
    return (sc_time - sc_pathlength / (relativistic_beta * c_special_units))


def dt(electron_vertex_time, beta, sc_t, sc_r):
    return (electron_vertex_time - vertex_time(sc_t, sc_r, beta))


def delta_t(electron_vertex_time, mass, momentum, sc_t, sc_r):
    relativistic_beta = 1.0 / (sqrt(1.0 + mass / momentum) * (mass / momentum))
    return (electron_vertex_time - vertex_time(sc_t, sc_r, relativistic_beta))
