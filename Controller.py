# This is a 1D CFD Thermal Solver using Practice B
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def converged(temperature1, temperature2, method, p, x):

    if method == 'max_diff':
        return np.absolute(temperature1 - temperature2).max()
    elif method == 'heat_flux_l23':
        q1 = p['k'](x[-2])*(temperature1[-2]-temperature1[-3])/(x[-2]-x[-3])
        q2 = p['k'](x[-2])*(temperature2[-2]-temperature2[-3])/(x[-2]-x[-3])
        return math.fabs(q1-q2)
    elif method == 'heat_flux_f0':
        q1 = p['k'](x[1]) * (temperature1[0] - temperature1[1]) / (x[1] - x[0])
        q2 = p['k'](x[1]) * (temperature2[0] - temperature2[1]) / (x[1] - x[0])
        return math.fabs(q1 - q2)


def init_guess(method, x, t):

    if method == 'constant':
        return np.ones(len(x), dtype=float) * t
    elif method == 'analytical':
        exact = []
        n = math.sqrt(4. * t[2] / (t[3](0.) * t[4]))
        for x_v in x:
            t_act = (t[0] - t[1]) * math.cosh(n * (x[-1] - x_v)) / math.cosh(n * x[-1]) + t[1]
            exact.append(t_act)
        return np.asarray(exact)


def hw4_iterator(prop1, alpha, c_method, g_method, g_value):

    num_cv_list = [10, 20, 40, 80]

    tolerance = 0.000001

    # geometry basics
    x_max = 0.02

    # This variable contains a set of properties for any given segment
    properties = [prop1]

    # This variable contains the boundary conditions
    bound_i = {'type': 'fixed', 'T': 400.}
    bound_c = {'type': None}
    boundaries = [bound_i, bound_c]

    for num_cv in num_cv_list:

        cv_boundaries = [it * x_max / float(num_cv) for it in xrange(num_cv + 1)]
        node_properties = [0] * (num_cv + 1)
        node_properties.append(0)

        x = solver.assemble_geometry(cv_boundaries)

        criterion = 1
        temperature_prev = init_guess(g_method, x, g_value)
        temp_totals = []

        counter = 0
        x_pos = None

        # iterate until temperature convergence is reached
        while criterion > tolerance:
            # calculate temperatures
            temperature_current, x_pos = solver.solve(cv_boundaries, boundaries, properties, node_properties,
                                                      temperature_prev, method='gs', alpha=alpha)

            criterion = converged(temperature_current, temperature_prev, c_method, prop1, x_pos)

            # set previous
            temperature_prev = temperature_current
            temp_totals.append(temperature_prev)

            counter += 1
        print counter

        # for t in temp_totals:
        #     plt.plot(x_pos, t)
        # plt.show()


def prob1():
    prop1 = {'k': lambda x: 401., 'h': 10., 'e': 0, 'T': 273., 'D': 0.003}
    c_method = 'heat_flux_f0'
    alpha = 1.
    hw4_iterator(prop1, alpha, c_method, 'constant', 400.)


def prob2():
    prop1 = {'k': lambda x: 401., 'h': 10., 'e': 1., 'T': 273., 'D': 0.003}
    c_method = 'heat_flux_f0'
    alpha_list = [1., 1.2]

    for alpha in alpha_list:
        hw4_iterator(prop1, alpha, c_method, 'constant', 400.)


def prob3():

    prop1 = {'k': lambda x: 401., 'h': 10., 'e': 1., 'T': 273., 'D': 0.003}
    c_method = 'heat_flux_f0'
    alpha_list = [1., 1.2]

    for alpha in alpha_list:
        hw4_iterator(prop1, alpha, c_method, 'analytical', [400., 273., prop1['h'], prop1['k'], prop1['D']])


if __name__ == '__main__':
    prob3()
