# This is a 1D CFD Thermal Solver using Practice B
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def exact_study(num_cv_list, x_dist_max, bc_list, p, exact, max_cv):

    max_temp = []

    def get_temps(num):

        cv_boundaries = [it * x_dist_max / float(num) for it in xrange(num + 1)]  # x boundary locations
        node_properties = [0 if y <= 0.030000000001 else 1 for y in cv_boundaries]
        node_properties.append(1)

        temps, x_dist = solver.solve(cv_boundaries, bc_list, p, node_properties, None)
        max_temp.append(temps.max())

    get_temps(num_cv_list[-1])
    num_cv_list.append(num_cv_list[-1]*2)
    get_temps(num_cv_list[-1])

    while num_cv_list[-1] < max_cv:

        num_cv_list.append(num_cv_list[-1]*2)
        get_temps(num_cv_list[-1])

    print num_cv_list
    print max_temp

    diff_temps = []
    exact_x = []
    exact_y = []

    for i in xrange(len(max_temp)):
        diff_temps.append((exact-max_temp[i])/exact)
        exact_x.append(num_cv_list[i])
        exact_y.append(exact)
    print diff_temps

    plt.semilogx(exact_x, exact_y, 'r--')
    plt.semilogx(num_cv_list, max_temp)
    plt.show()


def prob1():

    exact = 588.117139215064
    max_cv = 10240

    # properties dictionary
    prop1 = {'k': lambda x: 15, 'q': 4000000}
    prop2 = {'k': lambda x: 137 * math.exp(25 * x - 2), 'q': 0}

    # This variable contains a set of properties for any given segment
    properties = [prop1, prop2]

    x_max = 0.05

    # This variable contains the boundary conditions
    bound_i = {'type': 'insulated'}
    bound_c = {'type': 'convective', 'h': 1000, 'To': 300}
    boundaries = [bound_i, bound_c]

    # for grid convergence study
    number_control_volumes = [5]
    exact_study(number_control_volumes, x_max, boundaries, properties, exact, max_cv)


def prob2(number_control_volumes):

    tolerance = 0.0001

    # geometry basics
    x_max = 0.02

    # This variable contains a set of properties for any given segment
    prop1 = {'k': lambda x: 401, 'h': 1000, 'e': 0, 'T': 273, 'D': 0.003}
    properties = [prop1]

    cv_boundaries = [it * x_max / float(number_control_volumes) for it in xrange(number_control_volumes + 1)]
    node_properties = [0]*(number_control_volumes+1)
    node_properties.append(0)

    # This variable contains the boundary conditions
    bound_i = {'type': 'fixed', 'T': 400}
    bound_c = {'type': 'convective_radiative', 'h': 10, 'e': 0, 'T': 273}
    boundaries = [bound_i, bound_c]

    temperature_diffs = np.ones(number_control_volumes+2, dtype=float)
    temperature_prev = np.ones(number_control_volumes+2, dtype=float)*400.

    total_temps = [temperature_prev]
    x_pos = None

    # iterate until temperature convergence is reached
    while temperature_diffs.max() > tolerance:

        # calculate temperatures
        temperature_current, x_pos = solver.solve(cv_boundaries, boundaries, properties, node_properties,
                                                  temperature_prev)

        # get difference
        temperature_diffs = np.absolute(temperature_current-temperature_prev)

        # set previous
        temperature_prev = temperature_current

        total_temps.append(temperature_prev)

    exact = []
    t_inf = prop1['T']
    t_b = bound_i['T']
    n = math.sqrt(4.*prop1['h']/(prop1['k'](1.)*prop1['D']))

    for x_v in x_pos:
        t_act = (t_b-t_inf)*math.cosh(n*(x_max-x_v))/math.cosh(n*x_max)+t_inf
        exact.append(t_act)

    return temperature_prev, exact, x_pos


def prob2a():

    t_sim, t_exact, x = prob2(5)

    # for plotting final temperature vs exact temperature
    plt.plot(x, t_sim, 'b*')
    plt.plot(x, t_exact, 'r--')
    plt.show()


def prob2b():

    max_cvs = 10240
    ctr_vol = [5]
    diff_list = []

    def array_diff(v1, v2):
        max_diff = 0
        for i in xrange(len(v1)):
            diff = math.fabs(v1[i]-v2[i])
            if diff > max_diff:
                max_diff = diff
        return max_diff

    t_sim, t_exact, x = prob2(ctr_vol[-1])
    diff_list.append(array_diff(t_exact, t_sim))

    # loop until max cv
    while ctr_vol[-1] < max_cvs:

        ctr_vol.append(ctr_vol[-1]*2)
        t_sim, t_exact, x = prob2(ctr_vol[-1])
        diff_list.append(array_diff(t_exact, t_sim))

    # plot errors
    plt.semilogx(ctr_vol, diff_list)
    plt.show()


if __name__ == '__main__':

    prob2b()
