# This is a 1D CFD Thermal Solver using Practice B
import matplotlib.pyplot as plt
import math
import solver
import p3


def convergence_study(num_cv_list, x_dist_max, bc_list, p, crit):

    max_temp = []

    def get_temps(num):

        cv_boundaries = [it * x_dist_max / float(num) for it in xrange(num + 1)]  # x boundary locations
        node_properties = [0 if y <= 0.030000000001 else 1 for y in cv_boundaries]
        node_properties.append(1)

        temps, x_dist = solver.solve(cv_boundaries, bc_list, p, node_properties)
        max_temp.append(temps.max())
        print max_temp[-1]

    get_temps(num_cv_list[-1])
    num_cv_list.append(num_cv_list[-1]*2)
    get_temps(num_cv_list[-1])

    while abs(max_temp[-1]-max_temp[-2]) > crit:

        num_cv_list.append(num_cv_list[-1]*2)
        print num_cv_list[-1]
        get_temps(num_cv_list[-1])

    plt.plot(num_cv_list, max_temp)
    plt.show()


def prob3():
    p3.solve()


def prob4():

    # properties dictionary
    prop1 = {'k': lambda x: 15, 'q': 4000000}
    prop2 = {'k': lambda x: 60, 'q': 0}

    # This variable contains a set of properties for any given segment
    properties = [prop1, prop2]

    # use for equally spaced control volumes
    num_cvs = 20
    x_max = 0.05
    cv_boundaries = [it * x_max / float(num_cvs) for it in xrange(num_cvs + 1)]  # x boundary locations
    node_properties = [0 if y <= 0.030000000001 else 1 for y in cv_boundaries]
    node_properties.append(1)

    # This variable contains the boundary conditions
    bound_i = {'type': 'insulated'}
    bound_c = {'type': 'convective', 'h': 1000, 'To': 300}
    boundaries = [bound_i, bound_c]

    py, px = solver.solve(cv_boundaries, boundaries, properties, node_properties)

    plt.plot(px, py)
    plt.show()


def prob5ab():
    # properties dictionary
    prop1 = {'k': lambda x: 15, 'q': 4000000}
    prop2 = {'k': lambda x: 137 * math.exp(25 * x - 2), 'q': 0}

    # This variable contains a set of properties for any given segment
    properties = [prop1, prop2]

    # use for equally spaced control volumes
    num_cvs = 20
    x_max = 0.05
    cv_boundaries = [it * x_max / float(num_cvs) for it in xrange(num_cvs + 1)]  # x boundary locations
    node_properties = [0 if y <= 0.030000000001 else 1 for y in cv_boundaries]
    node_properties.append(1)

    # This variable contains the boundary conditions
    bound_i = {'type': 'insulated'}
    bound_c = {'type': 'convective', 'h': 1000, 'To': 300}
    boundaries = [bound_i, bound_c]

    py, px = solver.solve(cv_boundaries, boundaries, properties, node_properties)

    plt.plot(px, py)
    plt.show()


def prob5c():

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
    crt = 1
    convergence_study(number_control_volumes, x_max, boundaries, properties, crt)


if __name__ == '__main__':

    prob5c()
