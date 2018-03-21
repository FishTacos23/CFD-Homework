# This is a 1D CFD Thermal Solver using Practice A
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def prob_2_grid():

    """
    :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
    """

    bcs = ['fixed', 'adiabatic']

    # FIXED PROPERTIES
    rho = 2770.
    c = 875.
    k = 177.
    l = .01
    t_inf = 293.
    t_sur = 293.
    tb = 353.
    w = 1.
    b = 0.0003
    t_init = 293.
    crt = 0.0001
    h = 30.
    e = 1.

    lines = [':', '-.', '--', '-', '--', '-.', ':']

    # num grid
    it = 0
    max_div = 8
    for num_div in xrange(2, max_div+1):

        num_cells = math.pow(2, num_div)
        dx = l/num_cells
        t_stop = 10
        dt = 0.1

        alf = .5*float(num_cells)/math.pow(2., float(max_div)) + .5

        t, x, base, middle, tip, q_log, t_log = solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h,
                                                             t_inf, t_sur, True, crt)
        # plt.plot(x, t, color='k', alpha=alf, label=str(num_cells)+' Cells', linestyle=lines[it])
        plt.plot(t_log, q_log, color='k', alpha=alf, label=str(num_cells)+' Cells', linestyle=lines[it])

        it += 1

    plt.legend()
    plt.show()


def prob_2_time():

    """
    :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
    """

    bcs = ['fixed', 'adiabatic']

    # FIXED PROPERTIES
    rho = 2770.
    c = 875.
    k = 177.
    l = .01
    t_inf = 293.
    t_sur = 293.
    tb = 353.
    w = 1.
    b = 0.0003
    t_init = 293.
    crt = 0.0001
    h = 30.
    e = 1.

    lines = ['-', '--', '-.', ':']

    # num grid
    max_div = 4
    it = 0
    for num_div in xrange(1, max_div+1):

        num_cells = 32
        dx = l/num_cells
        t_stop = 10
        dt = 1./(math.pow(10., num_div))

        alf = .5*dt/math.pow(10., float(max_div)) + .5

        t, x, base, middle, tip, q_log, t_log = solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h,
                                                             t_inf, t_sur, True, crt)
        # plt.plot(x, t, color='k', alpha=alf, label=str(dt) + ' s', linestyle=lines[it])
        plt.plot(t_log, q_log, color='k', alpha=alf, label=str(dt)+' s', linestyle=lines[it])
        it += 1

    plt.legend()
    plt.show()


def compare_sample_solution():
    """
        :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
        """

    bcs = ['fixed', 'adiabatic']

    # FIXED PROPERTIES
    rho = 2770.
    c = 875.
    k = 177.
    l = .01
    t_inf = 293.
    t_sur = 293.
    tb = 353.
    w = 1.
    b = 0.0003
    t_init = 293.
    crt = 0.0001
    h = 5.
    e = 1.

    # num grid
    dx = l / 160
    t_stop = 10
    dt = .001

    t, x, base, middle, tip, q_log, t_log = solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h,
                                                         t_inf, t_sur, True, crt)

    plt.plot(t_log, base, 'k')
    plt.plot(t_log, middle, 'k--')
    plt.plot(t_log, tip, 'k-.')

    plt.ylim((290, 360))
    plt.show()


def analytical(x, h, k, w, b, l, t_inf, tb):
    """

    :param x: Position of nodes
    :type x: ndarray
    :param h: coefficient of heat transfer
    :type h: float
    :param k: thermal conductivity
    :type k: float
    :param w: unit width
    :type w: float
    :param b: thickness of plate
    :type b: float
    :param l: Domain Length
    :type l: float
    :param t_inf: temperature at infinity
    :type t_inf: float
    :param tb: boundary temperature
    :type tb: float
    :return: temperatures for nodes using the analytical solution
    """

    p = 2.*w + 2.*b
    a = w*b
    m = math.sqrt(h*p/(k*a))

    soln = np.cosh(m*(l-x))*(tb-t_inf)/math.cosh(m*l) + t_inf

    return soln


def compare_exact():
    """
            :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
            """

    bcs = ['fixed', 'adiabatic']

    # FIXED PROPERTIES
    rho = 2770.
    c = 875.
    k = 177.
    l = .01
    t_inf = 293.
    t_sur = 293.
    tb = 353.
    w = 1.
    b = 0.0003
    t_init = 293.
    crt = 0.0001
    h = 5.
    e = 0.

    # num grid
    dx = l / 32
    t_stop = 10
    dt = .01

    t, x, base, middle, tip, q_log, t_log = solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h,
                                                         t_inf, t_sur, True, crt)

    plt.plot(x, t, 'k--', label='h=5 numerical')

    t_exact = analytical(x, h, k, w, b, l, t_inf, tb)

    plt.plot(x, t_exact, 'k', label='h=5 exact')

    h = 30

    t, x, base, middle, tip, q_log, t_log = solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h,
                                                         t_inf, t_sur, True, crt)

    plt.plot(x, t, 'k--', label='h=30 numerical')

    t_exact = analytical(x, h, k, w, b, l, t_inf, tb)

    plt.plot(x, t_exact, 'k', label='h=30 exact')

    plt.legend()
    plt.show()


if __name__ == '__main__':
    compare_exact()
