# This is a 1D CFD Thermal Solver using Practice A
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def prob_2_4():

    """
    :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
    """

    hs = [5., 30.]
    eps = [0., 1.]

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
    t_init = 353.
    crt = 0.0001

    # num grid
    dx = l/50
    t_stop = 10
    dt = 0.1

    for h in hs:  # LOOP OVER EACH H
        print h

        for e in eps:  # LOOP OVER EACH EPSILON
            print e

            # DO GRID CONVERGENCE HERE

            # DO TIME STEP CONVERGENCE HERE

            solver.solve(l, dx, bcs, t_init, t_stop, dt, k, c, rho, tb, b, e, h, t_inf, t_sur, True, 0.0001)

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


if __name__ == '__main__':
    prob_2_4()
