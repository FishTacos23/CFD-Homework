# This is a 1D CFD Thermal Solver using Practice A
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def prob_1_3():

    """
    :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
    """

    schemes = ['central', 'upwind', 'hybrid', 'power_law']
    markers = ['o', '*', '+', 'v']
    us = [1., 20., 75.]

    boundaries = [1., 0.]
    rho = 1.
    gamma = 1.
    length = 1.
    x = np.linspace(0., length, 21)

    for u in us:
        print u
        phi_exact = analytical(x, rho, u, length, gamma)
        plt.figure()
        plt.plot(x, phi_exact)
        for i in xrange(len(schemes)):
            print schemes[i]
            f = rho / u
            phi_nm = solver.solve(x, boundaries, schemes[i], f, gamma)
            plt.plot(x, phi_nm, markers[i])
            err = 0
            for i in xrange(phi_nm.shape[0]):
                err += math.fabs(phi_nm[i]-phi_exact[i])
            print 'Phi Values:'
            print phi_nm
            print 'error = ' + str(err)
    plt.show()


def analytical(x, rho, u, l, gamma):
    """

    :param x: Position of nodes
    :type x: ndarray
    :param rho: Density
    :type rho: float
    :param u: Velocity
    :type u: float
    :param l: Domain Length
    :type l: float
    :param gamma:
    :type gamma: float
    :return: temperatures for nodes using the analytical solution
    """

    p = rho * u * l / gamma

    soln = 1. - (np.exp(p * x / l) - 1.) / (math.exp(p) - 1.)

    return soln


if __name__ == '__main__':
    prob_1_3()
