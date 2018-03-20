# This is a 1D CFD Thermal Solver using Practice A
import matplotlib.pyplot as plt
import math
import numpy as np
import solver


def prob_2_4():

    """
    :return: Prints a table for values of phi, and error for different schemes and velocities, and plots of phi
    """

    schemes = ['central', 'upwind', 'hybrid', 'power_law']
    markers = ['o', '*', '+', 'v']

    hs = [5., 30.]
    eps = [0., 1.]

    boundaries = [1., 0.]

    # FIXED PROPERTIES
    rho = 2770.
    c = 875.
    k = 177.
    length = .01
    t_inf = 293.
    t_surr = 293.
    tb = 353.
    w = 1.
    b = 0.0003

    # num grid
    num_grid = 50

    x = np.linspace(0., length, num_grid)

    t_init = 353.
    t = np.ones(x.shape)*t_init

    for h in hs:  # LOOP OVER EACH H
        print h

        for ep in eps:  # LOOP OVER EACH EPSILON
            print ep

            fig, lin = create_plotter(x, t)

            t_star = np.copy(t)
            t_old = np.copy(t)

        #     phi_exact = analytical(x, rho, h, length, gamma)
        #     plt.figure()
        #     plt.plot(x, phi_exact)
        #     for i in xrange(len(schemes)):
        #         print schemes[i]
        #         f = rho * h
        #         phi_nm = solver.solve(x, boundaries, schemes[i], f, gamma)
        #         plt.plot(x, phi_nm, markers[i])
        #         err = 0
        #         for j in xrange(phi_nm.shape[0]):
        #             err += math.fabs(phi_nm[j]-phi_exact[j])
        #         print 'Phi Values:'
        #         print phi_nm
        #         print 'error = ' + str(err)
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


def update_plots(fig, line, new_data):
    line.set_ydata(new_data)
    fig.canvas.draw()


def create_plotter(x, t):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line, = ax.plot(x, t)
    return fig, line


if __name__ == '__main__':
    prob_2_4()
