import numpy as np
import math


def assemble_cv_values(x, i, bc_list, gamma, f, scheme):
    """
    This function calculates all the necessary values for a given cv
    :return: Ap, Aw, Ae, b
    """

    dx_ww, dx_we, dx_ew, dx_ee, bc = get_cv_geometry(x, i, bc_list)

    if bc is not None:

        ap = 1
        aw = 0
        ae = 0
        b = bc

    else:

        dw = gamma / (dx_we + dx_ww)
        aw = dw*func_a(scheme, f/dw) + max(f, 0)

        de = gamma / (dx_ee + dx_ew)
        ae = de*func_a(scheme, f/de) + max(-f, 0)

        ap = aw + ae
        b = 0.

    return ap, aw, ae, b


def get_cv_geometry(x, i, bc_list):

    if i == 0:
        bc = bc_list[0]

        bdx = (x[i+1]-x[i])/2

        dx_ww = 0
        dx_we = 0
        dx_ew = bdx
        dx_ee = x[i + 1] - x[i] - bdx

    elif i == len(x)-1:
        bc = bc_list[1]

        bdx = (x[i]-x[i-1])/2

        dx_ww = x[i] - x[i-1] - bdx
        dx_we = bdx
        dx_ew = 0
        dx_ee = 0

    else:
        bc = None

        bdx = (x[i] - x[i - 1])

        dx_ww = x[i] - x[i-1] - bdx / 2
        dx_we = x[i] - x[i-1] - dx_ww
        dx_ee = x[i+1] - x[i] - bdx / 2
        dx_ew = x[i+1] - x[i] - dx_ee

    return [dx_ww, dx_we, dx_ew, dx_ee, bc]


def mat(x, bc_list, gamma, f, scheme):

    a = np.zeros((len(x), len(x)), dtype=float)
    b = np.zeros((len(x), 1), dtype=float)

    # loop over each node
    for i in xrange(len(x)):

        # get cv values
        ap, aw, ae, bp = assemble_cv_values(x, i, bc_list, gamma, f, scheme)

        # assemble into mat
        a[i, i] = ap
        b[i] = bp

        if aw != 0:
            a[i, i - 1] = -aw
        if ae != 0:
            a[i, i + 1] = -ae

    # solve mat equation
    a = np.matrix(a)
    b = np.matrix(b)

    return a.I * b


def func_a(scheme, p):
    if scheme == 'central':
        a = 1. - .5*math.fabs(p)
    elif scheme == 'upwind':
        a = 1.
    elif scheme == 'hybrid':
        a = max(0., 1. - .5*math.fabs(p))
    elif scheme == 'power_law':
        a = max(0., math.pow(1.-.1*math.fabs(p), 5.))
    else:
        raise ValueError("I'm sorry Dave, I can't do that.")
    return a


def solve(x, bc_list, t_init, t_stop, dt, criterion=0.0001):
    """
    This function assembles A and b then solves for T (AT = b)
    :param x: location of the nodes
    :type x: ndarray
    :param bc_list: [boundary condition 1, boundary condition 2]
    :type bc_list: [dict, dict]
    :param t_init: initial temperature
    :type t_init: float
    :param t_stop: final solution time
    :type t_stop: float
    :param dt: time step size
    :type dt: float
    :param criterion: convergence criterion
    :type criterion: float
    :return: iterable, T
    """

    # Set Initials
    t = np.ones(x)*t_init
    t_star = np.copy(t)
    t_old = np.copy(t)

    # Loop over time
    for i in xrange(int(t_stop/dt)):

        # Loop while
        while max(np.absolute(t_star-t)) > criterion:

            # calc k, sc, sp

            # calc ae, aw, ap, b

            # solve t
            t = mat(x, bc_list, gamma, f, scheme)

            # set t _ star
            t_star = np.copy(t)

        t_old = t

    return np.asarray(t)
