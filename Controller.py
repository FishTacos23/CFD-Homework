# This is a 1D CFD Thermal Solver using Practice B
import numpy as np
import matplotlib.pyplot as plt
import math


def assemble_geometry(cv_b):
    """
    This function handles all geometry set up
    :return: x_nodes, x_bounds
    """

    # get nodal locations
    x = [cv_b[0]]
    for i in xrange(len(cv_b)-1):
        x.append((cv_b[i]+cv_b[i+1])/2)
    x.append(cv_b[-1])

    return x


def assemble_cv_values(dx_ww, dx_we, dx_ew, dx_ee, bdx, kp, kw, ke, s, bc):
    """
    This function calculates all the necessary values for a given cv
    :return: Ap, Aw, Ae, b
    """

    # west node
    if dx_ww == 0 and dx_we == 0:
        aw = 0
    else:
        khw = (dx_ww + dx_we) * kw * kp / (kw * dx_we + kp * dx_ww)
        aw = khw / (dx_ww + dx_we)

    # east node
    if dx_ew == 0 and dx_ee == 0:
        ae = 0
    else:
        khe = (dx_ew + dx_ee) * ke * kp / (ke * dx_ew + kp * dx_ee)
        ae = khe / (dx_ew + dx_ee)

    ap = aw + ae
    b = s * bdx

    if bc is not None:
        if bc['type'] == 'convective':
            ap += bc['h']
            b += bc['h']*bc['To']

    return ap, aw, ae, b


def get_cv_geometry(x, cv_b, i, p, p_map, bc_list):

    kp = p[p_map[i]]['k'](x[i])
    s = p[p_map[i]]['q']

    if i == 0:
        bc = bc_list[0]

        bdx = 0

        dx_ww = 0
        dx_we = 0
        dx_ew = 0
        dx_ee = x[i + 1] - x[i] - bdx / 2

        kw = 0
        ke = p[p_map[i + 1]]['k'](x[i])

    elif i == len(x)-1:
        bc = bc_list[1]

        bdx = 0

        dx_ww = x[i] - x[i-1] - bdx / 2
        dx_we = 0
        dx_ew = 0
        dx_ee = 0

        kw = p[p_map[i - 1]]['k'](x[i])
        ke = 0

    else:
        bc = None

        bdx = (cv_b[i] - cv_b[i - 1])

        dx_ww = x[i] - x[i-1] - bdx / 2
        dx_we = x[i] - x[i-1] - dx_ww
        dx_ee = x[i+1] - x[i] - bdx / 2
        dx_ew = x[i+1] - x[i] - dx_ee

        kw = p[p_map[i-1]]['k'](x[i])
        ke = p[p_map[i+1]]['k'](x[i])

    return [dx_ww, dx_we, dx_ew, dx_ee, bdx, kp, kw, ke, s, bc]


def solve(cv_b, bc_list, p, p_map):
    """
    This function assembles A and b then solves for T (AT = b)
    :return: iterable, T
    """

    # get geometry
    x = assemble_geometry(cv_b)

    a = np.zeros((len(x), len(x)), dtype=float)
    b = np.zeros((len(x), 1), dtype=float)

    # loop over each node
    for i in xrange(len(x)):

        # get cv values
        cv_input = get_cv_geometry(x, cv_b, i, p, p_map, bc_list)
        ap, aw, ae, bp = assemble_cv_values(*cv_input)

        # assemble into mat
        a[i, i] = ap
        b[i] = bp

        if aw != 0:
            a[i, i-1] = -aw
        if ae != 0:
            a[i, i+1] = -ae

    # solve mat equation
    a = np.matrix(a)
    b = np.matrix(b)

    t = a.I*b

    return t, x


if __name__ == '__main__':

    # properties dictionary
    prop1 = {'k': lambda x: 15, 'q': 4000000}
    # prop2 = {'k': lambda x: 60, 'q': 0}
    prop2 = {'k': lambda x: 137 * math.exp(25 * x - 2), 'q': 0}

    # This variable contains a set of properties for any given segment
    properties = [prop1, prop2]

    # user defined boundary locations
    # cv_boundaries = [0, .01, .02, .03, .04, .05]
    # node_properties = [0, 0, 0, 0, 1, 1, 1]             # map cv to properties list

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

    py, px = solve(cv_boundaries, boundaries, properties, node_properties)

    plt.plot(px, py)
    plt.show()
