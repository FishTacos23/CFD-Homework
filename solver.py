import numpy as np
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
        if bc['type'] == 'fixed':
            ap = 1
            aw = 0
            ae = 0
            b = bc['T']

    return ap, aw, ae, b


def get_cv_geometry(x, cv_b, i, p, p_map, bc_list, t_prev):

    kp = p[p_map[i]]['k'](x[i])

    if t_prev is not None:
        boltz = 0.0000000567
        ts = t_prev[i]
        tinf = p[p_map[i]]['T']

        sp = (p[p_map[i]]['h'] + 4. * boltz * p[p_map[i]]['e'] * math.pow(ts, 3.)) * 4. / p[p_map[i]]['D']

        sc = (4./p[p_map[i]]['D']) * (p[p_map[i]]['h']*(ts-tinf) +
                                      boltz*p[p_map[i]]['e']*(math.pow(ts, 4.) - math.pow(tinf, 4.))) - sp*ts

        s = -(sc + sp*ts)
    else:
        s = p[p_map[i]]['q']

    if i == 0:
        bc = bc_list[0]

        bdx = 0

        dx_ww = 0
        dx_we = 0
        dx_ew = 0
        dx_ee = x[i + 1] - x[i] - bdx / 2

        kw = 0
        ke = p[p_map[i + 1]]['k'](x[i+1])

    elif i == len(x)-1:
        bc = bc_list[1]

        bdx = 0

        dx_ww = x[i] - x[i-1] - bdx / 2
        dx_we = 0
        dx_ew = 0
        dx_ee = 0

        kw = p[p_map[i - 1]]['k'](x[i-1])
        ke = 0

    else:
        bc = None

        bdx = (cv_b[i] - cv_b[i - 1])

        dx_ww = x[i] - x[i-1] - bdx / 2
        dx_we = x[i] - x[i-1] - dx_ww
        dx_ee = x[i+1] - x[i] - bdx / 2
        dx_ew = x[i+1] - x[i] - dx_ee

        kw = p[p_map[i-1]]['k'](x[i-1])
        ke = p[p_map[i+1]]['k'](x[i+1])

    return [dx_ww, dx_we, dx_ew, dx_ee, bdx, kp, kw, ke, s, bc]


def solve(cv_b, bc_list, p, p_map, t_prev):
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
        cv_input = get_cv_geometry(x, cv_b, i, p, p_map, bc_list, t_prev)
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

    return np.squeeze(np.asarray(t)), x
