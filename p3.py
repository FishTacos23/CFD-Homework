import scipy
import numpy as np
import math
import matplotlib.pyplot as plt


def solve():

    k_a = 15
    k_b = lambda x_i: 60
    q_a = 4000000
    h = 1000
    T_o = 300
    dx = 0.01
    bdx = 0.01
    x = []

    a_mat = np.zeros((7, 7), dtype=float)
    b_vec = np.zeros((7, 1), dtype=float)

    # Node 1 - LB
    x.append(0)
    a_mat[0, 0] = 2*k_a/dx
    a_mat[0, 1] = -2*k_a/dx

    # Node 2
    x.append(x[-1] + dx / 2)
    a_mat[1, 0] = -2*k_a/dx
    a_mat[1, 1] = k_a / dx + 2*k_a/dx
    a_mat[1, 2] = -k_a/dx
    b_vec[1] = q_a*bdx

    # Node 3
    x.append(x[-1] + dx)
    a_mat[2, 1] = -k_a / dx
    a_mat[2, 2] = k_a/dx + k_a/dx
    a_mat[2, 3] = -k_a/dx
    b_vec[2] = q_a * bdx

    # Node 4 - West of Surface
    x.append(x[-1] + dx)
    k_h = math.pow(((.5/k_a)+(.5/k_b(x[-1]+dx))), -1)
    a_mat[3, 2] = -k_a / dx
    a_mat[3, 3] = k_h/dx + k_a/dx
    a_mat[3, 4] = -k_h/dx
    b_vec[3] = q_a * bdx

    # Node 5 - East of Surface
    x.append(x[-1] + dx)
    k_h = math.pow(((.5 / k_a) + (.5 / k_b(x[-1]))), -1)
    a_mat[4, 3] = -k_h / dx
    k_h = math.pow(((.5 / k_b(x[-1])) + (.5 / k_b(x[-1]+dx))), -1)
    a_mat[4, 5] = -k_h/dx
    a_mat[4, 4] = -a_mat[4, 3] - a_mat[4, 5]

    # Node 6
    x.append(x[-1] + dx)
    k_h = math.pow(((.5 / k_b(x[-1]-dx)) + (.5 / k_b(x[-1]))), -1)
    a_mat[5, 4] = -k_h / dx
    k_h = dx*k_b(x[-1]+dx/2)*k_b(x[-1])/(k_b(x[-1]+dx/2)*dx/2)
    a_mat[5, 6] = -2 * k_h / dx
    a_mat[5, 5] = -a_mat[5, 4] - a_mat[5, 6]

    # Node 7 - RB
    x.append(x[-1] + dx/2)
    k_h = dx * k_b(x[-1] - dx / 2) * k_b(x[-1]) / (k_b(x[-1]) * dx / 2)
    a_mat[6, 5] = -2 * k_h / dx
    a_mat[6, 6] = h - a_mat[6,5]
    b_vec[6] = h*T_o

    a_mat = np.matrix(a_mat)
    b_vec = np.matrix(b_vec)

    t = a_mat.I*b_vec

    t = np.asarray(t)

    print t

    plt.plot(x, t)
    plt.show()


if __name__ == "__main__":
    solve()
