# This is a 1D CFD Thermal Solver


def assemble_geometry():
    """
    This function handles all geometry set up
    :return: x_nodes, x_bounds
    """
    pass

def assemble_cv_values():
    """
    This function calculates all the necessary values for a given cv
    :return: Ap, Aw, Ae, b
    """
    pass

def solve():
    """
    This function assembles A and b then solves for T (AT = b)
    :return: iterable, T
    """
    pass

if __name__ == '__main__':

    # properties dictionary
    prop1 = {k: 15, q: 4000000}

    # This variable contains a set of properties for any given segment
    properties = [prop1]

    # This variable contains the boundary conditions
    boundaries = []