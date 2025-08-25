import numpy as np
from .modelMCHS import modelMCHS

def parameters(m):
    """Returns number of decision variables, lower bounds, and upper bounds"""
    n   = 3
    # variables are: w_c, w_w, G_d
    lb  = np.array([1e-8, 1, 1e-12])
    ub  = np.array([1, 1e3, 1e-2])
    return n, lb, ub

def evaluate(P, m):
    """Evaluates a population for the DTLZ1_MINUS test problem"""

    if m != 2:
        raise Exception('RMCHS1 is only defined for 2 objectives.')

    N, n = np.shape(P)
    material = 'Cu'
    fluid = 'H2O'
    Q = 15.0  # W

    g = np.zeros((N, m))

    for ii in range(N):
        # design_var = [alpha, beta, G_d] = P[ii, :]

        # alpha = P[ii,0]
        # beta  = P[ii,1]
        # G_d   = P[ii,2]

        soo, (sgenht, sgenff), outpar = modelMCHS(
            design_var=P[ii,:], material=material, fluid=fluid, power=Q)

        # Objective 1: Minimize thermal resistance
        g[ii,0] = sgenht

        # Objective 2: Minimize pumping power
        g[ii,1] = sgenff

    return g
