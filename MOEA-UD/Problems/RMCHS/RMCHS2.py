import numpy as np
import microchannels as mchs
import materials as mats

def parameters(m):
    """Returns number of decision variables, lower bounds, and upper bounds"""
    n   = 3
    # variables are: w_c, w_w, G_d
    lb  = np.array([1e-8, 1, 1e-12])
    ub  = np.array([1, 1e8, 1e-2])
    return n, lb, ub

def evaluate(P, m):
    """Evaluates a population for the DTLZ1_MINUS test problem"""

    if m != 2:
        raise Exception('RMCHS1 is only defined for 2 objectives.')

    N, n = np.shape(P)

    base_material   = mats.Silicon()
    fluid_material  = mats.Air()

    mh = mchs.Microchannel(base_material, fluid_material)

    g = np.zeros((N, m))

    for ii in range(N):
        w_c = P[ii,0]
        w_w = P[ii,1]
        G_d = P[ii,2]

        mh.w_c = w_c
        mh.w_w = w_w
        mh.G_d = G_d

        # Objective 1: Minimize thermal resistance
        R_eq = mh.R_eq
        g[ii,0] = R_eq

        # Objective 2: Minimize pumping power
        DeltaP_total = mh.DeltaP_total()
        g[ii,1] = DeltaP_total

    return g
