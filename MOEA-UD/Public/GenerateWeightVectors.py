"""
Generate weight vectors.

I. Das and J. E. Dennis, "Normal-Boundary Intersection: A New Method for 
Generating the Pareto Surface in Nonlinear Multicriteria Optimization Problems," 
in SIAM Journal on Optimization, vol. 8, no. 3, pp. 631-657, 1998.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.

R. Y. Rubinstein, D. P. Kroese, "Simulation and the Monte Carlo Method," Third
Edition, John Wiley & Sons, 2016.
"""

import numpy as np
import itertools
import math

def nchoosek(v, k):
    """Returns a matrix with all combinations of v taken k at a time"""
    return np.array(list(itertools.combinations(v, k)))

def uniformWeightVectors(N, m):
    """Generates weight vectors uniformly distributed using the Simplex-Lattice Design"""
    H1 = 1
    while math.comb(H1+m, m-1) <= N:
        H1 += 1
    W = nchoosek(np.arange(1, H1+m), m-1)-np.tile(np.arange(0, m-1), (math.comb(H1+m-1, m-1), 1))-1
    W = (np.hstack((W, np.full((len(W), 1), H1)))-np.hstack((np.zeros((len(W), 1)), W)))/H1
    if H1 < m:
        H2 = 0
        while math.comb(H1+m-1, m-1)+math.comb(H2+m, m-1) <= N:
            H2 += 1
        if H2 > 0:
            W2 = nchoosek(np.arange(1, H2+m), m-1)-np.tile(np.arange(0, m-1), (math.comb(H2+m-1, m-1), 1))-1
            W2 = (np.hstack((W2, np.full((len(W2), 1), H2)))-np.hstack((np.zeros((len(W2), 1)), W2)))/H2
            W2 = W2/2+1/(2*m)
            W = np.vstack((W, W2))
    W[W<1e-6] = 1e-6
    return W

def randomWeightVectors(N, m):
    """Generates weight vectors uniformly distributed at random"""
    U = np.random.rand(N, m-1)
    Y = np.sort(U, axis=1)
    W = np.hstack((Y, np.ones((N, 1))))-np.hstack((np.zeros((N, 1)), Y))
    W[W<1e-6] = 1e-6
    return W
