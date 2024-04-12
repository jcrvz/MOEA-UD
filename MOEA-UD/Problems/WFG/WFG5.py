"""
WFG5 test problem.

S. Huband, L. Barone, L. While, and P. Hingston, "A Scalable Multi-objective 
Test Problem Toolkit," in Evolutionary Multi-Criterion Optimization, 
pp. 280-295, 2005.

Y. Tian, R. Cheng, X. Zhang, and Y. Jin, "PlatEMO: A MATLAB platform for 
evolutionary multi-objective optimization," in IEEE Computational Intelligence 
Magazine, vol. 12, no. 4, pp. 73-87, 2017.
"""

import numpy as np

def parameters(m):
    """Returns number of decision variables, lower bounds, and upper bounds"""
    k = m-1
    l = 24-(m-1)
    n = k+l
    lb = np.zeros(n)
    ub = np.arange(2.0, 2*n+1, 2)
    return n, lb, ub

def evaluate(P, m):
    """Evaluates a population for the WFG5 test problem"""
    N, n = np.shape(P)
    k = m-1
    D = 1
    A = np.ones(m-1)
    S = np.arange(2.0, 2*m+1, 2)
    
    z01 = P/np.tile(np.arange(2.0, 2*n+1, 2), (N, 1))
    
    t1 = np.zeros((N, n))
    t1 = s_decept(z01, 0.35, 0.001, 0.05)
    
    t2 = np.zeros((N, m))
    for i in range(0, m-1):
        t2[:,i] = r_sum(t1[:,i*k//(m-1):(i+1)*k//(m-1)], np.ones(k//(m-1)))
    t2[:,m-1] = r_sum(t1[:,k:], np.ones(n-k))
    
    x = np.zeros((N, m))
    for i in range(0, m-1):
        x[:,i] = np.maximum(t2[:,m-1], A[i])*(t2[:,i]-0.5)+0.5
    x[:,m-1] = t2[:,m-1]
    
    h = concave(x)
    
    return np.tile(D*x[:,m-1][:,np.newaxis], (1, m))+np.tile(S, (N, 1))*h

def s_decept(y, A, B, C):
    """Transformation function. Shift: Deceptive"""
    out = 1+(np.abs(y-A)-B)*(np.floor(y-A+B)*(1-C+(A-B)/B)/(A-B)+np.floor(A+B-y)*(1-C+(1-A-B)/B)/(1-A-B)+1/B)
    return out

def r_sum(y, w):
    """Transformation function. Reduction: Weighted Sum"""
    out = np.sum(y*np.tile(w, (len(y), 1)), axis=1)/np.sum(w)
    return out

def concave(x):
    """Shape function. Concave"""
    N = len(x)
    out = np.fliplr(np.cumprod(np.hstack((np.ones((N, 1)), np.sin(x[:,:-1]*np.pi/2))), axis=1))*np.hstack((np.ones((N, 1)), np.cos(x[:,-2::-1]*np.pi/2)))
    return out
