"""
MOEA-UD.
"""

import numpy as np
from scipy.spatial import distance

from Public.Population import population
from Public.UploadTestProblem import uploadTestProblem
from Public.GenerateWeightVectors import uniformWeightVectors, randomWeightVectors
from Public.RandomPopulation import randomPopulation
from Public.GenerateOffspring import generateOffspring
from Public.EfficientNonDominatedSort import efficientNonDominatedSort

def main(N, problem, m, max_generations):
    """Runs main framework of MOEA-UD"""
    parameters, evaluate = uploadTestProblem(problem)
    n, lb, ub = parameters(m)
    W = uniformWeightVectors(N, m)
    pc, nc = 1, 30
    pm, nm = 1/n, 20
    P = randomPopulation(N, n, m, lb, ub, evaluate)
    zmin = np.min(P.obj, axis=0)
    zrange = np.ones(m)
    generations = 0
    
    Fronts = efficientNonDominatedSort(P.obj)
    A = P.obj[Fronts[0]]
    Wu = np.copy(W)
    
    if len(Wu) < N:
        Wd = randomWeightVectors(N-len(Wu), m)
    else:
        Wd = np.empty((0, m))
    
    Data_gen = []
    
    while generations < max_generations:
        M = matingSelection(P, N)
        Q = generateOffspring(M, N, m, lb, ub, pc, nc, pm, nm, evaluate)
        zmin = np.min(np.vstack((zmin, Q.obj)), axis=0)
        R = population(np.vstack((P.dec, Q.dec)), np.vstack((P.obj, Q.obj)))
        P, zmax = survivalSelection(R, Wu, Wd, N, zmin, zrange)
        generations += 1
        
        if generations <= 0.9*max_generations:
            A, Appf, Wactive = archiveMaintenance(A, Q, W, N, zmin, zmax)
            if np.floor(generations % (0.05*max_generations)) == 0:
                Wu, Wd = adaptReferenceSets(A, Appf, Wactive, W, N, zmin, zmax)
                zrange = zmax-zmin
                zrange[zrange == 0] = 1e-12
        
        Data_gen.append(P.obj)
        
    return P, Data_gen

def matingSelection(P, N):
    """Selects random parent population"""
    O = N+1 if N%2 == 1 else N
    if len(P.dec) > 1:
        indexes = np.array([np.random.choice(len(P.dec), 2, replace=False) for i in range(0, O//2)])
    else:
        indexes = np.zeros((O//2, 2), dtype=int)
    M = P.dec[np.hstack((indexes[:,0], indexes[:,1]))]
    return M

def survivalSelection(R, Wu, Wd, N, zmin, zrange):
    """Returns population with best individuals"""
    unique = np.sort(np.unique(np.around(R.obj, 6), return_index=True, axis=0)[1])
    R.dec = R.dec[unique]
    R.obj = R.obj[unique]
    Fronts = efficientNonDominatedSort(R.obj)
    zmax = np.max(R.obj[Fronts[0]], axis=0)
    P, Fl = findCriticalFront(R, Fronts, N)
    if len(P.dec) < N:
        if len(P.dec) == 0:
            selected = nichingSelection(Fl.obj, Wu, zmin, zrange)[0]
            if len(Wd) > 0:
                candidates = np.setdiff1d(np.arange(0, len(Fl.obj)), selected)
                selected_d = nichingSelection(Fl.obj[candidates], Wd, zmin, zrange)[0]
                selected = np.hstack((selected, candidates[selected_d]))
            
            P = population(Fl.dec[selected], Fl.obj[selected])
            zmax = np.max(P.obj, axis=0)
            
            if len(P.dec) < N:
                candidates = np.setdiff1d(np.arange(0, len(Fl.dec)), selected)
                pruned = candidates[np.all(Fl.obj[candidates] <= zmax, axis=1)]
                if len(P.dec)+len(pruned) >= N and np.all(zmax-zmin > 1e-3):
                    C = population(Fl.dec[pruned], Fl.obj[pruned])
                    selected = subsetSelectionPPFs(P.obj, C.obj, N, zmin)
                    P.dec = np.vstack((P.dec, C.dec[selected]))
                    P.obj = np.vstack((P.obj, C.obj[selected]))
                else:
                    C = population(Fl.dec[candidates], Fl.obj[candidates])
                    selected = subsetSelectionPPFs(P.obj, C.obj, N, zmin)
                    P.dec = np.vstack((P.dec, C.dec[selected]))
                    P.obj = np.vstack((P.obj, C.obj[selected]))
                    zmax = np.max(P.obj, axis=0)
        else:
            selected = subsetSelectionPPFs(P.obj, Fl.obj, N, zmin)
            P.dec = np.vstack((P.dec, Fl.dec[selected]))
            P.obj = np.vstack((P.obj, Fl.obj[selected]))
    return P, zmax

def findCriticalFront(R, Fronts, N):
    """Returns population in best fronts and population in critical front"""
    P = population([], [])
    for front in Fronts:
        if len(P.dec)+len(front) > N:
            Fl = population(R.dec[front], R.obj[front])
            break
        else:
            P.dec = R.dec[front] if len(P.dec) == 0 else np.vstack((P.dec, R.dec[front]))
            P.obj = R.obj[front] if len(P.obj) == 0 else np.vstack((P.obj, R.obj[front]))
    return P, Fl

def nichingSelection(A, W, zmin, zrange):
    """Divides a given population into niches and returns best solution per niche and weight vectors representing each niche"""
    Aprime = (A-zmin)/zrange
    Dper = perpendicularDistance(Aprime, W)
    Daasf = augmentedAchievementScalarizingFunction(Aprime, W)
    pi = np.argmin(Dper, axis=1)
    d = Daasf[np.arange(len(A)),pi]
    selected = []
    J = np.unique(pi)
    for j in J:
        I = np.where(pi == j)[0]
        s = I[np.argmin(d[I])]
        selected.append(s)
    return selected, W[J]

def perpendicularDistance(A, W):
    """Returns fitness matrix using perpendicular distance"""
    U = np.tile(W, (len(A), 1))
    V = np.repeat(A, len(W), axis=0)
    Unorm = np.linalg.norm(U, axis=1)
    Proj_scalar = np.sum(V*U, axis=1)/Unorm
    Proj = Proj_scalar[:,None]*U/Unorm[:,None]
    return np.reshape(np.linalg.norm(Proj-V, axis=1), (len(A), len(W)))

def augmentedAchievementScalarizingFunction(A, W):
    """Returns fitness matrix using the augmented achievement scalarizing function"""
    Aext = np.repeat(A, len(W), axis=0)
    Wext = np.tile(W, (len(A), 1))
    Wext[Wext < 1e-6] = 1e-6
    return np.reshape(np.max(Aext/Wext, axis=1)+1e-4*np.sum(Aext/Wext, axis=1), (len(A), len(W)))

def archiveMaintenance(A, Q, W, N, zmin, zmax):
    """Updates archive based on niching and a given pair-potential energy function"""
    A = orderedNondominanceFilter(A, Q)
    unique = np.sort(np.unique(np.around(A, 6), return_index=True, axis=0)[1])
    A = A[unique]
    Appf = []
    Wactive = []
    if len(A) >= N:
        zrange = zmax-zmin
        zrange[zrange == 0] = 1e-12
        selected, Wactive = nichingSelection(A, W, zmin, zrange)
        Aniches = A[selected]
        if len(Aniches) == N:
            A = np.copy(Aniches)
        else:
            candidates = np.setdiff1d(np.arange(0, len(A)), selected)
            C = A[candidates][np.all(A[candidates] <= zmax, axis=1)]
            if len(Aniches)+len(C) >= N:
                selected = subsetSelectionPPFs(Aniches, C, N, zmin)
                Appf = C[selected]
                A = np.vstack((Aniches, Appf))
            else:
                A = np.vstack((Aniches, C))
    return A, Appf, Wactive

def orderedNondominanceFilter(A, Q):
    """Returns ordered non-dominated archive"""
    Fronts = efficientNonDominatedSort(Q.obj)
    O = Q.obj[Fronts[0]]
    Asurviving = np.ones(len(A), dtype='bool')
    Osurviving = np.ones(len(O), dtype='bool')
    Oelite = np.zeros(len(O), dtype='bool')
    for i in range(0, len(O)):
        for j in range(0, len(A)):
            if Asurviving[j]:
                if all(A[j] <= O[i]) and any(A[j] < O[i]):
                    Osurviving[i] = False
                    break
                elif all(O[i] <= A[j]) and any(O[i] < A[j]):
                    Asurviving[j] = False
                    Osurviving[i] = False
                    Oelite[i] = True
    return np.vstack((A[Asurviving], O[Oelite], O[Osurviving]))

def subsetSelectionPPFs(S, C, N, zmin):
    """Returns the selected subset using a given pair-potential function"""
    A = np.vstack((S, C))
    if len(A) <= N:
        selected = np.arange(0, len(C))
    else:
        zmax_arch = np.max(A, axis=0)
        denom = zmax_arch-zmin
        denom[denom == 0] = 1e-12
        Anorm = (A-zmin)/denom
        Aproj = projectionOnHyperplane(Anorm)
        Diss = dissimilarityMatrix(Aproj)
        selected = np.arange(0, N)
        candidates = np.arange(N, len(A))
        Memo = np.sum(Diss[:,selected][selected,:], axis=1)
        for candidate in candidates:
            Memo = Memo+Diss[selected,candidate]
            Memo = np.append(Memo, np.sum(Diss[candidate,selected]))
            selected = np.append(selected, candidate)
            I = np.arange(len(S), len(Memo))
            worst = I[np.argmax(Memo[I])]
            Memo = Memo-Diss[selected,selected[worst]]
            Memo = np.delete(Memo, worst)
            selected = np.delete(selected, worst)
        selected = selected[len(S):]-len(S)
    return selected

def projectionOnHyperplane(A):
    """Projects the given approximation set into the hyperplane using normal-boundary intersection"""
    m = np.shape(A)[1]
    s = np.sum(A, axis=1, keepdims=True)
    return A+np.tile((1-s)/m, (1, m))

def dissimilarityMatrix(A):
    """Returns dissimilarity matrix using Riesz s-energy"""
    m = np.shape(A)[1]
    s = 0.417531*m+4.655076
    d = distance.pdist(A, 'euclidean')
    denom = d**s
    denom[denom == 0] = 1e-12
    d = 1/denom
    return distance.squareform(d)

def adaptReferenceSets(A, Appf, Wactive, W, N, zmin, zmax):
    """Returns adapted reference sets of weight vectors for uniformity and diversity"""
    zmax_arch = np.max(A, axis=0)
    if len(A) == N and np.all(zmax_arch-zmin > 1e-3):
        Wu = np.copy(Wactive)
        if len(Appf) > 0:
            denom = zmax-zmin
            denom[denom == 0] = 1e-12
            V = (Appf-zmin)/denom
            denom = np.sum(V, axis=1)
            denom[denom == 0] = 1e-12
            Wd = V/denom[:,np.newaxis]
            Wd[Wd<1e-6] = 1e-6
        else:
            m = np.shape(W)[1]
            Wd = np.empty((0, m))
    else:
        Wu = np.copy(W)
        m = np.shape(W)[1]
        if len(Wu) < N:
            Wd = randomWeightVectors(N-len(Wu), m)
        else:
            Wd = np.empty((0, m))
    return Wu, Wd
