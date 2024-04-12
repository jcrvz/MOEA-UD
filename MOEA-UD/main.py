import sys
import numpy as np

from Public.ValidateArguments import validateArguments
from Public.ValidateObjectives import validateObjectives
from Public.UploadAlgorithm import uploadAlgorithm
from Public.SaveApproximationSet import saveApproximationSet
from Public.SaveGenerations import saveGenerations

if __name__ == '__main__':
    if (str(sys.argv[1]) == '--help'):
        f = open('../README.txt', 'r')
        contents = f.read()
        f.close()
        print(contents)
    else:
        if (len(sys.argv) != 7):
            sys.exit('Incorrect number of arguments. Use: main.py --help')
        algorithm = str(sys.argv[1])
        N = int(sys.argv[2])
        problem = str(sys.argv[3])
        m = int(sys.argv[4])
        max_generations = int(sys.argv[5])
        runs = int(sys.argv[6])
        
        validateArguments(algorithm, N, problem, m, max_generations, runs)
        m = validateObjectives(problem, m)
        
        seeds = np.random.choice(10000000, runs, replace=False)
        
        for run in range(1, runs+1):
            seed = seeds[run-1]
            np.random.seed(seed)
            print('Algorithm:', algorithm, '| Population size:', N, 
                  '| Problem:', problem, '| Objectives:', m, '| Generations:', 
                  max_generations, '| Run:', run)
            if algorithm == 'MOEA-UD':
                main = uploadAlgorithm(algorithm)
                P, Data_gen = main(N, problem, m, max_generations)
            saveApproximationSet(P.obj, algorithm, problem, run, seed, 'save_txt')
            saveGenerations(Data_gen, algorithm, problem, run, seed)
            del P, Data_gen, seed
