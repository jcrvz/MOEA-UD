import numpy as np
from modelMCHS import modelMCHS
import matplotlib.pyplot as plt

x_lb    = np.array([1e-3, 1.0, 1e-8])
x_ub    = np.array([1.0, 1e3, 5e-5])

q_lb    = np.array([1.0])
q_ub    = np.array([100.0])

Materials = ['Cu', 'Al']
Fluids   = ['H2O', 'Air']


def evaluate_models(design_var, power, material, fluid):
    # Evaluate the Python model
    soo_py, (sgenht_py, sgenff_py), outpar_py = modelMCHS(
        design_var=design_var, material=material, fluid=fluid, power=power)

    # Here you would call the MATLAB model, for example using MATLAB Engine API for Python
    # soo_mat, (sgenht_mat, sgenff_mat), outpar_mat = modelMCHS_matlab(
    #     design_var=design_var, material=material, fluid=fluid, power=power)

    return soo_py


if __name__ == '__main__':
    x1 = 1e-2
    x2 = 10
    x3 = np.linspace(x_lb[2], x_ub[2], 10000)
    q = 15.0

    plt.figure()

    for q in [10.0, 50.0, 100.0]:
        soo_vals = np.zeros(x3.shape)
        for i in range(len(x3)):
            design_var = np.array([x1, x2, x3[i]])
            soo_vals[i] = evaluate_models(
                design_var, q, Materials[0], Fluids[0])

        plt.plot(x3, soo_vals, label=f'q={q:.2f} W')

    plt.xscale('log')
    plt.xlabel(r'$x_3$ (m$^3$/s)')
    plt.ylabel(r'$S_{gen}$ (dimensionless)')
    plt.title('Soo vs x3 for fixed x1 and x2')
    plt.legend()
    plt.show()

