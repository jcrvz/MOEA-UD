
"""
Python translation of MATLAB function modelMCHS.m

Original signature:
    function [soo, moo] = modelMCHS(design_var, material, fluid, Q)

Returns:
    soo : float  # dimensionless entropy generation number S_gen * T_a / Q
    moo : tuple  # (R_eq, DeltaP_total)

Note:
- This is a line-by-line faithful port, preserving the same parameter names and formulas.
- Units follow the original comments (SI).
"""

from math import pi, sqrt, tanh, log
import numpy as np

def modelMCHS(design_var, material='Al', fluid='H2O', Q=15.0):
    # ---- Design variables ----
    alpha_c = float(design_var[0])
    beta    = float(design_var[1])
    G_d     = float(design_var[2])

    # ---- Variable parameters ----
    # Channel height (m)
    H_c = 1.7e-3

    # Half channel width (m)
    w_c = H_c * alpha_c / 2.0

    # Half fin width (m)
    w_p = w_c / beta

    # ---- Material properties ----
    if material == 'Al':       # Aluminium
        k_m   = 237.0
        rho_m = 2707.0
    elif material == 'Cu':     # Copper
        k_m   = 401.0
        rho_m = 8954.0
    elif material == 'SiC':    # Silicon carbide
        k_m   = 270.0
        rho_m = 3300.0
    elif material == 'HTCG':   # High thermal conductivity graphite (per original labels)
        k_m   = 1900.0
        rho_m = 1000.0
    elif material == 'kb1k':   # HTCG variant
        k_m   = 1000.0
        rho_m = 1000.0
    elif material == 'kb2k':   # HTCG variant
        k_m   = 2000.0
        rho_m = 1000.0
    elif material == 'Si':     # Silicon
        k_m   = 148.0
        rho_m = 2330.0
    else:
        raise ValueError('Material no definido')

    # ---- Fluid properties ----
    # We compute k_f, rho_f, nu (kinematic viscosity), c_p according to cases.
    if fluid == 'Air':
        k_f  = 0.0261
        rho_f = 1.1614
        nu    = 1.58e-5   # kinematic viscosity directly given (m^2/s)
        c_p   = 1007.0
    elif fluid == 'Air+10HR':
        k_f  = 0.02671
        rho_f = 1.17434
        nu    = 1.80045e-5 / rho_f  # dynamic? original writes nu = μ/ρ
        c_p   = 1008.45
    elif fluid == 'Air+50HR':
        k_f  = 0.028067
        rho_f = 1.16436
        nu    = 1.78858e-5 / rho_f
        c_p   = 1025.1
    elif fluid == 'Air+90HR':
        k_f  = 0.028583
        rho_f = 1.16182
        nu    = 1.83254e-5 / rho_f
        c_p   = 1034.8
    elif fluid == 'H2O':
        k_f  = 0.625
        rho_f = 994.2
        nu    = 7.25e-4 / rho_f  # μ/ρ
        c_p   = 4178.0

    # Nanofluids (water + nanoparticles)
    elif fluid == 'H2O+TiO2_01wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 8.4
        rho_np = 4250.0
        c_p_np = 710.0
        phi_nf = 0.01
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+TiO2_05wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 8.4
        rho_np = 4250.0
        c_p_np = 710.0
        phi_nf = 0.05
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+TiO2_09wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 8.4
        rho_np = 4250.0
        c_p_np = 710.0
        phi_nf = 0.09
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_01wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.01
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_05wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.05
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_09wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.09
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_1wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.01  # 1 wt% ~ 0.01 (as in the labels)
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_5wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.05
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    elif fluid == 'H2O+Al_9wt':
        k_fb  = 0.625
        rho_fb = 994.2
        mu_fb  = 7.25e-4
        c_p_fb = 4178.0
        k_np   = 237.0
        rho_np = 2700.0
        c_p_np = 897.0
        phi_nf = 0.09
        n      = 3.0
        k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np)) * k_fb / (k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np))
        rho_f = (1 + phi_nf)*rho_fb + phi_nf*rho_np
        mu_nf = mu_fb / (1 - phi_nf)**2.5
        nu    = mu_nf / rho_f
        c_p   = ((1 - phi_nf)*rho_fb*c_p_fb + phi_nf*rho_np*c_p_np) / rho_f

    else:
        raise ValueError('Fluido no definido')

    # Prandtl number
    Pr = nu * rho_f * c_p / k_f

    # ---- Fixed geometry & interface parameters ----
    H_b = 0.1e-3
    W_d = 51e-3
    L_d = 51e-3
    W_i = W_d
    L_i = L_d
    # q  = 5e4  # not used in original
    T_a = 300.0
    R_iA = 2.75e-4

    # ---- Geometry parameters ----
    N = (W_d/2.0 - w_p) / (w_c + w_p)

    A_c_ = 2.0 * w_c * H_c
    P_   = 2.0 * (H_c + 2.0*w_c)
    D_h  = 4.0 * A_c_ / P_

    # ---- Flow & heat parameters ----
    mdot  = rho_f * G_d / (2.0 * N)
    U_avg = mdot / (rho_f * w_c * H_c)
    Re_Dh = D_h * U_avg / nu

    A_b = W_d * L_d
    A_i = W_i * L_i

    # ---- Nusselt number & friction factor (laminar correlations) ----
    Nu_Dh = 2.253 + 8.164 * (1.0/(alpha_c + 1.0))**1.5
    f = np.sqrt( (3.2*(Re_Dh*D_h/L_d)**0.57)**2 + (4.7 + 19.64*(alpha_c**2 + 1.0)/(alpha_c + 1.0)**2)**2 ) / Re_Dh

    # ---- Thermal resistances ----
    R_i = R_iA / A_i
    h_avg = k_f * Nu_Dh / D_h

    mHc = np.sqrt( 2.0*(2.0*w_p + L_d)*h_avg / (k_m*2.0*w_p*L_d) ) * H_c
    eta_p = np.tanh(mHc) / mHc if mHc != 0 else 1.0

    A_eff = 2.0 * N * (eta_p*H_c + w_c) * L_d

    R_conv  = 1.0 / (h_avg * A_eff)
    R_f     = 1.0 / (rho_f * G_d * c_p)
    R_const = ((1.0 + beta) / (pi * k_m * A_b * beta)) * np.log(1.0 / np.sin(0.5*pi/(1.0 + beta))) * H_c * alpha_c

    a = np.sqrt(A_i/pi)
    b = np.sqrt(A_b/pi)

    tau     = H_b / b
    epsilon = a / b

    R_o = R_conv + R_f   # + R_const (left out in original R_o definition)

    Bi = 1.0 / (pi * k_m * b * R_o)

    lambda_c = pi + 1.0 / (np.sqrt(pi) * epsilon)
    Phi_c = (np.tanh(lambda_c*tau) + lambda_c/Bi) / (1.0 + lambda_c*np.tanh(lambda_c*tau)/Bi)
    Psi_avg = epsilon*tau/np.sqrt(pi) + 0.5*Phi_c*(1.0 - epsilon)**1.5

    R_b = Psi_avg / (np.sqrt(pi) * k_m * a)

    R_eq = R_o   # original uses R_o only (R_b and R_i commented out)

    # ---- Interface temperature ----
    T_i = T_a + R_eq * Q

    # ---- Pressure drop ----
    k_ce   = 1.79 - 2.32*(beta/(1.0 + beta)) + 0.53*(beta/(1.0 + beta))**2
    DeltaP = 0.5 * rho_f * U_avg**2 * (f * (L_d / D_h) + k_ce)

    # Tube losses (defined but excluded from total in original)
    L_tu = 0.5
    D_tu = 1.9e-2
    U_avgtub = 4.0 * G_d / (pi * D_tu**2)
    f_tu = 0.317 / (Re_Dh**0.25) if Re_Dh != 0 else 0.0  # note: placeholder; original ties to Re_Dh
    A_tuhs = pi * D_tu**2 / (2.0 * N * 2.0 * w_c * H_c)
    DeltaP_tu = 0.5 * rho_f * U_avgtub**2 * (0.42 + (1 - A_tuhs**2)**2 + 0.42*(1 - A_tuhs**2) + 1 + 2*f_tu*L_tu/D_tu)

    DeltaP_total = DeltaP  # + DeltaP_tu (omitted per original)

    Phi = G_d * DeltaP_total

    # ---- Entropy generation rates ----
    S_genht = Q**2 * R_eq / (T_a * T_i)
    S_genff = G_d * DeltaP_total / T_a
    S_gen   = S_genht + S_genff

    # ---- Heat sink mass ----
    Md = rho_m * L_d * (2.0*W_d*H_b + 2.0*w_p*H_c*(N + 1.0))

    # ---- Outputs ----
    soo = S_gen * T_a / Q
    moo = (R_eq, float(DeltaP_total))

    # Additional outputs (for debugging / analysis parity)
    outpar = dict(S_gen=S_gen, S_genht=S_genht, S_genff=S_genff, U_avg=U_avg, Phi=Phi, Md=Md,
                  h_avg=h_avg, Re_Dh=Re_Dh, R_i=R_i, R_b=R_b, R_const=R_const, R_conv=R_conv,
                  R_f=R_f, N=N, T_i=T_i, D_h=D_h, A_eff=A_eff)

    return soo, moo, outpar
