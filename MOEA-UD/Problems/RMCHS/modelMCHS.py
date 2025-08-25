"""
Python translation of MATLAB function modelMCHS.m

Original signature:
    function [soo, moo] = modelMCHS(design_var, material, fluid, Q)

Returns:
    soo : float  # dimensionless entropy generation number S_gen * T_a / power
    moo : tuple  # (R_eq, DeltaP_total)

Note:
- This is a line-by-line faithful port, preserving the same parameter names and formulas.
- Units follow the original comments (SI).
"""
import numpy as np
from math import pi

_MATERIALS = {
    'Al': { # Aluminium
        'k_m': 237.0,
        'rho_m': 2707.0,
        'c_p': 897.0
    },
    'Cu': { # Copper
        'k_m': 401.0,
        'rho_m': 8954.0,
        'c_p': 385.0
    },
    'SiC': { # Silicon carbide
        'k_m': 270.0,
        'rho_m': 3300.0,
        'c_p': 750.0
    },
    'HTCG': { # High Thermal Conductivity Graphite (HTCG)
        'k_m': 1900.0,
        'rho_m': 1000.0,
        'c_p': 710.0
    },
    'Si': { # Silicon
        'k_m': 148.0,
        'rho_m': 2330.0,
        'c_p': 700.0
    },
    'TiO2': { # Titanium dioxide
        'k_m': 8.4,
        'rho_m': 4250.0,
        'c_p': 710.0
    },
}

_FLUIDS = {
    'Air': {  # Air at 300K
        'k_f': 0.0261,
        'rho_f': 1.1614,
        'nu': 1.58e-5,  # kinematic viscosity directly given (m^2/s)
        'c_p': 1007.0
    },
    'Air+10HR': { # Air at 300K + 10% relative humidity
        'k_f': 0.02671,
        'rho_f': 1.17434,
        'nu': 1.80045e-5 / 1.80045e-5,  # dynamic? original writes nu = μ/ρ
        'c_p': 1008.45
    },
    'Air+50HR': { # Air at 300K + 50% relative humidity
        'k_f': 0.028067,
        'rho_f': 1.16436,
        'nu': 1.78858e-5 / 1.16436,  # μ/ρ
        'c_p': 1025.1
    },
    'Air+90HR': { # Air at 300K + 90% relative humidity
        'k_f': 0.028583,
        'rho_f': 1.16182,
        'nu': 1.83254e-5 / 1.16182,  # μ/ρ
        'c_p': 1034.8
    },
    'H2O': { # Water at 300K
        'k_f': 0.625,
        'rho_f': 994.2,
        'nu': 7.25e-4 / 994.2,  # μ/ρ
        'c_p': 4178.0
    }
}

_NANOFLUID_KEYS = [
        'H2O+TiO2_01wt',
        'H2O+TiO2_05wt',
        'H2O+TiO2_09wt',
        'H2O+Al_01wt',
        'H2O+Al_05wt',
        'H2O+Al_09wt',
        'H2O+Al_10wt',
        'H2O+Al_50wt',
        'H2O+Al_90wt'
    ]

def _get_nanofluid_properties(fluid_key):
    # Parses fluid_key and returns nanofluid properties.
    # Example fluid_key: 'H2O+TiO2_01wt'
    try:
        base_fluid, rest = fluid_key.split('+')
        nanoparticle_material, weight_str = rest.split('_')
        if not weight_str.endswith('wt'):
            raise ValueError
        weight_fraction = float(weight_str[:-2]) / 1000.0 # convert from wt% to fraction
    except Exception as e:
        raise ValueError(f"Invalid fluid key format: {fluid_key}") from e

    k_fb    = _FLUIDS[base_fluid]['k_f']
    rho_fb  = _FLUIDS[base_fluid]['rho_f']
    mu_fb   = _FLUIDS[base_fluid]['nu'] * rho_fb  # dynamic viscosity
    c_p_fb  = _FLUIDS[base_fluid]['c_p']

    k_np    = _MATERIALS[nanoparticle_material]['k_m']
    rho_np  = _MATERIALS[nanoparticle_material]['rho_m']
    c_p_np  = _MATERIALS[nanoparticle_material]['c_p']

    n       = 3.0  # shape factor for spherical particles

    k_f = ((k_np + (n - 1) * k_fb - (n - 1) * weight_fraction
           * (k_fb - k_np)) * k_fb
           / (k_np + (n - 1) * k_fb + weight_fraction * (k_fb - k_np)))
    rho_f = (1 + weight_fraction) * rho_fb + weight_fraction * rho_np
    mu_nf = mu_fb / (1 - weight_fraction) ** 2.5
    nu = mu_nf / rho_f
    c_p = ((1 - weight_fraction) * rho_fb * c_p_fb
           + weight_fraction * rho_np * c_p_np) / rho_f

    output = dict(
        k_f=k_f, rho_f=rho_f, mu_nf=mu_nf, nu=nu, c_p=c_p)

    return output


def modelMCHS(design_var,
              material='Al', fluid='H2O', power=15.0):
    # Define the design variables
    # design_var = [alpha, beta, flow]
    alpha   = float(design_var[0])      # := w_c / (H_c/2)
    beta    = float(design_var[1])      # := w_c / w_p
    flow    = float(design_var[2])      # := G_d [m^3/s]

    # Variable parameters
    H_c     = 1.7e-3                    # Channel height (m)
    w_c     = H_c * alpha / 2.0         # Half channel width (m)
    w_p     = w_c / beta                # Half fin width (m)

    # Material properties
    k_m = _MATERIALS[material]['k_m']    # Thermal conductivity of the base material (W/m-K)
    rho_m = _MATERIALS[material]['rho_m']# Density of the base material (kg/m^3)

    # Fluid properties
    if fluid in _NANOFLUID_KEYS:
        # Get nanofluid properties
        props = _get_nanofluid_properties(fluid)

        # Thermal conductivity of the nanofluid (W/m-K)
        k_f    = props['k_f']
        # Density of the nanofluid (kg/m^3)
        rho_f  = props['rho_f']
        # Kinematic viscosity of the nanofluid (m^2/s)
        nu     = props['nu']
        # Specific heat capacity of the nanofluid (J/kg-K)
        c_p    = props['c_p']
    else:
        # Thermal conductivity of the fluid (W/m-K)
        k_f    = _FLUIDS[fluid]['k_f']
        # Density of the fluid (kg/m^3)
        rho_f  = _FLUIDS[fluid]['rho_f']
        # Kinematic viscosity of the fluid (m^2/s)
        nu     = _FLUIDS[fluid]['nu']
        # Specific heat capacity of the fluid (J/kg-K)
        c_p    = _FLUIDS[fluid]['c_p']

    # Prandtl number
    Pr = nu * rho_f * c_p / k_f

    # Fixed geometry & interface parameters
    H_b = 0.1e-3        # Base thickness (m)
    W_d = 51e-3         # Device width (m)
    L_d = 51e-3         # Device length (m)
    W_i = W_d           # Interface width (m)
    L_i = L_d           # Interface length (m)
    # q  = 5e4          # not used in original
    T_a = 300.0         # Ambient temperature (K)
    R_iA = 2.75e-4      # Interface thermal resistance (m^2-K/W)

    # Geometry parameters
    # -------------------
    # Number of channels
    N       = (W_d/2.0 - w_p) // (w_c + w_p)
    # Cross-sectional area of one channel (m^2)
    A_c_    = 2.0 * w_c * H_c
    # Wetted perimeter of one channel (m)
    P_      = 2.0 * (H_c + 2.0*w_c)
    # Hydraulic diameter of one channel (m)
    D_h     = 4.0 * A_c_ / P_
    # Base area (m^2)
    A_b     = W_d * L_d
    # Interface area (m^2)
    A_i     = W_i * L_i

    # Flow & heat parameters
    # ----------------------
    # Mass flow rate per channel (kg/s)
    mdot  = rho_f * flow / (2.0 * N)
    # Average velocity in one channel (m/s)
    U_avg = mdot / (rho_f * w_c * H_c)
    # Reynolds number based on hydraulic diameter
    Re_Dh = D_h * U_avg / nu

    # Nusselt number (Nu) & friction factor (f)
    if Re_Dh < 2300: # Laminar flow
        # Nusselt number
        Nu_Dh   = 2.253 + 8.164 * (1.0 / (alpha + 1.0)) ** 1.5
        # Friction factor
        f       = np.sqrt((3.2*(Re_Dh*D_h/L_d)**0.57) ** 2
            + (4.7 + 19.64 * (alpha ** 2 + 1.0) / (alpha + 1.0) ** 2) ** 2) / Re_Dh
    else: # Turbulent flow
        # Frinction factor
        f      = 1 / (0.79 * np.log(Re_Dh) - 1.64) ** 2
        # Nusselt number
        Nu_Dh   = (f/2.0) * (Re_Dh - 1000.0) * Pr / (
                1.0 + 12.7 * np.sqrt(f/2.0) * (Pr**(2.0/3.0) - 1.0))

    # Thermal resistances
    # -------------------
    # Interface thermal resistance (K/W)
    R_i     = R_iA / A_i

    # Average convective heat transfer coefficient (W/m^2-K)
    h_avg   = k_f * Nu_Dh / D_h
    # Parameter for fin efficiency
    mHc     = np.sqrt(2.0 * (2.0 * w_p + L_d) * h_avg / (k_m * 2.0 * w_p * L_d)) * H_c
    # Fin efficiency
    eta_p   = np.tanh(mHc) / mHc
    # Effective heat transfer area (m^2)
    A_eff   = 2.0 * N * (eta_p*H_c + w_c) * L_d

    # Thermal resistance due to convection (K/W)
    R_conv  = 1.0 / (h_avg * A_eff)
    # Thermal resistance due to fluid flow (K/W)
    R_f     = 1.0 / (rho_f * flow * c_p)
    # Thermal resistance due to conduction (constriction) in the base (K/W)
    R_const = ((1.0 + beta) / (pi * k_m * A_b * beta)) * np.log(1.0 / np.sin(0.5*pi/(1.0 + beta))) * H_c * alpha

    # Geometrical parameters and symetry factors for thermal resistance network
    a       = np.sqrt(A_i/pi)
    b       = np.sqrt(A_b/pi)
    tau     = H_b / b
    epsilon = a / b

    # Thermal resistances between ambient and interface (K/W)
    R_o     = R_conv + R_f   # + R_const

    # Nondimensional Biot number
    Bi      = 1.0 / (pi * k_m * b * R_o)

    # Diemnsionless temperature parameters
    lambda_c    = pi + 1.0 / (np.sqrt(pi) * epsilon)
    Phi_c       = (np.tanh(lambda_c*tau) + lambda_c/Bi) / (1.0 + lambda_c*np.tanh(lambda_c*tau)/Bi)
    Psi_avg     = epsilon*tau/np.sqrt(pi) + 0.5*Phi_c*(1.0 - epsilon)**1.5

    # Thermal resistance due to diffusion in the base (K/W)
    R_b         = Psi_avg / (np.sqrt(pi) * k_m * a)

    # Equivalent thermal resistance (K/W)
    R_eq        = R_o  # + R_b + R_i

    # Interface temperature (K)
    T_i         = T_a + R_eq * power

    # Pressure Drops
    # ----------------
    # Entrance & exit losses
    k_ce   = 1.79 - 2.32*(beta/(1.0 + beta)) + 0.53*(beta/(1.0 + beta))**2
    # Pressure drop inside the channels
    DeltaP = 0.5 * rho_f * U_avg**2 * (f * (L_d / D_h) + k_ce)

    # Supply tube length & diameter
    L_tu = 0.5
    D_tu = 1.9e-2
    # Average velocity in the supply tube (m/s)
    U_avgtub    = 4.0 * flow / (pi * D_tu ** 2)
    # Friction factor in the supply tube (f_tu)
    f_tu        = (4 * (0.0929 + 1.01612 / (L_tu/D_tu)) *
                   Re_Dh ** (-0.268 - 0.3193 / (L_tu/D_tu)))
    # Relative area of the tube to the total channel area
    A_tuhs      = 0.25 * pi * D_tu**2 / (W_d * H_c)

    # Pressure drop in the supply tube (Pa)
    DeltaP_tu   = 0.5 * rho_f * U_avgtub**2 * (
            0.42 + (1 - A_tuhs**2)**2
            + 0.42*(1 - A_tuhs**2) + 1 + 2*f_tu*L_tu/D_tu)

    # Total pressure drop (Pa)
    DeltaP_total = DeltaP  # + DeltaP_tu

    # Pumping power (W)
    Phi = flow * DeltaP_total

    # Entropy generation rates (W/K)
    # ------------------------------
    # Due to heat transfer (thermal irreversibility)
    S_genht     = power**2 * R_eq / (T_a * T_i)
    # Due to fluid flow (fluid friction irreversibility)
    S_genff     = flow * DeltaP_total / T_a
    # Total entropy generation rate (W/K)
    S_gen       = S_genht + S_genff

    # Heat sink mass (kg)
    Md  = rho_m * L_d * (2.0*W_d*H_b + 2.0*w_p*H_c*(N + 1.0))

    # Outputs
    # -------
    # Single objective (dimensionless entropy generation number)
    soo = S_gen * T_a / power

    # Multiple objectives (R_eq, DeltaP_total)
    moo = (S_genht * T_a / power, S_genff * T_a / power)
    # moo = (R_eq, Phi)

    # Additional outputs (for debugging / analysis parity)
    outpar = dict(
        S_gen=S_gen, S_genht=S_genht, S_genff=S_genff,
        R_eq=R_eq, DeltaP=DeltaP_total,
        U_avg=U_avg, Phi=Phi, Md=Md,
        h_avg=h_avg, Re_Dh=Re_Dh, R_i=R_i, R_b=R_b, R_const=R_const, R_conv=R_conv, R_f=R_f, N=N, T_i=T_i, D_h=D_h, A_eff=A_eff)

    return soo, moo, outpar
