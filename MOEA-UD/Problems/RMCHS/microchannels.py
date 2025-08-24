import numpy as np
import matplotlib.pyplot as plt
from . import materials


class Microchannel():
    W_d = 51e-3  # Heat sink's width (m)
    L_d = 51e-3  # Heat sink's length (m)
    H_c = 1.7e-3  # Channel's height (m)
    H_b = 0.1e-3  # Channel's base (m)
    H_d = H_c + H_b  # Heat sink height (m)
    W_i = 51e-3  # Interface's width (m)
    L_i = 51e-3  # Interface's length (m)

    T_a = 300  # Ambient temperature (K)
    R_i_by_A = 2.75e-4  # (K.m2/W)

    L_tu = 0.5  # Supply tube length (m)
    D_tu = 1e-2  # Supply tube diameter (m)

    def __init__(self, base=materials.Copper(), coolant=materials.Water()):
        # Read materials
        self._bmat = base
        self._cmat = coolant

        # Fluid flow rate (m3/s)
        self.G_d = 7e-3

        # Heat flux density (W/m2)
        self.q = 15e4

        # Geometry properties (store as backing fields to avoid recursive properties)
        self._w_c = (250e-6 / 2)   # channel half-width (m)
        self._w_w = (140e-6 / 2)   # wall half-width (m)

        # Derived ratios (initialised from backing fields)
        self._alpha = 2 * self._w_c / self.H_c
        self._beta = self._w_c / self._w_w

    @property
    def Pr(self):
        # Prandtl number
        return self._cmat.nu * self._cmat.rho * self._cmat.c_p / self._cmat.k

    @property
    def N(self):
        # Number of channels (continuous, per MATLAB model)
        pitch = self.w_c + self.w_w
        usable = max(self.W_d / 2.0 - self.w_w, 0.0)
        N = usable // max(pitch, 1e-12)
        return max(float(N), 1e-12)

    # --- Geometry: widths and ratios (consistent, non-recursive) ---
    @property
    def w_c(self):
        """Channel half-width (m)."""
        return self._w_c

    @w_c.setter
    def w_c(self, value):
        self._w_c = float(max(value, 1e-12))
        # keep ratios in sync
        self._alpha = 2 * self._w_c / self.H_c
        self._beta = self._w_c / max(self._w_w, 1e-12)

    @property
    def w_w(self):
        """Wall half-width (m)."""
        return self._w_w

    @w_w.setter
    def w_w(self, value):
        self._w_w = float(max(value, 1e-12))
        # keep ratios in sync
        self._beta = self._w_c / self._w_w

    @property
    def alpha(self):
        """Channel aspect ratio = 2*w_c/H_c."""
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = float(max(value, 1e-12))
        # update w_c from alpha, keep beta consistent
        self._w_c = self._alpha * self.H_c / 2.0
        self._beta = self._w_c / max(self._w_w, 1e-12)

    @property
    def beta(self):
        """Channel-wall ratio = w_c/w_w."""
        return self._beta

    @beta.setter
    def beta(self, value):
        self._beta = float(max(value, 1e-12))
        # update w_w from beta, keep alpha consistent
        self._w_w = self._w_c / self._beta
        self._alpha = 2 * self._w_c / self.H_c

    @property
    def A_b(self):
        return self.W_d * self.L_d

    @property
    def A_i(self):
        return self.W_i * self.L_i

    @property
    def Q(self):
        return self.q * self.A_i

    # %% Thermal resistance calculations
    @property
    def R_f(self):
        # Fluid thermal resistance [K/W]
        return 1.0 / max(self._cmat.rho * self.G_d * self._cmat.c_p, 1e-12)

    @property
    def T_f(self):
        # Fluid bulk temperature at outlet [K] relative to ambient rise
        return self.T_a + self.R_f * self.Q

    @property
    def h_avg(self):
        # Film (convection) coefficient [W/K.m2]
        return self._cmat.k * self.Nu_Dh / self.D_h

    @property
    def eta_p(self):
        # Fin efficiency
        mHc = np.sqrt(2.0 * (2.0 * self.w_w + self.L_d) * self.h_avg / (
                self._bmat.k * 2.0 * self.w_w * self.L_d)) * self.H_c
        return np.tanh(mHc) / mHc

    @property
    def R_conv(self):
        # Convection thermal resistance [K/W]
        A_eff = 2.0 * self.N * (self.eta_p * self.H_c + self.w_c) * self.L_d
        return 1.0 / max(self.h_avg * max(A_eff, 1e-12), 1e-12)

    @property
    def T_s(self):
        # Solid–fluid surface temperature at the channel interface [K]
        # T_s = T_a + (R_f + R_conv) * Q (per MATLAB model structure)
        return self.T_a + (self.R_f + self.R_conv) * self.Q

    @property
    def R_o(self):
        # Thermal resistance of the microchannel heat sink [K/W]
        return self.R_conv + self.R_f

    @property
    def Bi(self):
        # Biot's number
        return 1.0 / (np.pi * self._bmat.k * self.R_o)

    @property
    def R_b(self):
        # Thermal resistance due contraction/dispersion
        sqrtpi = np.sqrt(np.pi)
        a = np.sqrt(self.A_i / np.pi)
        b = np.sqrt(self.A_b / np.pi)
        tau = self.H_b / b
        epsilon = a / b
        lambda_c = np.pi + 1.0 / (sqrtpi * epsilon)
        Phi_c = (np.tanh(lambda_c * tau) + lambda_c / self.Bi ) / (
                1.0 + lambda_c * np.tanh(lambda_c * tau) / self.Bi)
        Psi_avg = epsilon * tau / sqrtpi + 0.5 * Phi_c * \
                  (1.0 - epsilon) ** (3 / 2)
        return Psi_avg / (sqrtpi * self._bmat.k * a)

    @property
    def T_b(self):
        # Base temperature just below the interface [K]
        # R_b is omitted in R_eq per MATLAB; keep T_b equal to T_s unless enabled
        return self.T_s  # + self.R_b * self.Q

    @property
    def R_i(self):
        # Thermal interface resistance [K/W]
        return self.R_i_by_A / self.A_i

    @property
    def T_i(self):
        # Interfacial (chip) temperature [K]
        # MATLAB: T_i = T_a + R_eq * Q
        return self.T_a + self.R_eq * self.Q

    @property
    def R_eq(self):
        # Equivalent thermal resistance [K/W]
        # MATLAB model uses only R_o (R_conv + R_f) for R_eq; R_b and R_i are optional
        return self.R_o  # + self.R_b + self.R_i if you decide to include them

    # %% Sgen_ff calculations
    @property
    def DeltaP_mh(self):
        # Pressure drop inside microchannels
        k_ce = 1.79 - 2.32 * (self.beta / (1.0 + self.beta)) + 0.53 * (self.beta / (1.0 + self.beta)) ** 2
        return 0.5 * self._cmat.rho * self.U_avg ** 2 * (self.f * (self.L_d / max(self.D_h, 1e-12)) + k_ce)

    @property
    def U_avg_tu(self):
        # Average fluid velocity inside the supply tube
        return 4.0 * self.G_d / (np.pi * self.D_tu ** 2)

    @property
    def Re_D_tu(self):
        # Reynolds' number inside the supply tube
        return self.U_avg_tu * self.D_tu / self._cmat.nu

    @property
    def f_tu(self):
        # Friction factor inside the supply tube
        return 4.0 * (0.09290 + 1.01612 / (self.L_tu / self.D_tu)) * \
               self.Re_D_tu ** (-0.268 - 0.3193 / (self.L_tu / self.D_tu))

    @property
    def A_tu_hs(self):
        # Ratio between areas of tube's transversal section and heat sink
        return 0.25 * np.pi * self.D_tu ** 2 / max(self.W_d * self.H_c, 1e-12)

    @property
    def DeltaP_tu(self):
        # Pressure drop inside the supply tube
        return 0.5 * self._cmat.rho * self.U_avg_tu ** 2 * (
                0.42 + (1 - self.A_tu_hs ** 2) ** 2 + 0.42 * (
                1 - self.A_tu_hs ** 2) + 1. + 2 * self.f_tu * self.L_tu / self.D_tu)

    @property
    def DeltaP_total(self):
        """Total Pressure Drop."""
        return self.DeltaP_mh  # + self.DeltaP_tu

    @property
    def Phi(self):
        # Pumping power
        return self.G_d * self.DeltaP_total()

    # %% Entropy-related calculations
    @property
    def sgen_ht(self):
        return self.Q ** 2 * self.R_eq / (self.T_a * max(self.T_i, 1e-12))

    @property
    def sgen_ff(self):
        return self.G_d * self.DeltaP_total() / self.T_a

    @property
    def sgen(self):
        return self.sgen_ht + self.sgen_ff

    @property
    def Be(self):
        # Bejan's number
        return self.sgen_ht / self.sgen

    # %%
    # @property
    # def U_avg_sound(self):
    #     return (self.gamma * self.Rgas * self.T_a) ** (1 / 2)

    @property
    def D_h(self):
        # Hydraulic diameter for the channel unit cell (per MATLAB):
        # A_c = 2*w_c*H_c, P = 2*(H_c + 2*w_c), D_h = 4*A_c/P
        A_c = 2.0 * self.w_c * self.H_c
        P = 2.0 * (self.H_c + 2.0 * self.w_c)
        return 4.0 * A_c / max(P, 1e-12)

    @property
    def mdot(self):
        # Mass flow per channel (flow split into two halves with N channels each)
        return self._cmat.rho * self.G_d / max(2.0 * self.N, 1e-12)

    @property
    def U_avg(self):
        # Average velocity inside one microchannel [m/s]
        area = max(self.w_c * self.H_c, 1e-12)
        return self.mdot / (self._cmat.rho * area)

    # @property
    # def Ma(self):
    #     # Mach number
    #     return self.U_avg / self.U_avg_sound

    @property
    def Re_D(self):
        return self.U_avg * self.D_h / self._cmat.nu

    @property
    def flow_regime(self):
        return 'laminar' if self.Re_D < 2300 else 'turbulent'

    @property
    # Friction coefficient
    def f(self):
        if self.flow_regime == 'laminar':
            ab_constant = (self.alpha ** 2 + 1.0) / ((self.alpha + 1.0) ** 2.0)
            val = max(self.Re_D * self.D_h / self.L_d, 1e-12)
            num = (3.2 * (val ** 0.57)) ** 2.0 + (4.7 + 19.64 * ab_constant) ** 2.0
            return np.sqrt(num) / max(self.Re_D, 1e-12)
        else:
            return 1.0 / (0.79 * np.log(max(self.Re_D, 1e-12)) - 1.64) ** 2

    @property
    # Nusselt number
    def Nu_Dh(self):
        if self.flow_regime == 'laminar':
            return 2.253 + 8.164 * (1.0 / (self.alpha + 1.0)) ** 1.5
        else:
            denom = 1.0 + 12.7 * np.sqrt(max(self.f, 1e-12) / 2.0) * (self.Pr ** (2.0 / 3.0) - 1.0)
            return (self.f / 2.0) * max(self.Re_D - 1000.0, 0.0) * self.Pr / max(denom, 1e-12)

    def plotc(self):
        z = 0
        sgenp = np.arange(10000)
        wcp = np.arange(10000)
        wwp = np.arange(10000)
        sgenp = sgenp.astype(float)
        wcp = wcp.astype(float)
        wwp = wwp.astype(float)
        for x in range(100):
            self.w_c = (250e-6 / 2) - (x * 1e-6)
            self.w_w = (250e-6 / 2)
            for y in range(100):
                self.w_w = (250e-6 / 2) - (y * 1e-6)
                sgenp[z] = self.sgen
                wcp[z] = self.w_c
                wwp[z] = self.w_w
                z = z + 1

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.xlabel('Wc')
        plt.ylabel('Ww')
        ax.set_zlabel(r'$\dot{S}_{gen}$')
        ax.scatter(wcp, wwp, sgenp)

        plt.show()

        # wc = np.linspace(125e-6, 0.0004, 100)
        # sgenp = []
        #wc = []

        # for wc_i in wc:
        #     self.wc = wc_i
        #     sgenp.append(self.tegr())

        # # ax = plt.axes(projection='3d')
        # plt.xlabel(r'$w_c$')
        # # plt.ylabel(r'$\beta$')
        # plt.ylabel(r'$\dot{S}_{gen}$')
        # plt.plot(np.array(wc), np.array(sgenp))
        # plt.show()

#class Micromixer():


if __name__ == '__main__':
    from . import materials

    mc = Microchannel()
    mc.plotc()
   # mm = Micromixer()
