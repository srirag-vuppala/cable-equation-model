import numpy as np

class HodgkinHuxley():
    # membrane capacitance, in uF/cm^2
    C_m = 1.0
    # Sodium (Na) maximum conductances, in mS/cm^2
    g_Na = 120.0
    # Postassium (K) maximum conductances, in mS/cm^2
    g_K = 36.0
    # Leak maximum conductances, in mS/cm^2
    g_L = 0.3
    # Sodium (Na) Diffusion potentials, in mV
    V_Na = 115.0
    # Postassium (K) Diffusion potentials, in mV
    V_K  = -12.0
    # Leak current Diffusion potentials, in mV
    V_L  = 10.6

    counter = 0 

    def alpha_m(self, V):
        # Channel gating kinetics. Functions of membrane voltage
        return 0.1*(25.0 - V)/(np.exp(2.5 - (0.1 * V)) - 1.0)
    def beta_m(self, V):
        # Channel gating kinetics. Functions of membrane voltage
        return 4.0*np.exp(-(V / 18.0))

    def alpha_h(self, V):
        # Channel gating kinetics. Functions of membrane voltage
        return 0.07*np.exp(-(V / 20.0))
    def beta_h(self, V):
        # Channel gating kinetics. Functions of membrane voltage
        return 1.0/(1.0 + np.exp(3.0 - (0.1* V)))

    def alpha_n(self, V):
        # Channel gating kinetics. Functions of membrane voltage"""
        return (0.01*(10.0 - V))/(np.exp(1.0 - (0.1*V)) - 1.0)
    def beta_n(self, V):
        # Channel gating kinetics. Functions of membrane voltage
        return 0.125*np.exp(-(V / 80.0))
    
    def n_inf(self, Vm):
        # Inflection point potassium conductance to easily write gK
        final = []
        for v in Vm:
            final.append(self.alpha_n(v) / (self.alpha_n(v) + self.beta_n(v)))
        return final
    def m_inf(self, Vm):
        # Sodium activation variable 
        final = []
        for i in Vm:
            final.append(self.alpha_m(i) / (self.alpha_m(i) + self.beta_m(i)))
        return final
    def h_inf(self, Vm):
        # Sodium inactivation variable 
        final = []
        for i in Vm:
            final.append(self.alpha_h(i) / (self.alpha_h(i) + self.beta_h(i)))
        return final

    def main(self, V_old):
        n=self.n_inf(V_old)
        m=self.m_inf(V_old)
        h=self.h_inf(V_old)

        GK = []
        GNa = []
        GL = []
        for i in range(len(V_old)):
            GK.append((self.g_K / self.C_m) * np.power(n[i], 4.0))
            #GNa.append((self.g_Na / self.C_m) * np.power(m[i], 3.0) * h[i])
            GNa.append((self.g_Na / self.C_m) * np.power(m[i], 3.0) * 0.596)
            GL.append(self.g_L / self.C_m) 
        I = []
        for i in range(len(V_old)):
            I.append((GK[i] * (V_old[i] - self.V_K)) + (GNa[i] * (V_old[i] - self.V_Na)) + (GL[i] * (V_old[i] - self.V_L)) - 0) 

        return I

if __name__ == '__main__':
    HodgkinHuxley().main()