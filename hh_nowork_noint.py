import numpy as np
from cable import Input_stimuli
from scipy.integrate import odeint 
from scipy.integrate import quad 
from scipy.integrate import solve_ivp 

class HodgkinHuxley():
    C_m = 1
    """membrane capacitance, in uF/cm^2"""
    g_Na = 120
    """Sodium (Na) maximum conductances, in mS/cm^2"""
    g_K = 36
    """Postassium (K) maximum conductances, in mS/cm^2"""
    g_L = 0.3
    """Leak maximum conductances, in mS/cm^2"""
    V_Na = 115
    """Sodium (Na) Diffusion potentials, in mV"""
    V_K  = -12
    """Postassium (K) Diffusion potentials, in mV"""
    V_L  = -11
    """Leak current Diffusion potentials, in mV"""
    counter = 0 

    def alpha_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(0.1*(25 - V)/(np.exp((25-V) / 10) - 1), 2)
        except:
            print("hi")
            return 0
    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(4.0*np.exp(-(V / 18.0)), 2)
        except:
            print("hi")
            return 0

    def alpha_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(0.07*np.exp(-(V / 20.0)), 2)
        except:
            print("hi")
            return 0
    def beta_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(1.0/(1.0 + np.exp((30 - V) / 10.0)), 2)
        except:
            print("hi")
            return 0

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(0.01*(10 - V)/(np.exp((10 - V) / 10.0) - 1), 2)
        except:
            print("hi")
            return 0
    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        try:
            return round(0.125*np.exp(-(V / 80.0)), 2)
        except:
            print("hi")
            return 0
    
    # Add round here as well
    def n_inf(self, Vm=0.0):
        """ Inflection point potassium conductance to easily write gK"""
        try:
            return round(self.alpha_n(Vm) / (self.alpha_n(Vm) + self.beta_n(Vm)), 2)
        except:
            return 0
    def m_inf(self, Vm=0.0):
        """ Sodium activation variable """
        try:
            return round(self.alpha_m(Vm) / (self.alpha_m(Vm) + self.beta_m(Vm)), 2)
        except:
            return 0
    def h_inf(self, Vm=0.0):
        """ Sodium inactivation variable """
        try:
            return round(self.alpha_h(Vm) / (self.alpha_h(Vm) + self.beta_h(Vm)), 2)
        except:
            return 0

    logger = []

    def derivatives(self, y, t, V):
        dy = [0]*3
        n = y[0]
        m = y[1]
        h = y[2]
        #print(t)
        self.logger.append(t)
            

        # dn/dt
        dy[0] = (self.alpha_n(V[round(self.counter) ]) * (1.0 - n)) - (self.beta_n(V[round(self.counter) ]) * n)
        # dm/dt
        dy[1] = (self.alpha_m(V[round(self.counter) ]) * (1.0 - m)) - (self.beta_m(V[round(self.counter) ]) * m)
        # dh/dt
        dy[2] = (self.alpha_h(V[round(self.counter) ]) * (1.0 - h)) - (self.beta_h(V[round(self.counter) ]) * h)
        self.counter += 1 
        print(self.counter)
        #I need to figure out a way to control the number of time the integration occurs
        #self.counter = min(self.counter, 99)
        # if self.counter == 100:
        #     self.counter = 0
        return dy

    def main(self, t, V_old, J):
            
        Y = np.array([self.n_inf(V_old[0]), self.m_inf(V_old[0]), self.h_inf(V_old[0])])
        self.counter = 0 
        # the present implementation works becouse of rtol which essentially limits the actual integration dt within the odeint function.
        self.counter = 0 
        dn = [self.n_inf(V_old[0])]
        dm = [self.m_inf(V_old[0])]
        dh = [self.h_inf(V_old[0])]
        for i in range(1, len(V_old)):
            try:
                temp_dn = round(round((self.alpha_n(V_old[i]) * (1.0 - round(dn[i-1], 2))),2) - round((self.beta_n(V_old[i]) * round(dn[i-1],2)),2), 2)
                temp_dm = round(round((self.alpha_m(V_old[i]) * (1.0 - round(dm[i-1], 2))),2) - round((self.beta_m(V_old[i]) * round(dm[i-1],2)),2), 2)
                temp_dh =round(round((self.alpha_h(V_old[i]) * (1.0 - round(dh[i-1], 2))),2) - round((self.beta_h(V_old[i]) * round(dh[i-1],2)),2), 2)
            except:
                temp_dh =0
                temp_dm =0
                temp_dn =0
                print("hi")

            dn.append(temp_dn)
            dh.append(temp_dh)
            dm.append(temp_dm)
        # print(dn)
        # print(dh)
        # print(dm)

        n = dn 
        m = dm
        h = dh 

        GK = []
        GNa = []
        GL = []
        for i in range(len(V_old)):
            GK.append(round((self.g_K / self.C_m) * np.power(n[i], 4.0), 2))
            GNa.append(round((self.g_Na / self.C_m) * np.power(m[i], 3.0) * h[i], 2))
            GL.append(round(self.g_L / self.C_m, 2)) 
        I = []
        for i in range(len(V_old)):
            I.append(round(-1*(Input_stimuli(t[i]) / self.C_m) + (GK[i] * (V_old[i] - self.V_K)) + (GNa[i] * (V_old[i] - self.V_Na)) + (GL[i] * (V_old[i] - self.V_L)), 2))

        #pad 
        for i in range(J-len(V_old)):
            I.append(0)
            
        return I

if __name__ == '__main__':
    HodgkinHuxley().main()