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
        return round(0.1*(25 - V)/(np.exp((25-V) / 10) - 1), 2)
    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return round(4.0*np.exp(-(V / 18.0)), 2)

    def alpha_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return round(0.07*np.exp(-(V / 20.0)), 2)
    def beta_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return round(1.0/(1.0 + np.exp((30 - V) / 10.0)), 2)

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return round(0.01*(10 - V)/(np.exp((10 - V) / 10.0) - 1), 2)
    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return round(0.125*np.exp(-(V / 80.0)), 2)
    
    # Add round here as well
    def n_inf(self, Vm=0.0):
        """ Inflection point potassium conductance to easily write gK"""
        return self.alpha_n(Vm) / (self.alpha_n(Vm) + self.beta_n(Vm))
    def m_inf(self, Vm=0.0):
        """ Sodium activation variable """
        return self.alpha_m(Vm) / (self.alpha_m(Vm) + self.beta_m(Vm))
    def h_inf(self, Vm=0.0):
        """ Sodium inactivation variable """
        return self.alpha_h(Vm) / (self.alpha_h(Vm) + self.beta_h(Vm))

    logger = []

    def derivatives(self, y, V):
        dy = [0]*3
        n = y[0]
        m = y[1]
        h = y[2]
        #print(t)

        # dn/dt
        dy[0] = (self.alpha_n(V) * (1.0 - n)) - (self.beta_n(V) * n)
        # dm/dt
        dy[1] = (self.alpha_m(V) * (1.0 - m)) - (self.beta_m(V) * m)
        # dh/dt
        dy[2] = (self.alpha_h(V) * (1.0 - h)) - (self.beta_h(V) * h)
        self.counter += 1 
        #print(self.counter)
        #I need to figure out a way to control the number of time the integration occurs
        #self.counter = min(self.counter, 99)
        # if self.counter == 100:
        #     self.counter = 0
        return dy

    def main(self, t, V_old, J):
        init_len = len(V_old)
        """
        pad =  [0 for i in range(J-init_len)]
        print(pad)
        print(V_old)
        np.concatenate((V_old,pad))
        """
        for i in range(J-init_len):
            V_old = np.append(V_old, 0)
            
        Y = np.array([self.n_inf(V_old[0]), self.m_inf(V_old[0]), self.h_inf(V_old[0])])
        self.counter = 0 
        # the present implementation works becouse of rtol which essentially limits the actual integration dt within the odeint function.
        Vy = odeint(self.derivatives, Y,V_old, rtol = 1e0)
        self.counter = 0 
        #print(self.logger)
        #use only dt only for the integrate
        #don't calculate future values
        n = Vy[:,0]
        m = Vy[:,1]
        h = Vy[:,2]

        GK = []
        GNa = []
        GL = []
        for i in range(len(V_old)):
            GK.append((self.g_K / self.C_m) * np.power(n[i], 4.0))
            GNa.append((self.g_Na / self.C_m) * np.power(m[i], 3.0) * h[i])
            GL.append(self.g_L / self.C_m) 
        I = []
        for i in range(len(V_old)):
            I.append(-1*(Input_stimuli(t[i]) / self.C_m) + (GK[i] * (V_old[i] - self.V_K)) + (GNa[i] * (V_old[i] - self.V_Na)) + (GL[i] * (V_old[i] - self.V_L)))
        #padding I with zeros for values unknown to us
        # This is the value of J
        # shape = np.shape(I)
        # shape_arr = np.zeros(J)
        return I

if __name__ == '__main__':
    HodgkinHuxley().main()