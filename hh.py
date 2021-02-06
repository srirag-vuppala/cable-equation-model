from os import times
import numpy as np
import cable
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
        return 0.1*(25 - V)/(np.exp((25-V) / 10) - 1)
    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 4.0*np.exp(-(V / 18.0))

    def alpha_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.07*np.exp(-(V / 20.0))
    def beta_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 1.0/(1.0 + np.exp((30 - V) / 10.0))

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.01*(10 - V)/(np.exp((10 - V) / 10.0) - 1)
    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.125*np.exp(-(V / 80.0))
    
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

    def derivatives(self, y, t, V):
        dy = [0]*3
        n = y[0]
        m = y[1]
        h = y[2]
        #print(t)
            
        v_temp = V[self.counter] 

        # dn/dt
        dy[0] = (self.alpha_n(v_temp) * (1.0 - n)) - (self.beta_n(v_temp) * n)
        # dm/dt
        dy[1] = (self.alpha_m(v_temp) * (1.0 - m)) - (self.beta_m(v_temp) * m)
        # dh/dt
        dy[2] = (self.alpha_h(v_temp) * (1.0 - h)) - (self.beta_h(v_temp) * h)
        self.counter += 1 
        #print(self.counter)
        self.counter = min(self.counter, cable.J - 1)
        return dy

    def main(self, t, V_old, J, timestep):
        #build t arr upto timestep
        # to make sure to not use future values
        t = np.array([n*cable.dt for n in range(timestep)])
        
        Y = np.array([self.n_inf(V_old[timestep-1]), self.m_inf(V_old[timestep-1]), self.h_inf(V_old[timestep-1])])
        # Y = np.array([self.n_inf(V_old[0]), self.m_inf(V_old[0]), self.h_inf(V_old[0])])

        # #prepad the values of V_old
        # for i in range(timestep):
        #     V_old.insert(0, 0)
        # print(V_old)

        # #postpad the values of V_old
        # init_len = len(V_old)
        # for i in range(J-init_len):
        #     V_old = np.append(V_old, 0)
        #print(V_old)
            
        self.counter = 0 
        #Vy = odeint(self.derivatives, Y,t,  args=(V_old,), rtol = 1e9)
        Vy = odeint(self.derivatives, Y,t,  args=(V_old,))
        self.counter = 0 
        n = Vy[:,0]
        m = Vy[:,1]
        h = Vy[:,2]

        #pad the values of V_old
        init_len = len(n)
        for i in range(J-init_len):
            n = np.append(n, 0)
            m = np.append(m, 0)
            h = np.append(h, 0)

        GK = []
        GNa = []
        GL = []
        #for i in range(len(V_old)):
        for i in range(len(V_old)):
            GK.append((self.g_K / self.C_m) * np.power(n[i], 4.0))
            GNa.append((self.g_Na / self.C_m) * np.power(m[i], 3.0) * h[i])
            GL.append(self.g_L / self.C_m) 

        I = []
        # Going upto that point
        for i in range(timestep+1):
            I.append(-1*(1.0/ self.C_m) + (GK[i] * (V_old[i] - self.V_K)) + (GNa[i] * (V_old[i] - self.V_Na)) + (GL[i] * (V_old[i] - self.V_L)))

            #I.append(-1*(1.0/ self.C_m) + (GK[i] * (V_old[timestep] - self.V_K)) + (GNa[i] * (V_old[timestep] - self.V_Na)) + (GL[i] * (V_old[timestep] - self.V_L)))
        #I.append(-1*(1.0/ self.C_m) + (GK[timestep] * (V_old[timestep] - self.V_K)) + (GNa[timestep] * (V_old[timestep] - self.V_Na)) + (GL[timestep] * (V_old[timestep] - self.V_L)))

        #prepad the values of V_old
        # for i in range(timestep):
        #     I.insert(0, 0)

        # #postpad the values of V_old
        init_len = len(I)
        for i in range(J-init_len):
            I = np.append(I, 0)
        print("this is the Ion current")
        print(I)
        return I

if __name__ == '__main__':
    HodgkinHuxley().main()