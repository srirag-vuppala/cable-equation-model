from os import times
import numpy as np
import cable
from scipy.integrate import odeint 
from scipy.integrate import quad 
from scipy.integrate import solve_ivp 

class HodgkinHuxley():
    C_m = 1.0
    """membrane capacitance, in uF/cm^2"""
    g_Na = 120.0
    """Sodium (Na) maximum conductances, in mS/cm^2"""
    g_K = 36.0
    """Postassium (K) maximum conductances, in mS/cm^2"""
    g_L = 0.3
    """Leak maximum conductances, in mS/cm^2"""
    V_Na = 115.0
    """Sodium (Na) Diffusion potentials, in mV"""
    V_K  = -12.0
    """Postassium (K) Diffusion potentials, in mV"""
    V_L  = -11.0
    """Leak current Diffusion potentials, in mV"""
    counter = 0 

    def arr_alpha_m(self, V):
        fin = []
        for v in V:
            fin.append(self.alpha_m(v))
        return fin

    def arr_alpha_n(self, V):
        fin = []
        for v in V:
            fin.append(self.alpha_n(v))
        return fin
        
    def arr_alpha_h(self, V):
        fin = []
        for v in V:
            fin.append(self.alpha_h(v))
        return fin
 
    def arr_beta_n(self, V):
        fin = []
        for v in V:
            fin.append(self.beta_n(v))
        return fin
 
    def arr_beta_m(self, V):
        fin = []
        for v in V:
            fin.append(self.beta_m(v))
        return fin
            
    def arr_beta_h(self, V):
        fin = []
        for v in V:
            fin.append(self.beta_h(v))
        return fin
 
    def alpha_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        # THis is 25 because apparently the exp fucntion doesn't like np.exp(0)
        if V == 25:
            return 0.0
            #print(0.1*(25.0 - V)/(np.exp(((25.0-V) / 10.0)) - 1.0))
        return 0.1*(25.0 - V)/(np.exp(2.5 - (0.1 * V)) - 1.0)
    def beta_m(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 4.0*np.exp(-(V / 18.0))

    def alpha_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.07*np.exp(-(V / 20.0))
    def beta_h(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 1.0/(1.0 + np.exp(3.0 - (0.1* V)))

    def alpha_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        # THis is 25 because apparently the exp fucntion doesn't like np.exp(0)
        if V == 10:
            return 0
        return (0.01*(10.0 - V))/(np.exp(1.0 - (0.1*V)) - 1.0)
    def beta_n(self, V):
        """Channel gating kinetics. Functions of membrane voltage"""
        return 0.125*np.exp(-(V / 80.0))
    
    # Add round here as well
    def n_inf(self, Vm):
        """ Inflection point potassium conductance to easily write gK"""
        final = []
        for v in Vm:
            #print(v)
            final.append(self.alpha_n(v) / (self.alpha_n(v) + self.beta_n(v)))
        return final
    def m_inf(self, Vm):
        """ Sodium activation variable """
        final = []
        for i in Vm:
            final.append(self.alpha_m(i) / (self.alpha_m(i) + self.beta_m(i)))
        return final
    def h_inf(self, Vm):
        """ Sodium inactivation variable """
        final = []
        for i in Vm:
            final.append(self.alpha_h(i) / (self.alpha_h(i) + self.beta_h(i)))
        return final

    logger = []

    # def derivatives(self, y, t, V):
    #     dy = [0]*3
    #     n = y[0]
    #     m = y[1]
    #     h = y[2]
    #     #print(t)
            
    #     v_temp = V[self.counter] 

    #     # dn/dt
    #     dy[0] = (self.alpha_n(v_temp) * (1.0 - n)) - (self.beta_n(v_temp) * n)
    #     # dm/dt
    #     dy[1] = (self.alpha_m(v_temp) * (1.0 - m)) - (self.beta_m(v_temp) * m)
    #     # dh/dt
    #     dy[2] = (self.alpha_h(v_temp) * (1.0 - h)) - (self.beta_h(v_temp) * h)
    #     self.counter += 1 
    #     #print(self.counter)
    #     self.counter = min(self.counter, cable.J - 1)
    #     return dy
    # def derivatives_n(self, t,y, V):
    #     n = [0]
    #     v_temp = V[self.counter] 
    #     # dn/dt
    #     n = (self.alpha_n(v_temp) * (1.0 - y)) - (self.beta_n(v_temp) * y)
    #     self.counter += 1 
    #     #print(self.counter)
    #     self.counter = min(self.counter, cable.J - 1)
    #     return n 
        
    # def derivatives_m(self, t,y, V):
    #     m = [0]
    #     v_temp = V[self.counter]
    #     # dn/dt
    #     m = (self.alpha_m(v_temp) * (1.0 - y)) - (self.beta_m(v_temp) * y)
    #     self.counter += 1 
    #     self.counter = min(self.counter, cable.J - 1)
    #     return m 

    # def derivatives_h(self, t,y, V):
    #     h = [0]
    #     v_temp = V[self.counter]  
    #     # dn/dt
    #     h = (self.alpha_h(v_temp) * (1.0 - y)) - (self.beta_h(v_temp) * y)
    #     self.counter += 1 
    #     self.counter = min(self.counter, cable.J - 1)
    #     return h 

    def main(self, V_old):
        # print("this is V_old")
        #print(V_old)
        #build t arr upto timestep
        # to make sure to not use future values
        #t = np.array([n*cable.dt for n in range(timestep)])
        #Y = np.array([self.n_inf(V_old[timestep-1]), self.m_inf(V_old[timestep-1]), self.h_inf(V_old[timestep-1])])
        #Y = np.array([self.n_inf(V_old), self.m_inf(V_old), self.h_inf(V_old)])

        # #prepad the values of V_old
        # for i in range(timestep):
        #     V_old.insert(0, 0)
        # print(V_old)

        #postpad the values of V_old
        # init_len = len(V_old)
        # for i in range(J-init_len):
        #     V_old = np.append(V_old, 0)
        # print(V_old)
            
        # self.counter = 0 
        # #Vy = odeint(self.derivatives, Y,t,  args=(V_old,), rtol = 1e90)
        # Vy = odeint(self.derivatives, Y,t,  args=(V_old,))
        # self.counter = 0 
        # n = Vy[:,0]
        # m = Vy[:,1]
        # h = Vy[:,2]

        # print("hi")
        # self.counter = 0 
        # Vyn = solve_ivp(self.derivatives_n, (0, cable.dt), self.n_inf(V_old),  args=(V_old,))
        # print("hi1")
        # self.counter = 0 
        # Vym = solve_ivp(self.derivatives_m, (0, cable.dt), self.m_inf(V_old), args=(V_old,))
        # print("hi2")
        # self.counter = 0 
        # Vyh = solve_ivp(self.derivatives_h, (0, cable.dt), self.h_inf(V_old),  args=(V_old,))
        # print("hi3")
        # self.counter = 0 
        # n = Vyn.y
        # m = Vym.y
        # h = Vyh.y

        n = []
        m = []
        h = []
        dn = [self.n_inf(V_old)]
        dm = [self.m_inf(V_old)]
        dh = [self.h_inf(V_old)]
        #print(dn[0][0])
        
        a_n = []
        for i in range(1, len(V_old)):
            temp_dn = []
            temp_dm = []
            temp_dh = []
            for j in range(len(V_old)):
                temp_dn.append(self.alpha_n(V_old[j]) * (1.0 - dn[i-1][j]) - (self.beta_n(V_old[j]) * dn[i-1][j]))
                temp_dm.append(self.alpha_m(V_old[j]) * (1.0 - dm[i-1][j]) - (self.beta_m(V_old[j]) * dm[i-1][j]))
                temp_dh.append(self.alpha_h(V_old[j]) * (1.0 - dh[i-1][j]) - (self.beta_h(V_old[j]) * dh[i-1][j]))

            dn.append(temp_dn)
            dh.append(temp_dh)
            dm.append(temp_dm)
        n = dn
        m = dm
        h = dh

        #return [n,m,h]
      
        #pad the values of V_old
        # init_len = len(n)
        # for i in range(J-init_len):
        #     n = np.append(n, 0)
        #     m = np.append(m, 0)
        #     h = np.append(h, 0)

        GK = []
        GNa = []
        GL = []
        for i in range(len(V_old)):
            k_temp = []
            na_temp = []
            l_temp = []
            for j in range(len(V_old)):
                #k_temp.append(self.g_K * np.power(n[i][j], 4.0))
                k_temp.append(self.g_K * n[i][j]*n[i][j]*n[i][j]*n[i][j])
                #na_temp.append(self.g_Na * np.power(m[i][j], 3.0) * h[i][j])
                na_temp.append(self.g_Na * m[i][j]*m[i][j]*m[i][j] * h[i][j])
                l_temp.append(self.g_L) 
            GK.append(k_temp)
            GNa.append(na_temp)
            GL.append(l_temp)

        I = []
        # Going upto that point
        #for i in range(timestep+1):
        for i in range(len(V_old)):
            dI = []
            for j in range(len(V_old)):
                dI.append((GK[i][j] * (V_old[i] - self.V_K)) + (GNa[i][j] * (V_old[i] - self.V_Na)) + (GL[i][j] * (V_old[i] - self.V_L)) - cable.Input_stimuli(i*cable.dt))
            I.append(dI)
        
        return [n, m, h, I]

            #I.append(-1*(1.0/ self.C_m) + (GK[i] * (V_old[i] - self.V_K)) + (GNa[i] * (V_old[i] - self.V_Na)) + (GL[i] * (V_old[i] - self.V_L)))
            #I.append(-1*(1.0/ self.C_m) + (GK[i] * (V_old[timestep] - self.V_K)) + (GNa[i] * (V_old[timestep] - self.V_Na)) + (GL[i] * (V_old[timestep] - self.V_L)))
        #I.append(-1*(1.0/ self.C_m) + (GK[timestep] * (V_old[timestep] - self.V_K)) + (GNa[timestep] * (V_old[timestep] - self.V_Na)) + (GL[timestep] * (V_old[timestep] - self.V_L)))

        #prepad the values of V_old
        # for i in range(timestep):
        #     I.insert(0, 0)

        # #postpad the values of V_old
        # init_len = len(I)
        # for i in range(J-init_len):
        #     I = np.append(I, 0)
        # print("this is the Ion current")
        # print(I)
        return I

if __name__ == '__main__':
    HodgkinHuxley().main()