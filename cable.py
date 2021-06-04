import numpy as np
import hh as hh
import math
import matplotlib.pyplot as plt
import os 


np.set_printoptions(precision=2)

T = 1.5
N = 1000
dt = float(T)/float(N)

t = np.array([n*dt for n in range(N)])

""" Cell length , in cm"""
Cell_len = 0.01
""" Cell radius is width/2 here width = Cell length/6 (assumption) , in cm """
Cell_rad = .00241/2#Cell_len/12
""" Cell perimeter = perimeter of circle , in cm""" 
Cell_peri = 2*math.pi*Cell_rad
   
n_cells = 10 # number of nodes per cell
J = n_cells*10 # total number of nodes
dx = float(Cell_len*n_cells)/float(J)
x_grid = np.array([j*dx for j in range(J)])

def main():
    """ An implementation of the CableEquation where the Hodgkin Huxley equations
    describe the nature of Ionic current """ 

    """Adding in the constants for the equation """
    """membrane capacitance, in uF/cm^2"""
    C_m = 1
    """ resistance per unit len intracellular space, in omega-cm  """
    r_i = 1/(math.pi*Cell_rad*Cell_rad*6.7)
    """ resistance per unit len extracellular space, in omega-cm """
    r_e = 1/(math.pi*(Cell_rad+0.00001)*(Cell_rad+0.00001)*20) 
    
    r_i = 178330.0
    r_e = 888740.0

    #resting potential = -70
    V_old= [-70 for j in range(J)]

    """ This is the D*delta(t)/delta(x)^2 """
    #D_v = dt/(2*C_m*Cell_peri*(r_i + r_e)*dx*dx)
    D_v = dt/((r_i + r_e)*dx*dx*Cell_peri)
    
    A_v = np.diagflat([-D_v/2 for i in range(J-1)], -1) + np.diagflat([1. + D_v for i in range(J)]) + np.diagflat([-1*D_v/2 for i in range(J-1)], 1)
    A_v[0][0] = 1+D_v/2
    A_v[-1][-1] = 1+D_v/2

    B_v = np.diagflat([ D_v/2 for i in range(J-1)], -1) + np.diagflat([1. - D_v for i in range(J)]) + np.diagflat([D_v/2 for i in range(J-1)], 1)
    B_v[0][0] = 1-D_v/2
    B_v[-1][-1] = 1-D_v/2

    print(A_v)
    print(B_v)
   
    #change this value for the actual nth printing of the graph
    nplot = 10
    c=0

    val = 30
    V_old[0] = -70+val 
    V_old[1] = -70+val 
    V_old[2] = -70+val 
    V_old[-1] =  -70-val 
    V_old[-2] =  -70-val 
    V_old[-3] =  -70-val 

    dn = hh.HodgkinHuxley().n_inf([j+70.0 for j in V_old])
    dm =hh.HodgkinHuxley().m_inf([j+70.0 for j in V_old])
    dh = hh.HodgkinHuxley().h_inf([j+70.0 for j in V_old])

    print(V_old)
   
    # Actually making
    for i in range(1,N):
        # Here i is the ith time step in the entire simulation
                    
        hh_arr= hh.HodgkinHuxley().main([j+70.0 for j in V_old], dn, dm , dh)
        I_ion_old = hh_arr[3]
        
        B_part = np.matmul(B_v, V_old) - dt*np.asarray(I_ion_old)
        V_new = np.linalg.solve(A_v, B_part)
        
        # turn off I ion 
        #I_ion_old = [0]*J
        dn = hh_arr[0]
        dm = hh_arr[1]
        dh = hh_arr[2]
        
        
        if(i%nplot==0): #plot results every nplot timesteps
            plt.plot(x_grid,V_new,linewidth=2)
            plt.ylim([-100, 100])
            #plt.xlim([-2,35])
            filename = 'foo' + str(c+1).zfill(3) + '.jpg'
            plt.xlabel("x in cm")
            plt.ylabel("Transmembrane pot in mv")
            plt.axhline(0, color='black')
            plt.axvline(-0, color='black')
            plt.title("t = %2.2f"%(dt*(i+1)))
            plt.savefig(filename)
            plt.clf()
            c += 1
        
        V_old = V_new
        if i < 200:
            V_old[0] = -70+val 
            V_old[1] = -70+val 
            V_old[2] = -70+val 
            V_old[-1] =  -70-val 
            V_old[-2] =  -70-val 
            V_old[-3] =  -70-val 
        
os.system("ffmpeg -y -i 'foo%03d.jpg' cable_eqn.m4v")
os.system("rm -f *.jpg")

if __name__ == '__main__':
    main()