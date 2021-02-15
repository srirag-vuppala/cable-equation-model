import numpy as np
import hh as hh
import math
import matplotlib.pyplot as plt
import os, sys


def Input_stimuli(t):
    """ Current applied to create stimulus which is dependent on time, in milli Ampere(A)/cm^2 """
    # if 0 < t < 10:
    #     return 5.0
    # elif 5.0 < t < 10.0:
    #     return 1.0
    return 0.0 

np.set_printoptions(precision=2)

#T = 1000 
T = 6 
#N = 2000
N = 100
dt = float(T)/float(N)
t = np.array([n*dt for n in range(N)])

""" Creating a one dimensional domain with unit length and J=100 equally spaced grid points """
L = 30
#J = 1000
J = 100
dx = float(L)/float(J)
x_grid = np.array([j*dx for j in range(J)])

def main():
    """ An implementation of the CableEquation where the Hodgkin Huxley equations
    describe the nature of Ionic current """ 

    """Adding in the constants for the equation """
    """membrane capacitance, in uF/cm^2"""
    C_m = 1
    """ Cell length , in cm"""
    Cell_len = 0.01
    """ Cell radius is width/2 here width = Cell length/6 (assumption) , in cm """
    Cell_rad = Cell_len/12
    """ Cell perimeter = perimeter of circle , in cm""" 
    Cell_peri = 2*math.pi*Cell_rad
    """ resistance per unit len intracellular space, in omega-cm  """
    r_i = 1/(math.pi*Cell_rad*Cell_rad*6.7)
    #r_i = 150/(math.pi*Cell_rad*Cell_rad)
    """ resistance per unit len extracellular space, in omega-cm """
    r_e = 1/(math.pi*(Cell_rad+0.00001)*(Cell_rad+0.00001)*20) 

    #resting potential = -70
    V_old= [-70 for j in range(J)]

    """ This is the D*delta(t)/delta(x)^2 """
    #D_v = dt/(2*C_m*Cell_peri*(r_i + r_e)*dx*dx)
    D_v = dt/((r_i + r_e)*dx*dx*Cell_peri)
    print(D_v)
    
    A_v = np.diagflat([-D_v/2 for i in range(J-1)], -1) + np.diagflat([1. + D_v for i in range(J)]) + np.diagflat([-1*D_v/2 for i in range(J-1)], 1)
    A_v[0][0] = 1+D_v/2
    A_v[-1][-1] = 1+D_v/2

    B_v = np.diagflat([ D_v/2 for i in range(J-1)], -1) + np.diagflat([1. - D_v for i in range(J)]) + np.diagflat([D_v/2 for i in range(J-1)], 1)
    B_v[0][0] = 1-D_v/2
    B_v[-1][-1] = 1-D_v/2


    # for i in A_v[0]:
    #     print(i)

    print(A_v)
    print(B_v)
   
    #change this value for the actual nth printing of the graph
    nplot = 1
    c=0

    #the Zeroth time step
    # 0 time step index
    # I_ion_old = hh.HodgkinHuxley().main(t, [V_old[0]], J, 0)
    # V_old= [-70 for j in range(J)]
    #print(V_old)
    val = 5
    V_old[0] -= val 
    V_old[1] -= val 
    V_old[2] -= val 
    V_old[-1] +=  val 
    V_old[-2] +=  val 
    V_old[-3] +=  val 

    dn = hh.HodgkinHuxley().n_inf(V_old)
    dm =hh.HodgkinHuxley().m_inf(V_old)
    dh = hh.HodgkinHuxley().h_inf(V_old)
    hh_arr= hh.HodgkinHuxley().main(V_old, dn, dm , dh)
    I_ion_old = hh_arr[3]
    # for now I_ion_old is [n, m, h]

    #TODO Dr.Lin notes:
    # If you notice the graphs the graph of n and m spike down at the ending this is due to the alpha_n and alpha_m returning 0 to prevent nan
    # fig = plt.figure()
    # #plt.plot(t, hh_arr[0], 'b', label='n')
    # #plt.plot(t, hh_arr[1], 'r', label='m')
    # #plt.plot(t, hh_arr[2], 'g', label='h')
    # plt.plot(t, hh_arr[3], '', label='I ion')
    # plt.ylabel('Current mA')
    # plt.xlabel('Time (ms)')
    # plt.legend()
    # plt.show()

    # Actually making
    for i in range(1,N):
        # Here i is the ith time step in the entire simulation
                    
        # turn off I ion 
        #I_ion_old = [0]*J
        dn = hh_arr[0]
        dm = hh_arr[1]
        dh = hh_arr[2]
        B_part = np.matmul(B_v, V_old) - I_ion_old
        V_new = np.linalg.solve(A_v, B_part)
        if(i%nplot==0): #plot results every nplot timesteps
            plt.plot(x_grid,V_new,linewidth=2)
            plt.ylim([-1000, 100])
            plt.xlim([-2,35])
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
        if i < 9:
            val = 5 
            V_old[0] -= val 
            V_old[1] -= val 
            V_old[2] -= val 
            V_old[-1] +=  val 
            V_old[-2] +=  val 
            V_old[-3] +=  val 
        print(V_new)
        #hh_arr = hh.HodgkinHuxley().main(t, V_old)
        hh_arr = hh.HodgkinHuxley().main(V_old, dn, dm, dh)
        I_ion_new = hh_arr[3]
        dn = hh_arr[0]
        dm = hh_arr[1]
        dh = hh_arr[2]
        I_ion_old = I_ion_new

os.system("ffmpeg -y -i 'foo%03d.jpg' cable_eqn.m4v")
os.system("rm -f *.jpg")

if __name__ == '__main__':
    main()