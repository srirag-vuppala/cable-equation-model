import numpy as np
import hh_read as hh
import math
import matplotlib.pyplot as plt
import os, sys
from scipy.linalg import solve_banded


np.set_printoptions(precision=3)

def Input_stimuli(t):
    """ Current applied to create stimulus which is dependent on time, in milli Ampere(A)/cm^2 """
    # if 120.0 < t < 132.0:
    #     return 150.0
    # elif 5.0 < t < 10.0:
    #     return 1.0
    return 0.0

T = 10 
N = 100
dt = float(T)/float(N)
t = np.array([n*dt for n in range(N)])

""" Creating a one dimensional domain with unit length and J=100 equally spaced grid points """
L = 10
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
    D_v = dt/((r_i + r_e)*dx*dx)
    
    A_v = np.diagflat([-D_v/2 for i in range(J-1)], -1) + np.diagflat([1.+D_v]+[1. + D_v for i in range(J-2)]+[1.+D_v]) + np.diagflat([-1*D_v/2 for i in range(J-1)], 1)
    B_v = np.diagflat([ D_v/2 for i in range(J-1)], -1) + np.diagflat([1.-D_v]+[1. - D_v for i in range(J-2)]+[1.-D_v]) + np.diagflat([D_v/2 for i in range(J-1)], 1)

    Identity = np.identity(J)

    #change this value for the actual nth printing of the graph
    nplot = 1
    c=0

    #the Zeroth time step
    # 0 time step index
    #I_ion_old = hh.HodgkinHuxley().main(t, [V_old[0]], 0)
    I_ion_old = hh.HodgkinHuxley().main(t, V_old, 0)
    #I_ion_old = [i*(dt/C_m) for i in I_ion_old] 

    for i in range(1,N):
        # Here i is the ith time step in the entire simulation
        # if i < 10:
        #     val = 5 
        #     V_old[0] -= val 
        #     V_old[1] -= val 
        #     V_old[2] -= val 
        #     V_old[-1] +=  val 
        #     V_old[-2] +=  val 
        #     V_old[-3] +=  val 
        # turn off I ion 
        #I_ion_old = [0]*J
        B_part = np.matmul(B_v, V_old) - I_ion_old
        V_new = np.linalg.solve(A_v, B_part)
        if(i%nplot==0): #plot results every nplot timesteps
            plt.plot(x_grid,V_new,linewidth=2)
            plt.ylim([-200, 100])
            plt.xlim([-2,15])
            filename = 'foo' + str(c+1).zfill(3) + '.jpg'
            plt.xlabel("x in cm")
            plt.ylabel("Transmembrane pot in mv")
            plt.axhline(0, color='black')
            plt.axvline(-0, color='black')
            plt.title("t = %2.2f"%(dt*(i+1)))
            plt.savefig(filename)
            plt.clf()
            c += 1
        I_ion_new = hh.HodgkinHuxley().main(t, V_old, i)
        V_old = V_new
       # I_ion_new = [i*(dt/C_m) for i in I_ion_new]
        I_ion_old = I_ion_new

os.system("ffmpeg -y -i 'foo%03d.jpg' cable_eqn.m4v")
os.system("rm -f *.jpg")

if __name__ == '__main__':
    main()