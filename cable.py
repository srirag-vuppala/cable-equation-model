# imports
import numpy as np
import math
import matplotlib.pyplot as plt
import os 
import hh 


def main():
    """ An implementation of the CableEquation where the Hodgkin Huxley equations
    describe the nature of Ionic current """ 

    """ Defining our constants and assumptions """
    # Total duration of the simulation (in sec)
    T = 1.5
    # The number of time steps that we track for 
    N = 1000
    # The time that passes by between each time step 
    dt = float(T)/float(N)


    # Cell length (in cm)
    Cell_len = 0.01
    # Cell radius is width/2 (in cm)
    # width = Cell length/6 [assumption: the cell would be roughly 6:1 | length:width] 
    Cell_rad = .00241/2 
    # Cell perimeter = perimeter of circle (in cm)
    Cell_peri = 2*math.pi*Cell_rad
    
    # A node is point in space where we observe the action potential
    n_cells = 10 # number of nodes per cell
    J = n_cells*10 # total number of nodes
    # The seperation between two nodes (in cm)
    dx = float(Cell_len*n_cells)/float(J)
    # Every node that we observe in an array
    x_grid = np.array([j*dx for j in range(J)])


    # membrane capacitance, (in uF/cm^2)
    C_m = 1

    # resistance per unit len intracellular space, (in omega-cm)
    # r_i = 1/(math.pi*Cell_rad*Cell_rad*6.7)
    # resistance per unit len extracellular space, (in omega-cm) 
    # r_e = 1/(math.pi*(Cell_rad+0.00001)*(Cell_rad+0.00001)*20) 
    
    # We use the constant values below for the simulation but they are governed as shown above
    r_i = 178330.0
    r_e = 888740.0

    # A array of transmembrane potentials that exists at each node (in mV)
    # resting potential = -70
    V_old= [-70 for _ in range(J)]

    # Laplacian
    # The D_v is the constant that gets multiplied to the laplacian matrix
    D_v = dt/((r_i + r_e)*dx*dx*Cell_peri*C_m)
    
    # Constructing the Laplacian for the (n+1)th time step
    # The + and - cover the identity matrix that is expressed in our final equation
    A_v = np.diagflat([-D_v/2 for _ in range(J-1)], -1) + np.diagflat([1. + D_v for _ in range(J)]) + np.diagflat([-1*D_v/2 for _ in range(J-1)], 1)
    # Defining the boundary conditions
    # Note there are only two boundaries for the cable equation because there are only two ends for a cable
    A_v[0][0] = 1+D_v/2
    A_v[-1][-1] = 1+D_v/2

    # Constructing the Laplacian for the (n)th time step
    B_v = np.diagflat([ D_v/2 for _ in range(J-1)], -1) + np.diagflat([1. - D_v for _ in range(J)]) + np.diagflat([D_v/2 for _ in range(J-1)], 1)
    # Defining the boundary conditions
    B_v[0][0] = 1-D_v/2
    B_v[-1][-1] = 1-D_v/2

    #change this value for the actual nth printing of the graph
    nplot = 10
    c = 0

    # injected current
    val = 30

    # Actually making
    for i in range(N):
        # Continually inject current for the first 200 timesteps
        # This will serve as our excitation that we artificially give to spark the propagation of the action potential
        if i < 200:
            V_old[0]  = -70+val 
            V_old[1]  = -70+val 
            V_old[2]  = -70+val 
            V_old[-1] = -70-val 
            V_old[-2] = -70-val 
            V_old[-3] = -70-val
        # Here i is the ith time step in the entire simulation
        I_ion_old = hh.HodgkinHuxley().main([j+70.0 for j in V_old])
        B_part = np.matmul(B_v, V_old) - dt*np.asarray(I_ion_old)
        V_new = np.linalg.solve(A_v, B_part)
        
        if i % nplot == 0: #plot results every nplot timesteps
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
        # setup for the next iteration  
        V_old = V_new
         
# snapshots the operating system takes of the plotted graphs to make it into a video
os.system("ffmpeg -y -i 'foo%03d.jpg' cable_eqn.m4v")
os.system("rm -f *.jpg")

if __name__ == '__main__':
    main()