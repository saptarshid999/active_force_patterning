"""Dedalus script for solving the system of equations (Eq. 1a, 1b, 1c: as written down in the manuscript).

To run the script, install the latest version of Dedalus in the system (details can be found at https://dedalus-project.org/)

To compile and run this code:

$ conda activate dedalus3

$ python3 stresslets_segregation.py


"""

import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
import os
logger = logging.getLogger(__name__)



# Parameters
Lx = 35
#Nx = 384
Nx = 500
dealias = 1
#stop_sim_time = 40
stop_sim_time = 10
timestepper = d3.SBDF2 # d3.RK443
timestep = 1e-4
dtype = np.float64
b=10
c=1
d_0=1
k1=1
k2=0
k5=0
k3=1
zeta1=3
zeta2=2
k4=2
chi=1
chip=0.01
chipp=0.001
#chippp=-0.001
chippp=0.001
eta = 0.001
a = 1

# Bases
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=dtype)
xbasis = d3.RealFourier(xcoord, size=Nx, bounds=(0, Lx), dealias=dealias)

# Fields
u = dist.Field(name='u', bases=xbasis)
v = dist.Field(name='v', bases=xbasis)
r = dist.Field(name='r', bases=xbasis)
p = dist.Field(name='p', bases=xbasis)
e = dist.Field(name='e', bases=xbasis)

# Substitutions
dx = lambda A: d3.Differentiate(A, xcoord)

# Problem
problem = d3.IVP([u, v, r, p, e], namespace=locals())


#system of equations

problem.add_equation("dt(u) - eta*dx(dx((dt(u)))) - dx(2*chi*(zeta1*r + zeta2*p))  =  dx((b - (c*c/a)- 2*chip*(c/a)*(zeta1*r + zeta2*p))*dx(u)) + dx((chipp*((c/a)**2)*(zeta1*r + zeta2*p))*((dx(u))**2)) + dx(((-chippp/3)*((c/a)**3)*(zeta1*r + zeta2*p))*((dx(u))**3))")
problem.add_equation("v - dt(u) = 0") 
problem.add_equation("dt(r) - d_0*dx(dx((r))) = - dx(v*r) + 1-(c/a)*dx(u) - (k1 + k3*dx(u))*r  - (k2 + k4*dx(u))*p")              
problem.add_equation("dt(p) - d_0*dx(dx((p))) = - dx(v*p) + k5*(1-(c/a)*dx(u)) -  (k1 + k3*dx(u))*p  - (k2 + k4*dx(u))*r")
problem.add_equation("e - dx(u) = 0 ")


# Initial conditions
x = dist.local_grid(xbasis)
#n = 20
n = 40
u['g'] = 0
v['g'] = 0
e['g'] = 0
r['g'] = 1
np.random.seed(5)
p['g'] = 0.2*np.random.random((500,))
#plt.plot(x,p['g'])
#plt.show()
# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time

#array1= np.zeros((500,2))
array1= np.zeros((500,3))

# Main loop
u.change_scales(1)
r.change_scales(1)
p.change_scales(1)
u_list = [np.copy(u['g'])]
r_list = [np.copy(r['g'])]
p_list = [np.copy(p['g'])]
t_list = [solver.sim_time]


ff = open("max_phi.txt", "w")


while solver.proceed:
    solver.step(timestep)
    try:
        os.stat('segre')
    except:
        os.mkdir('segre')
    
    ttime = round(solver.sim_time,7)
    print(ttime,"\t",np.amax(p['g']), file=ff)
    
    if solver.iteration % 100 == 0:
        logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
    if solver.iteration % 100 == 0:
        p.change_scales(1)
        p_list.append(np.copy(p['g']))
        r.change_scales(1)
        r_list.append(np.copy(r['g']))
        t_list.append(solver.sim_time)
        
        array1[:,0]= np.ones(500)*(solver.sim_time)
        array1[:,1] = x.ravel()
        array1[:,2] = p['g']

        #for printing snapshots of \phi and \epsilon at different time points
        
        if solver.iteration >= 40500:
        	np.savetxt("phi%.7d.txt"%solver.iteration, array1[320:420,:], delimiter = '\t')
        
        #This section is for plotting. This code generates images which are stiched together to form movies. 
        	
        figure,axis = plt.subplots(2, 1)
        axis[0].plot(x.ravel(),p['g'])
        axis[0].set_title("phi profile")
        axis[1].plot(x.ravel(),e['g'])
        axis[1].set_title("strain profile")
        plt.savefig("segre/phi%0.7d.png"%solver.iteration, dpi = 200)
        plt.clf()
        
        


ff.close()

# Plot

plt.figure(figsize=(6, 4))
plt.pcolormesh(x.ravel(), np.array(t_list), np.array(p_list), cmap='RdBu_r', shading='gouraud', rasterized=True)
plt.xlim(0, Lx)
plt.ylim(0, stop_sim_time)
plt.colorbar()
#plt.clf()
plt.xlabel('x')
plt.ylabel('t')
#plt.title(f'KdV-Burgers, (a,b)=({a},{b})')
plt.tight_layout()
plt.savefig('1dsegre.pdf')
plt.savefig('1dsegre.png', dpi=200)
