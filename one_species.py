#To run the code, first activate the dedalus environment by typing 'conda activate dedalus3' in the terminal followed by typing 'python3 one_species.py'.
#The code will run, generating HDF5 files which will be saved in the folder 'snapshots'.
#The ipynb notebook analysis.ipynb can be used to analyse the contents of the folder snapshots and make movies.

import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
import os
logger = logging.getLogger(__name__)


# Parameters
Lx = 128
Nx = 1024
dealias = 2
stop_sim_time = 1000
timestepper = d3.SBDF2
timestep = 1e-6
step= timestep
dtype = np.float64
b=6
c=1
d_0=1.3
zeta1=3.5
chi=1
chip=0.5
chipp=0.05
chippp=-2
a = 1
alpha=2.5
ku=1
eta=0.1
k=1

# Bases
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=dtype)
xbasis = d3.RealFourier(xcoord, size=Nx, bounds=(0, Lx), dealias=dealias)

# Fields
u = dist.Field(name='u', bases=xbasis)
v = dist.Field(name='v', bases=xbasis)
r = dist.Field(name='r', bases=xbasis)


# Substitutions
dx = lambda A: d3.Differentiate(A, xcoord)

ux= dx(u)
uxx=dx(ux)
k_u= ku*np.exp(alpha*ux)
regular= dx(dx(k*uxx))

# Problem
problem = d3.IVP([u, v, r], namespace=locals())



problem.add_equation("dt(u) - eta*dx(dx((dt(u)))) - dx(chi*(zeta1*r)) + regular = dx((b - (c*c/a)- chip*(c/a)*(zeta1*r))*dx(u)) + dx(  (chipp*0.5*((c/a)**2)*(zeta1*r)*((dx(u))**2))) + dx(((-chippp/6)*((c/a)**3)*(zeta1*r))*((dx(u))**3))")
problem.add_equation("v - dt(u) = 0") 
problem.add_equation("dt(r) - d_0*dx(dx((r))) = - dx(v*r) + 1-(c/a)*dx(u) - k_u*r")


# Initial conditions
x = dist.local_grid(xbasis)
#n = 20
n = 40
u['g'] = 0
v['g'] = 0
r.fill_random('g', seed=234, distribution='normal', scale=1e-01) # Random noise
r['g']+=1


# Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time



# Main loop
u.change_scales(1)
r.change_scales(1)
v.change_scales(1)




# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1, max_writes=10000000)
snapshots.add_task(r, name='rho')
snapshots.add_task(dx(u), name='strain')
snapshots.add_task(u, name='u')
snapshots.add_task(v, name='v')


try:
    logger.info('Starting main loop')
    while solver.proceed:
        timestep = step
        solver.step(timestep)
        if (solver.iteration % 5000 == 0):
            logger.info('Iteration=%i, Time=%e, dt=%e' %(solver.iteration, solver.sim_time, timestep))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    solver.log_stats()


