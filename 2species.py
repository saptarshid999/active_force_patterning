import numpy as np
import matplotlib.pyplot as plt
import dedalus.public as d3
import logging
import os
logger = logging.getLogger(__name__)



# Parameters
Lx= 10
timestepper = d3.RK443
dtype = np.float64
dealias=4/2
N=256




# Create bases and domain
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=dtype)
xbasis = d3.RealFourier(xcoord, size=N, bounds=(0, Lx), dealias=dealias)


#Fields
u = dist.Field(name = 'u',bases=xbasis)
v = dist.Field(name = 'v',bases=xbasis)
r = dist.Field(name='r', bases=xbasis)
p = dist.Field(name='p', bases=xbasis)


dx = lambda A: d3.Differentiate(A, xcoord)



# 2D hydrodynamics
problem = d3.IVP([u,v,r,p],namespace=locals())
#parameters


b=6
c=1
zavg=2
zrel=2


d_0=1.0
eta=1.0
k=0.1

ku0=1
kbavg=1
kbrel=0
k3=2
k4=2
alpha1=k3+k4
alpha2=k3-k4




chi0=1.0
chi1=1.0
chi2=0.5
chi3=10



eps= dx(u)


sigma_lin=  (b-c**2)*eps

sigma_nonlin=2*chi0*(zavg*r+zrel*p) - 2*chi1*eps*(zavg*r+zrel*p)   +  chi2*((eps)**2)*(zavg*r+zrel*p) +  (chi3/3)*((eps)**3)*(zavg*r+zrel*p)


sigma=sigma_lin + sigma_nonlin

force_lin=dx(sigma_lin)
force_bs=dx(2*chi0*(zavg*r+zrel*p))
force_nonlin=dx(sigma_nonlin-2*chi0*(zavg*r+zrel*p))
force_total=dx(sigma)


problem.add_equation("dt(u)  - dx(sigma_lin) +k*dx(dx(dx(eps))) - dx(eta*dt(eps))  =   dx(sigma_nonlin) ")

problem.add_equation("v - dt(u) = 0")

problem.add_equation("dt(r) - d_0*dx(dx(r))   =   - dx(v*r)  +kbavg*(1-c*eps) - ((ku0*np.exp(alpha1*eps) + ku0*np.exp(alpha2*eps))/2)*r - ((ku0*np.exp(alpha1*eps) - ku0*np.exp(alpha2*eps))/2)*p ")

problem.add_equation("dt(p) - d_0*dx(dx(p))   =   - dx(v*p)  +kbrel*(1-c*eps) - ((ku0*np.exp(alpha1*eps) + ku0*np.exp(alpha2*eps))/2)*p - ((ku0*np.exp(alpha1*eps) - ku0*np.exp(alpha2*eps))/2)*r ")






x = dist.local_grid(xbasis)

# Build solver
solver = problem.build_solver(timestepper)


# Initial conditions

u['g'] = 0.0+ 0.0*(np.random.random(N)-0.5)
r['g'] = 1+0.0*(np.random.random(N)-0.5)
p['g'] = 0+0.1*(np.random.random(N)-0.5)







dt= 0.001

step=dt

# Integration parameters
solver.stop_sim_time = 10000
solver.stop_wall_time = np.inf
solver.stop_iteration = np.inf



# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.0001, max_writes=10000000)
snapshots.add_task(r+p, name='rho1')
snapshots.add_task(r-p, name='rho2')
snapshots.add_task(dx(u), name='strain')
snapshots.add_task(dx(u)*dx(u), name='strain2')
snapshots.add_task(u, name='u')
snapshots.add_task(v, name='v')
snapshots.add_task(force_total, name='force')



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
