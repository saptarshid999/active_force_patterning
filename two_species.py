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
timestepper = d3.SBDF1
timestep = 1e-6
step= 1e-6
dtype = np.float64
b=14
c=1
d_0=1.2 #diffusion
k1=1

#zeta1_zeta2
zeta1_num=3.25
zeta2_num=1

#zeta_avg and zeta_rel
zeta1= (zeta1_num+zeta2_num)/2
zeta2=(zeta1_num-zeta2_num)/2


chi=1
chip=0.5
chipp=0.05
chippp=-2
eta = 0.1 #viscosity
a = 1
#alpha1=3
alpha1=7
alpha2=5
ku0=1


# Bases
xcoord = d3.Coordinate('x')
dist = d3.Distributor(xcoord, dtype=dtype)
xbasis = d3.RealFourier(xcoord, size=Nx, bounds=(0, Lx), dealias=dealias)

# Fields
u = dist.Field(name='u', bases=xbasis)
v = dist.Field(name='v', bases=xbasis)
r = dist.Field(name='r', bases=xbasis)
p = dist.Field(name='p', bases=xbasis)


# Substitutions
dx = lambda A: d3.Differentiate(A, xcoord)

ux= dx(u)
uxx= dx(ux)

rx= dx(r)
rxx= dx(rx)
rxxx=dx(rxx)
rxxxx=dx(rxxx)

px= dx(p)
pxx= dx(px)
pxxx=dx(pxx)
pxxxx=dx(pxxx)

vx= dx(v)
vxx= dx(vx)

motiff = ((zeta1*r) + (zeta2*p))
d_motiff= ((zeta1*rx) + (zeta2*px))

T1_lin= 2*chi*d_motiff
T2_lin= (b-1)*uxx
T2_nonlin= (-2*chip*motiff*uxx) + (-2*chip*d_motiff*ux)
T3_nonlin= (chipp*motiff*2*ux*uxx) + (chipp*d_motiff*ux*ux)
T4_nonlin= (-chippp*motiff*ux*ux*uxx) + ((-chippp/3)*d_motiff*ux*ux*ux)
T5_lin= eta*vxx

#unbinding rates
ku_avg= (np.exp(alpha1*ux) +  np.exp(alpha2*ux))*ku0/2
ku_rel= (np.exp(alpha1*ux) -  np.exp(alpha2*ux))*ku0/2


full_stress= (2*chi*motiff) + ((b-1-(2*chip*motiff))*ux) + (chipp*motiff*(ux**2)) + ((-chippp/3)*motiff*(ux**3)) + (eta*vx)


k=1

regular= dx(dx(k*uxx))


safety= (1+np.tanh(100*(r-0.2)))/2


# Problem
problem = d3.IVP([u, v, r, p], namespace=locals())


#system of equations

problem.add_equation("dt(u) - T1_lin - T2_lin - T5_lin + regular = T2_nonlin + T3_nonlin + T4_nonlin ")
problem.add_equation("v - dt(u) = 0") 
problem.add_equation("dt(r) - d_0*rxx = - v*rx - vx*r + (safety*((1-ux) - (ku_avg)*r  - (ku_rel)*p))")              
problem.add_equation("dt(p) - d_0*pxx = - v*px - vx*p + (safety* (-  (ku_avg)*p  - (ku_rel)*r))")




# Initial conditions

x = dist.local_grid(xbasis)
#n = 40
u['g'] = 0
v['g'] = 0
r['g'] = 1
p.fill_random('g', seed=234, distribution='normal', scale=1e-01) # Random noise


#Solver
solver = problem.build_solver(timestepper)
solver.stop_sim_time = stop_sim_time


# Main loop
u.change_scales(1)
r.change_scales(1)
p.change_scales(1)
v.change_scales(1)



# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.1, max_writes=10000000)
snapshots.add_task(r+p, name='rho1')
snapshots.add_task(r-p, name='rho2')
snapshots.add_task(dx(u), name='strain')
snapshots.add_task(u, name='u')
snapshots.add_task(full_stress, name='stress')
snapshots.add_task((dx(dx(u))), name='strain_dx')
snapshots.add_task(dx((dx(dx(u)))), name='strain_ddx')



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
