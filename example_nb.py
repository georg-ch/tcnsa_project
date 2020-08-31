from Synapse_nb import synapse
import matplotlib.pyplot as plt
import numpy as np

### parameters ###
t_end = 120
dt = 0.0001

tau = 150
rho_star = 0.5
sigma = 70#2.8284

gamma_p = 321.808
gamma_d = 200
theta_p = 1.3
theta_d = 1

tau_ca = 0.02
c_pre = 1.3
c_post = 2
D = 13.7

### simulation ###
c, rho = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, c_pre, c_post, D, '1hz_pre', 1, pre_freq=60, post_freq=120)

times = np.arange(0, t_end, dt)

### plotting ###
fig, ax = plt.subplots(1, 1, dpi=200)
ax.plot(times, rho)
ax.axhline(1, ls='--', color='k')
ax.axhline(rho_star, ls='--', color='gray')
ax.axhline(0, ls='--', color='k')

ax.set_xlabel('time (s)')
ax.set_ylabel('synaptic efficacly $\\rho$')

plt.show()