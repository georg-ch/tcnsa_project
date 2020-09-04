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
#c_pre = 1.3
c_post = 2
D = 20

### simulation ###
c_111, rho_111 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 0)
c_112, rho_112 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 0)
c_113, rho_113 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 0)
c_11avg, rho_11avg = synapse(t_end, dt, tau, rho_star, 0, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 0)

c_121, rho_121 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 1)
c_122, rho_122 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 1)
c_123, rho_123 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 1)
c_12avg, rho_12avg = synapse(t_end, dt, tau, rho_star, 0, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1hz_pre', 1)


c_211, rho_211 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 0)
c_212, rho_212 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 0)
c_213, rho_213 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 0)
c_21avg, rho_21avg = synapse(t_end, dt, tau, rho_star, 0, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 0)

c_221, rho_221 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 1)
c_222, rho_222 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 1)
c_223, rho_223 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 1)
c_22avg, rho_22avg = synapse(t_end, dt, tau, rho_star, 0, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1hz_pre', 1)

times = np.arange(0, t_end, dt)

### plotting ###
fig, axes = plt.subplots(1, 2, sharey=True, sharex=True, dpi=200)

#axes[0].axhline(theta_p, ls='--', color='yellow')
axes[0].plot(times, rho_111, color='k', alpha=0.5)
axes[0].plot(times, rho_112, color='k', alpha=0.5)
axes[0].plot(times, rho_113, color='k', alpha=0.5)
axes[0].plot(times, rho_11avg, color='blue')

axes[0].plot(times, rho_121, color='k', alpha=0.5)
axes[0].plot(times, rho_122, color='k', alpha=0.5)
axes[0].plot(times, rho_123, color='k', alpha=0.5)
axes[0].plot(times, rho_12avg, color='blue')

#axes[0].axhline(theta_d, ls='--', color='cyan')
axes[0].set_xlabel('time (s)')
axes[0].set_ylabel('synaptic efficacy $\\rho$')
axes[0].set_title('$C_{pre} = 1.3$')
axes[1].set_title('$C_{pre} = 3$')

axes[1].plot(times, rho_211, color='k', alpha=0.5)
axes[1].plot(times, rho_212, color='k', alpha=0.5)
axes[1].plot(times, rho_213, color='k', alpha=0.5)
axes[1].plot(times, rho_21avg, color='red')

axes[1].plot(times, rho_221, color='k', alpha=0.5)
axes[1].plot(times, rho_222, color='k', alpha=0.5)
axes[1].plot(times, rho_223, color='k', alpha=0.5)
axes[1].plot(times, rho_22avg, color='red')

axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)

axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['left'].set_visible(False)
axes[1].tick_params(left=False)

axes[1].set_xlabel('time (s)')

plt.tight_layout()

plt.show()