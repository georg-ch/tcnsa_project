from Synapse_nb import synapse
import matplotlib.pyplot as plt
import numpy as np

### parameters ###
t_end = 0.1
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
c_1, rho_1 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 1.3, c_post, D, '1pre', 1)
c_2, rho_2 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, 3, c_post, D, '1pre', 1)

times = np.arange(0, t_end, dt)

### plotting ###
fig, axes = plt.subplots(1, 2, sharey=True, sharex=True, dpi=200)

axes[0].axhline(theta_p, ls='--', color='yellow')
axes[0].plot(times*1000, c_1, color='blue')
axes[0].axhline(theta_d, ls='--', color='cyan')
axes[0].set_xlabel('times (ms)')
axes[0].set_ylabel('calcium')
axes[0].set_title('$C_{pre} = 1.3$')

axes[1].axhline(theta_p, ls='--', color='yellow')
axes[1].axhline(theta_d, ls='--', color='cyan')
axes[1].plot(times*1000, c_2, color='red')
axes[1].text(t_end * 1100, theta_p, '$\\theta_p$', color='yellow')
axes[1].text(t_end * 1100, theta_d, '$\\theta_d$', color='cyan')
axes[1].set_xlabel('times (ms)')
axes[1].set_title('$C_{pre} = 3$')

axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)

axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
axes[1].spines['left'].set_visible(False)
axes[1].tick_params(left=False)

#plt.tight_layout()

plt.show()