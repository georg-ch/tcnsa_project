from Synapse_nb import synapse
import matplotlib.pyplot as plt
import numpy as np

### parameters ###
t_end = 5
dt = 0.0001

tau = 346.3615
rho_star = 0.5
sigma = 3.3501

gamma_p = 725.085
gamma_d = 331.909
theta_p = 1.3
theta_d = 1

tau_ca = 0.0226936
c_pre = 0.5617539
c_post = 1.23964

D = 4.6098

w0 = 1
w1 = 5
rho_init = 0.5

c_pre_range = np.linspace(0, c_pre, 50)
c_post_range = np.linspace(0, c_post * 1.5, 50)

rho_mat = np.zeros((c_pre_range.shape[0], c_post_range.shape[0]))

for i, c_pre in enumerate(c_pre_range):
    for j, c_post in enumerate(c_post_range):
        c, rho = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, c_pre, c_post, D,
                         'poisson',
                         rho_init, pre_freq=30, post_freq=30, dt_spike=None, time_start=-0.05)
        rho_mat[i, j] = rho[-1]

fig, axes = plt.subplots(1, 2, dpi=200)
a = axes[1].imshow((w0 + rho_mat * (w1 - w0))/(w0 + rho_init*(w1 - w0)), origin='lower left', cmap='seismic', extent=[0, 1, 0, 1.5], aspect=1/1.5)
plt.colorbar(a)
axes[0].plot((w0 + rho_mat[int(c_pre_range.shape[0]*2/3), :] * (w1 - w0))/(w0 + rho_init*(w1 - w0)))
plt.show()