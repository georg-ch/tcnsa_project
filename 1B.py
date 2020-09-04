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
c_pre = 1.3
c_post = 2
D = 20

c_pre = np.linspace(0.001, 4, 100)

T = t_end

print(T)
print(t_end)

alpha_p = tau_ca * (np.log(c_pre/theta_p) > 0) * np.log(c_pre/theta_p)/T
alpha_d = tau_ca * (np.log(c_pre/theta_d) > 0) * np.log(c_pre/theta_d)/T

gp = gamma_p * alpha_p
gd = gamma_d * alpha_d
rho_bar = gp/(gp + gd)

fig, ax1 = plt.subplots(dpi=200)

ax1.plot(alpha_p)
ax1.plot(alpha_d)

ax1.set_ylim(-0.3, 0.3)

ax2 = ax1.twinx()

ax2.plot(rho_bar, color='k')

plt.show()