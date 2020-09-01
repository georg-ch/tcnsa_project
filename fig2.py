from Synapse_nb import synapse
import matplotlib.pyplot as plt
import numpy as np

### parameters ###
t_end = 0.1
dt = 0.00001

tau = 150 #sec
rho_star = 0.5
sigma = 2.8284

gamma_p = 321.808
gamma_d = 200
theta_p = 1.3
theta_d = 1

tau_ca = 0.02 #in sec
c_pre = 1
c_post = 2
D = 13.7 #3 im ms

def fill_above_threshoold(times, c, thresh_p, thresh_d, ax):
    c_p= c.copy()
    c_p[c_p < thresh_p] = 0
    c_d= c.copy()
    c_d[c_d >= thresh_p] = 0
    c_d[c_d < thresh_d] = 0
    above_thresh_d_idx = np.where(c >= thresh_d)[0]
    ax.fill_between(times, c_d, [0]*len(c_d), color="b", alpha=0.25)
    ax.fill_between(times, c_p, [0]*len(c_p), color="orange", alpha=0.25)

def fig_A(time_to_check = 20):
    ### simulation ###
    c, rho = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d,
                     theta_p, theta_d, tau_ca, c_pre, c_post, D, '', 1,
                     dt_spike = -time_to_check)
    print("mid")
    c2, rho2 = synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d,
                     theta_p, theta_d, tau_ca, c_pre, c_post, D, '', 1,
                     dt_spike = time_to_check)

    times = np.arange(-0.05, t_end, dt)
    print(times.shape, c.shape)
    ### plotting ###
    fig, ax = plt.subplots(1, 2, dpi=200)
    ax[0].plot(times, c, label="-20", color="b")
    ax[1].plot(times, c2, label="+20", color="r")
    ax[0].axhline(theta_p, ls='--', color='orange')
    ax[0].axhline(theta_d, ls='--', color='blue')
    ax[1].axhline(theta_p, ls='--', color='orange')
    ax[1].axhline(theta_d, ls='--', color='blue')
    fill_above_threshoold(times, c, theta_p, theta_d, ax[0])
    fill_above_threshoold(times, c2, theta_p, theta_d, ax[1])


    ax[0].set_xlabel('time (s)')
    ax[1].set_xlabel('time (s)')
    ax[0].set_ylabel('calcium')
    ax[0].set_title("dt = -"+str(time_to_check)+" ms")
    ax[1].set_title("dt = "+str(time_to_check)+" ms")
    ax[0].set_xticks(np.arange(-0.05, t_end+dt, 0.05))
    ax[0].set_xticklabels(np.arange(-50, t_end*1000+1, 50).astype(int))
    ax[1].set_xticks(np.arange(-0.05, t_end+dt, 0.05))
    ax[1].set_xticklabels(np.arange(-50, t_end*1000+1, 50).astype(int))
    plt.savefig("fig2a")
    plt.close()

def fig_B(d_time_leg=2):
    area_d = []
    area_p = []
    min_P = []
    time_legs = np.arange(-100, 100, d_time_leg)
    theta_p=1.3
    dt=0.00001
    for dt_spike in time_legs:
        c, rho = synapse(max((abs(dt_spike)+50)/1000.0, 0.5), dt, tau, rho_star, sigma, gamma_p, gamma_d,
                         theta_p, theta_d, tau_ca, c_pre, c_post, D, '', 1,
                         pre_freq=0, post_freq=0, dt_spike=int(dt_spike), time_start=-max((abs(dt_spike)+50)/1000.0, 0.5))
        area_d.append((c>theta_d).sum()*dt)
        area_p.append((c>theta_p).sum()*dt)
        min_P.append(area_p[-1]*gamma_p/ (area_d[-1]*gamma_d + area_p[-1]*gamma_p))
    fig, ax1 = plt.subplots()
    ax1.plot(time_legs, area_d, color="b", label="a_d")
    ax1.plot(time_legs, area_p, color="orange", label="a_p")
    ax1.axvline(0, ls='--', color='grey')
    ax1.axhline(0, ls='--', color='grey')
    ax1.legend()
    ax2 = ax1.twinx()
    ax2.plot(time_legs, min_P, color="K", label="P_hat")
    ax2.set_ylim([0.4, 0.65])
    ax1.set_ylim([-0.02, 0.03])
    ax2.legend(loc='upper left')
    ax1.set_ylabel("fraction of time above threshold")
    ax1.set_xlabel("time leg delta t (ms)")
    ax2.set_ylabel("position of minimum P")

    fig.tight_layout()
    plt.savefig("fig2b")
    plt.close()

def fig_STDP_curves(d_time_leg=1):

    change_in_p = []
    min_P = []

    time_legs = np.arange(-100, 100, d_time_leg)
    dt=0.00001
    # elif type(dt_spike) == int and num_of_spike_pers > 0 and not pre_freq == None:

    for dt_spike in time_legs:
        c, rho = synapse(max((abs(dt_spike)+50)/1000.0, 0.5), dt, tau, rho_star, sigma, gamma_p, gamma_d,
                         theta_p, theta_d, tau_ca, c_pre, c_post, D, '', rho_init=1,
                         pre_freq=1, dt_spike=int(dt_spike), time_start=-max((abs(dt_spike)+50)/1000.0, 0.5), num_of_spike_pers=60)
        change_in_p.append(rho[-1]/rho[0])
        R_d = (c > theta_d).sum() * dt
        R_p =(c > theta_p).sum() * dt
        min_P.append(R_p/ (R_d+ R_p))
    fig, ax1 = plt.subplots()
    ax1.axvline(0, ls='--', color='grey')
    ax1.axhline(1, ls='--', color='grey')
    ax1.plot(time_legs, change_in_p, color="c", label="simulation")
    ax1.plot(time_legs, min_P, color="m", label="P_hat")
    # plt.savefig("fig2c")
    # plt.close()
    plt.show()
# fig_A(20)
# fig_B(1)
fig_STDP_curves()
