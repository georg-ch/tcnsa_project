import numpy as np
import numba as nb

@nb.jit(nopython=True)
def hs(x):
    '''
    numba-compatible implementation of heavyside function
    :param x: argument (float)
    :return: Theta(x)
    '''
    return x > 0

@nb.jit(nopython=True)
def synapse_wrapped(times, c, rho, eta, rho_star, gamma_p, gamma_d, theta_p, theta_d, theta_min, sigma, dt, tau, tau_sq, tau_ca, c_pre, spikes_pre, spikes_post, D_ind, c_post):
    '''
    numbified differential equation-solver loop. This function should only be called by synapse, and therefore there is
        no need to explain its parameters. It returns c(t), the calcium concentration and rho(t), the efficacy of
        the synapse
    '''
    for k, t in enumerate(times[:-1]):
        c[k + 1] = c[k] - (dt / tau_ca) * c[k] + c_pre * spikes_pre[k - D_ind] + c_post * spikes_post[k]
        rho[k + 1] = rho[k] + (dt / tau) * (-rho[k] * (1 - rho[k]) * (rho_star - rho[k])
                            + gamma_p * (1 - rho[k]) * hs(c[k] - theta_p)
                            - gamma_d * rho[k] * hs(c[k] - theta_d)
                            + sigma * tau_sq * np.sqrt(hs(c[k] - theta_d) + hs(c[k] - theta_p)) * eta[k])

    return c, rho

def synapse(t_end, dt, tau, rho_star, sigma, gamma_p, gamma_d, theta_p, theta_d, tau_ca, c_pre, c_post, D, spikes,
            rho_init, pre_freq=None, post_freq=None, dt_spike = None, time_start=-0.05, num_of_spike_pers = 0):
    '''
    Temporal evolution of synaptic efficacy and calcium concentration. This function creates some containers and
        calculates a few values that are passed to synapse_wrapped, which executes the loop within which
        the differential equations associated with the synapse are solved.
    :param t_end: ending time of simulation in seconds (int or float)
    :param dt: stepping increment of simulation (float)
    :param tau: time constant of synapse (int or float)
    :param rho_star: position of second well (int or float)
    :param sigma: scale of noise (float)
    :param gamma_p: potentiation constant (float)
    :param gamma_d: depression constant (float)
    :param theta_p: potentiation threshold (float)
    :param theha_d: depression threshold (float)
    :param tau_ca: time constant of calcium (float)
    :param c_pre: change in calcium in response to presynaptic spike (float)
    :param c_post: change in calcium in response to postsynaptic spike (float)
    :param D: delay of change in calcium in response to presynaptic spike
    :param spikes: way in which the pre- and postsynaptic spikes are generated. If '1hz_pre', the postsynaptic neuron
        is silent and the presynaptic neuron fires regularly with 1 Hz for 60 seconds. If 'poisson', pre_freq and
        post_freq must be set to values other than None. The neurons will fire poisson-like with the respective
        frequencies for the duration of the simulation
    :param rho_init: initial value of the synapse (float)
    :param pre_freq: frequency of presynaptic neuron id spikes is set to 'poisson'
    :param post_freq: frequency of postsynaptic neuron id spikes is set to 'poisson'
    :return: Tuple (c, rho) where c is the time-evolution of calcium and rho is the time-evolution of the synapse
    '''

    ### initializing containers and parameters ###
    D_ind = int(1 / 1000 * (D / dt))

    times = np.arange(0, t_end, dt)
    rho = np.zeros(times.shape[0])
    rho[0] = rho_init
    c = np.zeros(times.shape[0])

    tau_sq = np.sqrt(tau)
    theta_min = np.min([theta_p, theta_d])
    eta = np.random.normal(size=times.shape[0])
    spikes_pre = np.zeros(times.shape[0])
    spikes_post = np.zeros(times.shape[0])
    ### initializing pre- and postsynaptic spikes ###
    if spikes == '1hz_pre':

        pres_arr = np.zeros(int(60 * 1 / dt))
        for k, s in enumerate(pres_arr):
            if k % 10000 == 0:
                pres_arr[k] = 1
        spikes_pre[:int(60 * 1 / dt)] = pres_arr
    elif spikes == 'poisson' and pre_freq and post_freq:

        pre_prob = pre_freq * dt
        post_prob = post_freq * dt

        spikes_pre = np.random.choice([0, 1], size=(times.shape[0]), p=[1 - pre_prob, pre_prob])
        spikes_post = np.random.choice([0, 1], size=(times.shape[0]), p=[1 - post_prob, post_prob])

    elif type(dt_spike) == int and num_of_spike_pers > 0 and not pre_freq == None:
        times = np.arange(0, 1.0/pre_freq * num_of_spike_pers, dt)  # add 50 ms to simulation start
        rho = np.zeros(times.shape[0])
        rho[0] = rho_init
        c = np.zeros(times.shape[0])

        tau_sq = np.sqrt(tau)
        theta_min = np.min([theta_p, theta_d])
        eta = np.random.normal(size=times.shape[0])

        spikes_pre = np.zeros(times.shape[0])
        spikes_post = np.zeros(times.shape[0])
        delay_idx = dt_spike / 1000 / dt
        idx_freq = int(1.0 / pre_freq / dt) * 2
        spikes_pre[np.arange(0, len(spikes_pre), idx_freq).astype(int)] = 1
        spikes_post[np.arange(delay_idx, len(spikes_post), idx_freq).astype(int)] = 1

    elif type(dt_spike) == int:
        if dt_spike/1000.0 > t_end:
            return 0, 0 #this is a problem
        times = np.arange(time_start, t_end , dt) #add 50 ms to simulation start
        rho = np.zeros(times.shape[0])
        rho[0] = rho_init
        c = np.zeros(times.shape[0])

        tau_sq = np.sqrt(tau)
        theta_min = np.min([theta_p, theta_d])
        eta = np.random.normal(size=times.shape[0])

        spikes_pre = np.zeros(times.shape[0])
        spikes_post = np.zeros(times.shape[0])
        pre_spike_idx = int(abs(time_start) / dt) # pre time is in zero as convention
        post_spike_idx = int(abs(time_start) / dt + dt_spike/1000.0/dt)
        spikes_pre[pre_spike_idx] = 1
        spikes_post[post_spike_idx] = 1
    elif spikes == '1pre':
        spikes_pre[int(0.001/dt)] = 1
    else:
        return 0, 0

    ### run simulation ###
    return synapse_wrapped(times, c, rho, eta, rho_star, gamma_p, gamma_d, theta_p, theta_d, theta_min, sigma, dt, tau, tau_sq, tau_ca, c_pre, spikes_pre, spikes_post, D_ind, c_post)