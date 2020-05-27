import linf
import plotting
import os
import math
from scipy import signal
import numpy as np


def calc_psd(N_e, N_i, nu_e, nu_i, mu_e, mu_i, sigma_e, sigma_i, tau_e, tau_i, g_l, E_e, E_i, c_m, store):
    V_avg = np.mean(np.asarray(store['voltage']))

    Ge_avg = np.mean(np.asarray(store['ge']))

    Gi_avg = np.mean(np.asarray(store['gi']))

    r_eff = 1/(g_l + Ge_avg + Gi_avg)
    tau_eff = r_eff * c_m

    spec = []
    for f in range(20, 1000, 20):
        spec.append((
        (r_eff**2/(1+(2*math.pi*f*tau_eff)**2))*
        (((2*N_e*nu_e*(tau_e**2)*(mu_e**2 + sigma_e**2)*(E_e - V_avg)**2)/(1+(2*math.pi*f*tau_e)**2))+
        ((2*N_i*nu_i*(tau_i**2)*(mu_i**2 + sigma_i**2)*(E_i - V_avg)**2)/(1+(2*math.pi*f*tau_i)**2)))
        ) * 1000)
    return spec

def calc_psd_comp(N_e, N_i, nu_e, nu_i, mu_e, mu_i, sigma_e, sigma_i, tau_e, tau_i, g_l, E_e, E_i, c_m, store):
    V_avg = np.mean(np.asarray(store['voltage']))

    Ge_avg = np.mean(np.asarray(store['ge']))

    Gi_avg = np.mean(np.asarray(store['gi']))

    r_eff = 1/(g_l + Ge_avg + Gi_avg)
    tau_eff = r_eff * c_m

    spec = []
    excit = []
    inhib = []
    for f in range(20, 1000, 20):
        spec.append((
        (r_eff**2/(1+(2*math.pi*f*tau_eff)**2))*
        (((2*N_e*nu_e*(tau_e**2)*(mu_e**2 + sigma_e**2)*(E_e - V_avg)**2)/(1+(2*math.pi*f*tau_e)**2))+
        ((2*N_i*nu_i*(tau_i**2)*(mu_i**2 + sigma_i**2)*(E_i - V_avg)**2)/(1+(2*math.pi*f*tau_i)**2)))
        ) * 1000)

        excit.append((r_eff**2/(1+(2*math.pi*f*tau_eff)**2))*
        ((2*N_e*nu_e*(tau_e**2)*(mu_e**2 + sigma_e**2)*(E_e - V_avg)**2)/(1+(2*math.pi*f*tau_e)**2)) * 1000)

        inhib.append((r_eff**2/(1+(2*math.pi*f*tau_eff)**2))*
        ((2*N_i*nu_i*(tau_i**2)*(mu_i**2 + sigma_i**2)*(E_i - V_avg)**2)/(1+(2*math.pi*f*tau_i)**2)) * 1000)
    return spec, excit, inhib

def spectrogram(trace, sampling):
    freqs, times, spectrogram = signal.spectrogram(trace, fs=sampling)
    i = 0
    while os.path.exists("graphs/spectro_%s.png" % i):
        i += 1
    plotting.plot(spectrogram, None, "Spectrogram", "Time Window", "Frequency Band", "graphs/spectro_" + str(i) + ".png" )

def psd(trace, sampling, bar):
    #window = range(0, 1000, 1)
    freqs, psd = signal.welch(trace, fs=sampling, scaling="spectrum")
    freqs1 = freqs
    psd1 = psd
    freqs = freqs[freqs <= bar]
    psd = psd[:freqs.size]

    i = 0
    while os.path.exists("graphs/psd_%s.png" % i):
        i += 1
    plotting.log_plot(freqs1, psd1, "Power Spectrum", "Log Frequency/ Hz", "Log Power/ $mV^2Hz^{-1}$", "graphs/psd_" + str(i) + ".png" )
    #plotting.plot_width(freqs, psd, "Power Spectrum", "Frequency", "Power", "graphs/psd_nonlog_" + str(i) + ".png", 1.5 )

def run(plot, *args):
    #MILLI
    time = 1000
    vars = 'default'

    if(vars == 'default'):
        #MILLI volts
        E_e = 0
        E_i = -62
        E_l = -90

        tau_m = 25
        tau_e = 3
        tau_i = 10

        #In  KHz
        exc_rate = 0.08
        inh_rate = 0.02
        if(len(args) == 2):
            print(str(len(args)))
            exc_rate = args[0]
            inh_rate = args[1]

        exc = 800
        inh = 200

        #MICRO
        mu_e = 0.00075
        mu_i = 0.00075
        sigma_e = 0.00075
        sigma_i = 0.00075
        #MEGA
        r_m = 100
        g_l = 1/r_m
        #NANO
        c_m = 0.09
    elif(vars == 'paulo'):
        E_e = 0
        E_i = -75
        E_l = -62

        tau_e = 3
        tau_i = 10

        exc_rate = 8
        inh_rate = 1.5
        exc = 1
        inh = 1

        mu_e = 0.06
        mu_i = mu_e
        sigma_e = 0.06
        sigma_i = sigma_e

        g_l = 0.005
        r_m = 1/g_l
        c_m = 0.1
        tau_m = c_m * r_m


    #MILLI
    delta_t = 0.05
    sampling = (1/delta_t) *1000

    timesteps = int(time/delta_t)
    # Create a LINF neuron object with parameters
    neuron = linf.Linf_neuron(exc_rate, inh_rate, exc, inh, mu_e, mu_i,
                        sigma_e, sigma_i, r_m, c_m, delta_t, time, E_e,
                        E_i, E_l, tau_m, tau_e, tau_i)
    # Provide debug information on amplitude distribution
    neuron.test_dist()
    # For each timestep in run, step forward
    for i in range(timesteps):
        neuron.step()
    #spike_train = neuron.get_spikes()

    #for i, spike in enumerate(spike_train):
    #    if spike > time:
    #        del(spike_train[i])

    # Store contains variables pertinent to run, e.g trace, average exc conductance
    store = neuron.get_store()

    i = 0
    if not os.path.exists("graphs"):
        os.makedirs("graphs")
    while os.path.exists("graphs/mem_v_%s.png" % i):
        i += 1
    if(plot):
        plotting.plot_ticks(store['voltage'], None, "Membrane Voltage", "time/ms", "Voltage/mV", "graphs/mem_v_" + str(i) + ".png")
        plotting.plot_ticks(store['ge'], None, "Excitatory Conductance", "time/ms", "Conductance/nS", "graphs/g_e_" + str(i) + ".png")
        plotting.plot_ticks(store['gi'], None, "Inhibitory Conductance", "time/ms", "Conductance/nS", "graphs/g_i_" + str(i) + ".png")
        #plotting.event(spike_train, "Sample Spike Train", "timestep", "Events", "graphs/spike_train_" + str(i) + ".png")
        #spectrogram(np.asarray(store['voltage']), sampling)
        psd(np.asarray(store['voltage']), sampling, 1000)

        nu_e = exc_rate
        nu_i = inh_rate

        N_e = exc
        N_i = inh
        mu_e = 0.00075
        mu_i = 0.00075

        #freqs, psda = signal.welch(store['voltage'], fs=sampling, scaling="spectrum")
        #spec = calc_psd(N_e, N_i, nu_e, nu_i, mu_e, mu_i, sigma_e, sigma_i, tau_e, tau_i, g_l, E_e, E_i, c_m, store)
        #spec, excit, inhib = calc_psd_comp(N_e, N_i, nu_e, nu_i, mu_e, mu_i, sigma_e, sigma_i, tau_e, tau_i, g_l, E_e, E_i, c_m, store)
        #plotting.paulo_valid_plot(freqs, psda, spec, excit, inhib)

        #plotting.log_plot(range(20,1000,20), spec, "Calculated Power Spectrum", "Log Frequency / Hz", "Log Power / mV", "graphs/psd_calc" + str(i) + ".png" )
        #plotting.plot_width(range(20,1000,20), spec, "Calculated Power Spectrum", "Frequency", "Power", "graphs/psd_calc_nonlog_" + str(i) + ".png", 1.5 )


    return store

# No plotting version that accepts arguments for CMA-ES
def run_eval(params, args):

    # Support for variable number of optimisation parameters
    exc_rate = params[0]
    inh_rate = params[1]
    if(params.size >= 4):
        mu_e     = params[2]
        mu_i     = params[3]
        if(params.size >= 6):
            sigma_e  = params[4]
            sigma_i  = params[5]
            if(params.size >= 8):
                tau_e = params[6]
                tau_i = params[7]

    #power = args[0]
    #tau_eff = args[1]
    r_m   = args[4]
    #V_avg   = args[5]
    E_e     = args[6]
    E_i     = args[7]
    exc     = args[8]
    inh     = args[9]
    #ge_avg  = args[10]
    #gi_avg  = args[11]
    g_l      = args[12]
    E_l     = args[13]
    sampling = args[14]
    tau_m = args[15]
    if(params.size <= 6):
        tau_e   = args[2]
        tau_i   = args[3]
        if(params.size <= 4):
            sigma_e = args[18]
            sigma_i = args[19]
            if(params.size <= 2):
                mu_e = args[16]
                mu_i = args[17]
    #MILLI
    time = 1000
    #Micro
    mu_e = mu_e
    mu_i = mu_i

    c_m = 0.09
    #MILLI
    delta_t = 1000 / sampling
    timesteps = int(time/delta_t)
    neuron = linf.Linf_neuron(exc_rate, inh_rate, exc, inh, mu_e, mu_i,
                            sigma_e, sigma_i, r_m, c_m, delta_t, time,
                            E_e, E_i, E_l, tau_m, tau_e, tau_i)
    for i in range(timesteps):
        neuron.step()
    store = neuron.get_store()
    # The voltage store is the only thing we care about in optimisation
    # So scrap rest to save memory
    return store['voltage']

if __name__ == "__main__":
    run(True)
