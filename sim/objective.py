import math
import numpy as np
import simulation as sim
from scipy import signal
from time import perf_counter

#Params is a numpy array of size (n,)
#[0:10] to true range
def recode_new(params):
    params[0] = max(0,(2 - 0) * params[0]/10)
    params[1] = max(0,(2 - 0) * params[1]/10)
    if(params.size >= 4):
        params[2] = max(0, 0.00025 + (0.00125 - 0.00025) * params[2]/10)
        params[3] = max(0, 0.00025 + (0.00125 - 0.00025) * params[3]/10)
        if(params.size >= 6):
            params[4] = max(0, 0.00025 + (0.00125 - 0.00025) * params[4]/10)
            params[5] = max(0, 0.00025 + (0.00125 - 0.00025) * params[5]/10)
            if(params.size >= 8):
                params[6] = max(0, 0.5 + (15 - 0.5) * params[6]/10)
                params[7] = max(0, 0.5 + (15 - 0.5) * params[7]/10)
    return params


#Params is a numpy array of size (n,)
#[0:10] to true range
def recode(params):
    new_params = np.copy(params)
    new_params[0] = 0 + (2 - 0) * (1 - math.cos(math.pi - new_params[0]/10))/2
    new_params[1] = 0 + (2 - 0) * (1 - math.cos(math.pi - new_params[1]/10))/2
    if(new_params.size >= 4):
        new_params[2] = 0.00025 + (0.00125 - 0.00025) * (1 - math.cos(math.pi - new_params[2]/10))/2
        new_params[3] = 0.00025 + (0.00125 - 0.00025) * (1 - math.cos(math.pi - new_params[3]/10))/2
        if(new_params.size >= 6):
            new_params[4] = 0.00025 + (0.00125 - 0.00025) * (1 - math.cos(math.pi - new_params[4]/10))/2
            new_params[5] = 0.00025 + (0.00125 - 0.00025) * (1 - math.cos(math.pi - new_params[5]/10))/2
            if(new_params.size >= 8):
                new_params[6] = 0.5 + (15 - 0.5) * (1 - math.cos(math.pi - new_params[6]/10))/2
                new_params[7] = 0.5 + (15 - 0.5) * (1 - math.cos(math.pi - new_params[7]/10))/2
    return new_params

# v_trace is voltage trace for simulation run with params, args
# true values in power dictionary
def compute_error(params, args, v_trace):
    power = args[0]
    sampling = args[14] #

    freqs, psd = signal.welch(np.asarray(v_trace), fs=sampling, scaling="spectrum")

    assert freqs.size == power['freqs'].size

    psd_errors = np.zeros(freqs.size)

    for i in range(freqs.size):
        assert power['freqs'][i] == freqs[i]
        #take log to standardise error over magnitude changes
        error = (math.log(psd[i]) - power['log_psd'][i])**2
        psd_errors[i] = error

    # Scale psd error down to match other terms
    mean_psd_error = np.mean(psd_errors)

    v_avg_calc = sum(v_trace)/len(v_trace)

    v_error = (power['v_avg'] - v_avg_calc)**2

    # Variance of simulation signal
    #var_calc = sum((i - v_avg_calc) ** 2 for i in v_trace) / len(v_trace)
    # Use integral since thats what is used for calculating the value in `power['std_true']`
    var_calc = np.trapz(psd, freqs)
    std_calc = math.sqrt(var_calc)
    std_error = (power['std_true']- std_calc)**2

    error = mean_psd_error + v_error + std_error

    # Error Proportions
    p_x = mean_psd_error / error
    p_y = v_error / error
    p_z = std_error / error

    #print(str(p_x), str(p_y), str(p_z))


    return error, mean_psd_error, v_error, std_error

def get_trace_psd(params, args):
    params = recode(params)

    v_trace = sim.run_eval(params, args)
    sampling = args[14]
    freqs, psd = signal.welch(np.asarray(v_trace), fs=sampling, scaling="spectrum")

    return v_trace, freqs, psd

def sim_eval(params, *args):
    # recode to be out of [0:10]
    #print(params)
    params = recode(params)
    #print(params)

    v_trace = sim.run_eval(params, args)
    error = compute_error(params, args, v_trace)[0]


    return error

def sim_eval_two(params, args):
    # recode to be out of [0:10]
    #print(params)
    #params = recode(params)
    #print(params)

    v_trace = sim.run_eval(params, args)
    error = compute_error(params, args, v_trace)[0]


    return error
