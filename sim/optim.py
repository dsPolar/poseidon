import numpy as np
from scipy import signal
import simulation

import math
import plotting

import objective as ob
import multiprocessing as mp
from itertools import product
from itertools import repeat
from os import getpid
import os
import time
import random
import pickle

import cma as es
from cma.fitness_transformations import EvalParallel
def cma_mp_optimize(args, x0):
    sigma = 2.0
    options = {'verb_filenameprefix':'outcmaes/sixparam/', 'maxiter':1000}
    print("begin optimisation")
    evol = es.CMAEvolutionStrategy(x0, sigma, options)
    with EvalParallel(evol.popsize+1) as eval_all:
        while not evol.stop():
            X = evol.ask()
            evol.tell(X, eval_all(ob.sim_eval, X, args=args))
            evol.disp()
            evol.logger.add()
        evol.result_pretty()
    #evol.plot()
    return evol.result


def error_starmap_worker(e, second):
    i = second[0]
    args = second[1]
    x0 = np.array([e,i])
    x0 = ob.recode(x0)
    v_trace = simulation.run_eval(x0, args)
    error = ob.compute_error(x0, args, v_trace)[1:4]
    #print(str(getpid()) + " reporting e" + str(e) + " i" + str(i) + " error : " + str(error))
    return error

def error_starmap_worker_one(e, second):
    i = second[0]
    args = second[1]
    x0 = np.array([e,i])
    v_trace = simulation.run_eval(x0, args)
    error = ob.compute_error(x0, args, v_trace)[0]
    print(str(getpid()) + " reporting e" + str(e) + " i" + str(i) + " error : " + str(error))
    return error

def error_plot_starmap(args):
    pool = mp.Pool(7)

    error = pool.starmap(error_starmap_worker_one, product(map(lambda x: x / 10, range(0,301,5)),zip(map(lambda x: x / 10, range(0,301,5)), repeat(args))))
    np.save("starmap_error_test_300.npy", np.asarray(error))

    return error

def error_ratio_starmap(ratio, args):
    pool = mp.Pool()
    error = pool.starmap(error_starmap_worker, zip(map(lambda x: x * (ratio/250), range(1,501,1)), zip(map(lambda x: x * (1/250), range(1,501,1)), repeat(args))))
    np.save("ratio_starmap.npy", np.asarray(error))
    return error

def load_error():
    starmap_error = np.load("starmap_error_test_300.npy")

    res_array = np.reshape(starmap_error, (61,61))
    #res_array = np.delete(res_array, 0, axis=0)
    #res_array = np.delete(res_array, 0, axis=1)
    return res_array

def low_error(error_array):
    lowest = 100
    e_low = 0
    i_low = 0
    for e in range(0,11):
        for i in range(0,11):
            if error_array[e][i] < lowest:
                lowest = error_array[e][i]
                e_low = e
                i_low = i
                print(str(e_low), str(i_low), str(error_array[e][i]))
    print(str(e_low), str(i_low))

def plot_error(error_array):
    plotting.plot(error_array, None, "Error plane for nu", "nu_e coded", "nu_i coded", "nu_error_space.png")
    plotting.show(error_array, "nu_error_image_test_300.png")

def error_plane(args):
    #start = time.time()
    #error_plot_starmap(args)
    #print(str(time.time() - start))
    error_array = load_error()
    plot_error(error_array)

    #low_error(error_array)

def error_ratio(result, args):
    ratio = result[0] / result[1]
    error = error_ratio_starmap(ratio, args)
    psd_error = [a[0] for a in error]
    v_error = [a[1] for a in error]
    std_error = [a[2] for a in error]
    i = 0
    while os.path.exists("error_comps_%s.png" % i):
        i += 1
    plotting.error_comps(psd_error, v_error, std_error, "error_comps_" + str(i) + ".png")

#map(lambda x: x * (ratio/250), range(1,501,1)
def graph_recode():
    plotting.plot(list(map(lambda x: (1-math.cos(math.pi - x/10)), range(0,31,1))), None, "Recoding Function", "Optimiser Coding", "Firing Rate/kHz", "recoded.png")
    return

def graph_cmaes_error(path, name):
    errors = []
    f = open(path, "r+")
    for line in f.readlines():
        if line.find("iteration") < 0:
            #need indices of 5th and 6th space
            start = line.find(" ",line.find(" ",line.find(" ",line.find(" ", line.find(" ") + 1) +1) +1) +1)
            end = line.find(" ", start + 1)
            errors.append(float(line[start+1:end]))
    plotting.plot(errors, None, "Error over iterations", "iteration count", "Error value", name + "_iteration_error.png")

if __name__ == "__main__":
    jump = 2
    #graph_recode()
    #quit()


    if(jump == 1):
        for i in range(0,50):
            exc_rate = 0.08
            inh_rate = 0.02

            store = simulation.run(False, (exc_rate,inh_rate))

            tau_e = store['vars'][0]
            tau_i = store['vars'][1]
            E_exc = store['vars'][2]
            E_inh = store['vars'][3]
            r_m   = store['vars'][4]
            tau_m = store['vars'][5]
            c_m   = store['vars'][6]
            g_l   = store['vars'][7]
            N_e   = store['vars'][8]
            N_i   = store['vars'][9]
            g_l   = store['vars'][10]
            E_res = store['vars'][11]
            sampling = store['vars'][12]

            mu_e = store['vars'][13]
            mu_i = store['vars'][14]
            sigma_e = store['vars'][15]
            sigma_i = store['vars'][16]

            nu_e = store['vars'][17]
            nu_i = store['vars'][18]

            freqs, psd = signal.welch(np.asarray(store['voltage']), fs=sampling, scaling="spectrum")
            power = {'freqs':[],'psd':[]}
            power['freqs'] = freqs
            power['psd'] = psd
            V_avg = np.mean(np.asarray(store['voltage']))
            var_true = np.trapz(power['psd'], power['freqs'])
            std_true = math.sqrt(var_true)
            power['v_avg'] = V_avg
            power['std_true'] = std_true
            # We precompute the log PSD for the target here to save ops
            power['log_psd'] = list(map(math.log, power['psd']))

            Ge_avg = np.mean(np.asarray(store['ge']))
            Gi_avg = np.mean(np.asarray(store['gi']))
            r_eff = 1/(g_l + Ge_avg + Gi_avg)

            tau_eff = r_eff * c_m

            args=(power, tau_eff, tau_e, tau_i, r_m, V_avg, E_exc, E_inh,
                  N_e, N_i, Ge_avg, Gi_avg, g_l, E_res, sampling, tau_m,
                  mu_e, mu_i, sigma_e, sigma_i)

            #error_plane(args)
            x0 = ([20.0, 26.0])#, 10.0, 10.0, 10.0, 10.0])
            result = cma_mp_optimize(args, x0)
            f = open("six_error_out.txt", "a+")
            print(str(exc_rate), str(inh_rate), file=f)
            print(ob.recode(result[0]), file=f)
            f.close()
            f2 = open("six_error_out_error.txt", "a+")
            print(str(result[1]), file = f2)
            f2.close()

            #error_ratio(ob.recode(result[0]), args)

    elif(jump == 2):
        # Naughty traces with high (>100) error after minimisation
        bad_files = ["JP36_002.npy","JP37_001.npy","JP45_005.npy"]
        worst_files = ["JP37_001.npy"]
        best_files = ["JP37_005.npy"]
        files = ["JP37_003.npy","JP37_005.npy","JP38_003.npy","JP38_005.npy","JP41_006.npy","JP41_014.npy","JP43_005.npy"]
        for file in worst_files:
            print(file)
            # Set store to match Jon's experiment
            # Using JP37_005 atm
            # Trace available in poseidon/arrays/JP37_005.npy
            trace = np.load("../arrays/" + file)
            # Is a (N,1) array rather than (N,) so reshape here
            trace = trace.reshape((trace.shape[0],))
            start_cut = int(trace.size * 0.05)
            end_cut = int(trace.size * 0.8)
            # Trim trace to remove first 5% and final 20%
            trace = trace[start_cut:end_cut]


            # Need to set each of these
            tau_e = 3
            tau_i = 10
            E_exc = 0
            E_inh = -62
            r_m   = 100
            tau_m = 25
            c_m   = 0.09
            g_l   = 0.01
            N_e   = 800
            N_i   = 200
            #DUPLICATION OF GL IN STORE, CONCERNING?
            #g_l   = store['vars'][10]
            E_res = -90
            sampling = 20000

            # And these for two-param
            mu_e = 0.00075
            mu_i = 0.00075
            sigma_e = 0.00075
            sigma_i = 0.00075


            freqs,psd = signal.welch(trace, fs=sampling, scaling="spectrum")
            power = {'freqs':[], 'psd':[]}
            power['freqs'] = freqs
            power['psd'] = psd
            # We precompute the log PSD for the target here to save ops
            power['log_psd'] = list(map(math.log, power['psd']))
            V_avg = np.mean(trace)
            var_true = np.trapz(power['psd'], power['freqs'])
            std_true = math.sqrt(var_true)
            power['v_avg'] = V_avg
            power['std_true'] = std_true

            # Can't estimate easily from trace
            # Legacy use anyways since not using Paulo's equation
            Ge_avg = None
            Gi_avg = None
            tau_eff = None

            args=(power, tau_eff, tau_e, tau_i, r_m, V_avg, E_exc, E_inh,
                  N_e, N_i, Ge_avg, Gi_avg, g_l, E_res, sampling, tau_m,
                  mu_e, mu_i, sigma_e, sigma_i)

            x0 = [33.0,20.0]#, 1.09,17.8,-0.3,5.5, 15.0, 7.5]

            f = open("jon_error/" + file[:8] + "_two_out.txt", "a+")
            print("",file = f)
            f.close()
            for i in range(1):
                result = cma_mp_optimize(args,x0)
                f = open("jon_error/" + file[:8] + "_two_out.txt", "a+")
                #print(str(N_e) + "," + str(N_i), file = f)
                print(ob.recode(result[0]), file=f)
                print(str(result[1]))
                #print(ob.recode(result[5]), file=f)
                print(str(result[1]), file=f)
                #graph_cmaes_error("outcmaes/sixparam/fit.dat", "iteration_error/" + file[:8] + "_two_" + str(i) + "_")
                #print(args,file=f)
                fav = np.copy(result[0])
                trace, e_freqs, e_psd = ob.get_trace_psd(np.copy(fav), args)
                print(fav)
                fav_temp = np.copy(fav)


                plotting.overlap_psd(power['psd'], e_psd, freqs, file[:8] + " Target and Fitted Power Spectrum","Log Frequency/ Hz", "Log Power/ $mV^2Hz^{-1}$", file[:8] + "_two_overlap_bad.png")
                f.close()


                #retests = []
                #for y in range(500):
                #    print(str(y))
                #    retests.append(ob.sim_eval_two(ob.recode(np.copy(fav)), args))
                #mean = sum(retests) /len(retests)
                #variance = sum([((x - mean) ** 2) for x in retests]) / len(retests)
                #print(str(mean) + " " + str(variance))
                #f = open("jon_error/repeat_variance.txt", "w+")
                #for row in retests:
                #    print(row, file=f)
                #f.close()


    elif(jump == 3):
        file = "JP37_005.npy"
        print(file)
        # Set store to match Jon's experiment
        # Using JP37_005 atm
        # Trace available in poseidon/arrays/JP37_005.npy
        trace = np.load("../arrays/" + file)
        # Is a (N,1) array rather than (N,) so reshape here
        trace = trace.reshape((trace.shape[0],))
        start_cut = int(trace.size * 0.05)
        end_cut = int(trace.size * 0.8)
        # Trim trace to remove first 5% and final 20%
        trace = trace[start_cut:end_cut]


        # Need to set each of these
        tau_e = 3
        tau_i = 10
        E_exc = 0
        E_inh = -62
        r_m   = 100
        tau_m = 25
        c_m   = 0.09
        g_l   = 0.01
        N_e   = 800
        N_i   = 200
        #DUPLICATION OF GL IN STORE, CONCERNING?
        #g_l   = store['vars'][10]
        E_res = -90
        sampling = 20000

        # And these for two-param
        mu_e = 0.00075
        mu_i = 0.00075
        sigma_e = 0.00075
        sigma_i = 0.00075


        freqs,psd = signal.welch(trace, fs=sampling, scaling="spectrum")
        power = {'freqs':[], 'psd':[]}
        power['freqs'] = freqs
        power['psd'] = psd
        # We precompute the log PSD for the target here to save ops
        power['log_psd'] = list(map(math.log, power['psd']))
        V_avg = np.mean(trace)
        var_true = np.trapz(power['psd'], power['freqs'])
        std_true = math.sqrt(var_true)
        power['v_avg'] = V_avg
        power['std_true'] = std_true

        # Can't estimate easily from trace
        # Legacy use anyways since not using Paulo's equation
        Ge_avg = None
        Gi_avg = None
        tau_eff = None

        args=(power, tau_eff, tau_e, tau_i, r_m, V_avg, E_exc, E_inh,
              N_e, N_i, Ge_avg, Gi_avg, g_l, E_res, sampling, tau_m,
              mu_e, mu_i, sigma_e, sigma_i)
        error_plane(args)

    else:
        exc_rate = 0.4
        inh_rate = 0.28

        store = simulation.run(False, (exc_rate,inh_rate))

        tau_e = store['vars'][0]
        tau_i = store['vars'][1]
        E_exc = store['vars'][2]
        E_inh = store['vars'][3]
        r_m   = store['vars'][4]
        tau_m = store['vars'][5]
        c_m   = store['vars'][6]
        g_l   = store['vars'][7]
        N_e   = store['vars'][8]
        N_i   = store['vars'][9]
        g_l   = store['vars'][10]
        E_res = store['vars'][11]
        sampling = store['vars'][12]

        mu_e = store['vars'][13]
        mu_i = store['vars'][14]
        sigma_e = store['vars'][15]
        sigma_i = store['vars'][16]

        nu_e = store['vars'][17]
        nu_i = store['vars'][18]

        freqs, psd = signal.welch(np.asarray(store['voltage']), fs=sampling, scaling="spectrum")
        power = {'freqs':[],'psd':[]}
        power['freqs'] = freqs
        power['psd'] = psd
        V_avg = np.mean(np.asarray(store['voltage']))
        var_true = np.trapz(power['psd'], power['freqs'])
        std_true = math.sqrt(var_true)
        power['v_avg'] = V_avg
        power['std_true'] = std_true
        # We precompute the log PSD for the target here to save ops
        power['log_psd'] = list(map(math.log, power['psd']))

        Ge_avg = np.mean(np.asarray(store['ge']))
        Gi_avg = np.mean(np.asarray(store['gi']))
        r_eff = 1/(g_l + Ge_avg + Gi_avg)

        tau_eff = r_eff * c_m

        args=(power, tau_eff, tau_e, tau_i, r_m, V_avg, E_exc, E_inh,
              N_e, N_i, Ge_avg, Gi_avg, g_l, E_res, sampling, tau_m,
              mu_e, mu_i, sigma_e, sigma_i)

        with open("two_args.pickle", "wb") as f:
            pickle.dump(args,f)

        #error_plane(args)
        for y in range(2,9,1):
            x0 = ([y,y])
            result = cma_mp_optimize(args, x0)
            f = open("two_same_out.txt", "a+")
            print(str(exc_rate), str(inh_rate), file=f)
            print(ob.recode(result[0]), file=f)
            f.close()

        #error_ratio(ob.recode(result[0]), args)
