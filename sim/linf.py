import numpy as np
import lognormal
import random
import math
import plotting

# R = V/I
# G = I/V
# G = 1/R


class Linf_neuron:

    def __init__(self, exc_rate, inh_rate, exc, inh, mu_e, mu_i, sigma_e, sigma_i, r_m, c_m, delta_t, time, E_e, E_i, E_l, tau_m, tau_e, tau_i):
        #Initialiser

        #MILLI
        self.E_exc = E_e
        self.E_inh = E_i
        self.E_res = E_l

        self.r_m = r_m
        self.c_m = c_m

        #MILLI
        self.tau_m = tau_m

        self.delta_t = delta_t
        self.sampling = (1/delta_t) *1000

        self.timestep = 0

        self.exc_inputs = []
        self.inh_inputs = []

        self.exc_rate = exc_rate
        self.inh_rate = inh_rate

        self.exc_num = exc
        self.inh_num = inh

        self.g_e = 0
        self.g_i = 0
        self.g_l = 1/self.r_m

        #MILLI
        self.tau_e = tau_e
        self.tau_i = tau_i

        self.mem_v = -55

        self.time = time

        self.weight_dist_e = lognormal.Lognormal(mu_e, sigma_e)
        self.weight_dist_i = lognormal.Lognormal(mu_i, sigma_i)


        self.store = {'voltage':[],'ge':[],'gi':[], 'vars':()}

        self.store['vars'] = (self.tau_e, self.tau_i, self.E_exc, self.E_inh, self.r_m, self.tau_m, self.c_m, self.g_l, self.exc_num, self.inh_num, self.g_l, self.E_res, self.sampling, mu_e, mu_i, sigma_e, sigma_i, self.exc_rate, self.inh_rate)

        self.exc_weights = []
        self.exc_nexts = []

        self.inh_weights = []
        self.inh_nexts = []

        # Assign weights to each synapse
        # Assign a next firing time to each synapse
        # Model maintains a list of next firing times for each synapse
        for e in range(self.exc_num):
            self.exc_weights.append(self.weight_dist_e.gen_num())
            self.exc_nexts.append(random.expovariate(self.exc_rate))
        for i in range(self.inh_num):
            self.inh_weights.append(self.weight_dist_i.gen_num())
            self.inh_nexts.append(random.expovariate(self.inh_rate))

    def step(self):
        self.timestep += self.delta_t

        for e in range(self.exc_num):
            if (self.timestep >= self.exc_nexts[e]):
                self.g_e += self.exc_weights[e]
                self.exc_nexts[e] = self.timestep + random.expovariate(self.exc_rate)

        for i in range(self.inh_num):
            if (self.timestep >= self.inh_nexts[i]):
                self.g_i += self.inh_weights[i]
                self.inh_nexts[i] = self.timestep + random.expovariate(self.inh_rate)

        self.g_e += -self.g_e / self.tau_e
        self.g_i += -self.g_i / self.tau_i

        self.mem_v += ((
        (self.g_l * (self.E_res - self.mem_v))+
        (self.g_e * (self.E_exc - self.mem_v))+
        (self.g_i * (self.E_inh - self.mem_v))
        ) / self.c_m) * self.delta_t

        if self.mem_v > 0:
            self.mem_v = 0.0

        self.store['voltage'].append(self.mem_v)
        self.store['ge'].append(self.g_e)
        self.store['gi'].append(self.g_i)

    def get_mem_v(self):
        return self.mem_v

    def get_store(self):
        return self.store

    def get_timestep(self):
        return self.timestep


    def test_dist(self):
        lowest = 20.0
        highest = -20.0
        mean = 0
        runs = 100000
        tests = []
        for x in range(runs):
            y = self.weight_dist_e.gen_num()
            if(y > highest):
                highest = y
            elif(y < lowest):
                lowest = y
            tests.append(y)
        mean = sum(tests) / len(tests)
        variance = sum([((x - mean) ** 2) for x in tests]) / len(tests)
        std_dev = math.sqrt(variance)
        #plotting.hist(tests, "lognormal.png")

        print("Highest : " + str(highest) + " | Lowest : " + str(lowest) +" | Mean :" + str(mean) + " | StdDev :" + str(std_dev))
