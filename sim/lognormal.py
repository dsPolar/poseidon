import scipy.stats
import math
import plotting
import numpy as np

# Container class for a lognormal distribution
class Lognormal:
    def __init__(self, mu, sigma):
        # Lognorm Y, Normal X, exp(X) = Y
        # Mean of normal distribution
        self.mu = mu
        # Std Dev of normal distribution
        self.sigma = sigma
        #print(str(sigma))

    # generate
    def gen_num(self):
        res = (np.random.lognormal(self.mu, self.sigma, 1)[0] - 1)
        while(res < 0):
            res = (np.random.lognormal(self.mu, self.sigma, 1)[0] - 1)
        return res
