import time
import poisson
import lognormal
import random
rate = 0.04

start = time.time()
weight_dist = lognormal.Lognormal(0.00075, 0.00075)
weights = []
nextfire = []

for i in range(10000):
    weights.append(weight_dist.gen_num())
    nextfire.append(random.expovariate(rate))
end = time.time()

print("Time for objectlow : " + str(end-start))

start = time.time()
weight_dist_a = lognormal.Lognormal(0.00075, 0.00075)
excs = []

for i in range(10000):
    excs.append(poisson.Poisson_spikes(rate,weight_dist_a.gen_num(),0.05, 100))

end = time.time()

print("Time for object array : " + str(end-start))
