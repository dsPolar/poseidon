import math
import random

# DEPRECATED, TOO SLOW FOR LARGE SCALE POPULATIONS OR LARGE NUMBERS OF RUNS


# Class that handles generating a poisson spike train
# Interacted with simulataneously with time step increments
class Poisson_spikes:
    def __init__(self, rate, weight, delta_t, time):
        self.rate = rate
        self.weight = weight
        self.delta_t = delta_t
        self.nextfire = 0
        self.timestep = 0
        self.spikes = []
        self.time = time
        self.new_next()



    # Set a new next spike firing timestep
    def new_next(self):
        delay = random.expovariate(self.rate)
        self.nextfire = self.timestep + delay
        self.spikes.append(self.nextfire)

    def generate(self):
          temp_time = 0
          while(temp_time < self.time):
              temp_time += random.expovariate(self.rate)
              if(temp_time < self.time):
                  self.spikes.append(temp_time)

    def old_step(self):
        self.timestep += self.delta_t
        spikes = 0
        if(len(self.spikes) > 0):
            while(self.spikes[self.index] < self.timestep):
                spikes += 1
                if(self.index + 1 < len(self.spikes)):
                    self.index += 1
        return (spikes * self.weight)

    # Getter for next spike firing timestep
    def get_next(self):
        return self.nextfire

    # Increment the timestep
    # If spike at timestep after increment, return weighted spike
    def step(self):
        self.timestep += self.delta_t
        spike = 0
        while(self.timestep >= self.get_next()):
            spike += 1
            self.new_next()
        return (spike * self.weight)

    def get_spikes(self):
        return self.spikes
