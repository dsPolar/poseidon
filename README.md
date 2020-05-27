# poseidon

This repository contains the code utilised in my Master's thesis "Estimating Neural Input Statistics from Whole-Cell electrophysiology"

This project explores how the population statistics of excitatory and inhibitory inputs to a neuron can
be inferred via numerical optimisation of the power spectrum of a recording of membrane potential from
that neuron.

The optimisation algorithm used is called CMA-ES.
Traces have not been included in this repository in their original smr format or the pickled npy format
used by the optimiser code in optim.py

If you have any questions about my results or my thesis then email me  david.sharp266 at gmail.com and I will try to get back to you

Navigation

Data contains code used to import traces using neo from smr to python and then save as numpy files
It also plots the traces

sim contains the majority of the actual code for the project across a number of files

The key files for recreating any results found in my thesis are

simulation.py contains the simulation wrapper code and plots synthetic traces
    it uses linf.py which contains the neuron simulation code
    which in turn uses lognormal.py and poisson.py for handling synaptic strength distribution and poisson processes for input
        note. poisson.py is deprecated as it was too slow for the optimiser but has been left in here as it is clearer than the array based
        implementation

optim.py contains code for handling the optimiser as well as code for generating some of the data behind the weirder plots like the error plane
    objective.py is used for recoding parameters from CMA-ES representation to real-world representation
    
the workhorse part of optim lies in main in a switch statement though this will likely be changed as it's ugly as hell
the first case handles running cma-es on synthetic traces
the second loaded npy files into arrays and then ran cma-es on those though any file that could be read as an array would work
the third handles mapping an error plane for a given trace file

the code is inherently set up for multithreading using pycma's EvalParallel wrapper to multiprocessing's Pool
non cma code also uses Pool.starmap for multithreading
