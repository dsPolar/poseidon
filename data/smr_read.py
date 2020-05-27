import neo
import numpy as np


def smr_read(path):
    #set up a reader for .smr files
    reader = neo.io.Spike2IO(filename=path)
    #read the block from the smr
    data = reader.read(lazy=False)[0]

    length = 0
    for i, asig in enumerate(data.segments[0].analogsignals):
        length += 1

    ch_list = list( {} for i in range(length))

    for i, asig in enumerate(data.segments[0].analogsignals):
        times = asig.times.rescale('s').magnitude
        ch = asig.annotations['channel_id']
        fs = float(asig.sampling_rate)

        #Add to channel dictionary
        ch_list[ch]['times'] = times
        ch_list[ch]['signal'] = np.array(asig)
        ch_list[ch]['fs'] = fs
    #print(path)
    #print(ch_list)
    return ch_list
